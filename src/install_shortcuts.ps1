# Makes shortcuts, each with its own icon, for the .bat launchers: on the desktop and in the Start menu (folder "Brain anatomy tools").
# Run by INSTALL_APPS.bat. A shortcut runs cmd.exe with the .bat, which (unlike the .bat itself) Windows lets you pin to the taskbar.
# Each shortcut also gets its program's app ID, the one the program declares when it starts (app_identity.py), so that Windows shows
# the running program on the pinned button instead of as a separate Python button.
# Running it again replaces the shortcuts, so after moving this folder just run it again.
param([string]$Desktop = [Environment]::GetFolderPath("Desktop"),
      [string]$StartMenu = (Join-Path ([Environment]::GetFolderPath("Programs")) "Brain anatomy tools"))

$ErrorActionPreference = "Stop"
$here = $PSScriptRoot                  # src, with the icons
$top = Split-Path -Parent $here        # the repository's top folder, with the .bat files
$apps = @(   # AppId and Icon as in APPS in app_identity.py
    @{ Name = "Leica xlef viewer"; Bat = "RUN_XLEF_VIEWER.bat";       Icon = "xlef_viewer.ico";       AppId = "JhouLab.BrainAnatomyTools.XlefViewer" },
    @{ Name = "CCF browser";       Bat = "RUN_CCF_BROWSER.bat";       Icon = "ccf_browser.ico";       AppId = "JhouLab.BrainAnatomyTools.CCFBrowser" },
    @{ Name = "ABC atlas browser"; Bat = "RUN_ABC_ATLAS_BROWSER.bat"; Icon = "abc_atlas_browser.ico"; AppId = "JhouLab.BrainAnatomyTools.ABCAtlasBrowser" }
)

# WScript.Shell can't set a shortcut's app ID (the System.AppUserModel.ID property), so that is done through the shell's property store.
Add-Type -TypeDefinition @"
using System;
using System.Runtime.InteropServices;
public static class ShortcutAppId {
    [ComImport, Guid("886D8EEB-8CF2-4446-8D02-CDBA1DBDCF99"), InterfaceType(ComInterfaceType.InterfaceIsIUnknown)]
    interface IPropertyStore {
        void GetCount(out uint count);
        void GetAt(uint index, out PropertyKey key);
        void GetValue(ref PropertyKey key, IntPtr value);
        void SetValue(ref PropertyKey key, ref PropVariant value);
        void Commit();
    }
    [StructLayout(LayoutKind.Sequential, Pack = 4)]
    struct PropertyKey { public Guid Format; public uint Id; }
    [StructLayout(LayoutKind.Explicit, Size = 24)]
    struct PropVariant { [FieldOffset(0)] public ushort Type; [FieldOffset(8)] public IntPtr Text; }
    [DllImport("shell32.dll", CharSet = CharSet.Unicode, PreserveSig = false)]
    static extern void SHGetPropertyStoreFromParsingName(string path, IntPtr bindContext, int flags, ref Guid iid,
        [MarshalAs(UnmanagedType.Interface)] out IPropertyStore store);
    public static void Set(string shortcut, string appId) {
        Guid iid = typeof(IPropertyStore).GUID;
        IPropertyStore store;
        SHGetPropertyStoreFromParsingName(shortcut, IntPtr.Zero, 2 /* GPS_READWRITE */, ref iid, out store);
        try {
            PropertyKey key = new PropertyKey { Format = new Guid("9F4C2855-9F79-4B39-A8D0-E1D42DE1D5F3"), Id = 5 };
            PropVariant value = new PropVariant { Type = 31 /* VT_LPWSTR */, Text = Marshal.StringToCoTaskMemUni(appId) };
            try { store.SetValue(ref key, ref value); store.Commit(); }
            finally { Marshal.FreeCoTaskMem(value.Text); }
        } finally { Marshal.ReleaseComObject(store); }
    }
}
"@

$cmd = Join-Path $env:SystemRoot "System32\cmd.exe"
$shell = New-Object -ComObject WScript.Shell

foreach ($folder in @($Desktop, $StartMenu)) {
    New-Item -ItemType Directory -Force $folder | Out-Null
    foreach ($app in $apps) {
        $path = Join-Path $folder ($app.Name + ".lnk")
        $link = $shell.CreateShortcut($path)
        $link.TargetPath = $cmd
        # "call" keeps the quotes intact when a file dropped on the shortcut is added after them
        $link.Arguments = '/c call "' + (Join-Path $top $app.Bat) + '"'
        $link.WorkingDirectory = $top
        $link.IconLocation = (Join-Path $here ("icons\" + $app.Icon)) + ",0"
        $link.Description = $app.Name
        $link.Save()
        [ShortcutAppId]::Set($path, $app.AppId)
    }
    Write-Host "Shortcuts made in $folder"
}

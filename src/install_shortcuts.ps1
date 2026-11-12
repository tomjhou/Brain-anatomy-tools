# Makes shortcuts, each with its own icon, for the .bat launchers: on the desktop and in the Start menu (folder "Brain anatomy tools").
# Run by INSTALL_APPS_WINDOWS.bat. A shortcut runs cmd.exe with the .bat, which (unlike the .bat itself) Windows lets you pin to the taskbar.
# Each shortcut also gets its program's app ID, the one the program declares when it starts (app_identity.py), so that Windows shows
# the running program on the pinned button instead of as a separate Python button.
# Running it again replaces the shortcuts, so after moving this folder, or when an icon has changed, just run it again.
param([string]$Desktop = [Environment]::GetFolderPath("Desktop"),
      [string]$StartMenu = (Join-Path ([Environment]::GetFolderPath("Programs")) "Brain anatomy tools"),
      [string]$IconCopies = (Join-Path $env:LOCALAPPDATA "JhouLab\icons"),
      [string]$Pinned = (Join-Path $env:APPDATA "Microsoft\Internet Explorer\Quick Launch\User Pinned\TaskBar"),
      [string]$UninstallCopy = (Join-Path $env:LOCALAPPDATA "JhouLab\uninstall"),
      [string]$RegistryKey = "HKCU:\Software\Microsoft\Windows\CurrentVersion\Uninstall\JhouLab.BrainAnatomyTools",
      [string]$Classes = "HKCU:\Software\Classes",          # (for tests: where register_file_types.ps1 writes)
      [switch]$NoStartMenuRestart)   # (for tests: leave the Start menu alone)

$ErrorActionPreference = "Stop"
$here = $PSScriptRoot                  # src, with the icons
$top = Split-Path -Parent $here        # the repository's top folder, with the .bat files
$apps = @(   # AppId and Icon as in APPS in app_identity.py
    @{ Name = "Leica xlef viewer"; Bat = "run-xlef-viewer.bat";       Icon = "xlef_viewer.ico";       AppId = "JhouLab.BrainAnatomyTools.XlefViewer" },
    @{ Name = "CCF browser";       Bat = "run-ccf-browser.bat";       Icon = "ccf_browser.ico";       AppId = "JhouLab.BrainAnatomyTools.CCFBrowser" },
    @{ Name = "ABC atlas browser"; Bat = "run-abc-atlas-browser.bat"; Icon = "abc_atlas_browser.ico"; AppId = "JhouLab.BrainAnatomyTools.ABCAtlasBrowser" },
    @{ Name = "Section assembler"; Bat = "run-assembler.bat";         Icon = "section_assembler.ico"; AppId = "JhouLab.BrainAnatomyTools.SectionAssembler" }
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
    [DllImport("shell32.dll", CharSet = CharSet.Unicode)]
    static extern void SHChangeNotify(int eventId, uint flags, string item1, IntPtr item2);
    // Tells Explorer and the taskbar that a shortcut changed (SHCNE_UPDATEITEM), or with no path that icons may have (SHCNE_ASSOCCHANGED).
    public static void Notify(string shortcut) {
        if (shortcut == null) SHChangeNotify(0x08000000, 0, null, IntPtr.Zero);
        else SHChangeNotify(0x00002000, 0x0005 /* SHCNF_PATHW */, shortcut, IntPtr.Zero);
    }
}
"@

$cmd = Join-Path $env:SystemRoot "System32\cmd.exe"
$shell = New-Object -ComObject WScript.Shell

function Write-Shortcut($path, $app) {
    $link = $shell.CreateShortcut($path)
    $link.TargetPath = $cmd
    # "call" keeps the quotes intact when a file dropped on the shortcut is added after them
    $link.Arguments = '/c call "' + (Join-Path $top $app.Bat) + '"'
    $link.WorkingDirectory = $top
    $link.IconLocation = $app.IconCopy + ",0"
    $link.Description = $app.Name
    $link.Save()
    [ShortcutAppId]::Set($path, $app.AppId)
}

# Windows caches a shortcut's icon by the icon file's path, and the taskbar goes on showing a pinned one's old picture after the file
# changes. So the shortcuts point at a copy of each icon named after its contents (e.g. %LOCALAPPDATA%\JhouLab\icons\ccf_browser-1a2b3c4d.ico):
# a changed icon gets a new path, which nothing has cached.
New-Item -ItemType Directory -Force $IconCopies | Out-Null
$newIcon = $false
foreach ($app in $apps) {
    $source = Join-Path $here ("icons\" + $app.Icon)
    $stem = [IO.Path]::GetFileNameWithoutExtension($app.Icon)
    $hash = (Get-FileHash $source -Algorithm SHA256).Hash.Substring(0, 8).ToLower()
    $app.IconCopy = Join-Path $IconCopies "$stem-$hash.ico"
    if (-not (Test-Path $app.IconCopy)) { Copy-Item $source $app.IconCopy; $newIcon = $true }
    Get-ChildItem $IconCopies -Filter "$stem-*.ico" | Where-Object FullName -ne $app.IconCopy | Remove-Item -ErrorAction SilentlyContinue
}

foreach ($folder in @($Desktop, $StartMenu)) {
    New-Item -ItemType Directory -Force $folder | Out-Null
    foreach ($app in $apps) {
        Write-Shortcut (Join-Path $folder ($app.Name + ".lnk")) $app
    }
    Write-Host "Shortcuts made in $folder"
}

# Shortcuts already pinned to the taskbar are Windows' own copies of ours: bring them up to date too (the icon, and the .bat they run,
# which may have been renamed or moved).
foreach ($app in $apps) {
    $path = Join-Path $Pinned ($app.Name + ".lnk")
    if (Test-Path $path) {
        Write-Shortcut $path $app
        [ShortcutAppId]::Notify($path)
        Write-Host "Updated the pinned $($app.Name)"
    }
}
[ShortcutAppId]::Notify($null)

# Right-click a Leica file (.xlef, .xlif, .lif) > Open with > Leica xlef viewer
& (Join-Path $here "register_file_types.ps1") -Bat (Join-Path $top "run-xlef-viewer.bat") -Icon $apps[0].IconCopy -Classes $Classes

# An entry in Settings > Apps (for this user; no administrator needed), whose Uninstall runs uninstall_apps.ps1. It runs a copy kept
# outside the repository, so that it still removes the shortcuts and the entry itself if the repository folder was deleted first.
New-Item -ItemType Directory -Force $UninstallCopy | Out-Null
$uninstaller = Join-Path $UninstallCopy "uninstall_apps.ps1"
Copy-Item (Join-Path $here "uninstall_apps.ps1") $uninstaller -Force
$size = 0
foreach ($f in [IO.Directory]::EnumerateFiles($top, "*", [IO.SearchOption]::AllDirectories)) { $size += ([IO.FileInfo]$f).Length }
New-Item -Path $RegistryKey -Force | Out-Null
$entry = @{
    DisplayName     = "Brain anatomy tools"
    DisplayIcon     = $apps[0].IconCopy
    Publisher       = "Jhou Lab"
    InstallLocation = $top
    URLInfoAbout    = "https://github.com/tomjhou/Brain-anatomy-tools"
    UninstallString = "powershell.exe -NoProfile -ExecutionPolicy Bypass -File `"$uninstaller`" -Repo `"$top`""
}
foreach ($name in $entry.Keys) { New-ItemProperty -Path $RegistryKey -Name $name -Value $entry[$name] -PropertyType String -Force | Out-Null }
foreach ($name in @("NoModify", "NoRepair")) { New-ItemProperty -Path $RegistryKey -Name $name -Value 1 -PropertyType DWord -Force | Out-Null }
New-ItemProperty -Path $RegistryKey -Name EstimatedSize -Value ([int]($size / 1KB)) -PropertyType DWord -Force | Out-Null   # (KB)
Write-Host "Listed in Settings > Apps as ""Brain anatomy tools"" (uninstall from there, or with UNINSTALL_APPS_WINDOWS.bat)"

# The Start menu keeps icons in a cache of its own, read when its process starts, so after a new icon it shows the old one until the next
# sign-in. Restarting that process (Windows starts it again at once; a Start menu open just then closes) makes it read them now.
if ($newIcon -and -not $NoStartMenuRestart) {
    Stop-Process -Name StartMenuExperienceHost -Force -ErrorAction SilentlyContinue
    Write-Host "Restarted the Start menu, so that it shows the new icons"
}

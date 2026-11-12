# Lists the xlef viewer under "Open with" for Leica files (.xlef, .xlif, .lif) when they are right-clicked in Explorer, for this user (no
# administrator needed): a program ID with the viewer's name and icon whose "open" runs run-xlef-viewer.bat with the file, and that ID added
# to each extension's OpenWithProgids. It does not make the viewer the default for those files (Windows leaves that to the user: Open with >
# Choose another app > Always), and it leaves any other program's entries for them alone.
# Run by install_shortcuts.ps1 (INSTALL_APPS_WINDOWS.bat), and by the viewer's Settings ("Open Leica files with this viewer"). -Remove takes
# it all away again (uninstall_apps.ps1 does the same itself: keep $ProgId and the extensions there in step). Running it again (after moving
# the folder) points it at the new place.
param([string]$Bat = (Join-Path (Split-Path -Parent $PSScriptRoot) "run-xlef-viewer.bat"),
      [string]$Icon = (Join-Path $PSScriptRoot "icons\xlef_viewer.ico"),
      [string]$Classes = "HKCU:\Software\Classes",       # (for tests: another key)
      [string]$ProgId = "JhouLab.BrainAnatomyTools.XlefViewer.File",
      [switch]$Remove,
      [switch]$NoNotify)                                  # (for tests: don't tell Explorer)

$ErrorActionPreference = "Stop"
$extensions = @(".xlef", ".xlif", ".lif")
$name = "Leica xlef viewer"
$key = Join-Path $Classes $ProgId

if ($Remove) {
    if (Test-Path $key) { Remove-Item $key -Recurse }
    foreach ($ext in $extensions) {
        $with = Join-Path $Classes "$ext\OpenWithProgids"
        if (Test-Path $with) { Remove-ItemProperty -Path $with -Name $ProgId -ErrorAction SilentlyContinue }
    }
    Write-Host "Removed the xlef viewer from Open with"
} else {
    if (-not (Test-Path $Bat)) { throw "There is no $Bat" }
    # cmd.exe runs the .bat, as the shortcuts do ("call" keeps the file's quotes intact); the .bat finds Python and starts the viewer
    $command = '"' + (Join-Path $env:SystemRoot "System32\cmd.exe") + '" /c call "' + $Bat + '" "%1"'
    $iconRef = $Icon + ",0"
    New-Item -Path $key -Force | Out-Null
    Set-Item -Path $key -Value "Leica experiment"
    New-ItemProperty -Path $key -Name FriendlyTypeName -Value "Leica experiment" -PropertyType String -Force | Out-Null
    New-Item -Path "$key\DefaultIcon" -Force | Out-Null
    Set-Item -Path "$key\DefaultIcon" -Value $iconRef
    # The name and icon shown in Open with: without them Windows would show cmd.exe's ("Windows Command Processor")
    New-Item -Path "$key\Application" -Force | Out-Null
    foreach ($p in @{ ApplicationName = $name; ApplicationIcon = $iconRef; ApplicationCompany = "Jhou Lab" }.GetEnumerator()) {
        New-ItemProperty -Path "$key\Application" -Name $p.Key -Value $p.Value -PropertyType String -Force | Out-Null
    }
    New-Item -Path "$key\shell\open\command" -Force | Out-Null
    New-ItemProperty -Path "$key\shell\open" -Name FriendlyAppName -Value $name -PropertyType String -Force | Out-Null
    Set-Item -Path "$key\shell\open\command" -Value $command
    foreach ($ext in $extensions) {
        $with = Join-Path $Classes "$ext\OpenWithProgids"
        if (-not (Test-Path $with)) { New-Item -Path $with -Force | Out-Null }
        New-ItemProperty -Path $with -Name $ProgId -Value "" -PropertyType String -Force | Out-Null
    }
    Write-Host "Leica files ($($extensions -join ', ')) can be opened with the xlef viewer: right-click one, Open with"
}

if (-not $NoNotify) {
    # Tells Explorer that file associations changed (SHCNE_ASSOCCHANGED), so that the menu shows it without a sign-out
    Add-Type -Namespace JhouLab -Name Assoc -MemberDefinition @"
[DllImport("shell32.dll")] public static extern void SHChangeNotify(int eventId, uint flags, System.IntPtr a, System.IntPtr b);
"@
    [JhouLab.Assoc]::SHChangeNotify(0x08000000, 0, [IntPtr]::Zero, [IntPtr]::Zero)
}

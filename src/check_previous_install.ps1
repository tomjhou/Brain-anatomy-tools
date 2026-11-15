# Run by INSTALL_APPS_WINDOWS.bat after its Downloads check, before it makes the shortcuts: are the tools installed already?
#   - from this same folder (the installer run again, by mistake or as a "repair"): says so, and what running it again does and does not
#     fix, and asks whether to go on (exit code 3: no, stop);
#   - from another folder that is still there: asks whether this folder is to be the installed copy instead (its shortcuts replacing
#     the others). Yes: leaves a note (~\.brain_anatomy_tools_switch.txt) for the programs, the first of which opened from here offers to
#     bring over the old copy's settings, atlas download and results (see previous_install.py). No: exit code 3. An old copy in the
#     Downloads folder is warned about, strongly;
#   - nowhere, or from a folder that is gone: nothing is asked.
# Where the tools are installed from: the Settings > Apps entry install_shortcuts.ps1 makes, else (made before there was one) the folder
# the desktop shortcuts run. "The same folder" is judged by the folder on disk, not by the spelling of its path (case, short names,
# links, a drive letter or the network name of the same share): a file made here is looked for there.
param([string]$Folder = (Split-Path -Parent $PSScriptRoot),       # the repository's top folder
      [string]$RegistryKey = "HKCU:\Software\Microsoft\Windows\CurrentVersion\Uninstall\JhouLab.BrainAnatomyTools",
      [string]$Desktop = [Environment]::GetFolderPath("Desktop"),
      [string]$Note = (Join-Path $HOME ".brain_anatomy_tools_switch.txt"),
      [string[]]$Downloads = @(),                                 # (for tests: the Downloads folders)
      [string]$Answer)                                            # (for tests: "yes" or "no" instead of asking)

$ErrorActionPreference = "Stop"
$Folder = [IO.Path]::GetFullPath($Folder).TrimEnd("\")

function Installed-Folder {
    try {
        $location = (Get-ItemProperty -Path $RegistryKey -Name InstallLocation -ErrorAction Stop).InstallLocation
        if ($location) { return $location }
    } catch { }
    $shell = New-Object -ComObject WScript.Shell
    foreach ($name in @("ABC atlas browser", "Leica xlef viewer", "CCF browser", "Section assembler")) {    # as in install_shortcuts.ps1
        $path = Join-Path $Desktop "$name.lnk"
        if (-not (Test-Path -LiteralPath $path)) { continue }
        if ($shell.CreateShortcut($path).Arguments -match '"([^"]+)\\run-[^"\\]+\.bat"') { return $Matches[1] }
    }
    return $null
}

function Same-Folder([string]$a, [string]$b) {
    $probe = ".install_probe_" + [guid]::NewGuid().ToString("N")
    try { New-Item -ItemType File -Path (Join-Path $a $probe) -ErrorAction Stop | Out-Null }
    catch { return ([IO.Path]::GetFullPath($a).TrimEnd("\") -ieq [IO.Path]::GetFullPath($b).TrimEnd("\")) }   # (can't write here: by name)
    try { return (Test-Path -LiteralPath (Join-Path $b $probe)) }
    finally { Remove-Item -LiteralPath (Join-Path $a $probe) -Force -ErrorAction SilentlyContinue }
}

function In-Downloads([string]$path) {
    $list = $Downloads
    if ($list.Count -eq 0) {
        try { $list += (New-Object -ComObject Shell.Application).NameSpace("shell:Downloads").Self.Path } catch { }
        $list += (Join-Path $env:USERPROFILE "Downloads")
    }
    $here = [IO.Path]::GetFullPath($path).TrimEnd("\") + "\"
    foreach ($d in $list) {
        if ($d -and $here.StartsWith(([IO.Path]::GetFullPath($d).TrimEnd("\") + "\"), [StringComparison]::OrdinalIgnoreCase)) { return $true }
    }
    return $false
}

function Ask([string]$title, [string]$text, [string]$icon, [string]$default) {
    if ($Answer) { return ($Answer -eq "yes") }
    Add-Type -AssemblyName System.Windows.Forms
    $owner = New-Object System.Windows.Forms.Form -Property @{ TopMost = $true }      # (in front of the console window)
    try {
        $button = if ($default -eq "no") { "Button2" } else { "Button1" }
        return ([System.Windows.Forms.MessageBox]::Show($owner, $text, $title, "YesNo", $icon, $button) -eq "Yes")
    } finally { $owner.Dispose() }
}

$old = Installed-Folder
if (-not $old -or -not (Test-Path -LiteralPath $old -PathType Container)) { exit 0 }      # (not installed, or from a folder now gone)
$old = [IO.Path]::GetFullPath($old).TrimEnd("\")

if (Same-Folder $Folder $old) {
    $text = "The Brain anatomy tools are already installed from this folder:`n    $Folder`n`n" +
            "Running the installer again only remakes the desktop and Start menu shortcuts, the taskbar pins, Open with for Leica files " +
            "and the entry in Settings > Apps. It does not repair the programs themselves. For those:`n" +
            "  - missing or damaged program files: the Update button in one of the programs (or download the tools again);`n" +
            "  - missing Python libraries: run python-install.bat in this folder;`n" +
            "  - a broken Python environment: delete the .venv folder here; it is rebuilt when a program next starts.`n`n" +
            "Remake the shortcuts anyway?"
    if (Ask "Already installed from this folder" $text "Information" "no") { exit 0 }
    Write-Host "Already installed from this folder: nothing was changed."
    exit 3
}

$text = "The Brain anatomy tools are already installed from another folder:`n    $old`n`nThis new folder is:`n    $Folder`n`n" +
        "Make this new folder the installed copy? Its shortcuts replace the old ones, and the first time you open one of the programs " +
        "from here, it offers to bring over the old copy's settings, its ABC atlas download and your results, so that nothing is " +
        "downloaded again. The old copy itself is left as it is."
$icon = "Question"
if (In-Downloads $old) {
    $text += "`n`nThe old copy is in your Downloads folder. Files there are often cleared out, by hand or automatically (Storage Sense), " +
             "and the tools, the ABC atlas download (up to about 60 GB) and your results kept beside them would go too. Switching to " +
             "this folder is strongly recommended, and so is moving the atlas download and your results out of Downloads (offered " +
             "when you open a program from here)."
    $icon = "Warning"
}
if (-not (Ask "Already installed from another folder" $text $icon "yes")) {
    Write-Host "Nothing was changed: the tools are still installed from $old"
    exit 3
}
Set-Content -LiteralPath $Note -Value @("previous=$old", "new=$Folder") -Encoding UTF8
Write-Host "This folder will be the installed copy (the old one, $old, is left as it is)."
exit 0

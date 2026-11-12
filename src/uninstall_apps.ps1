# Undoes INSTALL_APPS_WINDOWS.bat, and offers to remove what the programs made: run by UNINSTALL_APPS_WINDOWS.bat, or by Uninstall in Windows' Settings >
# Apps (which runs a copy of this file kept in %LOCALAPPDATA%\JhouLab\uninstall, so that it works even after the repository folder was deleted).
# Removed without asking: the desktop and Start menu shortcuts, the icon copies, the logs, the xlef viewer's thumbnails, the Settings entry
# and the xlef viewer's place under Open with for Leica files.
# Asked one by one (Enter = no): the .venv, the browsers' caches and settings in src, the xlef viewer's settings, the ABC atlas download.
# Never touched: projects, images, the ABC browser's output folders, the repository's own files, Python itself, taskbar pins (unpin those by hand).
param([string]$Repo,                                   # the repository's top folder (may be gone)
      [string]$Desktop = [Environment]::GetFolderPath("Desktop"),
      [string]$StartMenu = (Join-Path ([Environment]::GetFolderPath("Programs")) "Brain anatomy tools"),
      [string]$JhouLab = (Join-Path $env:LOCALAPPDATA "JhouLab"),
      [string]$XlefSettings = (Join-Path $HOME ".xlef_viewer_settings.json"),
      [string]$RegistryKey = "HKCU:\Software\Microsoft\Windows\CurrentVersion\Uninstall\JhouLab.BrainAnatomyTools",
      [string]$Classes = "HKCU:\Software\Classes",          # (for tests: another key)
      [string]$FileProgId = "JhouLab.BrainAnatomyTools.XlefViewer.File",     # as in register_file_types.ps1
      [string]$Answers)                               # (for tests: the answers to the questions, e.g. "ynyn", instead of asking)

$ErrorActionPreference = "Stop"
if ($Repo) { $Repo = [IO.Path]::GetFullPath($Repo).TrimEnd('\') }   # (UNINSTALL_APPS_WINDOWS.bat passes "...\Brain-anatomy-tools\.")
$Shortcuts =@("Leica xlef viewer", "CCF browser", "ABC atlas browser", "Section assembler")     # as in install_shortcuts.ps1
$asked = 0

function Size-Of($path) {
    if (-not (Test-Path -LiteralPath $path)) { return 0 }
    $item = Get-Item -LiteralPath $path -Force
    if (-not $item.PSIsContainer) { return $item.Length }
    $total = 0
    foreach ($f in [IO.Directory]::EnumerateFiles($path, "*", [IO.SearchOption]::AllDirectories)) { $total += ([IO.FileInfo]$f).Length }
    return $total
}

function Size-Text($bytes) {
    if ($bytes -ge 1GB) { return "{0:N1} GB" -f ($bytes / 1GB) }
    if ($bytes -ge 1MB) { return "{0:N0} MB" -f ($bytes / 1MB) }
    return "{0:N0} KB" -f [Math]::Ceiling($bytes / 1KB)
}

function Ask($question) {
    if ($Answers) {
        $answer = if ($script:asked -lt $Answers.Length) { $Answers[$script:asked] } else { "n" }
        $script:asked++
        Write-Host "$question [y/N] $answer"
    } else {
        $answer = Read-Host "$question [y/N]"
    }
    return "$answer".Trim().ToLower().StartsWith("y")
}

function Remove-Path($path, $what) {
    if (Test-Path -LiteralPath $path) {
        Remove-Item -LiteralPath $path -Recurse -Force
        Write-Host "  Removed $what"
    }
}

Write-Host "==== Uninstalling the Brain anatomy tools ===="
Write-Host ""
$haveRepo = $Repo -and (Test-Path -LiteralPath (Join-Path $Repo "src"))
if ($Repo -and -not $haveRepo) { Write-Host "(The folder $Repo is gone, so only what was outside it is removed.)"; Write-Host "" }

# Programs running from the .venv keep files in it (and in the caches) open: they must be closed first.
if ($haveRepo) {
    $prefix = (Join-Path $Repo ".venv").TrimEnd('\') + '\'
    $running = Get-CimInstance Win32_Process | Where-Object {
        $_.ProcessId -ne $PID -and (
            ($_.ExecutablePath -and $_.ExecutablePath.StartsWith($prefix, [StringComparison]::OrdinalIgnoreCase)) -or
            ($_.CommandLine -and $_.CommandLine.IndexOf($prefix, [StringComparison]::OrdinalIgnoreCase) -ge 0)) }
    if ($running) {
        Write-Host "Some of the tools are still running. Close them (and any window started by a run-*.bat file), then uninstall again."
        exit 1
    }
}

# --- Removed without asking: what INSTALL_APPS_WINDOWS.bat and the programs made outside the repository, all small or rebuilt as needed
Write-Host "Removing the shortcuts, icons, logs and thumbnails:"
foreach ($name in $Shortcuts) {
    $path = Join-Path $Desktop "$name.lnk"
    if (Test-Path -LiteralPath $path) {      # (only one of ours: it runs one of the run-*.bat files)
        $args_ = (New-Object -ComObject WScript.Shell).CreateShortcut($path).Arguments
        if ($args_ -match 'run[-_][a-z0-9_-]+\.bat') { Remove-Path $path "the desktop shortcut $name" }
    }
}
Remove-Path $StartMenu "the Start menu folder"
foreach ($sub in @("icons", "logs", "xlef_viewer_thumbs")) { Remove-Path (Join-Path $JhouLab $sub) "$JhouLab\$sub" }
if (Test-Path $RegistryKey) { Remove-Item $RegistryKey -Recurse; Write-Host "  Removed the entry in Settings > Apps" }
# The xlef viewer under Open with for Leica files (register_file_types.ps1: here itself, as that file may be gone with the repository)
$progKey = Join-Path $Classes $FileProgId
if (Test-Path $progKey) {
    Remove-Item $progKey -Recurse
    foreach ($ext in @(".xlef", ".xlif", ".lif")) {
        $with = Join-Path $Classes "$ext\OpenWithProgids"
        if (Test-Path $with) { Remove-ItemProperty -Path $with -Name $FileProgId -ErrorAction SilentlyContinue }
    }
    Write-Host "  Removed the xlef viewer from Open with for Leica files"
}

# --- Asked about: larger things, or ones that would take time to get back
Write-Host ""
if ($haveRepo) {
    $venv = Join-Path $Repo ".venv"
    if (Test-Path -LiteralPath $venv) {
        if (Ask "Remove the private Python environment ($venv, $(Size-Text (Size-Of $venv)))? A run-*.bat file rebuilds it (needs the internet)") {
            Remove-Path $venv "the .venv"
        }
    }
    # The ABC atlas download: where the browsers were told it is (or the default, beside the repository). Asked before the caches,
    # which hold its location.
    $atlas = Join-Path (Split-Path -Parent $Repo) "abc_atlas_cache"
    $where = Join-Path $Repo "src\cache_local\abc_atlas_cache_location.json"
    if (Test-Path -LiteralPath $where) {
        try { $atlas = (Get-Content -LiteralPath $where -Raw | ConvertFrom-Json).path } catch { }
    }
    if ($atlas -and (Test-Path -LiteralPath $atlas)) {
        if (Ask "Remove the ABC atlas download ($atlas, $(Size-Text (Size-Of $atlas)))? Getting it back means downloading it again") {
            Remove-Path $atlas "the ABC atlas download"
        }
    }
    $caches = @((Join-Path $Repo "src\cache_local"), (Join-Path $Repo "src\viewer_settings.json")) | Where-Object { Test-Path -LiteralPath $_ }
    if ($caches) {
        $size = ($caches | ForEach-Object { Size-Of $_ } | Measure-Object -Sum).Sum
        if (Ask "Remove the ABC and CCF browsers' caches and settings (src\cache_local and src\viewer_settings.json, $(Size-Text $size))? Rebuilt as needed, slowly") {
            foreach ($c in $caches) { Remove-Path $c (Split-Path -Leaf $c) }
        }
    }
}
if (Test-Path -LiteralPath $XlefSettings) {
    if (Ask "Remove the xlef viewer's settings (${XlefSettings}: its recent files and preferences)?") {
        Remove-Path $XlefSettings "the xlef viewer's settings"
    }
}

# Last: this script's own copy (PowerShell has read it already), and the JhouLab folder if nothing else is left in it (the packaged ABC
# browser keeps its own things there).
Remove-Path (Join-Path $JhouLab "uninstall") "$JhouLab\uninstall"
if ((Test-Path -LiteralPath $JhouLab) -and -not (Get-ChildItem -LiteralPath $JhouLab -Force)) { Remove-Item -LiteralPath $JhouLab }

Write-Host ""
Write-Host "Done. Still to do by hand, if you want them gone:"
Write-Host "  - buttons pinned to the taskbar: right-click each, then Unpin from taskbar"
if ($haveRepo) { Write-Host "  - the program files themselves: delete the folder $Repo" }
if (-not $Answers) { Write-Host ""; Read-Host "Press Enter to close" | Out-Null }

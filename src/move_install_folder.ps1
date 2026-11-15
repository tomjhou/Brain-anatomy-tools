# Moves the tools' folder out of Downloads, then runs the installer from its new place. Started by check_install_folder.ps1, from a copy
# of this file in %TEMP% (a folder can't be moved while anything runs from it), once the user chose to move it; it waits for the installer
# (-WaitPid, its cmd.exe) to have ended first. Moves:
#   - the folder (-Source to -Dest): on the same drive a rename, at once (tried again for a few seconds, in case something holds it for
#     a moment); on another drive a copy, checked, and only then the original deleted (a failed copy leaves the original as it was);
#   - -Extra: the atlas download and output folder that were beside it, into the folder beside -Dest (only ever on the same drive);
#   - and puts right the records that name them (the copy's own atlas location and recent output folders, and the atlas location kept
#     per user), so that nothing is looked for in Downloads again.
# Nothing is moved while a program from the folder is open. -Wrapper: the folder it was unzipped into, removed if left empty.
param([Parameter(Mandatory = $true)][string]$Source,
      [Parameter(Mandatory = $true)][string]$Dest,
      [int]$WaitPid = 0,
      [string[]]$Extra = @(),
      [string]$Wrapper = "",
      [string]$UserSettings = (Join-Path $HOME ".brain_anatomy_tools_settings.json"),      # (as previous_install.py's SETTINGS_PATH)
      [switch]$NoInstall)                                     # (for tests: no installer, Explorer or waiting for Enter after)

$ErrorActionPreference = "Stop"
try { $Host.UI.RawUI.WindowTitle = "Moving the Brain anatomy tools" } catch { }
$utf8 = New-Object Text.UTF8Encoding($false)                  # (no byte-order mark: Python's json.load won't read one)

function Fail([string]$message) {
    Write-Host ""
    Write-Host $message
    if (-not $NoInstall) {
        Write-Host ""; Read-Host "Press Enter to close" | Out-Null
        Remove-Item -LiteralPath $PSCommandPath -Force -ErrorAction SilentlyContinue     # (this copy, in %TEMP%)
    }
    exit 1
}

Write-Host "Moving the Brain anatomy tools out of Downloads, to"
Write-Host "    $Dest"
Write-Host ""
if ($WaitPid) { try { Wait-Process -Id $WaitPid -Timeout 60 -ErrorAction SilentlyContinue } catch { } }
# Not in the folder itself (Windows won't move a folder a program is in), whatever it was started in
Set-Location -LiteralPath ([IO.Path]::GetTempPath())
[IO.Directory]::SetCurrentDirectory([IO.Path]::GetTempPath())

$inUse = ""
try { $inUse = & (Join-Path $Source "src\venv_in_use.ps1") -Venv (Join-Path $Source ".venv") } catch { }
if ($inUse) {
    Fail ("A program from the folder is open, so it can't be moved. Close the Brain anatomy tools, then run INSTALL_APPS_WINDOWS.bat " +
          "again (it is still in $Source). Nothing was moved.")
}

New-Item -ItemType Directory -Force (Split-Path -Parent $Dest) | Out-Null
if ([IO.Path]::GetPathRoot($Source) -ieq [IO.Path]::GetPathRoot($Dest)) {
    $moved = $false
    for ($try = 0; $try -lt 20 -and -not $moved; $try++) {
        try { [IO.Directory]::Move($Source, $Dest); $moved = $true } catch { Start-Sleep -Milliseconds 500 }
    }
    if (-not $moved) {
        Fail ("The folder could not be moved: something has it, or a file in it, open (a program from it, a Command Prompt in it, a " +
              "file open in an editor, or a virus scan). Close that, then run INSTALL_APPS_WINDOWS.bat again (it is still in $Source). " +
              "Nothing was moved.")
    }
} else {
    Write-Host "(On another drive: copying it, then deleting the original.)"
    & robocopy $Source $Dest /E /COPY:DAT /R:2 /W:1 /NFL /NDL /NJH /NJS /NP | Out-Null
    if ($LASTEXITCODE -ge 8) {
        Remove-Item -LiteralPath $Dest -Recurse -Force -ErrorAction SilentlyContinue
        Fail "The folder could not be copied to $Dest (robocopy code $LASTEXITCODE). Nothing was moved: it is still in $Source."
    }
    try { Remove-Item -LiteralPath $Source -Recurse -Force }
    catch { Write-Host "Copied, but the old folder could not be deleted entirely ($($_.Exception.Message)): delete $Source yourself." }
}
Write-Host "  Moved the tools to $Dest"

$renames = [ordered]@{ $Source = $Dest }                      # old -> new, for the records that name them
$newParent = Split-Path -Parent $Dest
foreach ($path in $Extra) {
    if (-not $path) { continue }
    $target = Join-Path $newParent (Split-Path -Leaf $path)
    if (Test-Path -LiteralPath $target) { Write-Host "  Left $path where it is ($target is there already)"; continue }
    try {
        [IO.Directory]::Move($path, $target)
        $renames[$path] = $target
        Write-Host "  Moved $path to $target"
    } catch {
        Write-Host "  Could not move $path ($($_.Exception.Message)); it is still there (the ABC atlas browser offers it again)."
    }
}

function Renamed([string]$path) {
    if (-not $path) { return $path }
    foreach ($old in $renames.Keys) {
        if ($path -ieq $old) { return $renames[$old] }
        if ($path.StartsWith($old.TrimEnd("\") + "\", [StringComparison]::OrdinalIgnoreCase)) { return $renames[$old] + $path.Substring($old.TrimEnd("\").Length) }
    }
    return $path
}

$cache = Join-Path $Dest "src\cache_local"
$record = Join-Path $cache "abc_atlas_cache_location.json"
try {
    if (Test-Path -LiteralPath $record) {
        $saved = (Get-Content -Raw -LiteralPath $record | ConvertFrom-Json).path
        if ((Renamed $saved) -ne $saved) { [IO.File]::WriteAllText($record, (@{ path = (Renamed $saved) } | ConvertTo-Json), $utf8) }
    }
} catch { Write-Host "  (Could not update ${record}: $($_.Exception.Message))" }
$recent = Join-Path $cache "recent_output_folders.json"
try {
    if (Test-Path -LiteralPath $recent) {
        $parsed = Get-Content -Raw -LiteralPath $recent | ConvertFrom-Json     # (into a variable first: PowerShell 5 passes a JSON
        $folders = @($parsed | ForEach-Object { Renamed ([string]$_) })        #  array on as one object, which then read as one path)
        [IO.File]::WriteAllText($recent, (ConvertTo-Json -InputObject @($folders)), $utf8)
    }
} catch { Write-Host "  (Could not update ${recent}: $($_.Exception.Message))" }
try {
    if (Test-Path -LiteralPath $UserSettings) {
        $settings = Get-Content -Raw -LiteralPath $UserSettings | ConvertFrom-Json
        if ($settings.atlas_location -and (Renamed $settings.atlas_location) -ne $settings.atlas_location) {
            $settings | Add-Member -NotePropertyName atlas_location -NotePropertyValue (Renamed $settings.atlas_location) -Force
            [IO.File]::WriteAllText($UserSettings, ($settings | ConvertTo-Json), $utf8)
        }
    }
} catch { Write-Host "  (Could not update ${UserSettings}: $($_.Exception.Message))" }

if ($Wrapper -and (Test-Path -LiteralPath $Wrapper) -and -not (Get-ChildItem -LiteralPath $Wrapper -Force)) {
    Remove-Item -LiteralPath $Wrapper -Force -ErrorAction SilentlyContinue
}

Write-Host ""
Write-Host "Done. The installer now runs from the new place."
if (-not $NoInstall) {
    if (-not $env:JHOULAB_INSTALL_TEST) { Start-Process explorer.exe -ArgumentList "`"$Dest`"" }   # (the old Explorer window shows a folder now gone)
    Start-Process -FilePath (Join-Path $Dest "INSTALL_APPS_WINDOWS.bat") -WorkingDirectory $Dest
    Start-Sleep -Seconds 2
    Remove-Item -LiteralPath $PSCommandPath -Force -ErrorAction SilentlyContinue        # (this copy, in %TEMP%)
}
exit 0

# Run by INSTALL_APPS_WINDOWS.bat before it makes the shortcuts: is this folder in the Downloads folder (or a folder inside it)?
# Downloads is a poor home for the tools: it gets cleared out (by hand, or by Storage Sense), and the ABC atlas browser keeps its atlas
# download (up to about 60 GB) and its results next to this folder by default, so they would go too. Then it offers to move the folder to
# a permanent place (a BrainTools folder in the home folder, or one chosen) and install from there:
#   - Move: starts move_install_folder.ps1 (a copy of it in %TEMP%: nothing may run from this folder while it moves), which waits for the
#     installer to finish, moves the folder (and the atlas download and output folder beside it, when they are there, on the same drive:
#     an instant rename), and runs the installer again from the new place. Exit code 4: the .bat then ends at once, letting go of the folder;
#   - Cancel: says why it matters and what to do, and exit code 2 (the .bat stops, no shortcuts made).
# Exit code 0: not in Downloads. Downloads is where Windows says it is (it can be moved, e.g. to another drive), and %USERPROFILE%\Downloads too.
param([string]$Folder = (Split-Path -Parent $PSScriptRoot),   # the repository's top folder
      [string[]]$Downloads = @(),                             # (for tests: the Downloads folders to check against)
      [string]$HomeFolder = $env:USERPROFILE,                 # (for tests: where the BrainTools folder goes)
      [string]$Answer,                                        # (for tests: "move", "cancel", or a folder (as if chosen), instead of asking)
      [switch]$Inline)                                        # (for tests: move here and now, not in a window of its own; no installer after)

$ErrorActionPreference = "Stop"
$Folder = [IO.Path]::GetFullPath($Folder).TrimEnd("\")
if ($Downloads.Count -eq 0) {
    try { $Downloads += (New-Object -ComObject Shell.Application).NameSpace("shell:Downloads").Self.Path } catch { }
    $Downloads += (Join-Path $env:USERPROFILE "Downloads")
}

function WithSlash([string]$path) { return ([IO.Path]::GetFullPath($path)).TrimEnd("\") + "\" }
function In-Downloads([string]$path) {
    foreach ($d in $Downloads) { if ($d -and (WithSlash $path).StartsWith((WithSlash $d), [StringComparison]::OrdinalIgnoreCase)) { return $true } }
    return $false
}
if (-not (In-Downloads $Folder)) { exit 0 }

$name = Split-Path -Leaf $Folder
$parent = Split-Path -Parent $Folder
$Beside = @(@{ Name = "abc_atlas_cache"; What = "the ABC atlas download" }, @{ Name = "ABC_atlas_browser_output"; What = "your ABC browser results" })

function Free-Destination([string]$base) {
    # <base>\<this folder's name>, or with -2, -3, ... after it if that is taken (never into a folder that is there already)
    $dest = Join-Path $base $name
    for ($k = 2; Test-Path -LiteralPath $dest; $k++) { $dest = Join-Path $base "$name-$k" }
    return $dest
}

function Tell-Cancelled {
    Write-Host "This folder is in your Downloads folder:"
    Write-Host "    $Folder"
    Write-Host ""
    Write-Host "Files in Downloads are often cleared out, by hand or automatically by Windows, which would remove the tools"
    Write-Host "(and the ABC atlas download and your results, which are kept beside this folder). Please move the whole $name"
    Write-Host "folder somewhere permanent, for example $(Join-Path $HomeFolder 'BrainTools') (not a OneDrive or network folder), then run"
    Write-Host "INSTALL_APPS_WINDOWS.bat again from there."
    if (Test-Path -LiteralPath (Join-Path $parent "abc_atlas_cache")) {
        Write-Host ""
        Write-Host "Move the abc_atlas_cache folder next to it as well, into the same new folder, so the ABC atlas browser finds"
        Write-Host "its download there."
    }
    Write-Host ""
    Write-Host "No shortcuts were made."
}

function Ask-Where([string]$dest, [string[]]$moving, [string[]]$staying) {
    # "move", "choose" or "cancel"
    if ($Answer) { if ($Answer -eq "move" -or $Answer -eq "cancel") { return $Answer } else { return "choose" } }
    Add-Type -AssemblyName System.Windows.Forms
    # sharp on a scaled screen, and in the system's own look and font, not Windows 95's
    try { Add-Type -MemberDefinition '[DllImport("user32.dll")] public static extern bool SetProcessDPIAware();' -Name Dpi -Namespace JhouLab
          [void][JhouLab.Dpi]::SetProcessDPIAware() } catch { }
    [System.Windows.Forms.Application]::EnableVisualStyles()
    $form = New-Object System.Windows.Forms.Form -Property @{
        Text = "The tools are in your Downloads folder"; TopMost = $true; StartPosition = "CenterScreen"; FormBorderStyle = "FixedDialog"
        MaximizeBox = $false; MinimizeBox = $false; AutoSize = $true; AutoSizeMode = "GrowAndShrink"; Padding = New-Object System.Windows.Forms.Padding(16)
        Font = New-Object System.Drawing.Font("Segoe UI", 9); AutoScaleMode = "Dpi" }
    $text = "This folder is in your Downloads folder:`n    $Folder`n`nFiles in Downloads are often cleared out, by hand or automatically " +
            "(Storage Sense), which would remove the tools, and the ABC atlas download (up to about 60 GB) and your results, which are " +
            "kept beside them.`n`nMove the folder here, and install from there?`n    $dest"
    if ($moving) { $text += "`n`nMoved along with it (instant: nothing is copied):`n" + (($moving | ForEach-Object { "    $_" }) -join "`n") }
    if ($staying) { $text += "`n`nLeft where it is (offered again when you open the ABC atlas browser):`n" + (($staying | ForEach-Object { "    $_" }) -join "`n") }
    $scale = 1.0                                                                  # (the screen's scaling: the width below is in its pixels)
    try { $g = [System.Drawing.Graphics]::FromHwnd([IntPtr]::Zero); $scale = $g.DpiX / 96; $g.Dispose() } catch { }
    $label = New-Object System.Windows.Forms.Label -Property @{ Text = $text; AutoSize = $true; MaximumSize = New-Object System.Drawing.Size([int](700 * $scale), 0) }
    $buttons = New-Object System.Windows.Forms.FlowLayoutPanel -Property @{ AutoSize = $true; FlowDirection = "RightToLeft"; Anchor = "Right"
                                                                              Margin = New-Object System.Windows.Forms.Padding(0, 14, 0, 0) }
    $result = @{ Value = "cancel" }
    foreach ($b in @(@("Cancel", "cancel"), @("Choose another place...", "choose"), @("Move and install", "move"))) {
        $button = New-Object System.Windows.Forms.Button -Property @{ Text = $b[0]; AutoSize = $true; Padding = New-Object System.Windows.Forms.Padding(8, 2, 8, 2) }
        $value = $b[1]
        $button.Add_Click({ $result.Value = $value; $form.Close() }.GetNewClosure())
        $buttons.Controls.Add($button)
        if ($value -eq "move") { $form.AcceptButton = $button }
        if ($value -eq "cancel") { $form.CancelButton = $button }
    }
    $layout = New-Object System.Windows.Forms.FlowLayoutPanel -Property @{ AutoSize = $true; FlowDirection = "TopDown" }
    $layout.Controls.Add($label)
    $layout.Controls.Add($buttons)
    $form.Controls.Add($layout)
    [void]$form.ShowDialog()
    $form.Dispose()
    return $result.Value
}

function Choose-Base {
    if ($Answer -and $Answer -ne "move" -and $Answer -ne "cancel") { return $Answer }
    Add-Type -AssemblyName System.Windows.Forms
    $dialog = New-Object System.Windows.Forms.FolderBrowserDialog -Property @{
        Description = "Choose where to put the $name folder (not Downloads, OneDrive or a network folder)"; SelectedPath = $HomeFolder; ShowNewFolderButton = $true }
    $owner = New-Object System.Windows.Forms.Form -Property @{ TopMost = $true }
    try { if ($dialog.ShowDialog($owner) -eq "OK") { return $dialog.SelectedPath } else { return $null } } finally { $owner.Dispose() }
}

$base = Join-Path $HomeFolder "BrainTools"
while ($true) {
    $dest = Free-Destination $base
    $sameDrive = [IO.Path]::GetPathRoot($dest) -ieq [IO.Path]::GetPathRoot($Folder)
    $moving = @(); $staying = @(); $extra = @()
    foreach ($b in $Beside) {
        $path = Join-Path $parent $b.Name
        if (-not (Test-Path -LiteralPath $path -PathType Container)) { continue }
        if ($sameDrive -and -not (Test-Path -LiteralPath (Join-Path (Split-Path -Parent $dest) $b.Name))) {
            $moving += "$path  ($($b.What))"; $extra += $path
        } else { $staying += "$path  ($($b.What))" }
    }
    $choice = Ask-Where $dest $moving $staying
    if ($choice -eq "cancel") { Tell-Cancelled; exit 2 }
    if ($choice -eq "move") { break }
    $picked = Choose-Base
    if (-not $picked) { continue }                                    # (the folder picker cancelled: ask again)
    if ((In-Downloads $picked) -or (WithSlash $picked).StartsWith((WithSlash $Folder), [StringComparison]::OrdinalIgnoreCase)) {
        if ($Answer) { Tell-Cancelled; exit 2 }
        [void][System.Windows.Forms.MessageBox]::Show("That is in Downloads, or in this folder itself. Choose another place.", "Choose another place")
        continue
    }
    $base = $picked
    if ($Answer) { $Answer = "move" }                                 # (tests: then the move, to the folder given)
}

# (an unzipped folder is often inside a folder of the same name, in Downloads: that one goes too, after, if nothing else is left in it)
$wrapper = if ($Downloads | Where-Object { $_ -and (WithSlash $_) -ieq (WithSlash $parent) }) { "" } else { $parent }

$helper = Join-Path $PSScriptRoot "move_install_folder.ps1"
if ($Inline) {
    & $helper -Source $Folder -Dest $dest -WaitPid 0 -Extra $extra -Wrapper $wrapper -NoInstall
    exit $LASTEXITCODE
}
$installerPid = (Get-CimInstance Win32_Process -Filter "ProcessId = $PID").ParentProcessId      # (the .bat's cmd.exe, to wait for)
$copy = Join-Path $env:TEMP ("JhouLab_move_" + [guid]::NewGuid().ToString("N").Substring(0, 8) + ".ps1")
Copy-Item -LiteralPath $helper -Destination $copy
function Quote([string]$s) { return "'" + ($s -replace "'", "''") + "'" }
$command = "& $(Quote $copy) -Source $(Quote $Folder) -Dest $(Quote $dest) -WaitPid $installerPid -Wrapper $(Quote $wrapper)"
if ($extra) { $command += " -Extra @(" + (($extra | ForEach-Object { Quote $_ }) -join ",") + ")" }
$encoded = [Convert]::ToBase64String([Text.Encoding]::Unicode.GetBytes($command))     # (encoded: paths reach it exactly, whatever is in them)
# (started in %TEMP%: it would otherwise be in this folder, as this script is, and Windows won't move a folder a program is in)
Start-Process powershell -WorkingDirectory $env:TEMP -ArgumentList @("-NoProfile", "-ExecutionPolicy", "Bypass", "-EncodedCommand", $encoded)
Write-Host "Moving the folder to $dest (in a window of its own); the installer then runs from there."
exit 4

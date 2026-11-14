# Run by INSTALL_APPS_WINDOWS.bat before it makes the shortcuts: is this folder in the Downloads folder (or a folder inside it)? Then it
# says so, asks for the folder to be moved somewhere permanent and the installer run again from there, and exits with code 2 (the .bat
# then stops). Downloads is a poor home for the tools: it gets cleared out (by hand, or by Storage Sense), and the ABC atlas browser
# keeps its atlas download (up to about 60 GB) next to this folder by default, so it would go too.
# Downloads is where Windows says it is (it can be moved, e.g. to another drive), and also the usual %USERPROFILE%\Downloads.
param([string]$Folder = (Split-Path -Parent $PSScriptRoot),   # the repository's top folder
      [string[]]$Downloads = @())                             # (for tests: the Downloads folders to check against)

if ($Downloads.Count -eq 0) {
    try { $Downloads += (New-Object -ComObject Shell.Application).NameSpace("shell:Downloads").Self.Path } catch { }
    $Downloads += (Join-Path $env:USERPROFILE "Downloads")
}

function WithSlash([string]$path) { return ([IO.Path]::GetFullPath($path)).TrimEnd("\") + "\" }

$here = WithSlash $Folder
foreach ($d in $Downloads) {
    if (-not $d) { continue }
    if ($here.StartsWith((WithSlash $d), [StringComparison]::OrdinalIgnoreCase)) {
        $name = Split-Path -Leaf $Folder
        $parent = Split-Path -Parent $Folder
        $suggest = Join-Path $env:USERPROFILE "BrainTools"
        Write-Host "This folder is in your Downloads folder:"
        Write-Host "    $Folder"
        Write-Host ""
        Write-Host "Files in Downloads are often cleared out, by hand or automatically by Windows, which would remove the tools"
        Write-Host "(and the ABC atlas download, which is kept beside this folder). Please move the whole $name folder somewhere"
        Write-Host "permanent, for example $suggest (not a OneDrive or network folder), then run INSTALL_APPS_WINDOWS.bat"
        Write-Host "again from there."
        if (Test-Path -LiteralPath (Join-Path $parent "abc_atlas_cache")) {
            Write-Host ""
            Write-Host "Move the abc_atlas_cache folder next to it as well, into the same new folder, so the ABC atlas browser finds"
            Write-Host "its download there."
        }
        Write-Host ""
        Write-Host "No shortcuts were made."
        exit 2
    }
}
exit 0

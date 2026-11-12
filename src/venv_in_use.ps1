# Says which programs are running from a virtual environment: writes the process ids (separated by commas) of the programs whose
# executable is inside the folder given as -Venv, and nothing if there are none. Called by ensure_python.bat before it deletes a
# virtual environment to rebuild it: a program that is running from it (the atlas browser, the CCF browser, ...) would crash when its
# files went, and deleting the folder while files are open can leave it half-deleted.
# (A venv's python.exe is a small launcher that starts the real Python as a child; the launcher's own executable is in the venv, and
# that is the process found here. pythonw.exe and any other program installed there, pip.exe for one, are found the same way.)
param([Parameter(Mandatory = $true)][string]$Venv)

$ErrorActionPreference = "SilentlyContinue"
$folder = [IO.Path]::GetFullPath($Venv).TrimEnd('\') + '\'
$ids = @(Get-CimInstance Win32_Process | Where-Object {
    $_.ExecutablePath -and $_.ExecutablePath.StartsWith($folder, [StringComparison]::OrdinalIgnoreCase) -and $_.ProcessId -ne $PID
} | ForEach-Object { $_.ProcessId })
if ($ids.Count -gt 0) { $ids -join ',' }

@echo off
REM Makes desktop and Start menu shortcuts, each with its own icon, for the xlef viewer, CCF browser, ABC atlas browser
REM and section assembler
REM (double-click this file in Windows Explorer). Unlike the run-*.bat files, the shortcuts can be pinned to the taskbar.
REM Run it again after moving this folder, or when an icon has changed (it updates pinned shortcuts too). See src\install_shortcuts.ps1.
setlocal
cd /d "%~dp0"
title Install app shortcuts
echo ==== Making shortcuts for the Brain anatomy tools ====
echo.

REM A folder from a downloaded zip has every file marked as "from the internet", and PowerShell then refuses to run the unsigned .ps1 scripts below on
REM a computer whose policy is RemoteSigned ("is not digitally signed"). So first clear the mark from the whole folder (an inline command, which that
REM policy allows; not the .venv or .git folders, which are big and were never downloaded). Harmless where there is no mark.
set "UNBLOCK_ROOT=%~dp0"
powershell -NoProfile -Command "Get-ChildItem -LiteralPath $env:UNBLOCK_ROOT -Recurse -Force -File -ErrorAction SilentlyContinue | Where-Object { $_.FullName -notlike '*\.venv\*' -and $_.FullName -notlike '*\.git\*' } | Unblock-File -ErrorAction SilentlyContinue" >nul 2>&1

REM Not from the Downloads folder (it gets cleared out): there, it offers to move this folder to a permanent place and install from
REM there. See src\check_install_folder.ps1. Code 4: the move has started, in a window of its own, which waits for this one to end (so
REM that nothing holds the folder), then runs the installer from the new place: so end now. Code 2: cancelled (no shortcuts made). If
REM the check itself can't run, the shortcuts are made anyway. (4 is looked at first: "errorlevel 2" means 2 or more.)
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0src\check_install_folder.ps1"
if errorlevel 4 exit /b 0
if errorlevel 2 (
  echo.
  pause
  exit /b 1
)

REM Installed already? From this folder (run again): says what that does and does not fix, and asks whether to go on. From another
REM folder: asks whether this one is to be the installed copy. Its "no", code 3, stops here. See src\check_previous_install.ps1.
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0src\check_previous_install.ps1"
if errorlevel 3 (
  echo.
  pause
  exit /b 0
)

powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0src\install_shortcuts.ps1"
if errorlevel 1 (
  echo.
  echo The shortcuts could not be made ^(see above^).
  pause
  exit /b 1
)

echo.
echo Done. The desktop shortcuts can be opened directly or dragged to the taskbar. They are also in the Start menu,
echo under "Brain anatomy tools". To open a Leica file ^(.xlef, .xlif, .lif^) in the viewer, right-click it, then
echo Open with ^> Leica xlef viewer.
echo.
pause

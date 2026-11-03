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

powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0src\install_shortcuts.ps1"
if errorlevel 1 (
  echo.
  echo The shortcuts could not be made ^(see above^).
  pause
  exit /b 1
)

echo.
echo Done. The desktop shortcuts can be opened directly or dragged to the taskbar. They are also in the Start menu,
echo under "Brain anatomy tools".
echo.
pause

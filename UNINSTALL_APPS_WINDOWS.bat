@echo off
REM Removes what INSTALL_APPS_WINDOWS.bat made (shortcuts, icons, the entry in Settings > Apps), and offers to remove the .venv, caches,
REM settings and the ABC atlas download. Projects and outputs are never touched. Also run by Uninstall in Settings > Apps.
REM See src\uninstall_apps.ps1.
setlocal
title Uninstall the Brain anatomy tools
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0src\uninstall_apps.ps1" -Repo "%~dp0."

@echo off
REM Starts the CCF browser (double-click this file in Windows Explorer).
REM Finds Python (installing Python 3.13 first if there is none): see ensure_python.bat.
setlocal
cd /d "%~dp0src"
title CCF browser
echo ==== CCF browser ====
echo.

call "%~dp0src\ensure_python.bat"

if not defined PYTHON_CMD (
  echo Python could not be found or installed automatically.
  echo Please install it from https://www.python.org/downloads/ ^(tick "Add python.exe to PATH"^)
  echo and run this file again.
  echo.
  pause
  exit /b 1
)

REM Once everything is in place, start it with pythonw.exe (no console window), and close this window. Its output goes to a log file
REM (%LOCALAPPDATA%\JhouLab\logs\ccf_browser.log), and errors are shown in a message box: see app_identity.py.
set "PYTHONW_CMD=%PYTHON_CMD:python.exe=pythonw.exe%"
REM (Only when its libraries are installed: the browser installs any missing itself, which needs the console.)
"%PYTHON_CMD%" -c "import numpy, pandas, PIL" >nul 2>&1
if errorlevel 1 set "PYTHONW_CMD="
if defined PYTHONW_CMD if exist "%PYTHONW_CMD%" (
  start "" "%PYTHONW_CMD%" CCF_browser.py %*
  exit /b 0
)

"%PYTHON_CMD%" CCF_browser.py %*
if errorlevel 1 (
  echo.
  echo The CCF browser stopped with an error ^(see above^). If a library is missing, run python-install.bat.
  pause
)

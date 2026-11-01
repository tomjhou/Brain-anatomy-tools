@echo off
REM Starts the Leica .xlef viewer (double-click this file in Windows Explorer, or drop an .xlef file on it).
REM Finds Python (installing Python 3.13 first if there is none): see ensure_python.bat.
REM Missing libraries (numpy, tifffile, pillow) are installed automatically on first run.
setlocal
cd /d "%~dp0src"
title Leica .xlef viewer
echo ==== Leica .xlef viewer ====
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

"%PYTHON_CMD%" -c "import numpy, tifffile, PIL, imagecodecs" >nul 2>&1
if errorlevel 1 (
  echo Installing the libraries the viewer needs ^(numpy, tifffile, pillow, imagecodecs^)...
  echo.
  "%PYTHON_CMD%" -m pip install numpy tifffile pillow imagecodecs
  if errorlevel 1 (
    echo.
    echo The libraries could not be installed ^(see above^).
    pause
    exit /b 1
  )
  echo.
)

REM tkinterdnd2 (dragging .xlef files onto the window) is optional: the viewer works without it.
"%PYTHON_CMD%" -c "import tkinterdnd2" >nul 2>&1
if errorlevel 1 "%PYTHON_CMD%" -m pip install tkinterdnd2 >nul 2>&1

REM Once everything is in place, start it with pythonw.exe (no console window), and close this window. Its output goes to a log file
REM (%LOCALAPPDATA%\JhouLab\logs\xlef_viewer.log), and errors are shown in a message box: see app_identity.py.
set "PYTHONW_CMD=%PYTHON_CMD:python.exe=pythonw.exe%"
if exist "%PYTHONW_CMD%" (
  start "" "%PYTHONW_CMD%" xlef_viewer.py %*
  exit /b 0
)

"%PYTHON_CMD%" xlef_viewer.py %*
if errorlevel 1 (
  echo.
  echo The viewer stopped with an error ^(see above^).
  pause
)

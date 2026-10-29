@echo off
REM Starts the Leica .xlef viewer (double-click this file in Windows Explorer, or drop an .xlef file on it).
REM Finds Python (installing Python 3.13 first if there is none): see ensure_python.bat.
REM Missing libraries (numpy, tifffile, pillow) are installed automatically on first run.
setlocal
cd /d "%~dp0"
title Leica .xlef viewer
echo ==== Leica .xlef viewer ====
echo.

call "%~dp0ensure_python.bat"

if not defined PYTHON_CMD (
  echo Python could not be found or installed automatically.
  echo Please install it from https://www.python.org/downloads/ ^(tick "Add python.exe to PATH"^)
  echo and run this file again.
  echo.
  pause
  exit /b 1
)

"%PYTHON_CMD%" -c "import numpy, tifffile, PIL" >nul 2>&1
if errorlevel 1 (
  echo Installing the libraries the viewer needs ^(numpy, tifffile, pillow^)...
  echo.
  "%PYTHON_CMD%" -m pip install numpy tifffile pillow
  if errorlevel 1 (
    echo.
    echo The libraries could not be installed ^(see above^).
    pause
    exit /b 1
  )
  echo.
)

"%PYTHON_CMD%" xlef_viewer.py %*
if errorlevel 1 (
  echo.
  echo The viewer stopped with an error ^(see above^).
  pause
)

@echo off
REM Opens the window that adds the surgery viewer's injection sites to an ABC browser 3D page as white balls
REM (double-click this file in Windows Explorer). See src\add_injections_3d.py.
REM Finds Python (installing Python 3.13 first if there is none): see ensure_python.bat.
setlocal
cd /d "%~dp0..\src"
title Injection sites in a 3D page
echo ==== Injection sites in a 3D page ====
echo.

call "%~dp0..\src\ensure_python.bat"

if not defined PYTHON_CMD (
  echo Python could not be found or installed automatically.
  echo Please install it from https://www.python.org/downloads/ ^(tick "Add python.exe to PATH"^)
  echo and run this file again.
  echo.
  pause
  exit /b 1
)

"%PYTHON_CMD%" add_injections_3d.py %*
if errorlevel 1 (
  echo.
  echo It stopped with an error ^(see above^). If a library is missing, run python-install.bat.
  pause
)

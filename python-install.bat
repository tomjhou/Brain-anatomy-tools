@echo off
REM Installs the Python libraries used by ABC_atlas_browser.py and CCF_browser.py into the private virtual environment (the .venv folder)
REM that the run-*.bat files use; see ensure_python.bat, which finds Python and makes the environment first if needed.
setlocal
cd /d "%~dp0"
call "%~dp0src\ensure_python.bat"
if not defined PYTHON_CMD (
  echo Python could not be found or installed automatically.
  echo Please install it from https://www.python.org/downloads/ ^(tick "Add python.exe to PATH"^)
  echo and run this file again.
  echo.
  pause
  exit /b 1
)
echo Installing into %PYTHON_CMD%
echo.

REM Required
"%PYTHON_CMD%" "%~dp0src\keep_awake.py" install numpy
"%PYTHON_CMD%" "%~dp0src\keep_awake.py" install pandas
"%PYTHON_CMD%" "%~dp0src\keep_awake.py" install matplotlib
"%PYTHON_CMD%" "%~dp0src\keep_awake.py" install pillow
"%PYTHON_CMD%" "%~dp0src\keep_awake.py" install anndata
"%PYTHON_CMD%" "%~dp0src\keep_awake.py" install scanpy
"%PYTHON_CMD%" "%~dp0src\keep_awake.py" install "abc_atlas_access[notebooks] @ https://github.com/alleninstitute/abc_atlas_access/archive/refs/heads/main.zip"

REM Optional: Leiden clustering of new UMAP runs
"%PYTHON_CMD%" "%~dp0src\keep_awake.py" install igraph

REM Parquet file support
"%PYTHON_CMD%" "%~dp0src\keep_awake.py" install pyarrow

REM Interactive 3D view in the browser
"%PYTHON_CMD%" "%~dp0src\keep_awake.py" install plotly

REM Drag-and-drop of CSV files into the 3D viewer and CCF browser windows
"%PYTHON_CMD%" "%~dp0src\keep_awake.py" install tkinterdnd2

REM CCF browser "Generate meshes" (region surfaces for MBF / NeuroInfo atlases)
"%PYTHON_CMD%" "%~dp0src\keep_awake.py" install vtk

echo.
echo Done.
pause

@echo off
setlocal
cd /d %~dp0

set "PYTHON=C:\Users\yhu39\AppData\Local\anaconda3\envs\omicsone\python.exe"

if not exist "%PYTHON%" (
    echo Cannot find Python environment:
    echo %PYTHON%
    echo Update PYTHON in run.bat or create the omicsone conda environment.
    pause
    exit /b 1
)

echo Starting OmicsOne Streamlit on http://localhost:8501
"%PYTHON%" -m omicsone --server.port 8501

pause

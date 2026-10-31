@echo off
rem windowcalc -- entry point for the fluctuation-analysis skill (Windows)
rem Usage: windowcalc.cmd rate --input counts.csv --nts 3e8
setlocal
set "SKILL_DIR=%~dp0"
set "VENV=%USERPROFILE%\.local\share\fluctuation-analysis\venv"
set "PYTHONPATH=%SKILL_DIR%"
set "PYTHONUTF8=1"

if exist "%VENV%\Scripts\python.exe" goto usevenv
where py >nul 2>nul
if %errorlevel%==0 goto usepy
where python >nul 2>nul
if %errorlevel%==0 goto usepython
echo [windowcalc] Python 3 not found. Install Python 3.10+ or run install.ps1 1>&2
exit /b 2

:usevenv
"%VENV%\Scripts\python.exe" -m windowcalc %*
exit /b %errorlevel%

:usepy
py -3 -m windowcalc %*
exit /b %errorlevel%

:usepython
python -m windowcalc %*
exit /b %errorlevel%

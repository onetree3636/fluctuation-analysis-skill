@echo off
rem ============================================================
rem  run.cmd -- Windows entry point for the fluctuation-analysis skill
rem
rem  Finds a usable Python 3 (PATH -> py launcher -> Miniconda/Anaconda
rem  default install dirs -> python.org default dirs), then runs run.py,
rem  which installs numpy/pyyaml into the skill's own .venv on first use.
rem
rem  Usage, from anywhere:
rem      "C:\path\to\fluctuation-analysis\run.cmd" rate --input counts.csv --nts 3e8
rem  From the skill folder:
rem      run.cmd rate --input counts.csv --nts 3e8
rem
rem  Why not "py -3 run.py": py.exe ships with python.org only. Miniconda and
rem  Anaconda users usually have no py.exe, and a terminal that conda never
rem  initialised (e.g. the shell an agent spawns) has no conda python on PATH.
rem ============================================================
setlocal EnableExtensions
set "SKILL_DIR=%~dp0"
set "PYTHONUTF8=1"
set "PYEXE="
set "PYMODE=exe"

if defined WINDOWCALC_PYTHON (
    if exist "%WINDOWCALC_PYTHON%" (
        set "PYEXE=%WINDOWCALC_PYTHON%"
        goto :run
    )
    echo [windowcalc] WINDOWCALC_PYTHON does not exist: %WINDOWCALC_PYTHON% 1>&2
    exit /b 2
)

rem 1) python on PATH (conda activate'd terminal, or python.org install)
rem    Skip the Windows Store stub under %LOCALAPPDATA%\Microsoft\WindowsApps:
rem    it appears in "where python" on many machines but only opens the Store.
for /f "delims=" %%P in ('where python 2^>nul') do (
    if /i not "%%~dpP"=="%LOCALAPPDATA%\Microsoft\WindowsApps\" (
        set "PYEXE=%%P"
        goto :run
    )
)

rem 2) py launcher (python.org only)
where py >nul 2>nul
if not errorlevel 1 (
    set "PYMODE=py"
    goto :run
)

rem 3) common conda and python.org install locations
for %%D in (
    "%USERPROFILE%\miniconda3"
    "%USERPROFILE%\anaconda3"
    "%LOCALAPPDATA%\miniconda3"
    "%LOCALAPPDATA%\anaconda3"
    "%ProgramData%\miniconda3"
    "%ProgramData%\anaconda3"
    "%USERPROFILE%\AppData\Local\Programs\Python\Python313"
    "%USERPROFILE%\AppData\Local\Programs\Python\Python312"
    "%USERPROFILE%\AppData\Local\Programs\Python\Python311"
    "%USERPROFILE%\AppData\Local\Programs\Python\Python310"
) do (
    if exist "%%~D\python.exe" (
        set "PYEXE=%%~D\python.exe"
        goto :run
    )
)

echo. 1>&2
echo [windowcalc] No Python 3 found. 1>&2
echo   Miniconda / Anaconda users: open "Anaconda Prompt" and run this command there, 1>&2
echo   or point at the interpreter explicitly: 1>&2
echo     set WINDOWCALC_PYTHON=C:\Users\YOU\miniconda3\python.exe 1>&2
echo     run.cmd ... 1>&2
echo   No Python at all? https://www.python.org/downloads/ 1>&2
exit /b 2

:run
if "%PYMODE%"=="py" (
    py -3 "%SKILL_DIR%run.py" %*
) else (
    "%PYEXE%" "%SKILL_DIR%run.py" %*
)
exit /b %errorlevel%

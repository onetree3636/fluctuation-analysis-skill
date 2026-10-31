# windowcalc installer for Windows (creates a venv with numpy + pyyaml)
# Usage:  powershell -ExecutionPolicy Bypass -File install.ps1
$ErrorActionPreference = "Stop"
$SkillDir = $PSScriptRoot
$Venv = Join-Path $env:USERPROFILE ".local\share\fluctuation-analysis\venv"
$Py = Join-Path $Venv "Scripts\python.exe"

if (-not (Test-Path $Py)) {
    Write-Host "Creating venv: $Venv"
    if (Get-Command py -ErrorAction SilentlyContinue) {
        py -3 -m venv $Venv
    } elseif (Get-Command python -ErrorAction SilentlyContinue) {
        python -m venv $Venv
    } else {
        Write-Host "ERROR: Python 3 not found. Install Python 3.10+ (tick 'Add to PATH') and retry."
        exit 2
    }
}

& $Py -m pip install --quiet --upgrade pip
& $Py -m pip install --quiet numpy pyyaml

$Cmd = Join-Path $SkillDir "windowcalc.cmd"
& $Cmd --help | Out-Null
if ($LASTEXITCODE -eq 0) {
    Write-Host "OK: use `"$Cmd`" (or add the skill folder to PATH)"
    Write-Host "   e.g. & `"$Cmd`" rate --input counts.csv --nts 3e8"
} else {
    Write-Host "FAILED: $Cmd --help returned $LASTEXITCODE"
    exit 1
}

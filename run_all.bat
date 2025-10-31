@echo off
setlocal enabledelayedexpansion
set PYTHONIOENCODING=utf-8
set MPLBACKEND=Agg

REM ============================================================
REM RH_Operator_Model - Full pipeline runner (Windows, CMD)
REM Steps:
REM   [1] Create/activate venv
REM   [2] Install requirements
REM   [3] Main pipeline (SCF, eigenvalues, zeta-raw)
REM   [4] Rescaling (BK/OT)
REM   [5] Paper summary bundle (JSON/CSV)
REM   [6] Diagnostics (SCF conv., dense Ψ, Weyl, λ-sweep)
REM   [7] Publication plots (use diagnostics; overwrite)
REM   [8] Optional: standalone SCF plot
REM ============================================================

REM --- move to repo root (this script’s folder) ---
cd /d "%~dp0"

echo [1/8] Creating/activating venv...
if not exist ".venv" (
  py -3 -m venv .venv
  if errorlevel 1 goto :fail
)
if exist ".venv\Scripts\activate.bat" (
  call ".venv\Scripts\activate.bat"
) else (
  call ".venv\Scripts\activate"
)
if errorlevel 1 goto :fail

echo [2/8] Installing requirements...
python -m pip install --upgrade pip >nul
python -m pip install -r requirements.txt
if errorlevel 1 goto :fail

echo.
echo [3/8] Running main pipeline (SCF, spectral checks, zeta-raw)...
python run_pipeline.py
if errorlevel 1 goto :fail

echo.
echo [4/8] Rescaling (BK / optimal transport)...
python -m src.asymptotic_rescaling
if errorlevel 1 goto :fail

echo.
echo [5/8] Building paper summary bundle...
python -m src.paper_results_summary
if errorlevel 1 goto :fail

echo.
echo [6/8] Computing diagnostics (SCF convergence, dense Ψ, Weyl, λ-fit)...
python -m src.compute_diagnostics
if errorlevel 1 goto :fail

echo.
echo [7/8] Creating publication-ready plots (now uses diagnostics)...
python -m src.paper_plots --overwrite
if errorlevel 1 goto :fail

echo.
echo [8/8] (Optional) Plot SCF convergence standalone...
python -m src.plot_scf_convergence
if errorlevel 1 echo (Optional SCF plot step failed; continuing.)

echo.
echo ============================================================
echo ✅ All steps completed successfully.
echo  - Results:            data\results.json
echo  - Rescaled:           data\results_rescaled.json
echo  - Paper summary:      data\paper_summary.json  and  data\paper_summary.csv
echo  - Diagnostics JSONs:  data\diagnostics\scf_history.json
echo                         data\diagnostics\scf_convergence.json
echo                         data\diagnostics\bridge_dense.json
echo                         data\diagnostics\weyl_comparison.json
echo                         data\diagnostics\lambda_fit.json
echo  - Figures:            data\figures\psi_raw_vs_s.png
echo                         data\figures\psi_rescaled_vs_s.png
echo                         data\figures\psi_compare_dense_log.png
echo                         data\figures\scf_convergence.png
echo                         data\figures\lambda_sweep.png
echo ============================================================
goto :eof

:fail
echo × Something failed. Check the console output above.
pause
exit /b 1

@echo off
setlocal ENABLEDELAYEDEXPANSION

REM ===============================================
REM  Proof 18 – Asymptotic Spectral Density (Batch)
REM  1) Lädt Eigenwerte
REM  2) Schätzt asymptotische Dichte ρ(E)
REM  3) Vergleicht mit BK-Asymptotik
REM  4) Kopplungs-/Stabilitäts-Check
REM  5) Schreibt Summary + Plots
REM ===============================================

REM ---- Konfiguration -------------------------------------------------
set PYTHON=python
set EIGENS=eigenvalues_clean.txt
set OUTDIR=.
set SCRIPT=proof18_Asymptotic_density_deep.py
REM Optional: Bandbreite/Glättung, Bins, etc. (nur wenn dein Script Flags hat)
REM set FLAGS=--bins 300 --smooth 15

echo ===============================================================
echo [Proof 18] Asymptotic Spectral Density  --  Start
echo ---------------------------------------------------------------
echo Script  : %SCRIPT%
echo Eigen   : %EIGENS%
echo OutDir  : %OUTDIR%
echo Flags   : %FLAGS%
echo ===============================================================
echo.

REM ---- Vorab-Pruefungen --------------------------------------------
if not exist "%SCRIPT%" (
  echo [ERROR] %SCRIPT% nicht gefunden.
  exit /b 1
)

if not exist "%EIGENS%" (
  echo [WARN] %EIGENS% nicht gefunden. Das Script erzeugt evtl. Testdaten.
)

REM ---- Ausführung ---------------------------------------------------
echo [RUN] %PYTHON% %SCRIPT% --eigenfile "%EIGENS%" --outdir "%OUTDIR%" %FLAGS%
%PYTHON% "%SCRIPT%" --eigenfile "%EIGENS%" --outdir "%OUTDIR%" %FLAGS%
if errorlevel 1 (
  echo [ERROR] Proof 18 ist mit einem Fehler beendet worden.
  exit /b 2
)

REM ---- Artifacts ----------------------------------------------------
set SUMMARY=%OUTDIR%\proof18_summary.json
set PLOT1=%OUTDIR%\proof18_asymptotic_density.png
set PLOT2=%OUTDIR%\proof18_stable_final.png

echo.
echo ---------------------------------------------------------------
echo  Erzeugte Dateien (Artifacts)
echo ---------------------------------------------------------------
if exist "%SUMMARY%" echo  - %SUMMARY%
if exist "%PLOT1%"  echo  - %PLOT1%
if exist "%PLOT2%"  echo  - %PLOT2%

REM ---- Kurz-Zusammenfassung aus JSON -------------------------------
if exist "%SUMMARY%" (
  echo.
  echo ---------------------------------------------------------------
  echo  Kurz-Zusammenfassung aus proof18_summary.json
  echo ---------------------------------------------------------------
  %PYTHON% - <<PYCODE
import json, sys, pathlib
p = pathlib.Path(r"%SUMMARY%")
try:
    d = json.loads(p.read_text(encoding="utf-8"))
except Exception as e:
    print("[WARN] Summary konnte nicht gelesen werden:", e)
    sys.exit(0)

def get(k, default="n/a"):
    return d.get(k, default)

print(f"  λ_fit        : {get('lambda_fit','n/a')}")
print(f"  b_fit        : {get('b_fit','n/a')}")
print(f"  res_BK (L2)  : {get('residual_to_BK','n/a')}")
print(f"  res_fit (L2) : {get('residual_to_fit','n/a')}")
print(f"  slope_logE   : {get('slope_logE','n/a')}")
print(f"  R^2          : {get('R2','n/a')}")
print(f"  n_eigs       : {get('n_eigs','n/a')}")
PYCODE
) else (
  echo [WARN] %SUMMARY% nicht gefunden – keine Zusammenfassung moeglich.
)

echo.
echo ===============================================================
echo [OK] Proof 18 fertig. Plots liegen in: %OUTDIR%
echo ===============================================================

endlocal

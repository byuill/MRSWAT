@echo off
REM =========================================================================
REM  MRSWAT_Launcher.bat
REM  Automated preprocessing, model execution, and postprocessing pipeline
REM  for the MRSWAT Mississippi River Salt-Water Intrusion model.
REM =========================================================================
setlocal enabledelayedexpansion

REM Save the base working directory (location of this batch file)
set "BASEDIR=%~dp0"
cd /d "%BASEDIR%"

REM =========================================================================
REM  STEP 1 — Create simulation folder  MRSWAT_YYYYMMDD_##
REM =========================================================================

REM Get current date in YYYYMMDD format (PowerShell for locale independence)
for /f "tokens=*" %%a in ('powershell -NoProfile -Command "Get-Date -Format yyyyMMdd"') do set "DATESTR=%%a"

REM Find the next available two-digit suffix (00, 01, 02, ...)
set /a SUFFIX=0
:FIND_SUFFIX
set "PADDED=0!SUFFIX!"
set "PADDED=!PADDED:~-2!"
set "SIMFOLDER=MRSWAT_!DATESTR!_!PADDED!"
if exist "!SIMFOLDER!\" (
    set /a SUFFIX+=1
    goto FIND_SUFFIX
)

echo =========================================================================
echo  Creating simulation folder: !SIMFOLDER!
echo =========================================================================

REM =========================================================================
REM  STEP 2 — Copy Template into the new simulation folder (exclude .venv)
REM =========================================================================

REM robocopy exit codes 0-7 are success; >=8 is failure
robocopy "Template" "!SIMFOLDER!" /E /XD ".venv" /NFL /NDL /NJH /NJS /NC /NS /NP >nul
if !errorlevel! geq 8 (
    echo  ERROR: Failed to copy Template folder.
    goto :EOF
)

REM Copy Python scripts from "Code Storage Main" into the Code Storage subfolder
robocopy "Code Storage Main" "!SIMFOLDER!\Code Storage" /E /XD "__pycache__" /NFL /NDL /NJH /NJS /NC /NS /NP >nul
if !errorlevel! geq 8 (
    echo  ERROR: Failed to copy Code Storage Main folder.
    goto :EOF
)

echo  Template and Code Storage copied successfully.
echo.

REM =========================================================================
REM  Set portable Python environment (self-contained in python_env folder)
REM =========================================================================
set "PYTHON=%BASEDIR%python_env\python.exe"
set "PYTHONIOENCODING=utf-8"
if not exist "%PYTHON%" (
    echo  ERROR: Portable Python not found at "%PYTHON%".
    echo         Ensure the python_env folder is present alongside this batch file.
    goto :EOF
)
echo  Portable Python: %PYTHON%
echo.

REM =========================================================================
REM  STEP 3 — Run Create_BC.py  (boundary-condition preprocessor)
REM =========================================================================

echo =========================================================================
echo  Running Create_BC.py ...
echo =========================================================================

pushd "!SIMFOLDER!\Code Storage"
"%PYTHON%" Create_BC.py > Create_BC.log 2>&1
set "BC_ERR=!errorlevel!"
popd

if !BC_ERR! neq 0 (
    echo  ERROR: Create_BC.py returned exit code !BC_ERR!.
    echo         Check "!SIMFOLDER!\Code Storage\Create_BC.log" for details.
    goto :EOF
)
echo  Create_BC.py completed.  Log: "!SIMFOLDER!\Code Storage\Create_BC.log"
echo.

REM =========================================================================
REM  STEP 4 — Determine Hot Start vs Cold Start
REM =========================================================================
REM  A hotstart is used when a previous simulation folder exists whose date
REM  is within 3 days of today.  The most recent qualifying folder (by date
REM  then by sequence number) is selected.
REM =========================================================================

echo =========================================================================
echo  Checking for previous simulation (hotstart eligibility) ...
echo =========================================================================

set "HOTSTART=0"
set "PREV_FOLDER="

REM --- PowerShell one-liner to find the most recent qualifying folder ------
for /f "usebackq tokens=*" %%f in (`powershell -NoProfile -Command ^
    "$today = Get-Date;" ^
    "Get-ChildItem -Directory -Path '.' -Filter 'MRSWAT_*_*' |" ^
    "  Where-Object { $_.Name -ne '!SIMFOLDER!' } |" ^
    "  ForEach-Object {" ^
    "    if ($_.Name -match '^MRSWAT_(\d{8})_(\d{2})$') {" ^
    "      try {" ^
    "        $d = [datetime]::ParseExact($Matches[1],'yyyyMMdd',$null);" ^
    "        $diff = ($today - $d).TotalDays;" ^
    "        if ($diff -ge 0 -and $diff -le 3) {" ^
    "          [PSCustomObject]@{Name=$_.Name; Date=$d; Seq=[int]$Matches[2]}" ^
    "        }" ^
    "      } catch {}" ^
    "    }" ^
    "  } |" ^
    "  Sort-Object Date,Seq -Descending |" ^
    "  Select-Object -First 1 -ExpandProperty Name;"`) do (
    set "PREV_FOLDER=%%f"
    set "HOTSTART=1"
)

REM =========================================================================
REM  STEP 5 — Run the appropriate RunFile creator
REM =========================================================================

if "!HOTSTART!"=="1" (
    echo  Previous simulation found: !PREV_FOLDER!
    echo  ^>^> HOTSTART mode activated.
    echo.
    echo =========================================================================
    echo  Running Create_RunFile_Hotstart.py ...
    echo =========================================================================

    pushd "!SIMFOLDER!\Code Storage"
    "%PYTHON%" Create_RunFile_Hotstart.py > Create_RunFile.log 2>&1
    set "RF_ERR=!errorlevel!"
    popd

    if !RF_ERR! neq 0 (
        echo  ERROR: Create_RunFile_Hotstart.py returned exit code !RF_ERR!.
        echo         Check "!SIMFOLDER!\Code Storage\Create_RunFile.log" for details.
        goto :EOF
    )
    echo  Create_RunFile_Hotstart.py completed.  Log: "!SIMFOLDER!\Code Storage\Create_RunFile.log"
) else (
    echo  No qualifying previous simulation found within 3 days.
    echo  ^>^> COLDSTART mode activated.
    echo.
    echo =========================================================================
    echo  Running Create_RunFile_Coldstart.py ...
    echo =========================================================================

    pushd "!SIMFOLDER!\Code Storage"
    "%PYTHON%" Create_RunFile_Coldstart.py > Create_RunFile.log 2>&1
    set "RF_ERR=!errorlevel!"
    popd

    if !RF_ERR! neq 0 (
        echo  ERROR: Create_RunFile_Coldstart.py returned exit code !RF_ERR!.
        echo         Check "!SIMFOLDER!\Code Storage\Create_RunFile.log" for details.
        goto :EOF
    )
    echo  Create_RunFile_Coldstart.py completed.  Log: "!SIMFOLDER!\Code Storage\Create_RunFile.log"
)
echo.

REM =========================================================================
REM  STEP 6 — Run the MRSWAT Fortran model
REM =========================================================================
REM  The executable reads the input-file name from stdin.  The file
REM  "input.txt" contains the single line "main.inp".
REM =========================================================================

echo =========================================================================
echo  Running MRSWAT model  (MRSWAT-111825.exe) ...
echo =========================================================================

pushd "!SIMFOLDER!"
"MRSWAT-111825.exe" < input.txt > MRSWAT_run.log 2>&1
set "MODEL_ERR=!errorlevel!"
popd

if !MODEL_ERR! neq 0 (
    echo  ERROR: MRSWAT-111825.exe returned exit code !MODEL_ERR!.
    echo         Check "!SIMFOLDER!\MRSWAT_run.log" for details.
    goto :EOF
)
echo  MRSWAT model completed successfully.  Log: "!SIMFOLDER!\MRSWAT_run.log"
echo.

REM =========================================================================
REM  STEP 7 — Run Postprocess_Results.py
REM =========================================================================

echo =========================================================================
echo  Running Postprocess_Results.py ...
echo =========================================================================

pushd "!SIMFOLDER!\Code Storage"
"%PYTHON%" Postprocess_Results.py > Postprocess_Results.log 2>&1
set "PP_ERR=!errorlevel!"
popd

if !PP_ERR! neq 0 (
    echo  ERROR: Postprocess_Results.py returned exit code !PP_ERR!.
    echo         Check "!SIMFOLDER!\Code Storage\Postprocess_Results.log" for details.
    goto :EOF
)
echo  Postprocess_Results.py completed.
echo.

REM =========================================================================
REM  STEP 8 — Write Available Observations
REM =========================================================================
REM  Run Observations_check.py to check for new observational salinity
REM  measurements and write them to the MRSWAT_Data.xlsx spreadsheet.
REM =========================================================================

echo =========================================================================
echo  Running Observations_check.py (Write Available Observations) ...
echo =========================================================================

pushd "!SIMFOLDER!\Code Storage"
"%PYTHON%" Observations_check.py > Observations_check.log 2>&1
set "OBS_ERR=!errorlevel!"
popd

if !OBS_ERR! neq 0 (
    echo  WARNING: Observations_check.py returned exit code !OBS_ERR!.
    echo           Check "!SIMFOLDER!\Code Storage\Observations_check.log" for details.
    echo           Continuing pipeline...
)
if !OBS_ERR! equ 0 (
    echo  Observations_check.py completed.  Log: "!SIMFOLDER!\Code Storage\Observations_check.log"
)
echo.

REM =========================================================================
REM  STEP 9 — Hotstart QC check (only when this was a hotstart simulation)
REM =========================================================================

if "!HOTSTART!" neq "1" goto :PIPELINE_DONE

echo =========================================================================
echo  STEP 9: Running HS_QC.py (hotstart quality control) ...
echo =========================================================================

pushd "!SIMFOLDER!\Code Storage"
"%PYTHON%" HS_QC.py > HS_QC.log 2>&1
set "QC_ERR=!errorlevel!"
popd

REM Echo the saved log to the terminal so the operator sees HS_QC output
REM in real time as well as in the persisted log file.
echo  --- HS_QC.py output (log: "!SIMFOLDER!\Code Storage\HS_QC.log") ---
type "!SIMFOLDER!\Code Storage\HS_QC.log"
echo  --- end of HS_QC.py output ---
echo.

if !QC_ERR! equ 0 (
    echo =========================================================================
    echo  HS_QC RESULT: MODELED HOTSTART ACCEPTED
    echo    The observed 9-ppt toe location agrees with the modeled hotstart
    echo    within the 15 river-mile threshold.  No corrective run is needed.
    echo    Continuing with the current modeled-hotstart simulation results.
    echo =========================================================================
    goto :PIPELINE_DONE
)

if !QC_ERR! neq 2 (
    echo  WARNING: HS_QC.py returned unexpected exit code !QC_ERR!.
    echo           Check "!SIMFOLDER!\Code Storage\HS_QC.log" for details.
    goto :PIPELINE_DONE
)

REM =========================================================================
REM  STEP 10 — QC mismatch: re-run with synthetic hotstart using observed toe
REM =========================================================================

echo =========================================================================
echo  HS_QC RESULT: SYNTHETIC HOTSTART TRIGGERED
echo    The observed 9-ppt toe location differs from the modeled hotstart
echo    toe by more than 15 river miles.  A corrected simulation will be
echo    built using a SYNTHETIC hotstart based on the observed toe location.
echo    Script: Create_RunFile_SyntheticHotstart.py
echo =========================================================================
echo.

set "MOD_SIMFOLDER=!SIMFOLDER!_Mod"

echo =========================================================================
echo  Creating modified simulation folder: !MOD_SIMFOLDER!
echo =========================================================================

robocopy "Template" "!MOD_SIMFOLDER!" /E /XD ".venv" /NFL /NDL /NJH /NJS /NC /NS /NP >nul
if !errorlevel! geq 8 (
    echo  ERROR: Failed to copy Template to !MOD_SIMFOLDER!.
    goto :PIPELINE_DONE
)

robocopy "Code Storage Main" "!MOD_SIMFOLDER!\Code Storage" /E /XD "__pycache__" /NFL /NDL /NJH /NJS /NC /NS /NP >nul
if !errorlevel! geq 8 (
    echo  ERROR: Failed to copy Code Storage Main to !MOD_SIMFOLDER!\Code Storage.
    goto :PIPELINE_DONE
)

REM Copy hs_qc_result.txt so Create_RunFile_SyntheticHotstart.py can read the observed toe RM
copy /Y "!SIMFOLDER!\hs_qc_result.txt" "!MOD_SIMFOLDER!\hs_qc_result.txt" >nul

echo  Modified folder created: !MOD_SIMFOLDER!
echo.

echo =========================================================================
echo  Running Create_BC.py for modified simulation ...
echo =========================================================================

pushd "!MOD_SIMFOLDER!\Code Storage"
"%PYTHON%" Create_BC.py > Create_BC.log 2>&1
set "BC_ERR=!errorlevel!"
popd

if !BC_ERR! neq 0 (
    echo  ERROR: Create_BC.py [Mod] returned exit code !BC_ERR!.
    echo         Check "!MOD_SIMFOLDER!\Code Storage\Create_BC.log" for details.
    goto :PIPELINE_DONE
)
echo  Create_BC.py (Mod) completed.  Log: "!MOD_SIMFOLDER!\Code Storage\Create_BC.log"
echo.

echo =========================================================================
echo  Running Create_RunFile_SyntheticHotstart.py for modified simulation ...
echo =========================================================================

pushd "!MOD_SIMFOLDER!\Code Storage"
"%PYTHON%" Create_RunFile_SyntheticHotstart.py > Create_RunFile.log 2>&1
set "RF_ERR=!errorlevel!"
popd

if !RF_ERR! neq 0 (
    echo  ERROR: Create_RunFile_SyntheticHotstart.py returned exit code !RF_ERR!.
    echo         Check "!MOD_SIMFOLDER!\Code Storage\Create_RunFile.log" for details.
    goto :PIPELINE_DONE
)
echo  Create_RunFile_SyntheticHotstart.py (Mod) completed.  Log: "!MOD_SIMFOLDER!\Code Storage\Create_RunFile.log"
echo.

echo =========================================================================
echo  Running MRSWAT model for modified simulation ...
echo =========================================================================

pushd "!MOD_SIMFOLDER!"
"MRSWAT-111825.exe" < input.txt > MRSWAT_run.log 2>&1
set "MODEL_ERR=!errorlevel!"
popd

if !MODEL_ERR! neq 0 (
    echo  ERROR: MRSWAT-111825.exe [Mod] returned exit code !MODEL_ERR!.
    echo         Check "!MOD_SIMFOLDER!\MRSWAT_run.log" for details.
    goto :PIPELINE_DONE
)
echo  MRSWAT model (Mod) completed successfully.  Log: "!MOD_SIMFOLDER!\MRSWAT_run.log"
echo.

echo =========================================================================
echo  Running Postprocess_Results.py for modified simulation ...
echo =========================================================================

pushd "!MOD_SIMFOLDER!\Code Storage"
"%PYTHON%" Postprocess_Results.py > Postprocess_Results.log 2>&1
set "PP_ERR=!errorlevel!"
popd

if !PP_ERR! neq 0 (
    echo  ERROR: Postprocess_Results.py [Mod] returned exit code !PP_ERR!.
    echo         Check "!MOD_SIMFOLDER!\Code Storage\Postprocess_Results.log" for details.
    goto :PIPELINE_DONE
)
echo  Postprocess_Results.py (Mod) completed.
echo.

echo =========================================================================
echo  Modified simulation completed: !MOD_SIMFOLDER!
echo =========================================================================
echo.

:PIPELINE_DONE
echo =========================================================================
echo  MRSWAT pipeline finished: !SIMFOLDER!
echo =========================================================================

endlocal

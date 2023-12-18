@echo off
SET venv_name=llambit_env

REM Check if the virtual environment already exists using Conda
conda env list | findstr /C:"%venv_name%"
IF %ERRORLEVEL% == 0 (
    echo Virtual environment "%venv_name%" already exists.
) ELSE (
    echo Creating virtual environment "%venv_name%"...
    conda create --name %venv_name% python=3.9.7 --yes
    IF %ERRORLEVEL% NEQ 0 GOTO Error
)

REM Activate the virtual environment
call conda activate %venv_name%

REM Check if pip is already installed
conda list pip | findstr pip
IF %ERRORLEVEL% NEQ 0 (
    conda install -c conda-forge pip
    IF %ERRORLEVEL% NEQ 0 GOTO Error
)

REM Check if pipreqs is already installed
conda list pipreqs | findstr pipreqs
IF %ERRORLEVEL% NEQ 0 (
    conda install -c conda-forge pipreqs
    IF %ERRORLEVEL% NEQ 0 GOTO Error
)

pipreqs
IF %ERRORLEVEL% NEQ 0 GOTO Error

pip install -r requirements.txt
IF %ERRORLEVEL% NEQ 0 GOTO Error

pip install git+https://github.com/Clarivate-SAR/woslite_py_client.git
IF %ERRORLEVEL% NEQ 0 GOTO Error

(echo Success) > install_success.txt
GOTO End

:Error
(echo Failure) > install_success.txt

:End

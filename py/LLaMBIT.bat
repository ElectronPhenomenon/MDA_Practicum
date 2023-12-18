@echo off
SET "AnacondaPath=%USERPROFILE%\anaconda3"
CALL "%AnacondaPath%\Scripts\activate.bat" %AnacondaPath%
CALL conda activate llambit_env
python main.py

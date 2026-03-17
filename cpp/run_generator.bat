@echo off
setlocal

set CONFIG=%1
if "%CONFIG%"=="" set CONFIG=cpp/configs/directional_smoke.json

set COMMAND=%2
if "%COMMAND%"=="" set COMMAND=generate

set THREADS=%3
if "%THREADS%"=="" set THREADS=16

powershell -ExecutionPolicy Bypass -File "%~dp0run_generator.ps1" -Config "%CONFIG%" -Command "%COMMAND%" -Threads %THREADS%
exit /b %ERRORLEVEL%

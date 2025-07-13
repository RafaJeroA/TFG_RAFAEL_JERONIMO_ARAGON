@echo off
echo Lanzando Julia con hilos automáticos...
cd /d "%~dp0"
julia -t auto main.jl
pause

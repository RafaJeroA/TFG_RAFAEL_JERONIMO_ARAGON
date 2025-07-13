@echo off
REM ─────────────────────────────────────────────────────────────────────
REM run_model.bat — Abre PowerShell, limita recursos y ejecuta el model
REM ─────────────────────────────────────────────────────────────────────

powershell -NoProfile -ExecutionPolicy Bypass -Command ^
  "& {
    # Detectar cores
    $cores = [Environment]::ProcessorCount
    # Calcular 90% de cores para Julia
    $jthreads = [math]::Max(1, [math]::Floor($cores * 0.9))
    # Workers = hilos - 1 (reserva 1 para el master)
    $workers = [math]::Max(1, $jthreads - 1)
    # Fijar variable de entorno para threads de Julia
    $env:JULIA_NUM_THREADS = $jthreads
    Write-Host 'Threads=' $jthreads ', Workers=' $workers
    # Cambiar al directorio del script
    Set-Location -LiteralPath '%~dp0'
    # Ejecutar Julia con live log y guardarlo en model.log
    & julia --threads $jthreads -p $workers --heap-size-hint=80%% main.jl --parallel 2>&1 | Tee-Object -FilePath model.log
  }"

pause
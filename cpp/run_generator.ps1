param(
    [string]$Config = "cpp/configs/directional_local_max600_fourphase.json",
    [string]$BuildDir = "cpp/build_fftw_omp16",
    [string]$Command = "generate",
    [string]$Compiler = "g++",
    [int]$Threads = 16
)

$ErrorActionPreference = "Stop"

$repoRoot = Split-Path -Parent $PSScriptRoot
Set-Location $repoRoot

Write-Host "Configuring CMake in $BuildDir ..."
& cmake -S cpp -B $BuildDir -G Ninja "-DCMAKE_CXX_COMPILER=$Compiler" "-DMF12_ENABLE_OPENMP=ON" "-DMF12_ENABLE_FFTW=ON"
if ($LASTEXITCODE -ne 0) {
    throw "CMake configure failed."
}

Write-Host "Building generator ..."
& cmake --build $BuildDir
if ($LASTEXITCODE -ne 0) {
    throw "CMake build failed."
}

$exe = Join-Path $BuildDir "ow3d_directional_generator.exe"
if (-not (Test-Path $exe)) {
    throw "Generator executable not found: $exe"
}

if ($Threads -gt 0) {
    $env:OMP_NUM_THREADS = "$Threads"
    $env:OMP_DYNAMIC = "FALSE"
    $env:OMP_PROC_BIND = "TRUE"
    $env:OMP_PLACES = "cores"
    $env:OMP_WAIT_POLICY = "ACTIVE"
    Write-Host "OMP_NUM_THREADS set to $Threads"
    Write-Host "OMP_DYNAMIC set to FALSE"
    Write-Host "OMP_PROC_BIND set to TRUE"
    Write-Host "OMP_PLACES set to cores"
    Write-Host "OMP_WAIT_POLICY set to ACTIVE"
}

Write-Host "Running: $Command $Config"
& $exe $Command $Config
if ($LASTEXITCODE -ne 0) {
    throw "Generator command failed."
}

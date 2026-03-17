#!/usr/bin/env bash
set -euo pipefail

CONFIG="${1:-cpp/configs/directional_smoke.json}"
COMMAND="${2:-generate}"
THREADS="${3:-16}"
BUILD_DIR="${BUILD_DIR:-cpp/build_linux}"
CXX_COMPILER="${CXX_COMPILER:-g++}"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

need_cmd() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "Missing required command: $1" >&2
    exit 1
  fi
}

need_cmd cmake
need_cmd "$CXX_COMPILER"

if ! command -v ninja >/dev/null 2>&1; then
  echo "ninja not found; CMake will try the default generator." >&2
  GENERATOR_ARGS=()
else
  GENERATOR_ARGS=(-G Ninja)
fi

echo "Environment check"
echo "- compiler: $("$CXX_COMPILER" --version | head -n 1)"
echo "- cmake: $(cmake --version | head -n 1)"
if command -v pkg-config >/dev/null 2>&1 && pkg-config --exists fftw3; then
  echo "- FFTW: found via pkg-config ($(pkg-config --modversion fftw3))"
elif [[ -f cpp/third_party/fftw/include/fftw3.h ]]; then
  echo "- FFTW: found in cpp/third_party/fftw"
else
  echo "- FFTW: not detected in system paths or cpp/third_party/fftw"
fi

echo "Configuring CMake in $BUILD_DIR ..."
cmake -S cpp -B "$BUILD_DIR" "${GENERATOR_ARGS[@]}" \
  -DCMAKE_CXX_COMPILER="$CXX_COMPILER" \
  -DMF12_ENABLE_OPENMP=ON \
  -DMF12_ENABLE_FFTW=ON

echo "Building generator ..."
cmake --build "$BUILD_DIR"

EXE="$BUILD_DIR/ow3d_directional_generator"
if [[ ! -x "$EXE" ]]; then
  echo "Generator executable not found: $EXE" >&2
  exit 1
fi

if [[ "$THREADS" != "0" ]]; then
  export OMP_NUM_THREADS="$THREADS"
  export OMP_DYNAMIC=FALSE
  export OMP_WAIT_POLICY=ACTIVE
  echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
fi

echo "Running: $COMMAND $CONFIG"
"$EXE" "$COMMAND" "$CONFIG"

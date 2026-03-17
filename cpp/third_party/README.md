## Local Third-Party Dependencies

This folder is the repo-local place for optional native dependencies when you do not want to rely only on system installs.

### FFTW

If you want to keep FFTW alongside the project, place it under:

```text
cpp/third_party/fftw/
  include/fftw3.h
  lib/libfftw3.so        # Linux
  lib/libfftw3.dll.a     # MinGW on Windows
  bin/libfftw3-3.dll     # optional Windows runtime DLL
```

The CMake build will look here automatically when `-DMF12_ENABLE_FFTW=ON`.

### OpenMP

OpenMP is not a standalone library that we can vendor in the same way. It comes from the compiler toolchain and runtime, for example:

- GCC/G++: `libgomp`
- Clang: `libomp`

So for Ubuntu, the practical path is to install a compiler with OpenMP support rather than copying an `OpenMP` folder into the repo.

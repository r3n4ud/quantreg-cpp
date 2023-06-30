# Quantile Regression C++ Library

## Prerequisites

This project has been designed to be built using CMake and depends on Armadillo and HDF5.

## Build and Install

```bash
export Armadillo_ROOT=~/.local && export HDF5_ROOT=~/.local
git clone https://github.com/r3n4ud/quantreg-cpp.git
cd quantreg-cpp
cmake -S . -B build-shared -DBUILD_SHARED_LIBS=ON  -DCMAKE_BUILD_TYPE=Release
cmake -S . -B build-static -DBUILD_SHARED_LIBS=OFF -DCMAKE_BUILD_TYPE=Release
cmake --build build-shared
cmake --build build-static
cmake --install build-shared --prefix _install
cmake --install build-static --prefix _install
```

:warning: For source release (i.e. without any git history), you will need a `.release.cmake` file
at the source root. For example, for a `v1.0.4` release:

```cmake
set(QUANTREG_CPP_VERSION_MAJOR 1)
set(QUANTREG_CPP_VERSION_MINOR 0)
set(QUANTREG_CPP_VERSION_PATCH 4)
set(QUANTREG_CPP_VERSION "${QUANTREG_CPP_VERSION_MAJOR}.${QUANTREG_CPP_VERSION_MINOR}.${QUANTREG_CPP_VERSION_PATCH}")
```

## Test!

```bash
cd tests
cmake -S . -B build
cmake --build build --target test
```

If anything goes wrong, you may:

```bash
cd build
ctest --rerun-failed --output-on-failure # or --verbose
```

## Usage

With the following structure:
```
.
├── CMakeLists.txt
├── engel.csv
└── foo.cpp
```

`CMakeLists.txt`:
```cmake
cmake_minimum_required(VERSION 3.24)
project(foo CXX)

find_package(BLAS REQUIRED)
find_package(HDF5 1.14.0 REQUIRED)
find_package(Armadillo 10.5.3 REQUIRED NO_MODULE PATHS ${Armadillo_ROOT} NO_DEFAULT_PATH)
find_package(quantreg-cpp REQUIRED)

add_executable(foo foo.cpp)
target_link_libraries(foo PUBLIC quantreg-cpp::quantreg-cpp ${HDF5_LIBRARIES} ${ARMADILLO_LIBRARIES})
```

`foo.cpp`:
```cpp
#include <quantreg-cpp>

int main() {
  try {
    arma::mat data;
    data.load("../engel.csv", arma::csv_ascii);
    data.shed_cols(0, 0);
    data.shed_rows(0, 0);
    arma::vec y = data.col(1);
    arma::mat X = data(arma::span::all, arma::span(0, 0));

    arma::vec coefficients;
    arma::vec residues;
    std::tie(coefficients, residues) = quantreg::qr(X, y, 0.5);
    std::cout.precision(9);
    coefficients.t().raw_print("coef=");

  } catch (const std::runtime_error& e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "Error: unknown exception" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
```

`engel.csv` can be found in the `data` directory of this repository.

```sh
cmake \
    -Dquantreg-cpp_ROOT=/opt/quantreg-cpp \
	-DArmadillo_ROOT=/opt/armadillo \
	-DHDF5_ROOT=/opt/hdf5 \
	-DBUILD_SHARED_LIBS=ON \
	-S . -B build
cmake --build build
cd build
./foo
```
will return:

```
coef=
81.4822474 0.560180551
```

## Issues and Feature Requests

Report issues and request features by opening an issue
[here](https://github.com/r3n4ud/quantreg-cpp/issues).

## Credits

### quantreg pfn
This is a C++ implementation of the Frisch-Newton algorithm and the preprocessing phase of the rq
algorithm as described in Portnoy and Koenker, Statistical Science, (1997) 279-300.
Original version consists in the [Roger Koenker's quantreg R
package](https://github.com/cran/quantreg) rq.fit.pfn and rq.fit.fnb functions (quantreg.R) and the
rqfnb fortran function (rqfnb.f).

Translated and adapted to C++ / Rcpp / RcppArmadillo by Dajun Xu, under the instructions of
Stephen L. Portnoy, 2015. Repository: https://github.com/xdj701/STAT391-Quantreg.
Adapted to C++ & Armadillo by Renaud Aubin, 2023.

This project is released under the 3-Clause BSD License, although the original versions were
not. Permissions was given by Professor Roger Koenker and Dajun Xu to Renaud Aubin on May/June
2023 via electronic correspondences.

Original correspondences can be provided upon request.

### CMake packaging
The cmake structure is a derivative work from Alex Reinking, licensed under the MIT License.

```
MIT License

Copyright (c) 2021 Alex Reinking

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

## License
This project is licensed under the 3-Clause BSD License.

```
3-Clause BSD License

Copyright (c) 2023 Renaud AUBIN

Redistribution and use in source and binary forms, with or without modification, are permitted
provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions
   and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions
   and the following disclaimer in the documentation and/or other materials provided with the
   distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse
   or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

```

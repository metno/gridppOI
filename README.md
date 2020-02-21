# Optimal interpolation library

OI algorithm from gridpp in a standalone C++ library

## Required dependencies

You need to install swig:

```
sudo apt install swig doxygen
```

## Installation of library

Install the library and packages using cmake. Create a build directory and perform the
installation there:

```
mkdir build
cd build
cmake ../
```

To compile just the library, do this:
```
make gridppOI
```

This creates the file libgridppOI.so

## Installation and use of python module

To **install** the python module, you only need the following (i.e. you do not need to also
install the gridppOI library)

```
sudo make install-python
```

The installation above should install the python package in a central place, which should
automatically be accessible in python3.

If you instead just want to **build** the package, do this:

```
make build-python
```

The package is then available in SWIG/python/gridppOI.py. Here is an example of how to use the package:

```python
import gridppOI
import numpy as np
N = 300
Y = 150
X = 100

# Define input grid
input = np.random.randn(Y, X) * 3;
blats = np.random.randn(Y, X);
blons = np.random.randn(Y, X);
belevs = np.random.rand(Y, X) * 100;
blafs = np.random.rand(Y, X) * 1;

# Define observations
pobs = np.random.randn(N) * 3;
pci = np.random.rand(N) + 0.5;
plats = np.random.randn(N);
plons = np.random.randn(N);
pelevs = np.random.rand(N) * 100;
plafs = np.random.rand(N);

# OI parameters
minRho = 0.0013
hlength = 10000
vlength = 10000
wmin = 0
maxElevDiff = 100
landOnly=False
maxLocations = 20
elevGradient = -0.0065
epsilon = 0.5

status, output = gridppOI.optimal_interpolation(input, blats, blons, belevs, blafs, pobs, pci,
        plats, plons, pelevs, plafs, minRho, hlength, vlength, wmin, maxElevDiff, landOnly,
        maxLocations, elevGradient, epsilon)
print(output)
```

## Run the test suite

Tests are written in the tests/ directory. The tests rely on the python interface to titanlib being installed (see above). Use nosetests3 to run all tests. First install the python library  and the nose package (pip3 install nose) then run the following:

```bash
nosetests3
```

## Build the documentation

Inside the build directory, type

```bash
make doc_doxygen
```

The html documentation will appear in build/doxygen/html.

## Build a debug version of the library

Just use the following options when running cmake:

```bash
cmake .. -DCMAKE_BUILD_TYPE=DEBUG
```

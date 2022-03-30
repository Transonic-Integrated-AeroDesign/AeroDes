# AeroDesign

* [Build](#Build)
* [Executables](#Executables)

## <a name="Build"></a> Build Instructions

Shown below are the two steps required to compile this code.

```
cd AeroDes

mkdir build       (make a build directory)
cd build

cmake ../.
make              (compiles the code)
make install      (installs the executables to your system)
```

## <a name="Executables"></a> Executables

There are two methods of using this code in your design process. You may either use the direct executables in linear fashion where you run each program separately. Or you may import the shared library resource. Lets go through both so you know how to use either of these resources.

### 1. Direct Executables

Here are the executables for various design tasks. Each will of course require their own input parameters found in their respective examples directory.

```
canareq
cfmactu
prandtl
wake
```

### 2. Shared Library

#### Quick test program for shared library

Compile the sample C++ program like so.

```g++ -o test -laerolib main.cpp```

Here's a look at the sample C++ program named, main.cpp.

```
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream> // std
#include <iomanip>  // setw
#include <sstream>  // istream
#include <fstream>  // fopen, ifstream
#include <string>
#include <stdio.h>  // strcpy
#include "aerodes.hpp"

int main(int argc, char** argv) {
    aerodes *aero = new aerodes(argc, argv);
    
    // ... stuff ...
    
    delete aero;
    return 1;
}

```

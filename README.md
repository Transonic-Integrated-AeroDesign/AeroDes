# AeroDesign
The original version of this software was written by Dr. Jean-Jacques Chattot in Fortran.  This C++ implementation has been created with the permission of Dr. Chattot.

The work done herein is based solely on the vortex line method outlined by Professors J.J. Chattot and M.M. Hafez [1].
The point of contact for questions should be directed to Carlos Pereyra.

Follow the sections below for compiling and using this package.

* [Build](#Build)
* [How to Use](#Usage)
* [References](#References)
## <a name="Build"></a> Build Instructions

Below are the steps required to compile this project.

**step one**, clone the repository
```
git clone https://github.com/Transonic-Integrated-AeroDesign/AeroDes
```

**step two**, change directory into AeroDes
```
cd AeroDes
```

**step three**, make a new build directory and change into it
```
mkdir build 
cd build
```

**step four**, build code with cmake utility and install to your system
```
cmake ../.
make
make install
```

```make``` compiles the code. ```make install``` installs the executables to your system.

### Re-syncing your code with latest stuff
```
git fetch
git pull origin master
```

### Depedencies

In order to compile the executables and the shared libraries, you will need the following packages on your computer.

**Mac OS X**
* cmake https://cmake.org/install/
* g++ https://developer.apple.com/xcode/
  
**Windows**
* WSL 2
* cmake
* g++

## <a name="Usage"></a> How to Use

There are two methods of using this code in your design process. You may either use the direct executables in linear fashion where you run each program separately. Or you may import the shared library resource. Lets go through both so you know how to use either of these resources.

### Executables

Here are the executables for various design tasks. Each will of course require their own input parameters found in their 
respective examples directory. 

```
canareq
cfmactu
prandtl
wake
```

See the examples/ directory for more details in running each example with these
executables.

### Shared Library

Outlined here are instructions for writing your own C++ program 
utilizing the AeroDes library. Follow these steps for an understanding of importing and using certain functions within the AeroDes framework.

Include the AeroDes library like so.
```c++
#include "aerodes.hpp"  // includes definitions and objects
```

Create a main function. This function effectively does whatever is placed inside the parenthesis.
```c++
int main(int argc, char** argv) {
    ...
    return 1;
}
```

Inside the main function create the AeroDes object for calling canard wake program ('wk').
```asm
int main(int argc, char** argv) {
    aerodes *aero = new aerodes(argc, argv); // create aero design object
    
    // aero->some_sub_object->function()
    
    // canard ADwake
    aero->wk->readInputPolar("polarbl.dat"); // reads in 2d polar from xfoil
    
    // canard equilibium
    aero->canary->readInputPolar("polarbl.dat"); // reads in 2d polar from xfoil
    
    delete aero;    // free memory at the end
    return 1;
}
```
Ultimately all the public functions in the src/wake, src/canareq, src/prandtline directories are accessible thorugh the ```*aero``` pointer.
Here's a look at the full template C++ program.

main.cpp
```c++
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

    // canard ADwake
    aero->wk->readInputPolar("polarbl.dat");

    // canard equilibium
    aero->canary->readInputPolar("polarbl.dat");
    
    delete aero;
    return 1;
}

```

Finally once the file above is complete; compile the sample C++ program like so.

```g++ -o test -laerolib main.cpp```

You absolutely must have the ```-laerolib``` flag in the compile line in order for your computer to find the ```aerodes.hpp``` library.

## Contributors
* Jean-Jacques Chattot, University of California, Davis, Professor
* Carlos Pereyra, University of California, Davis, Grad. Student (czpereyra@ucdavis.edu)

## <a name="References"></a> References

[1] J.J. Chattot and M.M. Hafez. 2015. Theoretical and Applied Aerodynamics: And Related Numerical Methods.
Springer Netherlands, Dordrecht

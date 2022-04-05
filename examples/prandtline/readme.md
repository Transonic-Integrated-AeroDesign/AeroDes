# Prandtline Instructions

This directory contains all the **required files to have in your working directory** before running. These
files are listed below:

* ***prandtline.data***
* ***polarbl[n].dat***

Typically the prandtline input data file is named ***prandtline.data***; however you 
may change it using the ```-in_prandtline``` flag (shown below). 

If you wish to integrate wing + fuselage
ensemble then read in the 9 polars named, ***polarbl[n].dat***. These will then be read into the prandtline program.


## Example:
You may run this example using two methods. Using the executables or with the AeroDes.hpp library.

### Using PrandtLine Exectuable
automatically reads in parameter file *prandtline.data*:
```
./../../prandtl
```

specify the parameter file *prandtline.data* or whatever file you want:
```
./../../prandtl -in_prandtline prandtline.data
```

### Using AeroDes-PrandtLine Functions

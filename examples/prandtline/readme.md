# Prandtline Instructions

This directory contains all the **required files to have in your working directory** before running. These
files are listed below:

* **prandtline.data**
* **polarbl[n].dat**

Typically the prandtline input data file is named **prandtline.data**; however you 
may change it using the ```--prandtl_in``` flag (shown below). 

If you wish to integrate wing + fuselage
ensemble then read in the 9 polars named, **polarbl[n].dat**. These will then be read into the prandtline program.


## Example:
You may run this example using two methods. Using the executables or with the AeroDes.hpp library.

### Using PrandtLine Exectuable
By default the executable automatically reads in a parameter file named, *prandtline.data*, and a polar file named, *polarbl.dat*:
```
prandtl
```

specify the input **parameter** file to whatever file you want:
```
prandtl --prandtl_in prandtline.data
```

specify the input **polar** file to whatever file you want:
```
prandtl --prantl_polar polarbl.dat
```

you can even use these flags simultaneously:
```
prandtl --prandtl_in prandtline.data --prantl_polar polarbl.dat
```
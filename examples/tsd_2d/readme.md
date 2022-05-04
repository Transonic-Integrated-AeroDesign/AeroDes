# TSD Instructions

![qjimage](https://github.com/Transonic-Integrated-AeroDesign/AeroDes/blob/master/docs/tsd/qj.png)

This directory contains all the **required files to have in your working directory** before running. These
files are listed below:

* **tsd.data**
* **geoprofortsd.xde**
* **tsd.in**

Typically the input data file is named **tsd.data**; however you 
may change it using the ```--tsd_in``` flag (shown below).

Users may interactively set the number of solver iterations using the input
flag ```--tsd_it```. Alternatively, you may also set it within the input file.

## Example:
You may run this example using two methods. Using the executables or with the AeroDes.hpp library.

### Using the TSD Exectuable
By default the executable automatically reads in a parameter file named, *tsd.data*:
```
ADtsd
```

To specify the input file to whatever file you want:
```
ADtsd --tsd_in tsd.data
```

To specify the number of solver iterations:

```
ADtsd --tsd_it 1000
```

To use these flags simultaneously do the following:
```
ADtsd --tsd_in prandtline.data --tsd_it 1000
```

## Input File Commands:

First off, lets start with these options. Previously (in the fortran codes) these 
options were real time user input. However in the cpp revamp we have strayed 
away from this to streamline the design process.

```
IMESH   1           # do you want to write the geometry mesh? Y/N=1/0
IWRITE  1           # do you want to write the geometry? Y/N=1/0
IFLOW   0           # do you want to read-in the flow? Y/N=1/0
ITER    1           # solver iterations, itx=?
```

The rest of the options are relatively straight forward. 
However, these are of importance, since these dictate the resolution of
your mesh around the imported geometry.

```
DX0     0.01		# dx0 minimum mesh step in x (leading & trailing edges)
DY0     0.01		# dy0 minimum mesh step in y (wing tip)
DZ0     0.05		# dz0 minimum mesh step in z (wing surface)
```

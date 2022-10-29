# Prandtline Instructions

This directory contains all the **required files to have in your working directory** before running. These
files are listed below:

* **prandtline.data**
* **polarbl.dat**

Typically the prandtline input data file is named **prandtline.data**; however you 
may change which file to use, by using the ```--prandtl_in``` flag. There is an example of this below. 

If you wish to integrate the ****wing**** and ****fuselage**** ensemble then read in the 9 polars named, **polarbl[n].dat**. These files can be 
found in the ```individual_polars/``` directory.

#### Polar File Format

The required format for the input polar file is demonstrated like so.
Typically the number of header lines are 12. The fortran code only expects this many, however
the cpp code can recognize and screen all headerline out automatically.

```
    [alpha]  [cz]   [cx]    [dummyval]      [cq]
     ...     ...    ...        ...         ...
     ...     ...    ...        ...         ...
  
     ...     ...    ...        ...         ...
    [left break-point] [right break-point]
```

Alpha is in degrees, and the rest of the drag coefficients are acquired from XFoil. 

Users must manually add in the left and right break points. Keep in mind the gamma integro equation is solved for the domain
[-1,1]. Essentially your leftmost break point cannot exceed -1, and the right most break point
cannot exceed +1.

## Example:
You may run this example using two methods. Using the executables or with the AeroDes.hpp library.

### Using PrandtLine Exectuable
By default the executable automatically reads in a parameter file named, *prandtline.data*, and a polar file named, *polarbl.dat*:
```
ADprandtline
```

specify the input **parameter** file to whatever file you want:
```
ADprandtline --prandtl_in prandtline.data
```

specify the input **polar** file to whatever file you want:
```
ADprandtline --prandtl_polar polarbl.dat
```

you can even use these flags simultaneously:
```
ADprandtline --prandtl_in prandtline.data --prandtl_polar polarbl.dat
```
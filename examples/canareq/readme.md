# Canareq Instructions

This directory contains all the **required files to have in your working directory** before running. These
files are listed below:

* **canareq.data**
* **canarpolar.dat**

Typically the canarline input data file is named **canarline.data**; however you may change it using the ```--ceq_in``` flag (shown below). 

## Example:
You may run this example using two methods. Using the executables or with the AeroDes.hpp library.

### Using Canary Exectuable
By default the executable automatically reads in a parameter file named, *canareq.data*, and a polar file named, *canarpolar.dat*.
Here we are setting the canard angle to 1 degree.
```
canary --ceq_tcd 1
```

specify the input **parameter** file to whatever file you want:
```
canary --ceq_tcd 1 --ceq_in canareq.data
```

specify the input **polar** file to whatever file you want:
```
canary --ceq_tcd 1 --ceq_polar canarpolar.dat
```

you can even use these flags simultaneously:
```
canary --ceq_tcd 1 --ceq_in canareq.data --ceq_polar canarpolar.dat
```

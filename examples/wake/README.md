# Canard Wake

Please inspect the wake.data file for configurational settings. The following is a showcase of how to run the canard wake executable. 
**Required files to have in your working directory** are as follows:

* wake.data 
* polarbl.dat
* wing.yxlexte 
* geocanard.xzmses

***wake.data***: is your input parameter file

***polarbl.dat***: is your input polar file

***wing.yxlexte***: is your wing geometry

***geocanard.xzmses***: is your canard geometry

## Running the code:
Basic usage (assumes input filename is "wake.data")

```
./../../build/wake
```

Custom usage (user defined input filename is "custom.data")

```
./../../build/wake -in custom.data
```

## Results:

Using the current data files, a preview of the results are shown here. The following figures demonstrate the canard wake and distributions that should be produced.

### Canard wake

![This is an image](https://github.com/carlos-pereyra/AeroDes/blob/master/docs/wake/canarwake_xz.png)

### Distributions

![This is an image](https://github.com/carlos-pereyra/AeroDes/blob/master/docs/wake/prandtline_ygw.png)

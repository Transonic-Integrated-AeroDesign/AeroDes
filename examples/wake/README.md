# Canard Wake

Please inspect the wake.data file for configurational settings. The following is a showcase of how to run the canard wake executable. Essential files to have in your directory include the following.

* wake.data       (input parameter file)
* polarbl.dat
* wing.yxlexte
* geocanard.xzmses

## Running the code:

Basic usage (assumes input filename is "wake.data")

```
./wake
```

Custom usage (user defined input filename is "custom.data")

```
./wake -in custom.data
```

## Results:

Using the current data files, a preview of the results are shown here. The following figures demonstrate the canard wake and distributions that should be produced.

### Canard wake

![This is an image](https://github.com/carlos-pereyra/AeroDes/blob/master/docs/wake/canarwake_xz.png)

### Distributions

![This is an image](https://github.com/carlos-pereyra/AeroDes/blob/master/docs/wake/prandtline_ygw.png)

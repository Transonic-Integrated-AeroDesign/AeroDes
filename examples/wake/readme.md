# Canard Wake

Please inspect the wake.data file for configurational settings. The following is a showcase of how to run the canard wake executable. 
**Required files to have in your working directory** are as follows:

* **wake.data** 
* **polarbl.dat**
* **wing.yxlexte**
* **geocanard.xzmses**

By default the ```wake``` executable reads a parameter file named, **wake.data** in your current working directory. 
Similiarly, the default polar filename is **polarbl.dat**.

Both of these options may be user specified using commandline flags (shown below).

Currently the default filenames for the two geometry files are **wing.yxlexte** (your wing geometry) and **geocanard.xzmses** (your canard geometry).

## Running the code:
Basic usage which assumes input filename is "wake.data".

```
wake
```

Custom usage where user defined input filename is "custom.data".

```
wake --wk_in wake.data
```

Custom user defined input **polar** file is whatever file you want.

```
wake --wk_polar polarbl.dat
```
Then of course all of these input flags can be used simultaneously.
```
wake --wk_in wake.data --wk_polar polarbl.dat
```

## Results:

Using the current data files, a preview of the results are shown here. The following figures demonstrate the canard wake and distributions that should be produced.

### Canard wake

![This is an image](https://github.com/carlos-pereyra/AeroDes/blob/master/docs/wake/canarwake_xz.png)

### Distributions

![This is an image](https://github.com/carlos-pereyra/AeroDes/blob/master/docs/wake/prandtline_ygw.png)

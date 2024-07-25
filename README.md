## CBraid and Braiding

TODO.

## Installation

To compile the example programs, from the base folder run
```
cd programs; make clean; make
```
This will create the executable `braiding` in the `programs` folder, as well as the `libcbraid.a` library in the `libs` folder.

By default, `braiding` will use Artin braids. To change this and use `XXX` instead, compile instead with 
```
cd programs; make clean; make USE_FOR_BRAIDING=XXX
```
Currently supported: Artin classic Braid groups, classic (`XXX` = `ARTIN`) and dual (`XXX` = `BAND`) Garside structures, B-series Artin groups (dual Garside structure) (`XXX` = `OCTAHEDRAL`), I-series Artin groups (dual Garside structure) (`XXX` = `DIHEDRAL`), and complex reflection braid groups B(e, e, n), dual (`XXX` = `DUAL_COMPLEX`) and semi-classic (`XXX` = `STANDARD_COMPLEX`, NOT FULLY WORKING) Garside structures.

To compile just the library `libcbraid.a`, from the base folder run
```
cd lib; make clean; make
```

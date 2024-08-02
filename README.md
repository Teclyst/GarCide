# GarCide (prerelease)

## Overview

_Garcide_ is a C++ library and application for computations in Garside groups. It was written during a L3 internship, under the supervision of **[Bert Wiest](https://perso.univ-rennes1.fr/bertold.wiest/)**.

The project contains two parts: on one side, the _GarCide_ library, that provides functions for various computations in Garside groups, and is meant to allow users to easily define classes for their own Garside groups.

On the other side, the _Braiding 2_ executable, provides a shell interface to easily run some calculations..

This library is based on 2001 _CBraid_ by **[Jae Choon Cha](http://gt.postech.ac.kr/~jccha/)** and 2004 _Braiding_ by **[Juan Gonzalez-Meneses](http://personal.us.es/meneses/)**, with contributions from **[Maria Cumplido](https://personal.us.es/cumplido/)**.

These two projects, maintained by **[Jean-Luc Thiffeault](http://www.math.wisc.edu/~jeanluc)**, can be found at [the repository this is a fork of](https://github.com/jeanluct/cbraid).

Although names may suggest otherwise, _GarCide_ roughly corresponds to a merge of _CBraid_ and computation-oriented parts of the original _Braiding_, with _Braiding 2_ only encompassing shell-handling parts.

## Building the project

### The first time

0)  You will need _CMake_ to compile the project.

    If you don't have it, you can get it with:

    ```
    sudo apt install cmake
    ```

1)  Create a new directory called `build` where to run _CMake_ from.
    ```
    mkdir build; cd build
    ```
    Henceforth only use _CMake_ from that `build` directory.

2)  Create the _CMake_ binary tree with: 
    ```
    cmake ..
    ```

### Building it

1)  Set options (see next subsection) with:
    ```
    cmake [options] ..
    ```
    Skip if you don't need to change options.

2)  Build with
    ```
    cmake --build .
    ```

3)  The _GarCide_ library can then be found at `build/lib/libgarcide.a`, and the _Braiding 2_ executable at `build/src/braiding.exe`.

    To run the latter, use
    ```
    src/braiding.exe
    ```

### Options

_CMake_ does not change its cached variables between runs. Therefore a binding will remain until explicitly changed: e.g. after running
```
cmake -DUSE_PAR=FALSE ..
```
`USE_PAR` will be set to `FALSE` for all subsequent calls to anything _CMake_ related that does not explicitly change `USE_PAR`.

Now for the available options (in bold are the values before any change is ever done):

*   `USE_PAR` (possible values __`TRUE`__, `FALSE`) - Whether parallelism should be used. Notice that you need Intel's TBB for the parallel code to actually work, so this option will not do anything if you do not have TBB.

*   `RANDOMIZE_ON_WORDS` (possible values `TRUE`, __`FALSE`__) - Whether randomizing braids should be understood as taking a random word in the atoms of a given length. Notice that for many Garside groups it is hard to provide better polling methods.

    You can ignore that option if you only care about _Braiding 2_.

*   `USE_FOR_BRAIDING` (possible values __`ARTIN`__, `BAND`, `OCTAHEDRAL`, `DIHEDRAL`, `DUAL_COMPLEX`, `STANDARD_COMPLEX`, `EUCLIDEAN_LATTICE`) - Selects which group should be used for _Braiding 2_.

    Currently supported: 
    * Regular braid groups (a.k.a. A-series Artin groups), classic Garside structure (`ARTIN`).
    * Regular braid groups, dual Garside structure (`BAND`).
    * B-series Artin groups, dual Garside structure (`OCTAHEDRAL`)
    * I-series Artin groups, dual Garside structure (`DIHEDRAL`)
    * Complex reflection braid groups $B(e, e, n)$, dual Garside structure (`DUAL_COMPLEX`)
    * Complex reflection braid groups $B(e, e, n)$, semi-classic Garside structure (`STANDARD_COMPLEX`), NOT FULLY WORKING AS OF NOW.
    * Euclidean lattices $\mathbb Z^n$ (`EUCLIDEAN_LATTICE`).

## Implementing a Garside group using GarCide

See `doc/implementing_garside_groups.md`.

## Generating documentation

TODO
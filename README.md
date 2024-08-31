# _GarCide_

[![GitHub release (latest by date)](https://img.shields.io/github/v/release/Teclyst/GarCide)](https://github.com/Teclyst/GarCide/releases/latest)
[![GitHub](https://img.shields.io/github/license/Teclyst/GarCide)](https://github.com/Teclyst/GarCide/blob/master/LICENSE.md)

_Garcide_ is a C++ library and application for computations in Garside groups. It was written during a L3 internship, under the supervision of **[Bert Wiest](https://perso.univ-rennes1.fr/bertold.wiest/)**.

## Overview

The project contains two parts: on one side, the _GarCide_ library, that provides functions for various computations in Garside groups, and is meant to allow users to easily define classes for their own Garside groups.

On the other side, the _Braiding_ executable provides a shell interface to easily run computations.

This library is based on 2001 _CBraid_ by **[Jae Choon Cha](http://gt.postech.ac.kr/~jccha/)** and 2004 _Braiding_ by **[Juan González-Meneses](http://personal.us.es/meneses/)**, with contributions from **[Maria Cumplido](https://personal.us.es/cumplido/)**.

These two projects, maintained by **[Jean-Luc Thiffeault](http://www.math.wisc.edu/~jeanluc)**, can be found [here](https://github.com/jeanluct/cbraid).

Although names may suggest otherwise, _GarCide_ roughly corresponds to a merge of _CBraid_ and computation-oriented parts of the original _Braiding_, while _Braiding_ only encompasses shell interaction.

## Requirements

* [_**CMake**_](https://cmake.org/), to compile the project. On Linux you can get it with

    ```shell
    sudo apt install cmake
    ```

* Intel's [_**TBB**_](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onetbb.html), needed for multithreading. Not strictly necessary, although using it often results in significant speedup. On Linux you can get it with

    ```shell
    sudo apt install libtbb-dev
    ```

* [_**Doxygen**_](https://www.doxygen.nl/index.html), to generate documentation. If you don't care about it, you can ignore it, but if you plan to play around with the code it is probably a good idea to get it. On Linux, you can do this with

    ```shell
    sudo apt install doxygen
    ```

## Building the project

### The first time

1) Create a new directory called `build` where to run _CMake_ from.

    ```shell
    mkdir build; cd build
    ```

    Henceforth only use _CMake_ from that `build` directory.

2) Create the _CMake_ binary tree with:

    ```shell
    cmake ..
    ```

### Building it

1) Set options (see next subsection) with:

    ```shell
    cmake [options] ..
    ```

    To change an option, use `-D[option]=[value]`.

    Skip if you don't need to change options.

2) Build with

    ```shell
    cmake --build .
    ```

3) The _GarCide_ library can then be found at `build/lib/libgarcide.a`, and the _Braiding_ executable at `build/src/braiding.exe`.

    To run the latter, use

    ```shell
    src/braiding.exe
    ```

### Options

_CMake_ does not change its cached variables between runs. Therefore a binding will remain until explicitly changed: _e.g._ after running

```shell
cmake -DUSE_PAR=FALSE ..
```

`USE_PAR` will be set to `FALSE` for all subsequent calls to anything _CMake_ related that does not explicitly change `USE_PAR`.

Now for the available options (in bold are the values before any change is ever done):

* `USE_PAR` (possible values **`TRUE`**, `FALSE`) - Whether parallelism should be used. Notice that as you need Intel's _TBB_ for the parallel code to actually work, this option will not do anything if _TBB_ is not on your machine.

    Parallelism should significantly speed up computations (for super summit, ultra summit and sliding circuits sets and centralizers) in most case, but may have the opposite effect for small cases (_e.g._ very small number of strands and very short braids) and depending on architecture.

* `RANDOMIZE_ON_WORDS` (possible values `TRUE`, **`FALSE`**) - Whether randomizing braids should be understood as taking a random word in the atoms of a given length. Notice that for many Garside groups it is hard to provide better polling methods.

    You can ignore that option if you only care about _Braiding_.

* `USE_FOR_BRAIDING` (possible values **`ARTIN`**, `BAND`, `OCTAHEDRAL`, `DIHEDRAL`, `DUAL_COMPLEX`, `STANDARD_COMPLEX`, `EUCLIDEAN_LATTICE`) - Selects which group should be used for _Braiding_.

    Currently supported:
  * Regular braid groups (a.k.a. $\mathbf A$-series Artin groups), classic Garside structure (`ARTIN`).
  * Regular braid groups, dual Garside structure (`BAND`).
  * $\mathbf B$-series Artin groups, dual Garside structure (`OCTAHEDRAL`).
  * $\mathbf I$-series Artin groups, dual Garside structure (`DIHEDRAL`).
  * Complex reflection braid groups $\mathrm B(e, e, n)$, dual Garside structure (`DUAL_COMPLEX`).
  * Complex reflection braid groups $\mathrm B(e, e, n)$, semi-classic Garside structure (`STANDARD_COMPLEX`), NOT FULLY WORKING AS OF NOW.
  * Euclidean lattices $\mathbb Z^n$ (`EUCLIDEAN_LATTICE`).

* `GENERATE_DOC` (possible values **`TRUE`**, `FALSE`) - whether documentation should be generated when building the project.

* `DOXYGEN_WARNINGS` (possible values `YES`, **`NO`**) - whether _Doxygen_ should be allowed to output warnings.

* `DOXYGEN_QUIET` (possible values **`YES`**, `NO`) - whether _Doxygen_ should avoid outputing build details.

* `CMAKE_BUILD_TYPE` (possible values `Debug`, **`Release`**) - whether the project should be built in debug mode (debug symbols, no compiler optimizations, better for development) or release mode (compiler optimizations, no debug symbols).

    Build with the latter for benchmarking.

## Implementing a Garside group using _GarCide_

See `doc/implementing_garside_groups.md`.

## Generating documentation

If you have _Doxygen_, and assuming that `GENERATE_DOC` is set to `TRUE`, then documentation will be automatically generated when building the project.

To read it, open `build/html/index.html` with a browser. By the command line and with _Firefox_,

```shell
firefox build/html/index.html
```

The _Doxygen_ pages use the _**[Doxygen Awesome](https://github.com/jothepro/doxygen-awesome-css)**_ theme by **[jothepro](https://github.com/jothepro)**.

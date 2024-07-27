## GarCide (prerelease)

This project (written in <img src = https://raw.githubusercontent.com/devicons/devicon/blob/master/icons/cplusplus/cplusplus-plain.svg height = 16>) contains two parts: on one side, the _GarCide_ library, that allows various computations on Garside groups (notably super summit, ultra summit and sliding circuits set calculations). The idea behind this project is that since algorithms for Garside groups are totally generic once base factor operations have been taken care of, it should be easy and fast to implement a Garside group.

The goal of _GarCide_ is to give anyone (even someone that does not have extensive programming experience) who wishes to do calculations in a Garside group a tool to do so, with (relatively) little work having to be done; it is the author's opinion that such a tool should focus on efficiency rather than stark minimality, and this reflects on the way things are implemented.

For a guide on how to implement a Garside Group using _GarCide_, see `/doc`.

The other (much smaller) part is the _Braiding 2_ executable. It provides a shell interface to easily run some calculations, and is easy to specialize for Garside group implementation based on _GarCide_ (see, again, `/doc`).

This library is based on 2001 _CBraid_ by **[Jae Choon Cha](http://gt.postech.ac.kr/~jccha/)** and 2004 _Braiding_ by **[Juan Gonzalez-Meneses](http://personal.us.es/meneses/)**, with contributions from **[Maria Cumplido](https://personal.us.es/cumplido/)**.

These two projects, maintained by **[Jean-Luc Thiffeault](http://www.math.wisc.edu/~jeanluc)**, can be found at TODO [the repository this is a fork of](https://github.com/jeanluct/cbraid).

Although names may suggest otherwise, _GarCide_ roughly corresponds to a merge of _CBraid_ and computation-oriented parts of the original _Braiding_, with _Braiding 2_ only encompassing shell-handling parts.

## Installation

To compile the braiding executable, from the base folder run
```
cd programs; make clean; make
```
This will create the executable `braiding` in the `programs` folder, as well as the `libcbraid.a` library in the `libs` folder.

By default, `braiding` will use braids from the regular braid group, with its classic Garside structure as an Artin Group. To change this and use `XXX` instead, compile with
```
cd programs; make clean; make USE_FOR_BRAIDING=XXX
```
Currently supported are regular braid groups (a.k.a. A-series Artin groups), classic (`XXX` = `ARTIN`) and dual (`XXX` = `BAND`) Garside structures, B-series Artin groups (dual Garside structure) (`XXX` = `OCTAHEDRAL`), I-series Artin groups (dual Garside structure) (`XXX` = `DIHEDRAL`), and complex reflection braid groups B(e, e, n), dual (`XXX` = `DUAL_COMPLEX`) and semi-classic (`XXX` = `STANDARD_COMPLEX`, NOT FULLY WORKING AS OF NOW) Garside structures.

To compile just the library `libcbraid.a`, from the base folder run
```
cd lib; make clean; make
```
This will not compile any of the group implementations. You must opt in to compiling them (to reduce compilation time). To compile a given group implementation, use option `USE_XXX=1` (this is not needed for the `XXX` passed to `USE_FOR_BRAIDING` if compiling from `/programs`).

Note that the library uses parallel algorithms by default, using the Intel library Threading Building Blocks (TBB), to speed up calculations. Therefore compilation will fail if you do not have it. To fix this, you can either get it, or compile with option `USE_PAR=0` to disable parallelism altogether.

You may want to run tests with random braids (typically for benchmarking purposes), but it can be hard to implement uniform polling over factors. If you want to be able to get reasonably random braids of a given size no matter what, but do not care for uniformity over factors, then you can choose instead to randomize braids as random words in atoms. To do this, compile with option `RANDOMIZE_ON_ATOMS=1`. Note that you do not need to worry about randomization if you only care about getting _Braiding_ to work with your implementation.
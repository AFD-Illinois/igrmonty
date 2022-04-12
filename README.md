# igrmonty

A monte carlo code for relativistic radiative transport, 
including synchrotron emission and absorption, bremsstrahlung
emission (absorption is negligible), and compton scattering.

For publications, please cite this journal article:
Dolence, J.C., Gammie, C.F., Moscibrodzka, M., and Leung, P.K.: 2009, ApJS 184, 387.
https://ui.adsabs.harvard.edu/abs/2009ApJS..184..387D/abstract
And this repository:
https://github.com/AFD-Illinois/igrmonty

Substantial contributions made by
Ben Ryan
George Wong
Ricardo Castro-Yarza
Chi-kwan Chan

prepped for public release
cfg 18 Apr 2020


## quick-and-dirty overview

... can be found [here](https://github.com/AFD-Illinois/igrmonty/blob/master/docs/tutorial.pdf).

## tests

Tests of the bremsstrahlung emission functions can be built and run with cmake:

```bash
$ mkdir build && cd build
$ cmake ..
$ make
$ make test
```

This also builds the code.  YMMV using this build system vs the standard `make` solution.

## adding new electron distribution functions

1. add new names/parameters to src/model_radiation.h

2. add absorptivity to src/radiation.c:alpha_inv_abs(...)

3. add eDF in src/hotcross.c:dNdgammae_[X], modifying dNdgammae as appropriate

4a. modify eDF functions in src/compton.c
    - dfdgam(...)
    - fdist(...)  -- note overlap with dNdgammae in src/hotcross.c

4b. ensure src/compton.c:sample_beta_distr_num(...) produces reasonable results on edge cases

5a. implement jnu_[X] and int_jnu_[X] in src/jnu_mixed.c. it's likely the integral calls into a table

5b. modify switches in src/jnu_mixed.c  -- jnu(...), jnu_ratio_brems(...), int_jnu(...)



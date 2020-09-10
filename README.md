# Neutrino Tomography of Earth

## Description:
Set of codes to solve the complex diff. equation of neutrino propagation through vacuum/material
combined with a customized discretized Earth model with Radioactive Isotope distributon and Density models.

The Geoneutrino flux model is based on the previous work:

O. Sramek et al. Earth and Planetary Science Letters 361 (2013) 356-366.

The Neutrino Oscillation in vacuum/material was written by A. Rubbia and A. Bueno.

An example work based on this code: 

A. Bueno, M. Campanelli, A. Rubbia Nuclear Physics B 589 (2000) 577-608.


## Build instrunctions:
```
$ mkdir build
$ cd build
$ cmake3 ../nutomo
$ make
```

## Then try to run the various executables in the `Apps` folder

Example apps:
 * `Run_NuFlux_Simple.cpp` : fast calculation of antineutrino flux at a given surface position integrating over the whole planet
 * `Run_NuFlux_OscLongitude.cpp` : calculation of directed antineutrino flux at a given surface position from a plane at a given Longitude

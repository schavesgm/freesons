# Free mesons at non-zero temperature
This repository contains a code implementing the calculation of
mesonic correlation functions for different channels at non-zero
temperature. The meson correlators are calculated in the free theory,
therefore, no gauge interactions are present.

The lattice where the fermion fields are located is an hypercube of `4`
dimensions, three space directions and one time direction. The
boundary conditions in the time direction are anti-periodic while
being periodic in the spatial directions. The anti-periodicity in the
time direction puts the system in a thermal bath of temperature `T`.
The temperature is defined as `T = 1 / (Nt * at)`, where `Nt` is the
number of points in the time direction and `at` is the lattice spacing
in the time direction. The implemented lattice can be anysotropy,
therefore the lattice spacing in the time direction `as` and the lattice
spacing in the temporal direction `at` might not be equal. The factor
of anysotropy is defined as `X = as / at`. The code uses `at = 1`,
therefore `as = X * at`.

The code calculates correlation functions for different operators in a
single run. Defined the number of points in the spatial direction
`Ns`, the number of points in the temporal direction `Nt`, the quark
mass in lattice units `mq * Ns` and some other parameters; the code
generates `8` different files in the `out` directory. The file whose
name contains `c0` is the scalar operator. The file containing `c1` is
the pseudoscalar. The next three files are the temporal vector, the
spatial vector and the spatial + temporal vector. The last three files
are the temporal axial vector, the spatial axial vector and both the
temporal and the axial vectors.

The code can be compiled using `make`,

```bash
make
```

Information about how to run the code and the command line parameters
can be found by using,

```bash
./freesons --help
```

A standard run could be

```bash
./freesons -Ns 64 -Nt 24 -mq 0.64 -Xi 1.0
```

This run will produce a lattice of `64 x 64 x 64 x 24` sites. The
correlators will contain a mass of `0.01` lattice units and anysotropy
factor of `1.0`.

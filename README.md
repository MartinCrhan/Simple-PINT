# Simple-PINT
A simple FORTRAN code to run path integral simmulations of particles in external potentials. Uses normal modes and exact harmonic integration.

The program is compatible with OpenMP parallelism. Note, that during a parallel run a degraded performance is possible, when using low values of "stride" (such as stride = 1).

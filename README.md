# pes_bnplush2_sitespecific
A Fortran 90 interface for the site-specific potential energy surface for the H2 + BN dissociation reaction

#Currently implemented: Ortho site#
How to use?

the PES_h2bn.f90 file by default comes with a program to call a function energyget(X,position) which results in a 64-bit real number which is the energy interpolated on the PES, generating a value of the energy at a point.
We tested the OpenMP/MPI parallel execution of the PES driver routine as well, by running 4 MPI ranks at 4 OMP threads each, all generating correct results. 

How to compile? 
The necessary file for the neural network input generation schemes is pipnn.f90, which contains all the necessary subroutines to compute the symmetry functions of the PES. It's compiled as a Fortran modulefile and linked statically. We provide a Makefile which links the example program PES_h2bn.f90 against the module and compiles them together. If you're using it in a bigger program, include all other files BEFORE pipnn.f90, except the file which calls the function directly. REPLACE PES_h2bn.f90 with the name of your file

Compiler: intel oneapi ifx/mpiifx (both tested) 2025.1, intel MPI, intel MKL



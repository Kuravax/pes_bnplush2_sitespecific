# pes_bnplush2_sitespecific
A Fortran 90 interface for the site-specific potential energy surface for the H2 + BN dissociation reaction

#Currently implemented: Ortho site#
How to use?

the PES_h2bn.f90 file by default comes with a program to call a function energyget(X,position) which results in a 64-bit real number which is the energy interpolated on the PES, generating a value of the energy at a point.

Compiler: intel oneapi ifx 2025.1, intel MKL



# Hydrogen_polarizability_and_vdW_energy
The purpose of this code is to compute the polarizability of the Hydrogen atom and the van der Waals energy of a H-H dimer.

Second order perturbation theory is utilized to compuet these quantities.

The polarizability (and also the vdW energy) has a contribution from the bound-to-bound transitions and from the
bound-to-continuum transitions as well. These contributions are also considered separately and printed as output.

The computation of the bound and free state Coulomb wavefunction is based on:
"A Fortran program to calculate the matrix elements of the Coulomb interaction involving hydrogenic wave functions",
L. Sarkadi, Computer Physics Communications, 133, (2000), 119–127

Run compile.sh for compilation and for running:

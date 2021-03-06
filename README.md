# Hydrogen_polarizability_and_vdW_energy
The purpose of this code is to compute the static polarizability of the Hydrogen atom and the van der Waals energy of a H-H dimer. (The exact value of the static polariziabiity in atomic units is 9/2)

Second order perturbation theory is utilized to compute these quantities.

The polarizability (and also the vdW energy) has a contribution from the bound-to-bound transitions and from the
bound-to-continuum transitions as well. These contributions are also considered separately and printed as output.

The computation of the free state Coulomb wavefunction is based on:
"A Fortran program to calculate the matrix elements of the Coulomb interaction involving hydrogenic wave functions",
L. Sarkadi, Computer Physics Communications, 133, (2000), 119–127

Run ./compile.sh for compilation and for running:

./Hydrogen_vdW_version_3.x < input_1_okk_parameters.inp

The code can be easily modified to compute the frequency dependent polarizability too.

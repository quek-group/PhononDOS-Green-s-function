# PhononDOS-Green-s-function
How to use the code:
-------------------
1. INPUT: dyn_3Nx3N-bulk.dat and example-dyn_3Nx3N-connection.dat are from the folder "bulk"

format: 
3Nx3N matrix. N is the number of atoms; 3 refers to in x, y and z. Each element is FC_pi,qj/sqrt(M_p * M_q), where FC_pi,qj is the force constant between atom p (move along i) and atom q (move along j). i, j = x,y,z.  FC_pi,qj has the unit of Ry/bohr which is directly from .dyn file computed by QuantumESPRESSO.

dyn_3Nx3N-bulk.dat: dynamical matrix elements between atoms within one unitcell without peridic boundary. 
dyn_3Nx3N-connection: upper triangular matrix, containing interactions between atoms in two unitcell. Because at the bottom surface region contains some atoms from the bulk region, so this matrix is also used as the connection matrix betwene the surface region and the bulk region. So please make sure this is correct for your system. Otherwise you can modify the code to read in these two different connection matrices.


2. INPUT: dyn_3Rx3R-device.dat is from the folder "surface"
format: 3Rx3R matrix.

4. Some parameters need to be given in the code: 1> n_atm(number of atom in bulk and connection maxtrices, same); n_atm_dvc (number of atoms in device/surface); zero(positive infinitesimal parameter); omega, omega_max and omega_step (frequencies of the GF DOS); it_cycle(number of iterative cycle)
5. Module load numpy.
6. Run it.

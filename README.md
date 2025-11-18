
# VecStokesTest

Short test set up to test legacy vs vectorised stokes solver for idealised marine ice sheet in Elmer/Ice

## Note on running

Compile Elmer fortran:

elmerf90 FISOC_Elmer_geometries.f90 -o FISOC_Elmer_geometries.so

elmerf90 USF_CoV.F90 -o USF_CoV.so

Running on a laptop with mpirun:

mpirun -np 3 ElmerSolver AK_legacy.sif 

## Different sifs available

Three binary choices are present, resulting in 8 sifs.

Starting point for the sif was either Rupert's old FISOC example 5 (F5) or Ann-Kristin's recent Dibble sif (AK).

Stokes solver was either legacy or Vectorized Stokes (VS).

Run was either transient (default) or steady state (SS).

## Differences in outputs

F5 vs AK gives small differences that we're not worried about.

Transient vs steady gives no discernable difference.

Legacy vs VS gives huge differences, especially in the shelf. More specifically, the legacy solver looks ok but the VS solver gives massive velocities out of the domain through the calving front. The difference decreases toward the grounded part of the domain, where the differences are small. Could the two solvers be interpreting the calving BC differently somehow?


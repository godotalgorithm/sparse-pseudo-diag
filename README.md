# sparse-pseudo-diag

This repository contains a simple reference implementation of the sparse pseudo-diagonalization algorithm that is used by the fast MOZYME solver of [MOPAC](https://github.com/openmopac/mopac) along with two example matrices (for the Crambin protein and for a linear polyethylene polymer).

The current version replicates the simple two-pass structure of the MOZYME solver, whereby the Fock matrix is projected into the local molecular orbital (LMO) basis and that projected matrix is used to define Jacobi rotations applied to the LMO basis for the purpose of decoupling the occupied and virtual LMOs. However, there are additional, rather complicated mechanisms in MOZYME for damping rotations and pruning the list of transformed occupied-virtual pairs that are not yet replicated here. These mechanisms are not needed for the polymer example, but they will be needed to effectively maintain sparsity and converge the protein example. As such, this implementation is a reasonable proxy for the computational workload of MOZYME calculations, but it is not yet a fully reliable solver algorithm.

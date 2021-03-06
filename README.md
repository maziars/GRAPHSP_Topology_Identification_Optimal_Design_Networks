# Topology Identification and Optimal Design of Noisy Networks
This is a Memory and Cache Efficient C/C++ Implementation of Customized Algorithms for Topology Identification and Optimal Design of Noisy Networks. It uses the compiled Frotran BLAS binary files to improve the speed of the linear algebraic computations. So, if you want to use this code you need to build the BLAS package for the CPU architecture you use and link the binary files in the project.

See: http://people.ece.umn.edu/users/mihailo/software/graphsp/index.html for more details on the problem formulation and the methods. The code implemets three different methods for finding the topology identification and optimal design of undirected consensus networks with additive stochastic disturbances (see this paper for details: http://www-bcf.usc.edu/~mihailo/papers/mogjovTCNS17.pdf):

- the infeasible primal-dual IP method;
- the proximal gradient method;
- the proximal Newton method.

The proximal gradient updates the controller graph Laplacian via convenient use of the soft-thresholding operator. In the IP method, the Newton direction is obtained using an inexact iterative procedure based on the preconditioned conjugate gradients and, in the proximal Newton method it is computed using cyclic coordinate descent over the set of active variables. 
This C/C++ proximal gradient implementation has been used to solve the poblem for graphs with millions of edges on single computer in minutes (see: http://www-bcf.usc.edu/~mihailo/papers/mogjovTCNS17.pdf)

Explanation of the files:
- Matrix_save.m: this is a Matlab file that saves the info of the graph into the txt files such that it could be used by the C/C++ code

- Constants.h : includes the constants needed by algorithms (stopping criteria, max # iterations, ...) as well as debugging and using flags to compile stuff in or out

- main.c : runs the algorithms. For simplicity, the main assumes that the .txt folders containing the Edges (E), Laplacian (L), M and N (MN) and $\gamma$ (gamma) is in the same folder. See the matlab save_matrix.m file to see how the files are read into memory and should be formatted. 

- IO_manager files: read the graph information from the file

- D_comp: This file contains the main dense matrix computation routines that use BLAS routines

- PCG_Solver: Contains the diagonaly preconditioned conjugate gradient method for finding the Newtone step in IP

- IP_solver: Contains the routin for doing interior point method

- Proximal_solver: it contains two implementations of proximal gradien method with BB step size. The first one, PROX1 uses customized sparse implementation. The second one, PROX2 does not use the sparse implementations. The sparse implementation still uses BLAS routines and does not use sparse versions.

-SparseTools: A simple implementation of sparse class using iterators. Note that we still use the BLAS methods for dense matrices to do linear algebra for sparse vectors using customized wrappers.

- QuadCD: contains the implementations of Newton Proximal method with backtracking for step-size selection. It uses Coordinate descent with active set implementation to find the Newton direction.




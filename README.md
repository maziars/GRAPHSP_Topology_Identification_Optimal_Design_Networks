# Topology Identification and Optimal Design of Noisy Networks
This is a Memory and Cache Efficient C/C++ Implementation of Customized Algorithms for Topology Identification and Optimal Design of Noisy Networks. It uses the compiled Frotran BLAS binary files to improve the speed of the linear algebraic computations. So, if you want to use this code you need to build the BLAS package for the CPU architecture you use and link the binary files in the project.

See: http://people.ece.umn.edu/users/mihailo/software/graphsp/index.html for more details on the problem formulation and the methods. The code implemets three different methods for finding the topology identification and optimal design of undirected consensus networks with additive stochastic disturbances (see this paper for details: http://www-bcf.usc.edu/~mihailo/papers/mogjovTCNS17.pdf):

- the infeasible primal-dual IP method;
- the proximal gradient method;
- the proximal Newton method.

The proximal gradient updates the controller graph Laplacian via convenient use of the soft-thresholding operator. In the IP method, the Newton direction is obtained using an inexact iterative procedure based on the preconditioned conjugate gradients and, in the proximal Newton method it is computed using cyclic coordinate descent over the set of active variables. 

Explanation of the files:
- Constants.h : includes the constants needed by algorithms (stopping criteria, max # iterations, ...) as well as debugging and using flags to compile stuff in or out

- main.c : runs the algorithms. For simplicity, the main assumes that the .txt folders containing the Edges (E), Laplacian (L), M and N (MN) and $\gamma$ (gamma) is in the same folder. See the IO manager to see how the files are read into memory and should be formatted. 

- IO_manager files: read the graph information from the file

- D_comp: This file contains the main dense matrix computation routines that use BLAS routines

- PCG_Solver: Contains the diagonaly preconditioned conjugate gradient method for finding the Newtone step in IP

- IP_solver: Contains the routin for doing interior point method

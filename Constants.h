//
//  Constants.h
//  PCG
//
//  Created by Maziar Sanjabi on 12/27/15.
//  Copyright © 2015 Maziar Sanjabi. All rights reserved.
//
// This file contains the constants to be defined before running the algorithms


#ifndef Constants_h
#define Constants_h

// Constants for IP Solver
#define EPS_DG          1.e-6// The duality gap stoping criterion for IP solver.
#define EPS_RES         1.e-3// The residual stoping criterion for IP solver. The IP solver stops only if residual and duality gap is below these thresholds
#define MAXIT           100 // Max number of iterations for IP solver
#define MAXIT_PCG       1000// Max number of Conjugate gradient steps to perform in each step of IP solver
#define DELTA_PCG_AFF   0.3 // Constant for Preconditioned Conjugate gradient method for solving linear system of equation
#define DELTA_PCG       0.3// Constant for Preconditioned Conjugate gradient method for solving linear system of equation



// Constants for Proximal Gradient Method
#define EPS_ISTA        1.e-4 // The stoping criteria for Proximal solver. It stops when the size of the proximal step is below this threshold


// Constants for Proxiaml Newtone method that uses CD with active set to find NEwton direction
#define QUADCD_MAXIT    1000// Maximum number of QUADCD iterations
#define QUADCD_TOL3     0.001//Constant for stoping the inside CD for finding the Newton Direction
#define EPS_CD          0.000001 // Stoping criteria for QUADCD. Stops when the objective difference between two consecutive steps is below this threshold



// Debugging flags

//#define QUAD_DEBUG_ON //Enables debugging on QUAD_CD
//#ifndef DEBUG_ON
//#define DEBUG_ON // Enables Debugging in all files
//#endif
#ifdef DEBUG_ON
//#define DEBUG_ON_DEEP // Deep debugging print many intermediate steps
#endif

//#define INIT1 // Use alternative method for initalizing the IP solver method


// If defined for each method, it is compiled and used in the main file to solve the problem
#define IPSOLVER
#define PROX1
#define PROX2
#define QUAD_CD
#endif /* Constants_h */

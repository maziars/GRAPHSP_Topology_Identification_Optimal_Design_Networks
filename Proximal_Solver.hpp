//
//  Proximal_Solver.hpp
//  PCG
//
//  Created by Maziar Sanjabi on 1/3/16.
//  Copyright © 2016 Maziar Sanjabi. All rights reserved.
//

#ifndef Proximal_Solver_hpp
#define Proximal_Solver_hpp



#include <stdio.h>
#include "SparseTools.hpp"
#include "D_Comp.hpp"
#include "Constants.h"


//#define BETA        0.5
#define MAX_BT      20 // maximum number of backtracking steps (if used)

void Proximal_Solver(double *Lp,int* E, double* Q,double gamma, SparseVec *xs, double* y, int m, int n);

void Proximal_Solver2(double *Lp,int* E, double* Q,double gamma, double *xs, double* y, int m, int n);
double Compute_Obj2(double* Ginv,double* Q, double* xs, double gamma, int m , int n);
void Compute_grad(double* Y, double* grad, int* E, double gamma, int m, int n);

#endif /* Proximal_Solver_hpp */

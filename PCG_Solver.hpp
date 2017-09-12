//
//  PCG_Solver.hpp
//  PCG
//
//  Created by Maziar Sanjabi on 12/24/15.
//  Copyright © 2015 Maziar Sanjabi. All rights reserved.
//

#ifndef PCG_Solver_hpp
#define PCG_Solver_hpp

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "D_comp.hpp"
#include "Constants.h"




int diagonal_PCG(double* Y, int *E, double *Ginv, double*x, double* b, int m, int n, double e, double* D, double *xbar, double *ybar);
void Compute_grad(double* Y, double* grad, int* E, double gamma, int m, int n);
#endif /* PCG_Solver_hpp */

//
//  QuadCD.hpp
//  PCG
//
//  Created by Maziar Sanjabi on 3/30/16.
//  Copyright © 2016 Maziar Sanjabi. All rights reserved.
//

#ifndef QuadCD_hpp
#define QuadCD_hpp

#include <stdio.h>
#include "Constants.h"

#define BETA 0.5 // step-size reduction for backtrcking
#define ALPHA 0.3 // amount of acceptable reduction in objective value for back tracking
#define MAX_BACK_TRACK 20 // max number of backtrack steps


void Quad_Solver(double *Lp,int* E, double* Q,double gamma, double* x, double* y, int m, int n);




#endif /* QuadCD_hpp */

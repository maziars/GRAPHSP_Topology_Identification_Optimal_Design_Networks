//
//  IP_Solver.hpp
//  PCG
//
//  Created by Maziar Sanjabi on 12/24/15.
//  Copyright © 2015 Maziar Sanjabi. All rights reserved.
//

#ifndef IP_Solver_hpp
#define IP_Solver_hpp

#include <stdio.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include "D_Comp.hpp"
#include "PCG_Solver.hpp"
#include "Constants.h"





void IP_Solver(double *Lp,int* E, double* Q,double gamma, double *x, double* y, int m, int n);


#endif /* IP_Solver_hpp */

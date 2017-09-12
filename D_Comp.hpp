//
//  D_Comp.hpp
//  
//
//  Created by Maziar Sanjabi on 12/13/15.
//  Copyright © 2015 Maziar Sanjabi. All rights reserved.
//

#ifndef D_Comp_hpp
#define D_Comp_hpp

#include <stdio.h>
#include <iostream>
#include <cmath>
#include "Constants.h"


double innder_prod(double* a, double* b, int n);
void matrix_prod(double* A, double* x, double* y, int m, int n);
void vector_sub(double* x, double* y, double* z, int n);
void vector_add(double* x, double* y, double* z, int n);
void vector_deepcopy(double* x, double* y, int n);
void CG(double* A, double* x, double* b, int m, double e);
void vector_weight_sum(double* x, double* y, double* z, double alpha, int n);
void vector_print(double* x, int n);
void vector_prod(double* a, double* b, double*c, int n);
void PCG(double* A, double*x, double* b, int m, double e, double* D1, double* D2, double* D3);
//double *cholesky(double *A, int n);
void cholesky(double *A, int n);
void show_matrix(double *A, int n);
void compute_G(double *Lp, int* E, double* G, double *x, int m, int n);
void G_inv(double *G, double *Ginv, int n);
void show_vector(double *v, int n);
void Compute_Y(double *Ginv, double *Qp, double *Y, int n);
void Compute_rd(double *Y, double *y, int* E, double gamma, int m, int n, double* rd);
void Compute_Ap(double *Y, int *E, double *Ginv, double *p, int m, int n, double *Ap, double *x, double *y);
void Compute_diagA(double *Y, double *Ginv, double *x, double *y, int *E, int m, int n, double* D);
void Compute_rd_b(double *Y, double *y, int* E, double gamma, int m, int n, double* rd, double* b);
void Compute_b(double *Y, double *y, int* E, double gamma, int m, int n, double* b);
void Initialize_y_rd_b(double *Y, double *y, int* E, double gamma, int m, int n, double* rd, double* b);



void show_E(int* E, int n);


#endif /* D_Comp_hpp */

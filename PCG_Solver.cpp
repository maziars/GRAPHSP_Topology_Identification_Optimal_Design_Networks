//
//  PCG_Solver.cpp
//  PCG
//
//  Created by Maziar Sanjabi on 12/24/15.
//  Copyright © 2015 Maziar Sanjabi. All rights reserved.
//

#include "PCG_Solver.hpp"
//#include "D_Comp.cpp"






void diag_pre_cond(double *r, double *D, double* z, int m){// z = r./D
    for (int i = 0; i<m; i++) {
        z[i] = r[i]/D[i];
    }

}


int diagonal_PCG(double* Y, int *E, double *Ginv, double*x, double* b, int m, int n, double e, double* D, double *xbar, double *ybar){
    
    double* r = (double*) calloc(m, sizeof(double));
    double* z = (double*) malloc(m*sizeof(double));
    double* p = (double*) malloc(m*sizeof(double));
    double* Ap = (double*) malloc(m*sizeof(double));
    double res;
    double res_new;
    double alpha;
    int i = -1;
    

    // Setting initial x to zero
    std::fill_n(x, m, 0.0);
    
#ifdef DEBUG_ON_DEEP
    printf("initial x in PCG is: \n");
    show_vector(x, m);
    printf("\n\n");
    
    printf("initial Ap in PCG is: \n");
    show_vector(Ap, m);
    printf("\n\n");
#endif
    
    //Compute_Ap(Y, E, Ginv, x, m, n, r, xbar, ybar);// r = A*x // No need to do this as x=0
    // r = Ap = A*x = A * 0 = 0
    

    vector_sub(b, r, r, m);// r = b - A*x
    diag_pre_cond(r, D, z, m);// z = M^{-1}*r
    vector_deepcopy(z, p, m);// p = z
    res = innder_prod(z, r, m); // res = <r,z> = r'*M^{-1}*r
    
    //double res_stop = innder_prod(r, r, m); //res_stop is the norm of r to power two (used it to make the code similar to MATLAB version)
    
#ifdef DEBUG_ON_DEEP
    printf("res stop is: %f\n", res_stop);
#endif
    
    
//    if(res_stop>e){
    if(res>e){
        int maxit = ( m < MAXIT_PCG)? m:MAXIT_PCG;
        
        for (i = 0; i< maxit; i++) {
            
#ifdef DEBUG_ON_DEEP
            printf("p at PCG iteration %d is:\n", i);
            show_vector(p, m);
            printf("\n\n");
#endif
            Compute_Ap(Y, E, Ginv, p, m, n, Ap, xbar, ybar);// Ap = A*p
#ifdef DEBUG_ON_DEEP
            printf("Ap at PCG iteration %d is:\n", i);
            show_vector(Ap, m);
            printf("\n\n");
#endif

            alpha = res/innder_prod(p, Ap, m);// alpha = res/p'*A*p
            vector_weight_sum(x, p, x, alpha, m); // x = x + alpha*p;
            vector_weight_sum(r, Ap, r, -alpha, m);// r = r -alpha*Ap;
            diag_pre_cond(r, D, z, m);// z = M^{-1}*r
            res_new = innder_prod(z, r, m); // res_new = <r,z> = r'*M^{-1}*r
            //res_stop = innder_prod(r, r, m); //res_stop is the norm of r to power two (used it to make the code similar to MATLAB version)

#ifdef DEBUG_ON_DEEP
            printf("res stop is: %f\n", res_stop);
#endif
            //if(res_stop<e){
            if (res_new<e) {
                break;
            }
            else{
                vector_weight_sum(z, p, p, res_new/res, m);// p = z + (res_new/res)*p;
                res = res_new;
                
            }
            
            
        }
        
        
    }
    //printf("Used %d PCG iterations\n", i+1);
    
    free(r);
    free(z);
    free(p);
    free(Ap);
    
    return i;
    
    
}


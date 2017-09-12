//
//  QuadCD.cpp
//  PCG
//
//  Created by Maziar Sanjabi on 3/30/16.
//  Copyright © 2016 Maziar Sanjabi. All rights reserved.
//

#include "QuadCD.hpp"
#include "stdlib.h"
#include "D_Comp.hpp"
#include "Proximal_Solver.hpp"

#define defined_abs(x)  (x>=0)? x:-x
#define index_find(i,j,n)    i*n+j


unsigned int find_active(double* x, double* grad, unsigned int* active, int m){
    unsigned int size = 0;
    for (int i=0; i<m; i++) {
        if (grad[i]<0 || x[i] > 0) {
            active[size] = i;
            size++;
        }
    }
    return size;
}




void Compute_grad2(double* Y, double* grad, int* E, double gamma, int m, int n){
    int k;
    int j;
    int i = 0;
    if(m>16){
        int new_m = m - 16*(m/16);
        for(;i<new_m; i++){
            // 0
            j = (i<<1);
            k = E[j++];
            j = E[j];
            grad[i] = 2 + gamma - ( Y[k*n+k] + Y[j*n+j] - Y[k*n+j] - Y[j*n+k] );
            // 1
            i++;
            j = (i<<1);
            k = E[j++];
            j = E[j];
            grad[i] = 2 + gamma - ( Y[k*n+k] + Y[j*n+j] - Y[k*n+j] - Y[j*n+k] );
            // 2
            i++;
            j = (i<<1);
            k = E[j++];
            j = E[j];
            grad[i] = 2 + gamma - ( Y[k*n+k] + Y[j*n+j] - Y[k*n+j] - Y[j*n+k] );
            // 3
            i++;
            j = (i<<1);
            k = E[j++];
            j = E[j];
            grad[i] = 2 + gamma - ( Y[k*n+k] + Y[j*n+j] - Y[k*n+j] - Y[j*n+k] );
            // 4
            i++;
            j = (i<<1);
            k = E[j++];
            j = E[j];
            grad[i] = 2 + gamma - ( Y[k*n+k] + Y[j*n+j] - Y[k*n+j] - Y[j*n+k] );
            // 5
            i++;
            j = (i<<1);
            k = E[j++];
            j = E[j];
            grad[i] = 2 + gamma - ( Y[k*n+k] + Y[j*n+j] - Y[k*n+j] - Y[j*n+k] );
            // 6
            i++;
            j = (i<<1);
            k = E[j++];
            j = E[j];
            grad[i] = 2 + gamma - ( Y[k*n+k] + Y[j*n+j] - Y[k*n+j] - Y[j*n+k] );
            // 7
            i++;
            j = (i<<1);
            k = E[j++];
            j = E[j];
            grad[i] = 2 + gamma - ( Y[k*n+k] + Y[j*n+j] - Y[k*n+j] - Y[j*n+k] );
            // 8
            i++;
            j = (i<<1);
            k = E[j++];
            j = E[j];
            grad[i] = 2 + gamma - ( Y[k*n+k] + Y[j*n+j] - Y[k*n+j] - Y[j*n+k] );
            // 9
            i++;
            j = (i<<1);
            k = E[j++];
            j = E[j];
            grad[i] = 2 + gamma - ( Y[k*n+k] + Y[j*n+j] - Y[k*n+j] - Y[j*n+k] );
            // 10
            i++;
            j = (i<<1);
            k = E[j++];
            j = E[j];
            grad[i] = 2 + gamma - ( Y[k*n+k] + Y[j*n+j] - Y[k*n+j] - Y[j*n+k] );
            // 11
            i++;
            j = (i<<1);
            k = E[j++];
            j = E[j];
            grad[i] = 2 + gamma - ( Y[k*n+k] + Y[j*n+j] - Y[k*n+j] - Y[j*n+k] );
            // 12
            i++;
            j = (i<<1);
            k = E[j++];
            j = E[j];
            grad[i] = 2 + gamma - ( Y[k*n+k] + Y[j*n+j] - Y[k*n+j] - Y[j*n+k] );
            // 13
            i++;
            j = (i<<1);
            k = E[j++];
            j = E[j];
            grad[i] = 2 + gamma - ( Y[k*n+k] + Y[j*n+j] - Y[k*n+j] - Y[j*n+k] );
            // 14
            i++;
            j = (i<<1);
            k = E[j++];
            j = E[j];
            grad[i] = 2 + gamma - ( Y[k*n+k] + Y[j*n+j] - Y[k*n+j] - Y[j*n+k] );
            // 15
            i++;
            j = (i<<1);
            k = E[j++];
            j = E[j];
            grad[i] = 2 + gamma - ( Y[k*n+k] + Y[j*n+j] - Y[k*n+j] - Y[j*n+k] );
            
        }
    }
    for (; i<m; i++) {
        j = (i<<1);
        k = E[j++];
        j = E[j];
        grad[i] = 2 + gamma - ( Y[k*n+k] + Y[j*n+j] - Y[k*n+j] - Y[j*n+k] );
    }
}


void Compute_diagH(double *Y, double *Ginv, int *E, int m, int n, double* diagH){// D = diag(H1 . H2) (hadamard prod of H1 and H2) where H1 = E^T * Y * E and H2 = E^T * Ginv * E
    int j;
    int k;
    for (int i=0; i<m; i++) {
        j = (i<<1);
        k = E[j++];
        j = E[j];
        diagH[i] = (Y[k*n+k] + Y[j*n+j] - Y[k*n+j] - Y[j*n+k])*(Ginv[k*n+k] + Ginv[j*n+j] - Ginv[k*n+j] - Ginv[j*n+k]);
        //        rd[i] = (Y[k*n+k]-1) + (Y[j*n+j]-1) - Y[k*n+j] - Y[j*n+k];
    }
    
}


void multiplyE(double* y, double* x, int* E, int m, int n){// y = E^T*x (y in mx1 and x in nx1 and E is nxm)
    int j;
    int k;
    for (int i = 0; i<m; i++) {
        j = (i<<1);
        k = E[j++];
        j = E[j];
        y[i] = x[k]-x[j];
    }
}


void update_Hdx(double* Y, double* Ginv, int*E, double* Hdx, double mu, int m, int n, int i){
    double* temp1 = (double*) malloc(m*sizeof(double));
    double* temp2 = (double*) malloc(m*sizeof(double));
    double* temp3 = (double*) malloc(m*sizeof(double));
    double* temp4 = (double*) malloc(m*sizeof(double));
    int j;
    int k;
    j = (i<<1);
    k = E[j++];
    j = E[j];
    for (int ind=0; ind<m; ind++) {
        temp1[ind] = Y[n*k+ind] - Y[n*j+ind];
//        printf("\n n: %d\t, j: %d\t, k: %d\t, i: %d\t", n, j, k, i);
//        fflush(stdout);
        temp2[ind] = Ginv[n*k+ind] - Ginv[n*j+ind];
    }
    multiplyE(temp3, temp1, E, m, n);
    multiplyE(temp4, temp2, E, m, n);
    for (int ind=0; ind<m; ind++) {
        Hdx[ind] += mu*temp3[ind]*temp4[ind];
    }
    free(temp1);
    free(temp2);
    free(temp3);
    free(temp4);
}


void update_Hdx2(double* Y, double* Ginv, int*E, double* Hdx, double mu, int m, int n, int i){
    int el_i;
    int k_i;
    el_i = (i<<1);
    k_i = E[el_i++];
    el_i = E[el_i];
    for (int ind=0; ind<m; ind++) {
        int el_ind = (ind<<1);
        int k_ind = E[el_ind++];
        el_ind = E[el_ind];
        Hdx[ind] += mu*(Y[n*k_ind+k_i]-Y[n*k_ind+el_i]-Y[n*el_ind+k_i]+Y[n*el_ind+el_i])*(Ginv[n*k_ind+k_i]-Ginv[n*k_ind+el_i]-Ginv[n*el_ind+k_i]+Ginv[n*el_ind+el_i]);
    }
}




void Quad_Solver(double *Lp,int* E, double* Q,double gamma, double*x, double* y, int m, int n)
{// Outputs vectors are x and y. x is the primal solution.
    // R is assumed to be Identity
    // Only works for resistive networks

    double* Ginv = (double*) malloc(n*n*sizeof(double));
    double* Y = (double*) malloc(n*n*sizeof(double));
    double* grad = (double*) malloc(m*sizeof(double));
    double* G = (double*) malloc(n*n*sizeof(double));
    double* dx = (double*) malloc(m*sizeof(double));
    double* Hdx = (double*) malloc(m*sizeof(double));
    double* diagH = (double*) malloc(m*sizeof(double));
    double* xnew = (double*) malloc(m*sizeof(double));
    unsigned int* active = (unsigned int*) malloc(m*sizeof(unsigned int));
    
    
    
    
    
    
    // Initialization
    
    for (int i = 0; i<m; i++) {
        x[i] = 0;
    }
    
    
    compute_G(Lp, E, G, x, m, n); // computing G
    
    
    
#ifdef DEBUG_ON
    printf("G is: \n");
    show_matrix(G, n);
    printf("\n\n");
#endif
    
    
    G_inv(G, Ginv, n); // Computing the G inverse
    
    // If R is Identity, then qvec(i) = 2 for all i.
    
    Compute_Y(Ginv, Q, Y, n);
    
    
    double obj_old = Compute_Obj2(Ginv, Q, x, gamma, m, n);
    
    
    
    
#ifdef QUAD_DEBUG_ON
    printf("\nInitial G is\n");
    show_matrix(G, n);
    printf("\nInitial Q is\n");
    show_matrix(Q, n);
    printf("\nInitial Ginv is\n");
    show_matrix(Ginv, n);
    printf("\nInitial Y is\n");
    show_matrix(Y, n);
    printf("\nInitial objective is: %f\n", obj_old);
    
#endif // QUAD_DEBUG_ON
    
    //////////////////////////////////////////////////////////////
    
    
    
    for (int outer=0; outer<QUADCD_MAXIT; outer++) {
        
        ///////////////////
//#ifdef QUAD_DEBUG_ON
        printf("\nIteration: %d\n", outer+1);
//#endif // QUAD_DEBUG_ON
        
        

        Compute_grad(Y, grad, E, gamma, m, n);
        
        
        
#ifdef QUAD_DEBUG_ON
        printf("\n x at the beginning of iteration %d:\n", outer+1);
        for (int i = 0; i<m; i++) {
            if (x[i] > 0) {
                printf("(index, value): (%d,%f)\n", i+1, x[i]);
            }
        }
        
        printf("objective is: %f\n", obj_old);
        
        printf("\n Gradient is: \n");
        show_vector(grad, m);
        printf("\n\n");
#endif
        
        
        
        ///////////////////
        
        
        
        for (int i = 0; i<m; i++) {
            Hdx[i] = 0;
            dx[i] = 0;
        }
        
        Compute_diagH(Y, Ginv, E, m, n, diagH);

#ifdef QUAD_DEBUG_ON
        
        printf("\ndiagH is: \n");
        show_vector(diagH, m);
        printf("\n\n");
#endif //QUAD_DEBUG_ON
        
        
        // Finding the active set
        unsigned int active_size = find_active(x, grad, active, m);

#ifdef QUAD_DEBUG_ON
        
        printf("Active set is: \n");
        for(int i=0; i<active_size; i++)
            printf("%d\n", active[i]+1);
        
        
#endif //QUAD_DEBUG_ON
        
        //printf("right before the newton direction loop\n");
        // Finding the Newton direction with coordiante descent
        int ind;
        for (ind = 0; ind<m; ind++) {
            double musum = 0;
            double normdx = 0;
            for (int iin=0; iin<active_size; iin++) {
                unsigned int i = active[iin];
                double b_over_daigH = (Hdx[i]+grad[i])/diagH[i];
                double c = x[i] + dx[i];
                double mu = (c-b_over_daigH>=0) ? -b_over_daigH : -c;
                dx[i]+=mu;
                normdx += defined_abs(dx[i]);
                musum += defined_abs(mu);
                update_Hdx2(Y, Ginv, E, Hdx, mu, m, n, i);
            }
            if (musum<normdx*QUADCD_TOL3) {
#ifdef QUAD_DEBUG_ON
                //printf("\nwe broke here 1 !\n");
#endif // QUAD_DEBUG_ON
                break;
            }
            
        }
        printf("Number of CD rounds: %d\n", ind+1);
        
        // line search and backtracking
        double tau = 1;
        double obj_new = obj_old;
        double grad_dx_inner_prod = innder_prod(grad, dx, m);
        /////////////////////////////////////////////////////////
        if(grad_dx_inner_prod>0)
            break;
        /////////////////////////////////////////////////////////
        int j=0;
        for (j=0; j< MAX_BACK_TRACK; j++){
            
            int flag = 0;
            for (int i=0; i<m; i++) {
                xnew[i] = x[i] + tau*dx[i];
                if (xnew[i]<0) {
                    tau = tau*BETA;
                    flag = 1;
#ifdef QUAD_DEBUG_ON
                    //printf("\nwe broke here 2 !\n");
#endif // QUAD_DEBUG_ON
                    break;
                }
            }
            if (flag ==0) {
                compute_G(Lp, E, G, xnew, m, n); // computing G
                G_inv(G, Ginv, n); // Computing the G inverse
                Compute_Y(Ginv, Q, Y, n);
                obj_new = Compute_Obj2(Ginv, Q, xnew, gamma, m, n);
                if (obj_new<= obj_old + tau*ALPHA*grad_dx_inner_prod) {
#ifdef QUAD_DEBUG_ON
                    //printf("\nwe broke here 3 !\n");
#endif // QUAD_DEBUG_ON
                    break;
                }
                else{
                    tau = BETA*tau;
                }
            }
            
            
        }
        printf("Number of line searches: %d\n", j+1);
        
        if (obj_old-obj_new< EPS_CD) {
#ifdef QUAD_DEBUG_ON
            //printf("\nwe broke here 4 !\n");
#endif // QUAD_DEBUG_ON
            obj_old = obj_new;
            vector_deepcopy(xnew, x, m);
            break;
        }
        
        obj_old = obj_new;
        vector_deepcopy(xnew, x, m);
        
        
    }


    free(Ginv);
    free(Y);
    free(grad);
    free(G);
    free(dx);
    free(Hdx);
    free(diagH);
    free(xnew);
    free(active);


}

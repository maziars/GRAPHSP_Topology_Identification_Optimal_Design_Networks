//
//  D_Comp.cpp
//  
//
//  Created by Maziar Sanjabi on 12/13/15.
//  Copyright © 2015 Maziar Sanjabi. All rights reserved.
//
// This file contains the main dense matrix computation routines that use BLAS routines
// It also contains the implementation of Conjugate gradient and Preconditioned Conjugate gradient method to solve the lienar system of equations

#include "D_Comp.hpp"


// The BLAS routines imported
extern "C" void dgemm_(const char *TRANSA, const char *TRANSB, const int *M, const int *N, const int *K, double *ALPHA, double *A, const int *LDA, double *B, const int *LDB, double *BETA, double *C, const int *LDC);

extern "C" void DGETRF_(const int *M,const int *N,double *A, const int *LDA, int *IPIV, const int *INFO );

extern "C" void DGETRI_(const int *N, double *A, const int *LDA, int *IPIV, double* WORK, const int * LWORK, const int *INFO );

extern "C" double DDOT_(const int *n, double* x, const int *incx, double *y, const int *incy);


// The C wrapper to call thos routines much easier
double innder_prod(double* a, double* b, int n){ // out = <a, b>
    
    
    double prod;
    int inc = 1;
    prod = DDOT_(&n, a, &inc, b, &inc);
    return prod;
}


void vector_prod(double* a, double* b, double*c, int n){// c = a .* b
    for (int i = 0; i<n; i++) {
        c[i] = a[i]*b[i];
    }
}





void matrix_prod(double* A, double* x, double* y, int m, int n){// y = A * x
    
    // This one uses LAPACK implementation and is much faster //todo: change it to dgbmv
    char TRANSA = 'T';
    char TRANSB = 'N';
    //int M = m;
    int N = 1;
    //int K = n;
    double ALPHA = 1.0;
    //int LDA = n;
    //int LDB = n;
    double BETA = 0.0;
    //int LDC = m;
    
    dgemm_(&TRANSA, &TRANSB, &m, &N, &n, &ALPHA, A, &n, x, &n, &BETA, y, &m);

}


void vector_add(double* x, double* y, double* z, int n){// z = x + y
    for (int i=0; i<n; i++) {
        z[i] = x[i] + y[i];
    }
}

void vector_sub(double* x, double* y, double* z, int n){// z = x - y
    for (int i=0; i<n; i++) {
        z[i] = x[i] - y[i];
    }
}

void vector_deepcopy(double* x, double* y, int n){// y = x

    memcpy (y, x, n*sizeof(double));
}

void vector_weight_sum(double* x, double* y, double* z, double alpha, int n){ // z = x + alpha*y
    for (int i=0; i<n; i++) {
        z[i] = x[i] + alpha * y[i];
    }
}


void vector_print(double* x, int n){
    printf("***");
    double p;
    for (int i = 0; i<n; i++) {
        p = x[i];
        printf("%f\n", p);
    }
    printf("***");
}


void pre_cond(double* r, double* D1, double* D2, double* D3, double* z, int n){// z = M^{-1}*r
    
    double* temp = (double*) malloc((n/2)*sizeof(double));
    double* r1 = r;
    double* r2 = &r[n/2];
    double* z1 = z;
    double* z2 = &z[n/2];
    int m = n/2;
    
    vector_prod(r1 , D3, z1, m); // z1 = D3*r1
    vector_prod(r2 , D1, z2, m); // z2 = D1*r2
    vector_prod(z1 , D2, temp, m); // temp = D3*D2*r1
    vector_add(temp, z2, z2, m); // z2 = D3*D2*r1 + D1*r2
    vector_prod(r2 , D2, temp, m); // temp = D2*r2
    vector_prod(temp , D3, temp, m); // temp = D3*D2*r2
    vector_add(temp, z1, z1, m); // z1 = D3*r1 + D3*D2*r2
    vector_prod(temp , D2, temp, m); // temp = D2*D3*D2*r2
    vector_add(temp, z2, z2, m); // z2 = D3*D2*r1 + (D1 + D2*D3*D2)*r2
    
    free(temp);
}



// Conjugate gradient method
void CG(double* A, double*x, double* b, int m, double e){
    
    double* r = (double*) malloc(m*sizeof(double));
    double* p = (double*) malloc(m*sizeof(double));
    double* Ap = (double*) malloc(m*sizeof(double));
    double res;
    double res_new;
    double alpha;
    
    

    matrix_prod(A, x, r, m, m);// r = A*x
    vector_sub(b, r, r, m);// r = b - A*x
    vector_deepcopy(r, p, m);// p = r
    res = innder_prod(r, r, m); // res = || r ||^2
    
    
    
    if(res>e){
        for (int i = 0; i< m; i++) {
            matrix_prod(A, p, Ap, m, m);// Ap = A*p
            alpha = res/innder_prod(p, Ap, m);// alpha = res/p'*A*p
            vector_weight_sum(x, p, x, alpha, m); // x = x + alpha*p;
            vector_weight_sum(r, Ap, r, -alpha, m);// r = r -alpha*Ap;
            res_new = innder_prod(r, r, m); // res_new = || r ||^2
            if(res_new<e)
                break;
            else{
                vector_weight_sum(r, p, p, res_new/res, m);// p = r + (res_new/res)*p;
                res = res_new;
            }
        }
        
        
    }
    
    free(r);
    free(p);
    free(Ap);
    
    
}




void PCG(double* A, double*x, double* b, int m, double e, double* D1, double* D2, double* D3){
    
    double* r = (double*) malloc(m*sizeof(double));
    double* z = (double*) malloc(m*sizeof(double));
    double* p = (double*) malloc(m*sizeof(double));
    double* Ap = (double*) malloc(m*sizeof(double));
    double res;
    double res_new;
    double alpha;
    
    
    
    matrix_prod(A, x, r, m, m);// r = A*x
    vector_sub(b, r, r, m);// r = b - A*x
    pre_cond(r, D1, D2, D3, z, m); // z = M^{-1}*r
    vector_deepcopy(z, p, m);// p = r
    res = innder_prod(z, r, m); // res = <r,z> = r'*M^{-1}*r
    
    
    
    if(res>e){
        for (int i = 0; i< m; i++) {
            printf("i is %d\n", i);
            matrix_prod(A, p, Ap, m, m);// Ap = A*p
            alpha = res/innder_prod(p, Ap, m);// alpha = res/p'*A*p
            vector_weight_sum(x, p, x, alpha, m); // x = x + alpha*p;
            vector_weight_sum(r, Ap, r, -alpha, m);// r = r -alpha*Ap;
            pre_cond(r, D1, D2, D3, z, m); // z = M^{-1}*r
            res_new = innder_prod(z, r, m); // res_new = <r,z> = r'*M^{-1}*r
            if(res_new<e)
                break;
            else{
                vector_weight_sum(z, p, p, res_new/res, m);// p = z + (res_new/res)*p;
                res = res_new;
                
            }
            
        }
        
        
    }
    
    
    free(r);
    free(z);
    free(p);
    free(Ap);
    
    
}



void cholesky(double *A, int n) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < (i+1); j++) {
            double s = 0;
            for (int k = 0; k < j; k++)
                s += A[i * n + k] * A[j * n + k];
            A[i * n + j] = (i == j) ?
            sqrt(A[i * n + i] - s) :
            (1.0 / A[j * n + j] * (A[i * n + j] - s));
            A[j*n + i] = A[i*n + j];
        }
}



void show_matrix(double *A, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            printf("%2.5f ", A[i * n + j]);
        printf("\n");
    }
}



void show_vector(double *v, int n) {
    for (int i = 0; i < n; i++) {
        printf("%2.5f ", v[i]);
        printf("\n");
    }
}




void compute_G(double *Lp, int* E, double* G, double *x, int m, int n){
    memcpy(G, Lp, n*n*sizeof(double));
    int itemp;
    int jtemp;
    for (int i = 0; i<m; i++) {
        jtemp = i<<1;
        itemp = E[jtemp++];
        jtemp = E[jtemp];
        
        G[itemp*n+itemp] += x[i];
        G[itemp*n+jtemp] -= x[i];
        G[jtemp*n+itemp] -= x[i];
        G[jtemp*n+jtemp] += x[i];
        
    }

}



void G_inv(double *G, double *Ginv, int n){

    
    char UPLO = 'L';
    int INFO;

    int LWORK = n*n;
    double* WORK =(double*) malloc(LWORK*sizeof(double));
    int* IPIV =(int*) malloc(n*sizeof(int));
    DGETRF_(&n,&n,G, &n, IPIV, &INFO );
    DGETRI_(&n, G, &n, IPIV, WORK, &LWORK, &INFO );
    free(WORK);
    free(IPIV);
    int N = n*n;
    
    
    vector_deepcopy(G, Ginv, N); // todo: remove the G and go directly for Ginv

    
}







void Compute_Y(double *Ginv, double *Qp, double *Y, int n){

    double* temp = (double*) malloc(n*n*sizeof(double));
    char TRANSA = 'N';
    char TRANSB = 'T';
    double ALPHA = 1.0;
    double BETA = 0.0;
    
    dgemm_(&TRANSA, &TRANSA, &n, &n, &n, &ALPHA, Ginv, &n, Qp, &n, &BETA, temp, &n);
    dgemm_(&TRANSA, &TRANSB, &n, &n, &n, &ALPHA, temp, &n, Ginv, &n, &BETA, Y, &n);
    
    free(temp);
    
    
}


void Compute_rd(double *Y, double *y, int* E, double gamma, int m, int n, double* rd){
    int k;
    int j;
    for (int i=0; i<m; i++) {
        j = (i<<1);
        k = E[j++];
        j = E[j];
        rd[i] = gamma - ( (Y[k*n+k]-1) + (Y[j*n+j]-1) - Y[k*n+j] - Y[j*n+k] + y[i]);
    }
    
}



void Compute_rd_b(double *Y, double *y, int* E, double gamma, int m, int n, double* rd, double* b){
    int k;
    int j;
    for (int i=0; i<m; i++) {
        j = (i<<1);
        k = E[j++];
        j = E[j];
        rd[i] = gamma - ( (Y[k*n+k]-1) + (Y[j*n+j]-1) - Y[k*n+j] - Y[j*n+k] + y[i]);
        b[i] = - (rd[i]+y[i]);
    }
    
}

void Compute_b(double *Y, double *y, int* E, double gamma, int m, int n, double* b){
    int k;
    int j;
    for (int i=0; i<m; i++) {
        j = (i<<1);
        k = E[j++];
        j = E[j];
        b[i] = - gamma + ( (Y[k*n+k]-1) + (Y[j*n+j]-1) - Y[k*n+j] - Y[j*n+k] );
    }
    
}






void Compute_Ap(double *Y, int *E, double *Ginv, double *p, int m, int n, double *Ap, double *x, double *y){ // Ap = A*p
    double* temp = (double*) calloc(n*n, sizeof(double));
    double* temp2 = (double*) calloc(n*n, sizeof(double));
//    double* b = (double*) malloc(n*sizeof(double));
    double* M = (double*) malloc(n*n*sizeof(double));
    int k;
    int j;
    for (int i = 0; i<m; i++) {
        j = i<<1;
        k = E[j++];
        j = E[j];
        
        temp[j*n+j] += p[i];
        temp[j*n+k] -= p[i];
        temp[k*n+j] -= p[i];
        temp[k*n+k] += p[i];
    }
    

    char TRANS = 'N';
    char TRANSB = 'T';
    double ALPHA = 1.0;
    double BETA = 0.0;
    
    dgemm_(&TRANSB, &TRANSB, &n, &n, &n, &ALPHA, temp, &n, Ginv, &n, &BETA, temp2, &n);
    dgemm_(&TRANSB, &TRANSB, &n, &n, &n, &ALPHA, Y, &n, temp2, &n, &BETA, M, &n);
    
    free(temp2);
    free(temp);
//    free(b);
    
    

#ifdef DEBUG_ON_DEEP
    printf("M is: \n");
    show_matrix(M, n);
    printf("\n\n");
    
    
    printf("y is: \n");
    show_vector(y, m);
    printf("\n\n");
    
    printf("x is: \n");
    show_vector(x, m);
    printf("\n\n");
#endif
    
    
    for (int i=0; i<m; i++) {
        j = (i<<1);
        k = E[j++];
        j = E[j];
        Ap[i] = 2*(M[k*n+k] + M[j*n+j] - M[k*n+j] - M[j*n+k]) + p[i]*y[i]/x[i];
        //        rd[i] = (Y[k*n+k]-1) + (Y[j*n+j]-1) - Y[k*n+j] - Y[j*n+k];
    }
    free(M);
    
    
#ifdef DEBUG_ON_DEEP
    printf("Ap is: \n");
    show_vector(Ap, m);
    printf("\n\n");
#endif
    
    
    
    
}



void Compute_diagA(double *Y, double *Ginv, double *x, double *y, int *E, int m, int n, double* D){// D = 2* diag(H1 . H2) + y./x (hadamard prod of H1 and H2) where H1 = E^T * Y * E and H2 = E^T * Ginv * E
    int j;
    int k;
    for (int i=0; i<m; i++) {
        j = (i<<1);
        k = E[j++];
        j = E[j];
        D[i] = 2*(Y[k*n+k] + Y[j*n+j] - Y[k*n+j] - Y[j*n+k])*(Ginv[k*n+k] + Ginv[j*n+j] - Ginv[k*n+j] - Ginv[j*n+k])+ y[i]/x[i];
        //        rd[i] = (Y[k*n+k]-1) + (Y[j*n+j]-1) - Y[k*n+j] - Y[j*n+k];
    }
    
}







void Initialize_y_rd_b(double *Y, double *y, int* E, double gamma, int m, int n, double* rd, double* b){
    int k;
    int j;
    double temp;
//    double small = 0.1;
    for (int i=0; i<m; i++) {
        j = (i<<1);
        k = E[j++];
        j = E[j];
        
        temp = gamma - ( (Y[k*n+k]-1) + (Y[j*n+j]-1) - Y[k*n+j] - Y[j*n+k]);
#ifdef INIT1
        y[i] = 1.01*fabs(temp);
        rd[i] = temp - y[i];
#endif
#ifndef INIT1
        y[i] = (temp>0)? temp:0;
        rd[i] = (temp>0)? 0: temp;
#endif
        b[i] = - (rd[i]+y[i]);
    }
    
}


void show_E(int* E, int n){
    for (int i = 0; i < n; i++) {
        int j = i<<1;
        printf("%d, %d", E[j], E[j+1]);
        printf("\n");
    }
    
}





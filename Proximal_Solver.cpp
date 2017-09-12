//
//  Proximal_Solver.cpp
//  PCG
//
//  Created by Maziar Sanjabi on 1/3/16.
//  Copyright © 2016 Maziar Sanjabi. All rights reserved.
//
// This includes two implementations of the proximal solver. The first one uses a customized sparse representation and the second one uses simple non-sparse implementation

#include "Proximal_Solver.hpp"
//#include "D_Comp.cpp"

void compute_G_Sparse(double *Lp, int* E, double* G, SparseVec *xs, int m, int n){
    
    memcpy(G, Lp, n*n*sizeof(double));
    int itemp;
    int jtemp;
    int i;
    double x;
    NZiterator it = xs->NZ.begin();
    while(it != xs->NZ.end())
    {
        i = it->first;
        x = it->second;
        jtemp = i<<1;
        itemp = E[jtemp++];
        jtemp = E[jtemp];
        
        
        G[itemp*n+itemp] += x;
        G[itemp*n+jtemp] -= x;
        G[jtemp*n+itemp] -= x;
        G[jtemp*n+jtemp] += x;
        
        it++;
    }
    
    
}


double Compute_Obj(double* Ginv,double* Q, SparseVec* xs, double gamma, int n){
    double obj = (2+gamma)*(xs->sum())+innder_prod(Ginv, Q, n*n);
    return obj;
}



double Compute_Obj2(double* Ginv,double* Q, double* xs, double gamma, int m , int n){
    double obj = 0;
    for (int i = 0; i<m; i++) {
        obj = obj + xs[i];
    }
    obj = (2+gamma)*obj+innder_prod(Ginv, Q, n*n);
    return obj;
}

void Compute_grad(double* Y, double* grad, int* E, double gamma, int m, int n){
    int k;
    int j;
    for (int i=0; i<m; i++) {
        j = (i<<1);
        k = E[j++];
        j = E[j];
        grad[i] = 2 + gamma - ( Y[k*n+k] + Y[j*n+j] - Y[k*n+j] - Y[j*n+k] );
    }
}



void Compute_new(SparseVec* xs, double* prox_step, double* grad, double ssize, int m, int n){// prox_step = x_new-x_old
    double x_old;
    double temp;
    SparseVec xd(m);
    for (int i = 0; i<m; i++) {
        xs->getValue(i, &x_old);
        temp = x_old - ssize*grad[i];
        if (temp>0) {
            //xs->update(i, temp);
            xd.insert(i, temp);
            prox_step[i] = temp-x_old;
        }
        else{
            //xs->erase(i);
            prox_step[i] = -x_old;
        }
    }
    *xs = xd;
}

void Compute_new2(double* xs, double* prox_step, double* grad, double ssize, int m, int n, SparseVec* xpss){// prox_step = x_new-x_old
    double x_old;
    double temp;
    for (int i = 0; i<m; i++) {
        x_old = xs[i];
        temp = x_old - ssize*grad[i];
        if (temp>0) {
            //xs->update(i, temp);
            xs[i] = temp;
            xpss->insert(i, temp);
            prox_step[i] = temp-x_old;
        }
        else{
            //xs->erase(i);
            prox_step[i] = -x_old;
            xs[i] = 0.0;
        }
    }
}




   
void Proximal_Solver(double *Lp,int* E, double* Q,double gamma, SparseVec *xs, double* y, int m, int n)
{// Outputs vectors are x and y. x is the primal solution.
    // R is assumed to be Identity
    // Only works for resistive networks
    
    double* Ginv = (double*) malloc(n*n*sizeof(double));
    double* Y = (double*) malloc(n*n*sizeof(double));
    double* grad = (double*) malloc(m*sizeof(double));
    double* grad_old = (double*) malloc(m*sizeof(double));
    double* prox_step = (double*) malloc(m*sizeof(double));
    double* G = (double*) malloc(n*n*sizeof(double));
    
    
    double ssize = 1.0;
    double obj;
    double prox_size;
    double norm_prox_size;
    
    // Initialzations

    compute_G_Sparse(Lp, E, G, xs, m, n); // computing G
    
    
#ifdef DEBUG_ON
    printf("G is: \n");
    show_matrix(G, n);
    printf("\n\n");
#endif
    
    
    G_inv(G, Ginv, n); // Computing the G inverse
    
    // If R is Identity, then qvec(i) = 2 for all i.
    
    Compute_Y(Ginv, Q, Y, n);
    
    obj = Compute_Obj(Ginv, Q, xs, gamma, n);
    Compute_grad(Y, grad, E, gamma, m, n);
    
#ifdef DEBUG_ON
    printf("Initial objective is: %f\n", obj);
    
    printf("\nInitial Gradient is: \n");
    show_vector(grad, m);
    printf("\n\n");
#endif
    
    
    for (int iter = 0; iter< MAXIT; iter++) {
    //for (int iter = 0; iter< 2; iter++) {
    
        

        
        
        
//        for (int i = 0; i<MAX_BT; i++) { // Back-Tracking
//            <#statements#>
//        }
        
        
        
        
        
        // Update
        Compute_new(xs, prox_step, grad, ssize, m, n);// prox_step = x_new-x_old
        
        
        // BB step-size
        compute_G_Sparse(Lp, E, G, xs, m, n); // computing G
        
#ifdef DEBUG_ON
        printf("G is: \n");
        show_matrix(G, n);
        printf("\n\n");
#endif

        G_inv(G, Ginv, n); // Computing the G inverse
        vector_deepcopy(grad, grad_old, m); // grad_old = grad
        
        Compute_Y(Ginv, Q, Y, n);
        Compute_grad(Y, grad, E, gamma, m, n); // grad = new grad
        obj = Compute_Obj(Ginv, Q, xs, gamma, n);
        
#ifdef DEBUG_ON
        printf("objective is: %f\n", obj);
        
        printf("\nGradient is: \n");
        show_vector(grad, m);
        printf("\n\n");
#endif
        
        vector_sub(grad, grad_old, grad_old, m); // grad_old = grad - grad_old
        
        prox_size = innder_prod(prox_step, prox_step, m);
        norm_prox_size = prox_size/(ssize*ssize);
        
        printf("iteration: %d, objective: %f, norm_prox_size: %f\n", iter, obj, norm_prox_size);
        
        if (norm_prox_size<(EPS_ISTA*EPS_ISTA)) {
            break;
        }
        
        ssize = prox_size/innder_prod(prox_step, grad_old, m);
        if (ssize<0) {
            ssize = 0.0001;
        }
        else if (ssize == NAN)
            ssize = 1;
        
        
        
        
        
    }

    free(Y);
    free(grad);
    free(grad_old);
    free(prox_step);
    free(G);
   








}








//////////////////////////////////////////////////////////////////////////////////////////////////////////

void Proximal_Solver2(double *Lp,int* E, double* Q,double gamma, double *xs, double* y, int m, int n)
{// Outputs vectors are x and y. x is the primal solution.
    // R is assumed to be Identity
    // Only works for resistive networks
    
    double* Ginv = (double*) malloc(n*n*sizeof(double));
    double* Y = (double*) malloc(n*n*sizeof(double));
    double* grad = (double*) malloc(m*sizeof(double));
    double* grad_old = (double*) malloc(m*sizeof(double));
    double* prox_step = (double*) malloc(m*sizeof(double));
    double* G = (double*) malloc(n*n*sizeof(double));
    
    
    double ssize = 1.0;
    double obj;
    double prox_size;
    double norm_prox_size;
    
    
    
    // Initialzations
    for (int i = 0; i<m; i++) {
        xs[i] = 0;
    }
    compute_G(Lp, E, G, xs, m, n); // computing G
    
    
    
#ifdef DEBUG_ON
    printf("G is: \n");
    show_matrix(G, n);
    printf("\n\n");
#endif
    
    
    G_inv(G, Ginv, n); // Computing the G inverse
    
    // If R is Identity, then qvec(i) = 2 for all i.
    
    Compute_Y(Ginv, Q, Y, n);
    
    obj = Compute_Obj2(Ginv, Q, xs, gamma, m, n);
    Compute_grad(Y, grad, E, gamma, m, n);
    
#ifdef DEBUG_ON
    printf("Initial objective is: %f\n", obj);
    
    printf("\nInitial Gradient is: \n");
    show_vector(grad, m);
    printf("\n\n");
#endif
    
    
    for (int iter = 0; iter< MAXIT; iter++) {
        clock_t  start_iter;
        start_iter = std::clock();
        
        //for (int iter = 0; iter< 2; iter++) {
        
        
        
        
        
        
        //        for (int i = 0; i<MAX_BT; i++) { // Back-Tracking
        //            <#statements#>
        //        }
        
        
        
        
        SparseVec xss(m);
        SparseVec* xpss = &xss;
        // Update
        Compute_new2(xs, prox_step, grad, ssize, m, n, xpss);// prox_step = x_new-x_old
        
        
        // BB step-size
        compute_G_Sparse(Lp, E, G, xpss, m, n); // computing G
        
        
#ifdef DEBUG_ON
        printf("G is: \n");
        show_matrix(G, n);
        printf("\n\n");
#endif
       
#ifdef TIME_PROFILE
        clock_t start_Ginv;
        start_Ginv = std::clock();
#endif
        
        G_inv(G, Ginv, n); // Computing the G inverse
        
#ifdef TIME_PROFILE
        printf("\nGinv time is: %f ms\n", (std::clock() - start_Ginv) / (double)(CLOCKS_PER_SEC / 1000));
#endif
        
        vector_deepcopy(grad, grad_old, m); // grad_old = grad
        
#ifdef TIME_PROFILE
        clock_t start_Y;
        start_Y = std::clock();
#endif
        
        Compute_Y(Ginv, Q, Y, n);
        
#ifdef TIME_PROFILE
        printf("\nY time is: %f ms\n", (std::clock() - start_Y) / (double)(CLOCKS_PER_SEC / 1000));
#endif
        
        Compute_grad(Y, grad, E, gamma, m, n); // grad = new grad
        obj = Compute_Obj2(Ginv, Q, xs, gamma, m, n);
//        printf("\nGinv, Y, grad, obj time is: %f ms\n", (std::clock() - start_Ginv) / (double)(CLOCKS_PER_SEC / 1000));
        
#ifdef DEBUG_ON
        printf("objective is: %f\n", obj);
        
        printf("\nGradient is: \n");
        show_vector(grad, m);
        printf("\n\n");
#endif
        
        vector_sub(grad, grad_old, grad_old, m); // grad_old = grad - grad_old
        
        prox_size = innder_prod(prox_step, prox_step, m);
        norm_prox_size = prox_size/(ssize*ssize);
        
        printf("iteration: %d, objective: %f, norm_prox_size: %f\n", iter, obj, norm_prox_size);
        
        if (norm_prox_size<(EPS_ISTA*EPS_ISTA)) {
            break;
        }
        
        ssize = prox_size/innder_prod(prox_step, grad_old, m);
        if (ssize<0) {
            ssize = 0.0001;
        }
        else if (ssize == NAN)
            ssize = 1;
//#ifdef TIME_PROFILE
        printf("\nIteration time is: %f ms\n", (std::clock() - start_iter) / (double)(CLOCKS_PER_SEC / 1000));
//#endif
        
        
    }
    
    free(Y);
    free(grad);
    free(grad_old);
    free(prox_step);
    free(G);
    
    
    
    
    
    
    
    
    
}





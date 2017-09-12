//
//  IP_Solver.cpp
//  PCG
//
//  Created by Maziar Sanjabi on 12/24/15.
//  Copyright © 2015 Maziar Sanjabi. All rights reserved.
//

#include "IP_Solver.hpp"
//#include "D_Comp.cpp"







void IP_Solver(double *Lp,int* E, double* Q,double gamma, double *x, double* y, int m, int n){// Outputs vectors are x and y. x is the primal solution.
// R is assumed to be Identity
// Only works for resistive networks
    
    double* xt = (double*) malloc(m*sizeof(double));
    double* xxt = (double*) malloc(m*sizeof(double));
    double* yt = (double*) malloc(m*sizeof(double));
    double* D = (double*) malloc(m*sizeof(double));
    double* b = (double*) malloc(m*sizeof(double));
    double* rd = (double*) malloc(m*sizeof(double));
    double* Ginv = (double*) malloc(n*n*sizeof(double));
    double* Y = (double*) malloc(n*n*sizeof(double));
    
    double eta;
    double rd_size;
    double e;
    double mu;
    double sigma;
    double alpha_x;
    double alpha_y;
    int I1 = -1;
    int I2 = -1;
    
    // Initilize x
    for (int i = 0; i<m; i++) {
        x[i] = 0.05*500.0/n;//0.05;// + (rand()%10)/10.0;
    }
#ifdef DEBUG_ON
    printf("Initial x is: \n");
    show_vector(x, m);
    printf("\n\n");
#endif
    
    
    
    for (int iter = 0; iter< MAXIT; iter++) {
        // Doing the iterations
        
        double* G = (double*) malloc(n*n*sizeof(double));
        compute_G(Lp, E, G, x, m, n); // computing G
        
#ifdef DEBUG_ON
        printf("G is: \n");
        show_matrix(G, n);
        printf("\n\n");
#endif
        
        
        G_inv(G, Ginv, n); // Computing the G inverse
        free(G); // Don't need G anymore
        
        Compute_Y(Ginv, Q, Y, n);
        
        if (iter==0){
            Initialize_y_rd_b(Y, y, E, gamma, m, n, rd, b);
        
#ifdef DEBUG_ON
            printf("Initial Ginv is: \n");
            show_matrix(Ginv, n);
            printf("\n\n");
            
            printf("Initial Y is: \n");
            show_matrix(Y, n);
            printf("\n\n");
            
            
            printf("Initial y is: \n");
            show_vector(y, m);
            printf("\n\n");
            
            printf("Initial rd is: \n");
            show_vector(rd, m);
            printf("\n\n");
#endif
        }
        else
            Compute_rd_b(Y, y, E, gamma, m, n, rd, b); // we have rd and b now
        
        
#ifdef DEBUG_ON
        printf("rd is: \n");
        show_vector(rd, m);
        printf("\n\n");
        
        printf("b is: \n");
        show_vector(b, m);
        printf("\n\n");
#endif


        
        // Check the stoping criteria
        rd_size = innder_prod(rd, rd, m);
        rd_size = sqrt(rd_size);
        eta = innder_prod(x, y, m);

#ifdef DEBUG_ON
        printf("eta is: %f\n", eta);
#endif
        
        if (iter == 0) {
            printf("it# \t Res_Size\t Opt_Gap\t I1 \t I2\n");
            printf("%d,\t\t %f,\t %f,\t %d,\t\t %d\n", 0, rd_size, eta, I1+1, I2+1);
            //printf("Initial Optimality Condition are, Residual size: %2.5f, Optimality Gap: %2.5f\n", rd_size, eta);
        }
        else {
            //printf("Optimality Condition at iteration %d, Residual size: %2.5f, Optimality Gap: %2.5f\n", iter, rd_size, eta);
            printf("%d,\t\t %f,\t %f,\t %d,\t\t %d\n", iter, rd_size, eta, I1+1, I2+1);
        }
        
        if ((rd_size< EPS_RES) && (eta< EPS_DG)) {
            break;
        }
        
        // Solve the linear equations
        
        e = (DELTA_PCG_AFF*eta);
        e = e * e/(innder_prod(b, b, m));
        
        if (e>0.01) {
            e = 0.01;
        }
#ifdef DEBUG_ON
        printf("e is: %f \n", e);
        
#endif
        
        // Calculating the diagonal entries of linear approximation D= diagA
        Compute_diagA(Y, Ginv, x, y, E, m, n, D);

#ifdef DEBUG_ON
        printf("D is: \n");
        show_vector(D, m);
        printf("\n\n");
#endif
        
        // Finding  xt
        
        I1 = diagonal_PCG(Y, E, Ginv, xt, b, m, n, e, D, x, y);
        
#ifdef DEBUG_ON
        printf("xt is: \n");
        show_vector(xt, m);
        printf("\n\n");
#endif
        
        // Finding yt
        for (int i = 0; i<m; i++) {
            yt[i] = -y[i]*(xt[i]/x[i]+1);
        }
        
#ifdef DEBUG_ON
        printf("yt is: \n");
        show_vector(yt, m);
        printf("\n\n");
#endif
        
        
        // Finding the step size for the linear approximation
        alpha_x = 1.0;
        alpha_y = 1.0;
        
        for (int i=0; i<m; i++) {
            if (xt[i]<0) {
                double temp = fabs(x[i]/xt[i]);
                alpha_x = (temp < alpha_x)? temp : alpha_x;
            }
            if (yt[i]<0) {
                double temp = fabs(y[i]/yt[i]);
                alpha_y = (temp < alpha_y)? temp : alpha_y;
            }
        }

#ifdef DEBUG_ON
        printf("alpha_x and alpha_y are: %f, %f\n", alpha_x, alpha_y);
#endif
#ifdef DEBUG_ON_DEEP
        printf("x before computing sigma is: \n");
        show_vector(x, m);
        printf("\n\n");
        
        printf("y before computing sigma is: \n");
        show_vector(y, m);
        printf("\n\n");
#endif
        
        
        // Computing mu and sigma
        mu = eta/m;
        sigma = 0;
        for (int i = 0; i<m; i++) {
            sigma = sigma + (x[i] + alpha_x*xt[i] )*(y[i] + alpha_y*yt[i]);
        }
        sigma = sigma/ eta;
        sigma = sigma * sigma * sigma;
        
#ifdef DEBUG_ON
        printf("mu and sigma are: %f, %f\n", mu, sigma);
#endif
        
        // Update the path
        double mu_sigma = mu*sigma;
        for (int i = 0; i<m; i++) {
            b[i] = (mu_sigma - yt[i]*xt[i])/x[i];
        }
        
        
#ifdef DEBUG_ON
        printf("Updated b is: \n");
        show_vector(b, m);
        printf("\n\n");
#endif
        
        // Solve second linear equality
        
        //finding xx_t
        e = (DELTA_PCG*eta);
        e = e * e/(innder_prod(b, b, m));
        
        if (e>0.01) {
            e = 0.01;
        }
        
        I2 = diagonal_PCG(Y, E, Ginv, xxt, b, m, n, e, D, x, y);
        
        
        
        // The solution for the second linear equation:
#ifdef DEBUG_ON
        printf("xxt is: \n");
        show_vector(xxt, m);
        printf("\n\n");
#endif
        
        // Update the xt based on new xxt
        for (int i = 0; i<m; i++) {
            xt[i] = xt[i] + xxt[i];
        }
        
        // Finding yt
        for (int i = 0; i<m; i++) {
            yt[i] = b[i]-y[i]*(xt[i]/x[i]+1);
        }
        
#ifdef DEBUG_ON
        printf("xt is: \n");
        show_vector(xt, m);
        printf("\n\n");
        
        printf("yt is: \n");
        show_vector(yt, m);
        printf("\n\n");
#endif
        
        // Finding the step size for the linear approximation
        alpha_x = 1.0;
        alpha_y = 1.0;
        
        for (int i=0; i<m; i++) {
            if (xt[i]<0) {
                double temp = fabs(x[i]/xt[i]);
                alpha_x = (temp < alpha_x)? temp : alpha_x;
            }
            if (yt[i]<0) {
                double temp = fabs(y[i]/yt[i]);
                alpha_y = (temp < alpha_y)? temp : alpha_y;
            }
        }
        alpha_x = alpha_x*.99;
        alpha_y = alpha_y*.99;

#ifdef DEBUG_ON
        printf("alpha_x and alpha_y are: %f, %f\n", alpha_x, alpha_y);
#endif
        
        // Update x and y
        
        for (int i=0; i<m; i++) {
            x[i] = x[i] + alpha_x * xt[i];
            y[i] = y[i] + alpha_y * yt[i];
        }
        
#ifdef DEBUG_ON
        printf("x is: \n");
        show_vector(x, m);
        printf("\n\n");
        
        printf("y is: \n");
        show_vector(y, m);
        printf("\n\n");
#endif

        
        
    }
    

    
#ifdef DEBUG_ON
    printf("x is: \n");
    show_vector(x, m);
    printf("\n\n");
    
    printf("y is: \n");
    show_vector(y, m);
    printf("\n\n");
#endif
    
    
    
    free(xt);
    free(yt);
    free(D);
    free(b);
    free(rd);
    free(Ginv);
    free(Y);
    
    
    
}
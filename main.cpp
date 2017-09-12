//
//  main.cpp
//  PCG
//
//  Created by Maziar Sanjabi on 12/12/15.
//  Copyright © 2015 Maziar Sanjabi. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <ctime>
#include "D_Comp.hpp"
#include "PCG_Solver.hpp"
#include "IP_Solver.hpp"
#include "Constants.h"
#include "IO_MANAGER.hpp"
#include "Proximal_Solver.hpp"
#include "QuadCD.hpp"








int main(int argc, const char * argv[]) {

    int n; // Number of nodes
    int m; // Number of possible controller edges
    double gamma; // Sparsity promoting constant
    int readl;
    int MN[2];
    time_t now;
    time_t after;
    
    
//    FILE *file1 = fopen("/Users/Maziar/Downloads/graphsp/resistive/MN.txt", "r");
    FILE *file1 = fopen("MN.txt", "r");
    if (file1 == NULL){
        perror("ERROR:");
    }
    
    readl = fread_int(file1, MN);
    if (readl != 2) {
        printf("I/O error in MN.txt\n");
    }
    m = MN[0];
    n = MN[1];
    
    
    
    
    
    int* E = (int*) malloc(2*m*sizeof(int));
    double* Lp = (double*) malloc(n*n*sizeof(double));
    double* Qp = (double*) malloc(n*n*sizeof(double));
    double* x = (double*) malloc(m*sizeof(double));
    double* y = (double*) malloc(m*sizeof(double));
    
    file1 = fopen("gamma.txt", "r");
    if (file1 == NULL){
        perror("ERROR:");
    }
    readl = fread_float(file1, &gamma);
    if (readl != 1)
        printf("I/O error in gamma.txt\n");
    
    
    

    file1 = fopen("Lp.txt", "r");
    if (file1 == NULL){
        perror("ERROR:");
    }
    readl = fread_float(file1, Lp);
    if (readl != n*n)
        printf("I/O error in Lp.txt\n");
    
#ifdef DEBUG_ON
    printf("Lp that is read from the file is:\n");
    show_matrix(Lp, n);
    
    printf("\n\n");
#endif
    
    

    file1 = fopen("Qp.txt", "r");
    if (file1 == NULL){
        perror("ERROR:");
    }
    readl = fread_float(file1, Qp);
    if (readl != n*n)
        printf("I/O error in Qp.txt\n");
    

#ifdef DEBUG_ON
    printf("Qp that is read from the file is:\n");
    show_matrix(Qp, n);
    
    printf("\n\n");
#endif
    
    

    file1 = fopen("E.txt", "r");
    if (file1 == NULL){
        perror("ERROR:");
    }
    readl = fread_int(file1, E);
    if (readl != 2*m)
        printf("I/O error in E.txt\n");
    
    
#ifdef DEBUG_ON
    printf("E that is read from the file is:\n");
    show_E(E, m);
#endif
    
    

    
    
    
    
//    ////////////////////////////////////////////////////////////////////////////////////////////////
    
#ifdef IPSOLVER
    
    printf("\n------IP_Solver--------\n");
    now = time(0);
    
    IP_Solver(Lp,E, Qp, gamma, x, y, m, n);

    after = time(0);
    
    std::cout <<"\nElapsed time is: "<< after - now << " s"<< std::endl;
    
    
    printf("\nPrinting the sparse solution x:\n ");
    double small = 0.00001;
    for (int i = 0; i<m; i++) {
        if (x[i] > small) {
            printf("(index, value): (%d,%f)\n", i+1, x[i]);
        }
    }
#endif //IPSOLVER
//
//    /////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    
///////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef PROX1
    
    printf("\n------Prox with Sparse Computation--------\n");
    now = time(0);
    
    SparseVec xs(m);
    SparseVec* xps = &xs;
    Proximal_Solver(Lp,E, Qp, gamma, xps, y, m, n);
    
    
    after = time(0);
    
    std::cout <<"\nElapsed time is: "<< after - now << " s\n"<< std::endl;
    
    xs.print();
#endif // PROX1
///////////////////////////////////////////////////////////////////////////////////////////////////////
    
#ifdef PROX2
    printf("\n------Prox without Sparse Computation--------\n");
//    time_t now = time(0);
    now = time(0);
    
//    SparseVec xs(m);
//    SparseVec* xps = &xs;
    Proximal_Solver2(Lp,E, Qp, gamma, x, y, m, n);
    
    
//    time_t after = time(0);
    after = time(0);
    
    std::cout <<"\nElapsed time is: "<< after - now << " s\n"<< std::endl;
    
    printf("\nPrinting the sparse solution x:\n ");
    
    for (int i = 0; i<m; i++) {
        if (x[i] > 0) {
            printf("(index, value): (%d,%f)\n", i+1, x[i]);
        }
    }

    
    
#endif // PROX2
    
    
//////////////////////////////////////////////////////////
#ifdef QUAD_CD
    printf("\n------Quad CD--------\n");
    now = time(0);
    Quad_Solver(Lp, E, Qp, gamma, x, y, m, n);
    after = time(0);
    
    std::cout <<"\nElapsed time is: "<< after - now << " s\n"<< std::endl;
    for (int i = 0; i<m; i++) {
        if (x[i] > 0) {
            printf("(index, value): (%d,%f)\n", i+1, x[i]);
        }
    }
#endif // QUAD_CD
///////////////////////////////////////////////////
    
    
    
    
    
    return 0;
    
    
    
}

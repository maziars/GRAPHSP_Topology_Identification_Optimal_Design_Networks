//
//  IO_MANAGER.cpp
//  PCG
//
//  Created by Maziar Sanjabi on 12/25/15.
//  Copyright © 2015 Maziar Sanjabi. All rights reserved.
//

#include "IO_MANAGER.hpp"




int fread_int(FILE *file, int *integers){
    
    int i=0;
    int num;
    while(fscanf(file, "%d", &num) > 0) {
        integers[i] = num;
        i++;
    }
    fclose(file);
    return i;
    
}


int fread_float(FILE *file, double *floats){
    
    int i=0;
    double num;
    while(fscanf(file, "%lf", &num) > 0) {
        floats[i] = num;
        i++;
    }
    fclose(file);
    return i;
}




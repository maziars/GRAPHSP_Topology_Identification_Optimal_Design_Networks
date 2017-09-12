//
//  IO_MANAGER.hpp
//  PCG
//
//  Created by Maziar Sanjabi on 12/25/15.
//  Copyright © 2015 Maziar Sanjabi. All rights reserved.
//

#ifndef IO_MANAGER_hpp
#define IO_MANAGER_hpp

#include <stdio.h>



int fread_int(FILE *file, int *integers);
int fread_float(FILE *file, double *floats);

#endif /* IO_MANAGER_hpp */

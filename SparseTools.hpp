//
//  SparseVector.hpp
//  SparseTools
//
//  Created by Maziar Sanjabi on 1/1/16.
//  Copyright © 2016 Maziar Sanjabi. All rights reserved.
//

#ifndef SparseVector_hpp
#define SparseVector_hpp

#include <stdio.h>
#include <iostream>
#include <map>
#include <string>
#include <iterator>

//void SparseVec::insertNZ(<#int idx#>, <#double value#>);


typedef std::map<int, double>::iterator NZiterator;
class SparseVec;




using namespace std;





class SparseVec {
    
    
public:
    int size, nz;
    
    
    std::map<int, double> NZ;
    
    SparseVec(int n){
        size = n;
        nz = 0;
    };
    

    
    int insert(int idx, double value){
        
        // Check if insertion is successful or not
        int result = 0;
        if (idx>-1 && idx< size ) {
            if(NZ.insert(std::make_pair(idx, value)).second == false)
            {
                NZ[idx] = value;
                if (NZ[idx] != value) {
                    printf("(index, value) pair (%d, %f) not properly set", idx, value);
                    result = -1;
                }
                
            }
            else
                result = 1;
            if (result>-1) {
                nz = nz + result;
                
            }
        }
        else{
            printf("index %d out of bound !!!", idx );
            result = -1;
        }
        
        return result;
    };
    
    NZiterator getNZbegin(){
        return NZ.begin();
    };
    
    NZiterator getNZend(){
        return NZ.end();
    };
    
    int print(){
        NZiterator it = NZ.begin();
        while(it != NZ.end())
        {
            std::cout<<it->first<<" :: "<<it->second<<std::endl;
            it++;
        }
        return 1;
    };
    
    
    int update(int idx, double value){
        // Check if insertion is successful or not
        int result = -1;
        if (idx>-1 && idx<size) {
            NZiterator f = NZ.find(idx);
            if(f!= NZ.end())
            {
                result = 0;
            }
            else
                result = 1;
            NZ[idx] = value;
            if(NZ[idx] != value)
                result = -1;
            else
                nz = nz + result;
        }
        return result;
    };
    
    
    int getValue(int idx, double *value){
        int result = 1;
        if (idx>-1 && idx<size) {
            NZiterator f = NZ.find(idx);
            if(f != NZ.end())
                *value = f->second;
            else
            {
                result = 0;
                *value = 0.0;
            }
            
        }
        else{
            result = -1;
        }
        return result;
        
        
    };
    
    
    int erase(int idx){
        int result = 1;
        std::map<int, double>::iterator f;
        f = NZ.find(idx);
        if (f != NZ.end()) {
            NZ.erase(f);
        }
        else
            result = 0;
        
        nz = nz - result;
        return result;
    };
    
    
    double sum(){
        double sum = 0.0;
        NZiterator it = NZ.begin();
        while(it != NZ.end())
        {
            sum = sum + it->second;
            it++;
        }
        return sum;
    }
    
    
    
};













#endif /* SparseVector_hpp */

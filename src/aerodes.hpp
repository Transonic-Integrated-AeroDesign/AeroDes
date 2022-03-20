#ifndef AERODES_H
#define AERODES_H

#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream> // std
#include <iomanip>  // setw
#include <fstream>  // fopen
#include <string>
#include <math.h>

class aerodes{
    public:
        aerodes(int argc, char** argv);
        ~aerodes();
        
        // memory
        double* create_1d_double_array(int n1, double *array);
        double** create_2d_double_array(int n1, int n2, double **array);
        void delete_2d_double_array(double **array);
    
    private:
};

#endif

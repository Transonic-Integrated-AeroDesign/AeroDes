#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream> // std
#include <iomanip>  // setw
#include <sstream>  // istream
#include <fstream>  // fopen, ifstream
#include <string>
#include <stdio.h>  // strcpy
//#include <string.h>
#include "aerodes.hpp"

#ifndef DBG
#define DBG 0
#endif

using namespace std; // g++ aerodes.cpp -c

aerodes::aerodes(int argc, char** argv) {
}

aerodes::~aerodes() {
}

double *aerodes::create_1d_double_array(int n2, double *array) {
    // create a n2 x 1 matrix
    array = (double *) malloc(n2*sizeof(double));
    return array;
}

double **aerodes::create_2d_double_array(int n1, int n2, double **array) {
    // create a n1 x n2 matrix
    int n=0;
    double *data = (double *) malloc(n1*n2*sizeof(double));
    for (int i=0; i<n1; i++) {
        array[i] = &data[n];
        n += n2;
    }
    return array;
}

void aerodes::delete_2d_double_array(double **array) {
    free(array[0]);
    free(array);
}

int main(int argc, char** argv) {
}



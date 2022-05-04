/*
 * ©2022 The Regents of the University of California.  All rights reserved.
 */

#ifndef AD_MEMORY_H
#define AD_MEMORY_H

#include <cstdlib>
#include "ADvariables.hpp"

//#ifndef DBG
//#define DBG 0
//#endif

class ADmemory : virtual public ADvariables {
public:
    ADmemory(AD *adshr) : ADvariables(adshr) {
        // intialize pointers
        alr = NULL;
        ald = NULL;
        cl_al = NULL;
        cd_al = NULL;
        cq_al = NULL;
    };

    ~ADmemory() {
    }

    int **create_2d_int_array(int n1, int n2, int **&array) {
        //
        // create a n1 x n2 matrix
        //
        int n = 0;
        int *__restrict data = (int *) malloc(n1 * n2 * sizeof(int));
        array = (int **) malloc(sizeof(int *) * n1);
        for (int i = 0; i < n1; i++) {
            array[i] = &data[n];
            n += n2;
        }
        return array;
    }

    void delete_2d_int_array(int **__restrict &array) {
        delete[] array[0];
        delete[] array;
    }

    double *create_1d_double_array(int n1, double *&array) {
        // create a n1 x 1 matrix
        array = (double *) malloc(n1*sizeof(double));
        return array;
    }

    double **create_2d_double_array(int n1, int n2, double **&array) {
        //
        // create a n1 x n2 matrix
        //
        int n = 0;
        double *data = (double *) malloc(n1 * n2 * sizeof(double));
        array = (double **) malloc(sizeof(double *) * n1);
        for (int i = 0; i < n1; i++) {
            array[i] = &data[n];
            n += n2;
        }
        return array;
    }

    double ***create_3d_double_array(int n1, int n2, int n3, double ***&array) {
        //
        // create a n1 x n2 x n3 matrix
        //
        int n = 0, m = 0;
        double *data = (double *) malloc(n1 * n2 * n3 * sizeof(double));  // vector R^{n} ; n = n1*n2*n3
        double **plane = (double **) malloc(n1 * n2 * sizeof(double *)); // matrix R^{n1 x n2}
        array = (double ***) malloc(n1 * sizeof(double **)); // matrix R^{n1 x n2}
        for (int i = 0; i < n1; ++i) {
            m = i * n2;
            array[i] = &plane[m];
            for (int j = 0; j < n2; ++j) {
                plane[m + j] = &data[n];
                n += n3;
            }
        }
        return array;
    }

    void delete_1d_double_array(double *array) {
        free(array);
    }

    void delete_2d_double_array(double **array) {
        free(array[0]);
        free(array);
    }

    void delete_3d_double_array(double ***array) {
        free(array[0][0]);
        free(array[0]);
        free(array);
    }
};
#endif
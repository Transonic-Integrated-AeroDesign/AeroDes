#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream> // std
#include <iomanip>  // setw
#include <fstream>  // fopen
#include <string>
#include <math.h>   // pow

#include "wake.hpp"

/*
 * compile:
 * 	g++ -Xpreprocessor -fopenmp -framework Accelerate -o test main.cpp algorithmic.cpp
 * 	g++ -o test main.cpp wake.cpp
 *  ./test -in wake.data
 *	valgrind --leak-check=yes ./test
 */

using namespace std;

int main(int argc, char** argv) {
    wake *wk = new wake();

    wk->cmdInput(argc, argv);
    wk->readInputParams();
    wk->printInputParams();

    wk->readInputPolar("polarbl.dat");
    //wk->printXFoilMaxValues();

    //wk->setMesh();
    //wk->printGeomSummary();

    //wk->solveLiftingLine();
    //wk->printDistributions();

    delete wk;
}

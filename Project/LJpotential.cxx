/*
    Philip Krause
    118470776@umail.ucc.ie
    Project 9 Quantum Assignment
    file for generating non-symmetrical potentials
*/

// Header Files
#include <iostream>
#include <string>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <cmath>
using namespace std;
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include "gnuplot.cxx"

// Lenard-Jones Potential function
double *LJV(double *V, double *X, double starting_x_nz, double endpoint, int size, double sigma, double epsilon) {
    
    
    // function to calculate the Lenard-Jones Potential for a range of points

    /*
    L-J potential given by the following formula
    V(x) = 4e( (sigma/x)^(12) - (sigma/x)^(6) )

    */

    // decalre del_x is the spacial step such that you have "size" number of steps between endpoint and starting_x_nz
    double del_x = (endpoint - starting_x_nz)/size; 

    // x is the current position i.e. x = x_0 +i*del_x kind of thing
    double x  = starting_x_nz;

    // loop through each point from starting_x_nz
    for (int i = 0; i < size; i++) {
        // assign the current calculated the position value, x = (x_0 + i*del_x)
        X[i] = x;
        // calculate the potential at the point x
        V[i] = 4*epsilon*(pow((sigma/x),12)  -pow((sigma/x), 6));
        // calculate x for the next iteration
        x += del_x;
    }
    // return the potential function
    return V;
}


/*

// test the above functio
int main() {
    int N = 10000;
    double end = 5;
    double V[N], X[N];
    double x_0 = 0.9;
    double sig =1;
    double ep = 10;


    LJV(V, X, x_0, end, N, sig, ep);
            cout << "hi 2" << endl;

    gnuplot_one_function("LJV", "linespoints", "x-axis", "V(x)", X, V, N);
    //gnuplot_one_function_jpg("Lenard-Jones Potential", "linespoints", "X-axis", "V(x)", X, V, N, "Lenard-Jones Potential function.jpg");
    return 0;
}
*/
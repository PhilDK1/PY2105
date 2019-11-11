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
    double del_x = (endpoint - starting_x_nz)/size;   
    double x  = starting_x_nz;
    for (int i = 0; i < size; i++) {
        X[i] = x;
        V[i] = 4*epsilon*((sigma/x)*(sigma/x)*(sigma/x)*(sigma/x)*(sigma/x)*(sigma/x)*(sigma/x)*(sigma/x)*(sigma/x)*(sigma/x)*(sigma/x)*(sigma/x) -(sigma/x)*(sigma/x)*(sigma/x)*(sigma/x)*(sigma/x)*(sigma/x));
        x += del_x;
    }

    return V;
}
/*
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
    gnuplot_one_function_jpg("Lenard-Jones Potential", "linespoints", "X-axis", "V(x)", X, V, N, "Lenard-Jones Potential function.jpg");
    return 0;
}
*/
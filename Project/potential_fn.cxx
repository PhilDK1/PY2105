/*
    Philip Krause
    118470776@umail.ucc.ie
    Project 9 Quantum Assignment
    file for generating symmetrical potentials
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

/*  file to genereate potential function
    using pointers which are the easiest ways of returning an array

    for the arguements
        double* V, is the pointer to the memory address for an array of doubles which will serve as
            the potenitial at a given point (as an array of potential for various points)
        double* X, is the pointer to the memory addess for the x axis and is used to determine where
            the potential stepts are located
        int size, the number of points
        double lV_step, the potential step on the left hand side
        double rV_step, the potential step on the right hand side
        double graphing_distance, the length of the ranges of x axis that will be plotted +/- half
            of the value
        double l_step_at, the x coord where the potential step is on the left
        double r_step_at, the x coord where the potential step is on the right
*/
double* gen_v(  double* V,
                double* X, 
                int size, 
                double V_L,  
                double graphing_distance, 
                double L) {
    

    // check that the graph contains the potential step, and isn't just graphing the potential of 0 inside +/-L
    if (abs(2*L) > graphing_distance) cout << "check points not all in range" << endl;

    //calculate the space step del_x
    double del_x = graphing_distance/size;
    //calculate the start point (left)
    double starting = -0.5*graphing_distance;
    
    // Generate the function X axis and the potential steps
    for (int i = 0; i < size; i++) {
        // the start point and each step incriments by del_x
        X[i] = starting + (i*del_x);
        // check if the X value is greater than +L
        if (X[i] > L){
            // set the potential to the potential at the point
            V[i] = V_L;
        // check if the X value is less than -L
        } else if (X[i] < -L) { 
            // set the potential to the potential at the point
            V[i] = V_L;
        // else the x value is between +/- L and is 0
        } else {
            V[i] = 0.0;
        }

    }
    // return the 2 arrays (not sure if this is required but it works so i didn't mess around with it)
    return V, X;
}

/*

runing the file to see if it all worked

int main() {
    // define variables and assign values

    int N = 1001;
    double V_at_step_l = 4000;
    double V_at_step_r = 4000;
    double point = 1;
    double distance = 3;
    double del_x = distance/N;
    double X[N], V[N], Index[N];


    // use function
    gen_v(V, X, N, V_at_step_l, V_at_step_r, distance, -point, point);

    // plot function
    gnuplot_one_function("Test of potential generation", "linespoints", "x", "V", X, V, N);
    return 0;
}
*/
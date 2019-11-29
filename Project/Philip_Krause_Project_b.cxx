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
        double V_L, the potential step on either side
        double graphing_distance, the length of the ranges of x axis that will be plotted +/- half
            of the value
        double L, the =/-x coord where the potential step is on either side
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


double *mini_step(double* V,
                double* X, 
                int size, 
                double V_L,  
                double graphing_distance, 
                double L,
                double mV_L,
                double m_L){
    /*
    args as in "gen_v" but with addition of,

    double mV_L, the mini potential step (barrier)
    double m_L, the +/- x coord, between which the barrier is located
    */


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
        } else if ((X[i] > -m_L) && (X[i] < m_L)) { 
            V[i] = mV_L;
        } else {
            V[i] = 0.0;
        }

    }
    // return the 2 arrays (not sure if this is required but it works so i didn't mess around with it)
    return V, X;
}


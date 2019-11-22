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

// #include "gnuplot.cxx"
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

double *mini_step(double* V,
                double* X, 
                int size, 
                double V_L,  
                double graphing_distance, 
                double L,
                double mV_L,
                double m_L){
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


/*
// runing the file to see if it all worked

int main() {
    // define variables and assign values

    int N = 1001;
    double V_at_step_l = 4000;
    double V_at_step_r = 4000;
    double mVstep = 100;
    double m_l = 0.25;
    double point = 1;
    double distance = 3;
    double del_x = distance/N;
    double X[N], V[N], Index[N];


    // use function
    // gen_v(V, X, N, V_at_step_l, V_at_step_r, distance, -point, point);
    mini_step(V, X, N, V_at_step_l, distance, point, mVstep, m_l);

    // plot function
    gnuplot_one_function("Test of potential generation", "linespoints", "x", "V", X, V, N);
    return 0;
}

*/

/*
int main() {
    int N = 10000;
    double end = 3;
    double V[N], X[N];
    double x_0 = 0.95;
    double sig =1;
    double ep = 10;


    LJV(V, X, x_0, end, N, sig, ep);
            cout << "hi 2" << endl;

    gnuplot_one_function("LJV", "linespoints", "x-axis", "V(x)", X, V, N);
    gnuplot_one_function_jpg("Lenard-Jones Potential", "linespoints", "X-axis", "V(x)", X, V, N, "Lenard-Jones Potential function2.jpg");
    return 0;
}
*/
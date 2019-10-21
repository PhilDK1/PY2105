/*
Program to model a dampened driven pendulum for various parameters

Philip Krause
118470776
*/

//Header files
#include <iostream>
#include <string>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
using namespace std;
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include "gnuplot-1.cxx"



//Main Function
main() {
    //Declare Variables to be used in program

    // f_0 the magnituted of the driving force
    double f_0;

    // T the final time
    double T;

    // Capital Omega, _O_ , as a double
    double _O_;


    // Capital Omega of the driving force, _O_D , as a double
    double _O_D;

    // gamma the dampening
    double y;

    // N the number of time steps
    int N=300;

    // time array, t
    double t[N];

    // the angular velocity array, w
    double w[N];

    // angular position theta
    double theta[N];

    // the time step
    double dt;





    // --- Run 1 --- //

    // Define variables from run 1
    T = 10;
    f_0 = 0;
    _O_ = 1;
    theta[0] = 1;
    w[0] = 0;

    //Array for the gammas
    double gamma_array[4] = {0.0, 0.5, 1.0, 1.5};

    //so there's no variable undeclared/unassigned warning
    _O_D = 1;
    
    for (int n = 0; n<4;n++) {
        // Calculate the time step dt
        dt = T / N;
        y = gamma_array[n];

        //define variable name for the output file
        // make array of chars big enough to hold entire filename
        char filename[30];
        // sprintf allows for the formatting of strings so the file name can be generated on each run
        sprintf(filename, "Krause_Philip_PS4_a_%.1lf.jpg", gamma_array[n]);
    
        // Let initial time be 0
        t[0] = 0.0;

        for (int iter = 1; iter < N+1; iter++) {
            t[iter] = iter*dt;
            w[iter] = w[iter - 1] 
                    + dt*(

                    -1*(_O_*_O_*sin(theta[iter - 1])) 

                    - 2*(y*w[iter - 1]) 

                    + f_0*cos(_O_D*t[iter])

                    );
            theta[iter] = theta[iter - 1] + dt*(w[iter]);
        }

    gnuplot_one_function_jpg("Graph of θ vs time","linespoints","time (sec)", "θ (rad)", t, theta, N,filename);
    }






    // --- Run 2 --- //


    // Define variables from run 2
    T = 50;
    f_0 = 0.5;
    _O_ = 1;
    _O_D = 2/3;
    theta[0] = 0.2;
    w[0] = 0;

    //Array for the gammas
    double gamma_array_2[2] = {0.25, 1.0};
    
    for (int n = 0; n<2;n++) {
        // Calculate the time step dt
        dt = T / N;
        y = gamma_array_2[n];

        //define variable name for the output file
        // make array of chars big enough to hold entire filename
        char filename[30];
        // sprintf allows for the formatting of strings so the file name can be generated on each run
        sprintf(filename, "Krause_Philip_PS4_b_%.2lf.jpg", gamma_array_2[n]);
    
        // Let initial time be 0
        t[0] = 0.0;

        for (int iter = 1; iter < N+1; iter++) {
            t[iter] = iter*dt;
            w[iter] = w[iter - 1] 
                    + dt*(

                    -1*(_O_*_O_*sin(theta[iter - 1])) 

                    - 2*(y*w[iter - 1]) 

                    + f_0*cos(_O_D*t[iter])

                    );
            theta[iter] = theta[iter - 1] + dt*(w[iter]);
        }

    gnuplot_one_function_jpg("Graph of θ vs time","linespoints","time (sec)", "θ (rad)", t, theta, N,filename);
    }


}
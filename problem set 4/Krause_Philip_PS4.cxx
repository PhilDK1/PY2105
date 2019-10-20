/*
Program to model a dampened driven pendulum

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
    int N;

    // Prompt for N so the next variables are declarable
    cout << "Enter the number of time steps to be taken: ";
    cin >> N;

    // time array, t
    double t[N];

    // the angular velocity array, w
    double w[N];

    // angular position theta
    double theta[N];

    // the time step
    double dt;

    // Prompt for the constants
    cout << "Enter the driving force: ";
    cin >> f_0;

    cout << "Enter the final time: ";
    cin >> T;

    cout << "Enter the value for Capital Omega: ";
    cin >> _O_;

    cout << "Enter the value for gamma, the dampening: ";
    cin >> y;

    cout << "Enter the value for the initial angular position: ";
    cin >> theta[0];

    cout << "Enter the value for the initial angular velocity: ";
    cin >> w[0];

    cout << "Enter the value for Omega D:";
    cin >> _O_D;

    // Calculate the time step dt
    dt = T / N;
    
    // Let initial time be 0
    t[0] = 0.0;

    for (int iter = 1; iter < N+1; iter++) {
        t[iter] = iter*dt;
        w[iter] = w[iter - 1] 
                + dt*(

                -1*(_O_*sin(theta[iter - 1])) 

                - 2*(y*w[iter - 1]) 

                + f_0*cos(_O_D*t[iter])

                );
        theta[iter] = theta[iter - 1] + dt*(w[iter]);
    }

    gnuplot_one_function("graph of theta vs time","linespoints","time (sec)", "theta (rad)", t, theta, N);


    
}
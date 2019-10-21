/*
 * Philip Krause
 * 
 * 118470776
 *
 * PY2105
 *
 *
 * File used to approximate the velocity using the Euler Method of numerical differentitation
 * and plots the velocity at a given time t on the interval [0, T]. Given that x(t) = x(0) + v(0)*t + 0.5*a*t^2
 * where the velocity is constant so the above equation goes to x(t) = x(0) + v*t
 * 
 *
 */


//header files 
#include <iostream>
#include <string>
#include <stdio.h>
#include <unistd.h>
#include <math.h>


// namespace set up to prevent the need to type "std::" in each place where it would normally be required
using namespace std;

//More headers
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include "gnuplot-1.cxx"

//Main Function
int main(int argc, char **argv) {
	// Initialize variables to be used in the program
	
	// Velocity variable, v
	double v;
	
	// position array with 1000 points
	double x[1000];
	
	// final time T
	double T;
	
	// time step, d_t
	double d_t;
	
	// time array, t
	double t[1000];
	
	// Get user input for the velocity, v
	cout << "Enter the velocity of the particle: ";
	cin >> v;
	
	//Get user input of the initial point x[0]
	cout  << "Enter the initial point of the particle: ";
	cin >> x[0];
	
	// Get User input for the fianl time to be plotted
	cout << "Enter the final time to be plotted: ";
	cin >> T;
	
	// calculate the time step d_t
	d_t = T/1000;
	
	t[1] = 0;
	for (int i = 1; i <1000 ; i++) {
		t[i] = d_t*i;
		x[i] = x[i-1] + d_t*v;
		}
	
	gnuplot_one_function ("solving the differential equation dx(t)/dt = v","linespoints", "time", "position", t, x, 1000);
	
	return 0;
	}

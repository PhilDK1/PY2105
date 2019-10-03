/*
 * Krause_Philip_2_4.cxx
 * 
 * Copyright 2019 118470776 <118470776@PHYSPML200>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */



 //header files 
#include <iostream>
#include <string>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
//so std::cin and std:: cout isn't required
using namespace std;
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include "gnuplot-1.cxx"


// declare constants

// accelerationdue to gravity, g
const double g = 9.81;

//Main Function
int main(int argc, char **argv) {

    //initialise intial value variables, initial velocity, height above ground level and final time
    double v0, x0, tf;

    //number of plot points N
    int N;

    // Prompt user for the above variables

    cout << "enter initial velocity, v0: ";
    cin >> v0;
    cout << "enter initial height, x0: ";
    cin >> x0;
    cout << "final time, tf: ";
    cin >> tf;

    cout << "enter number of plot points, N: ";
    cin >> N;
    
    // initialise arrays for gnuplot to plot
    double t [N];
    double Xoft [N];
    double Voft [N];

    
    // Set up variables for gnuplot and declare them
    // x-axis (time, t) start point, t_min
    double t_min = 0.0;
    // x-axis (time, t) end point, tf
    
    //time step, delta_t
    double delta_t = (tf - t_min)/(N- 1.0);

    cout << "t_min = " << t_min << ", t_max= " << tf << ", delta_t = " << delta_t << "\n";


    for (int i = 0; i < N; i++) {
        t[i] = t_min + (delta_t * i);
        Voft[i] = g*t[i] + v0;
		double x_val = (-0.5*g*(t[i]*t[i]))+ (v0*t[i])+ x0;
		if (x_val <0) {
			Xoft[i] = 0.0;
		} else {
			Xoft[i] = x_val;
		}

    }

	//edit some more

    // Plot x and y arrays 
	cout << "Plot V(t) = g*t + v0 \n";  gnuplot_one_function ("Plot of v(t)","linespoints", "x", "y", t, Voft, N);
    // Plot x and y arrays
	cout << "Plot x(t) = -(1/2)gt^2 + v0*t + x0\n";  gnuplot_one_function ("Plot of x(t)","linespoints", "x", "y", t, Xoft, N);
	
	

    return 0;
}

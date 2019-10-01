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

/*
 * defining up as positive
 */


// accelerationdue to gravity, g
const double g = 9.81;

//Main Function
int main(int argc, char **argv) {

    /*initialise intial value variables, initial velocity, height above 
     * ground level and final time
     */ 
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
    /*time step, delta_t
     * (final time - initial time) / number of points
     */
    double delta_t = (tf - t_min)/(N- 1.0);

	// prininting out the the start time, max time and the time step
    cout << "t_min = " << t_min << ",\nt_max= " << tf << ",\ndelta_t = " << delta_t << "\n";



	/* for loop to iterate through each of the points and calculate the
	 * at each time step t[i]
	 */
    for (int i = 0; i < N; i++) {
		// calculate t at iteration i (start time + iteration * time step)
        t[i] = t_min + (delta_t * i);
		
		/*check the value for height and if it's less than 0 then set
		 * the height as 0 as it's hit the ground
		 * 
		 * else evaluate height and store that value in the array
		 *
		 * 
		 * as above if the position is less than 0 assume the particle 
		 * has hit the ground and has come to a full stop hence velocity 
		 * is 0
		 * 
		 * else evaluate velocity and store it in array
		 */
		
		
		Xoft[i] = (0.5*g*(t[i]*t[i]))+ (v0*t[i])+ x0;
		Voft[i] = g*t[i] + v0;

	}

    // Plot x and y arrays 
	cout << "Plot V(t) = g*t + v0 \n";  gnuplot_one_function ("Plot of v(t)","linespoints", "t", " velocity v(t)", t, Voft, N);
    // Plot x and y arrays
	cout << "Plot x(t) = -(1/2)gt^2 + v0*t + x0\n";  gnuplot_one_function ("Plot of x(t)","linespoints", "t", "displacement x(t)", t, Xoft, N);
	
	

    return 0;
}

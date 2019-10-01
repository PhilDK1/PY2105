/*
 * test.cxx
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


const int number_points = 20; // Numer of Points in plot


int main(int argc, char **argv)
{
	
	double x [number_points];
	double y [number_points];
	
	double x_min; // Minimal x value of plot
	double x_max; // Maximal x value of plot  
	double delta_x; // Stepsize 
	int nr;
	
	
	x_min = 0.0;  x_max = 1.0;  delta_x = (x_max - x_min) / (number_points-1.0);
	
	cout << "x_min = " << x_min << ", x_max= " << x_max << ", delta_x = " << delta_x << "\n";
	
	nr = 0;
	
	while (nr < number_points) {
		x [nr] = x_min + delta_x * nr;
		y [nr] = x[nr] * x[nr];      
		nr = nr + 1;    
	} 
	
	// Plot x and y arrays
	cout << "Plot f(x) = x*x!\n";  gnuplot_one_function ("First plot","linespoints", "x", "y", x, y, number_points);
	
	
	return 0;
}


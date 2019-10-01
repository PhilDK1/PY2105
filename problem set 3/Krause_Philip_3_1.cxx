/*
 * Krause_Philip_3_1.cxx
 * 
 * Copyright 2019 118470776 <118470776@TRAJAN>
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

/*
 * Philip Krause
 * 
 * File that takes in 2 vectors the Velocity of a particle and the
 * magnetic field vector and the charge on a particle and prints the
 * resulting vector for the force on the product
 * 
 * where f = q*v X B, where v is the veolcity vector B the magnetic
 * field vector and q the charge on the particle * denotes scalar 
 * multiplication and X denotes the vector cross product
 * 
 */


//header files 
#include <iostream>
#include <string>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

using namespace std;

double * Compute_cross_product(double velocity[3],double B_field[3]) {
	
	static double vect_product[3];
	
	vect_product[0] = ((velocity[1]*B_field[2])- (velocity[2]*B_field[1]));
	vect_product[1] = -((velocity[0]*B_field[2])- (velocity[2]*B_field[0]));
	vect_product[2] = ((velocity[0]*B_field[1])- (velocity[1]*B_field[0]));
	
	
	return vect_product;
	
	}



int main(int argc, char **argv)
{
	double vel[3] = {1, 2, 3};
	double b[3] = {4, 7, 9};
	double q = 2.0;
	double *result;
	double vect_product[3];
	
	for (int i =0; i<3; i++) {
		
		vect_product[i] = *(result + i) ;
		
	}
	
	result = Compute_cross_product(vel, b);
	
	for (int i =0; i<3; i++) {
		cout << q * vect_product[i] << endl;
	}
	return 0;
}


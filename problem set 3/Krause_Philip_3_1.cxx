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


// Main Function
int main(int argc, char **argv)
{
	// Initialise variables that will be used in the program

	//Charge 
	double q;
	//array for the velocity vector
	double vel[3];
	//array for the magnetic field, B
	double B[3];
	// array for the vector product of the velocity and and B field, along with the charge
	double product[3];

	// Get user Input for the velocity Vector
	cout << "Enter the x value of the velocity: ";
	cin >> vel[0];
	cout << "Enter the y value of the velocity: ";
	cin >> vel[1];
	cout << "Enter the z value of the velocity: ";
	cin >> vel[0];

	// Get user Input for the Charge
	cout << "Enter the charge of the particle: ";
	cin >> q;

	//Get user input for the Magnetic Field, B
	cout << "Enter the x value of the Magnetic Field: ";
	cin >> B[0];
	cout << "Enter the y value of the Magnetic Field: ";
	cin >> B[1];
	cout << "Enter the z value of the Magnetic Field: ";
	cin >> B[2];

	// Calculate the cross product

	/*
	*					| i  j  k  |
	* F = q*vel x B	= 	| vx vy vz | = q(((vy*Bz)-(vz*By))i - ((vx*Bz)-(vz*Bx))j + ((vx*By)-(vy*Bx))k)
	*					| Bx By Bz |
	*
	*/

	product[0] = q*((vel[1]*B[2]) - (vel[2]*B[1]));
	product[1] = q*((vel[0]*B[2]) - (vel[2]*B[0]));
	product[2] = q*((vel[0]*B[1]) - (vel[1]*B[0]));


	return 0;
}

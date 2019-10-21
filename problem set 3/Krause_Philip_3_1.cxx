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

// namespace set up to prevent the need to type "std::" in each place where it would normally be required
using namespace std;


// Main Function
int main(int argc, char **argv)
{
	// Initialise variables that will be used in the program

	//Charge 
	double q;
	//array for the velocity vector, v_
	double v_[3];
	//array for the magnetic field, B_
	double B_[3];
	// array for the vector product of the velocity and and B field, along with the charge force, F_
	double F_[3];
	
	//direction array
	char axis[3] = {'x', 'y', 'z'};

	// Get user Input for the velocity Vector, v_
	
	for (int i = 0; i < 3; i++) {
		cout << "Enter the " << axis[i] << " value of the velocity: ";
		cin >> v_[i];
		}

	// Get user Input for the Charge
	cout << "Enter the charge of the particle: ";
	cin >> q;

	//Get user input for the Magnetic Field, B
	for (int i = 0; i < 3; i++) {
		cout << "Enter the " << axis[i] << " value of the  Magnetic Field: ";
		cin >> B_[i];
		}

	// Calculate the cross product

	/*
	*					| i_  j_  k_  |
	* F_ = q*v_ x B_ = 	| v_x v_y v_z |
	*					| B_x B_y B_z |
	*
	*
	*	= q(((v_y*B_z)-(v_z*B_y))i_ - ((v_x*B_z)-(v_z*B_x))j_ + ((v_x*B_y)-(v_y*B_x))k_)
	*/

	F_[0] = q*((v_[1]*B_[2]) - (v_[2]*B_[1]));
	F_[1] = q*((v_[0]*B_[2]) - (v_[2]*B_[0]));
	F_[2] = q*((v_[0]*B_[1]) - (v_[1]*B_[0]));

	// print out the answers
	for (int i = 0; i < 3; i++) {
		cout << "F_" << axis[i] << " = " <<F_[i] << endl;
		}
	
	

	return 0;
}

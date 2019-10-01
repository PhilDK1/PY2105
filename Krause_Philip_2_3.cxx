/*
 * Krause_Philip_2_3.cxx
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

// Main function
int main(int argc, char **argv)
{
	//intialise 2 arrays of doubles length 3 and a variable to be used to store the anser (scalar_product)
	double vector_1[3];
	double vector_2[3];
	double scalar_product = 0.0;
		
	// prompt user for each element of the vector/array
	for (int i = 0; i < 3; i++) {
		cout << "Enter number " << i+1 << " of vector 1: ";
		//store value in the ith position of the array
		cin >> vector_1[i];
		}
	

	// prompt user for each element of the vector/array 2
	for (int i = 0; i < 3; i++) {
		cout << "Enter number " << i+1 << " of vector 2: ";
		//store value in the ith position of the array

		cin >> vector_2[i];
		}
	
	
	// loop throught the elements in the array
	// multiply the the arrays element wise
	for (int i = 0; i < 3; i++) {
		// add the product of any 2 corresponding vectors together
		scalar_product = scalar_product + (vector_1[i] * vector_2[i]);		
		}
	// print the result
	cout << scalar_product<< endl;

	
	
	return 0;
}


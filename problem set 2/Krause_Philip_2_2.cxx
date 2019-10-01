/*
 * Krause_Philip_2_2.cxx
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

/*   
 * Philip Krause
 * 118470776
 * Problem 2.2
 * 
 * Get 10 inputs (doubles) and puts them in a 
 */


//header files 
#include <iostream>
#include <string>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

//so std::cin and std:: cout isn't required
using namespace std;


int main(int argc, char **argv)
{
	// initialise array
	// double arr[10] = input_func();
	
	//cout << arr;
	//initialise array of doubles of length 10 called "array"
	double array[10];

	//intitialise "max" variable to be used in comparisons
	double max;
	

	// loop through each loaction in the array and prompt for a user to store in said location
	for (int i = 0; i < 10; i++) {
		// when prompted for the 1st entry that value is stored in the 0th position
		cout << "Enter number " << i+1 << " of array: ";
		// store inout in array
		cin >> array[i];
		}
	

	//assume the max value is the first value
	max = array[0];
	
	/*loop through the values starting as the 1th position (as opposed to the 1st)
	* and compare the selected value with the current max and if the selected value is greater
	* than the current max values assign the selected value to max and continue else continue
	* for all values in the array.
	*/
	for (int n = 1; n < 10; n++) {
	if (array[n] > max) {
			max = array[n];
		}
	
	}

	//print out the max value
	cout << max << endl;
	
	return 0;
}


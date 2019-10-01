/*
 * Krause_Philip_2_1.cxx
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
 * Problem 2.1
 * make sure user inputs a number between 1 and 99
 * then prints that number to the screen
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


// separate nction for geting user input cause i'm lazy
// function breaks if a non number is entered eg "777+"
int get_input(){
	//variable for the user input
	int n;
	
	//user prompt
	cout << "Enter a number between 1 and 99: ";
	
	//actually get number
	cin >> n;
	
	//return number to calling function
	return n;
}


//Main function
int main(int argc, char **argv)
{
	//initialise user input variable "num" and get user input by calling "get_input fn"
	int num = get_input();
	
	
	//while loop to ensure number is betwee 1 and 99 prompt for valid input if input not between sepecified numbers
	while ((num < 1) || (num > 99)) {
		cout << "Invalid input!!\n";
		num = get_input();
	}
	
	//pint the valid number
	cout << num << endl;
	
	
	return 0;
}


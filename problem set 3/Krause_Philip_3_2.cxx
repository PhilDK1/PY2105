/*
 * Philip Krause
 *
 * Krause_Philip_3_2.cxx
 * 
 * File that calculates the work done by a 3 dimensional force vector dotted
 * with a 3 dimensional displacement vector
 *
 * Work, W, Force, F_, displacement, s_
 * 
 * 
 * W = F_ . s_
 */


//header files 
#include <iostream>
#include <string>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

using namespace std;

// Main Function
int main(int argc, char **argv) {

    // Initialize variables to be used in Program

    // Work scalar, W
    double W;

    // Force vector, F_
    double F_[3];

    // displacement vector, s_
    double s_[3];
    
    //direction array
	char axis[3] = {'x', 'y', 'z'};


    // Prompt user for the values of x, y and z of Force
    
    for (int i = 0; i < 3; i++) {
		cout << "Enter the value of Force in the "<< axis[i] << " direction: ";
		cin >> F_[i];
		}

    // Prompt user for the values of x, y and z of the displacement
    for (int i = 0; i<3; i++) {
		cout << "Enter the value of displacement in the "<< axis[i] <<" direction: ";
		cin >> s_[i];
		}

	
    //Calculate the dot product
    W = 0.0;
    for (int i = 0; i < 3; i++) {
        W += (F_[i] * s_[i]);
    }

    // Print Answer#
    cout << W<< endl;

    return 0;
}

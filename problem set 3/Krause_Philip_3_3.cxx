
/*
 * Philip Krause
 * 
 * File that takes in a vector of doubles, U_ and a 3x3 Matrix, M and returns the matrix vector product
 */


//header files 
#include <iostream>
#include <string>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

// namespace set up to prevent the need to type "std::" in each place where it would normally be required
using namespace std;

//Main Function
int main(int argc, char **argv) {

    //Initialize variables to be used in the program

    // Vector, U_
    double U_[3];

    // Matrix, M
    double M[3][3];

    // Product vector, V_
    double V_[3];
    
    //direction array
	char axis[3] = {'x', 'y', 'z'};
	
	//position array
	char position[3][3]  = {
							{'1', 's', 't'},
							{'2', 'n', 'd'},
							{'3', 'r', 'd'}
							};


    // Get user input for the Vector U_
    for (int i = 0; i < 3; i++) {
        cout << "Enter the " << axis[i] << " value of the the Vector U_: ";
        cin >> U_[i];
    }

    // Get user input for the matrix M
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            cout << "Enter the " << position[j][0] << position[j][1] <<  position[j][2] << " value of the " << position[i][0] << position[i][1] <<  position[i][2] << " row: ";
            cin >> M[i][j];
        }
    }

    /* Calculate the product of the the Matrix M and the Vector v
    *
    * Where the product of the two is defined as below
    *
    *
    *               | a b c |   | x0 |      | a*x0 + b*x1 + c*x2 |
    *   M x U_ =    | d e f | x | x1 | =    | d*x0 + e*x1 + f*x2 |
    *               | g h i |   | x2 |      | g*x0 + h*x1 + i*x2 |
    *
    *
    */

    for (int i = 0; i < 3; i++) {
        V_[i] = 0.0;
        for (int j = 0; j < 3; j++) {
            V_[i] += M[i][j] * U_[j];
        }
    }
    
    // Print the Answer
    for (int i = 0; i < 3; i++) {
        cout << V_[i] << endl;
    }


    return 0;
}

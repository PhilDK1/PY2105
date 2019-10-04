/*
 * Philip Krause
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

    // Prompt user for the values of x, y and z of Force
    cout << "Enter the value of Force in the x direction: ";
    cin >> F_[0];
    cout << "Enter the value of Force in the y direction: ";
    cin >> F_[1];
    cout << "Enter the value of Force in the z direction: ";
    cin >> F_[2];


    // Prompt user for the values of x, y and z of the displacement
    cout << "Enter the value of displacement in the x direction: ";
    cin >> s_[0];
    cout << "Enter the value of displacement in the y direction: ";
    cin >> s_[1];
    cout << "Enter the value of displacement in the z direction: ";
    cin >> s_[2];

    //Calculate the dot product
    W = 0.0;
    for (int i = 0; i < 3; i++) {
        W += (F_[i] * s_[i]);
    }

    // Print Answer#
    cout << W<< endl;

    return 0;
}
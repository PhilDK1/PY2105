/*
    Philip Krause
    118470776@umail.ucc.ie
    Project 9 Quantum Assignment
    testing for wave fn

    File to generate wave functions for symmetrical potentials as well as generating the wave function 
    for odd and even parity solutions for a square well using shooting method
*/

// Header Files
#include <iostream>
#include <string>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <cmath>
using namespace std;
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
// #include "gnuplot.cxx"
// #include "Philip_Krause_Project_b.cxx"

//calculate fn, decides whether the function is diverging and adusts the cahgne in energy accordingly
double* calc(double *V, 
            double *X, 
            double *wavefn,
            int size, 
            double* last_diverge,
            double graphing_distance,
            double* E_0, 
            double* del_E, 
            double cutoff) {
    
    // printiting out the value of the diverging wave
    // cout << wavefn[size -5] << endl;

    // if-else if ladder to check various conditions and making decisions based on that
    if (last_diverge[0] == 0 && wavefn[size - 5] > 0) {
        // first iteration, diverging up, increase energy
        E_0[0] = E_0[0] + del_E[0];
    } else if (last_diverge[0] == 0 && wavefn[size - 5] < 0) {
        // first iteration, diverging down, decrease energy
        E_0[0] = E_0[0] - del_E[0];
    } else if (last_diverge[0] == 1 && wavefn[size - 5] < 0) {
        // previously diverging up, now divering down, overshot
        del_E[0] = del_E[0] /2 ;
        E_0[0] = E_0[0] - del_E[0];
    } else if (last_diverge[0] == 1 && wavefn[size - 5] > 0) {
        // Still divering up, energy needs to increase
        E_0[0] = E_0[0] + del_E[0];
    } else if (last_diverge[0] == -1 && wavefn[size - 5] < 0) {
        // still diverging downwards, energy needs to decrease
        E_0[0] = E_0[0] - del_E[0];
    } else if (last_diverge[0] == -1 && wavefn[size - 5] > 0) {
        // previously diverging down, now diverging up
        // Energy needs to increase
        del_E[0] = del_E[0]/2;
        E_0[0] = E_0[0] + del_E[0];
    }
    // check the sign of the last divergence
    // if lessa than 0 and the absolute value is gr
    if (wavefn[size - 5] < 0 /*&& abs(wavefn[size - 5]) > cutoff*/) {
        last_diverge[0] = -1;
    } else {
        last_diverge[0] = 1;
    }
    //return the calculated values
    return E_0, del_E, last_diverge;


}

// double integrate(double *wavefn,
//                  double *X,
//                  )


double *gen_phi_even(double *wavefn,
                double *V,
                double *X,
                int size,
                double graphing_distance, 
                double E,
                double cutoff) {
    // for a symmetric even parity wave the initial conditions are always
    // these conditions are 1, 1, corresponging to phi(x) and d/dx(phi(x))
    double phi_0 = 1;
    double phi_n1 = 1;

    // calculate the the spacial step
    double del_x = graphing_distance/size;
    
    // intialise variable to be used as temperary location for a given answer for the wave function,
    // to test if it's diverging
    double phi_temp;
    
    // variable to get the mid point of the function
    int mid;
    // check the number of points if even
    if (size % 2 == 0) {
        // the mid point is half the total number of points
        mid = (size/2);
    } else {
        // else take one from the number of points and half it to get the midpoint
        mid = size -1;
        mid = (mid/2);
    }
    
    // work out intial conditions for the wave function
    wavefn[mid] = 2*phi_0 - phi_n1 - 2*del_x*del_x*(E - V[mid-1])*phi_0;
    // work out the initial condition for one step right of 0
    wavefn[mid + 1] = 2*wavefn[mid] - phi_0 - 2*del_x*del_x*(E - V[mid-1])*wavefn[mid];
    // work out the initial condition for one step left of 0
    wavefn[mid - 1] = 2*wavefn[mid] - phi_0 - 2*del_x*del_x*(E - V[mid-1])*wavefn[mid];

    // iterate from the 2nd (greater than 0) to the Nth step and evalue at the wave function as it goes
    for (int i = 2; i < mid; i++) {
        // store the value for the wave function at that point to check for diverence
        phi_temp = 2*wavefn[mid + i - 1] - wavefn[mid + i - 2] - 2*del_x*del_x*(E - V[mid + i - 1])*wavefn[mid + i - 1];

        // check if the wave function is diverging if it's not
        if (abs(phi_temp) < cutoff) {

            // assign the value of that wave fn to the correct location
            wavefn[mid + i] = phi_temp;
        } else {
            // else it is diverging and check whether it is diverging up or down
            if (phi_temp > 0) {
                // if divering up
                for (int n = 0; n < (mid - i); n++)
                {
                    // assign every value greater than the ith as 2.5
                    wavefn[mid + i + n] = 2.5;
                }

                // break out the loop and move on to the next section
                break;
            // else it's diverging down
            } else if (phi_temp <  0) {
                // divering down
                for (int n = 0; n < (mid - i); n++)
                {
                    // assign every value less than the ith as -2.5
                    wavefn[mid + i + n] = -2.5;
                }
                // break out the loop and move on to the next section
                break;
            }
            
            
        }
        
    }


    // iterate from 2 steps less than 0 (mid-th value) to the 0th value 
    for (int i = 2; i < mid+2; i++) {
        // store the value for the wave function at that point to check for diverence
        phi_temp = 2*wavefn[mid - i + 1] - wavefn[mid - i + 2] - 2*del_x*del_x*(E - V[mid - i + 1])*wavefn[mid - i + 1];
        
        // check if the wave function is diverging if it's not
        if (abs(phi_temp) < cutoff) {
            // assign the value of that wave fn to the correct location
            wavefn[mid - i] = phi_temp;
        } else {
            // else it is diverging and check whether it is diverging up or down
            if (phi_temp > 0) {
                // if divering up
                for (int n = 0; n < (mid - i); n++)
                {
                    // assign every value less than the ith as 2.5
                    wavefn[mid - i - n] = 2.5;
                }
                break;
            // else it's diverging down
            } else if (phi_temp <  0) {
                // if divering down
                for (int n = 0; n < ( mid -i); n++)
                {
                    // assign every value less than the ith as -2.5
                    wavefn[mid - i - n] = -2.5;
                }

                // break out of the loop
                break;
            }
        }
        
    }


    // return the calculated wavefunction
    return wavefn;
}


// fn to calculate the wave fn for the odd parity ie phi(x) = -phi(-x)
double *gen_phi_odd(double *wavefn,
                double *V,
                double *X,
                int size,
                double graphing_distance, 
                double E,
                double cutoff) {
    
    // calculate 
    double del_x = graphing_distance/size;
    double phi_0 = 0;
    double phi_n1 = -del_x;


    double phi_temp;

    int mid;
    if (size % 2 == 0) {
        mid = (size/2);
    } else {
        mid = size -1;
        mid = (mid/2);
    }
    
    wavefn[mid] = 2*phi_0 - phi_n1 - 2*del_x*del_x*(E - V[mid-1])*phi_0;
    wavefn[mid - 1] = -2*wavefn[mid] + phi_0 + 2*del_x*del_x*(E - V[mid])*wavefn[mid];
    wavefn[mid + 1] = 2*wavefn[mid] - wavefn[mid-1] - 2*del_x*del_x*(E - V[mid])*wavefn[mid];
    
    
    
    for (int i = 2; i < mid; i++) {
        phi_temp = 2*wavefn[mid+i-1] - wavefn[mid+i-2] - 2*del_x*del_x*(E - V[mid+i-1])*wavefn[mid+i-1];
        if (abs(phi_temp) < cutoff) {
            wavefn[mid + i] = phi_temp;
        } else {
            if (phi_temp > 0) {
                // divering up
                for (int n = 0; n < (mid - i); n++)
                {
                wavefn[mid + i + n] = 2.5;
                }
                break;
            } else if (phi_temp <  0) {
                // divering down
                for (int n = 0; n < (mid - i); n++)
                {
                wavefn[mid + i + n] = -2.5;
                }
                break;
            }
            
            
        }
        
    }

    
    for (int i = 2; i < mid+2; i++) {
        phi_temp = 2*wavefn[mid - i + 1] - wavefn[mid - i + 2] - 2*del_x*del_x*(E - V[mid - i + 1])*wavefn[mid - i + 1];
        if (abs(phi_temp) < cutoff) {
            wavefn[mid - i] = phi_temp;
        } else {
            if (phi_temp > 0) {
                // divering up
                for (int n = 0; n < (mid - i); n++)
                {
                    wavefn[mid - i - n] = 2.5;
                }
                break;
            } else if (phi_temp <  0) {
                // divering down
                for (int n = 0; n < ( mid -i); n++)
                {
                wavefn[mid - i - n] = -2.5;
                }
                break;
            }
        }
        
    }



    return wavefn;
}

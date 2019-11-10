/*
    Philip Krause
    118470776@umail.ucc.ie
    Project 9 Quantum Assignment
    testing for wave fn
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
// #include "potential_fn.cxx"

//calculate fn
double* calc(double *V, 
            double *X, 
            double *wavefn,
            int size, 
            double* last_diverge,
            double graphing_distance,
            double* E_0, 
            double* del_E, 
            double cutoff) {
    
    cout << wavefn[size -5] << endl;

    if (last_diverge[0] == 0 && wavefn[size - 5] > 0) {
        // first iteration, diverging up, increase energy
        E_0[0] = E_0[0] + del_E[0];
    } else if (last_diverge[0] == 0 && wavefn[size - 5] < 0) {
        // first iteration, diverging down, decrease energy
        E_0[0] = E_0[0] - del_E[0];
    } else if (last_diverge[0] == 1 && wavefn[size - 5] < 0) {
        // previously diverging up, now divering down, overshot
        del_E[0] = 2*del_E[0] /3 ;
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
        del_E[0] = 2*del_E[0]/3;
        E_0[0] = E_0[0] + del_E[0];
    }

    if (wavefn[size - 5] < 0 && abs(wavefn[size - 5]) > cutoff) {
        last_diverge[0] = -1;
    } else {
        last_diverge[0] = 1;
    }

    return E_0, del_E, last_diverge;


}

// double integrate(double *wavefn,
//                  double *X,
//                  )


double *gen_phi_even(double *wavefn,
                double *V,
                double *X,
                int size,
                double phi_0,
                double phi_n1,
                double graphing_distance, 
                double E,
                double cutoff) {
    

    double del_x = graphing_distance/size;
    double phi_temp;

    int mid;
    if (size % 2 == 0) {
        mid = (size/2);
    } else {
        mid = size -1;
        mid = (mid/2);
    }
    
    wavefn[mid] = 2*phi_0 - phi_n1 - 2*del_x*del_x*(E - V[mid-1])*phi_0;
    wavefn[mid + 1] = 2*wavefn[mid] - phi_0 - 2*del_x*del_x*(E - V[mid-1])*wavefn[mid];
    wavefn[mid - 1] = 2*wavefn[mid] - phi_0 - 2*del_x*del_x*(E - V[mid-1])*wavefn[mid];
    for (int i = 2; i < mid; i++) {
        phi_temp = 2*wavefn[mid + i - 1] - wavefn[mid + i - 2] - 2*del_x*del_x*(E - V[mid + i - 1])*wavefn[mid + i - 1];
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



double *gen_phi_odd(double *wavefn,
                double *V,
                double *X,
                int size,
                double phi_0,
                double phi_n1,
                double graphing_distance, 
                double E,
                double cutoff) {
    

    double del_x = graphing_distance/size;
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

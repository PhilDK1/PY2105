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
#include "gnuplot.cxx"
#include "potential_fn.cxx"

//calculate fn
double calc(double *V, 
            double *X, 
            double *wavefn,
            int size, 
            double graphing_distance, 
            double phi_0, 
            double phi_1, 
            double E_0, 
            double del_E, 
            double cutoff) {
    
    

}

double *gen_phi(double *wavefn,
               double *V,
               double *X,
               int size,
               double graphing_distance, 
               double phi_0, 
               double phi_1,
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
    int r_count = 2;
    int l_count = 1;
    wavefn[mid] = 2*phi_1 - phi_0 - 2*del_x*del_x*(E - V[mid-1])*phi_1;
    wavefn[mid + 1] = 2*wavefn[mid] - phi_1 - 2*del_x*del_x*(E - V[mid-1])*wavefn[mid];
    wavefn[mid - 1] = 2*wavefn[mid] - wavefn[mid + 1] - 2*del_x*del_x*(E - V[mid-1])*wavefn[mid];

    for (int i = 2; i < mid; i++) {
        phi_temp = 2*wavefn[mid + i - 1] - wavefn[mid + i - 2] - 2*del_x*del_x*(E - V[mid + i - 1])*wavefn[mid + i - 1];
        if (abs(phi_temp) < cutoff) {
            wavefn[mid + i] = phi_temp;
            r_count++;
        } else {
            for (int n = 0; n < (mid - i); n++)
            {
                wavefn[mid + i + n] = cutoff;
            }
            break;
            
        }
    }

    for (int i = 2; i < mid; i++) {
        phi_temp = 2*wavefn[mid - i + 1] - wavefn[mid - i + 2] - 2*del_x*del_x*(E - V[mid - i + 1])*wavefn[mid - i + 1];
        if (abs(phi_temp) < cutoff) {
            wavefn[mid - i] = phi_temp;
            l_count++;
        } else {
            for (int n = 0; n < (mid - i); n++)
            {
                wavefn[mid - i - n] = cutoff;
            }
            break;
        }
    }



    return wavefn;
}

// Main fn
int main() {
    int N = 1001;
    double V[N], X[N], wavefn[N];
    double graphing_distance = 3;
    double phi_0, phi_1 = 1;
    double E_0 = 1;
    double del_E = 0.5;
    double cutoff = 2.5;
    double lVstep = 4000;
    double rVstep = 4000;
    double VStepPoint_l = -1;
    double VstepPoint_r = 1;

    gen_v(V, X, N, lVstep, rVstep, graphing_distance, VStepPoint_l, VstepPoint_r);

    return 0;
}
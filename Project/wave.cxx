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
    
    cout << cutoff <<endl;

    double del_x = graphing_distance/size;
    cout << del_x << endl;
    double phi_temp;

    int mid;
    if (size % 2 == 0) {
        mid = (size/2);
    } else {
        mid = size -1;
        mid = (mid/2);
    }
    
    cout << mid << endl;
    // wavefn[mid] = 2*phi_1 - phi_0 - 2*del_x*del_x*(E - V[mid-1])*phi_1;
    wavefn[mid] = 2 - 1 - 2*del_x*del_x*(E - V[mid-1]);
    wavefn[mid + 1] = 2*wavefn[mid] - 1 - 2*del_x*del_x*(E - V[mid-1])*wavefn[mid];
    wavefn[mid - 1] = 2*wavefn[mid] - wavefn[mid + 1] - 2*del_x*del_x*(E - V[mid-1])*wavefn[mid];
    cout << cutoff << endl;
    for (int i = 2; i < mid; i++) {
        phi_temp = 2*wavefn[mid + i - 1] - wavefn[mid + i - 2] - 2*del_x*del_x*(E - V[mid + i - 1])*wavefn[mid + i - 1];
        cout << cutoff << endl;
        if (abs(phi_temp) < cutoff) {
            wavefn[mid + i] = phi_temp;
        } else {
            for (int n = 0; n < (mid - i); n++)
            {
                wavefn[mid + i + n] = 2.5;
            }
            break;
            
        }
        
    }

    for (int i = 2; i < mid; i++) {
        phi_temp = 2*wavefn[mid - i + 1] - wavefn[mid - i + 2] - 2*del_x*del_x*(E - V[mid - i + 1])*wavefn[mid - i + 1];
        cout << cutoff << endl;
        if (abs(phi_temp) < cutoff) {
            wavefn[mid - i] = phi_temp;
        } else {
            for (int n = 0; n < (mid - i); n++)
            {
                wavefn[mid - i - n] = 2.5;
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
    double phi_0 = 1;
    double phi_1 = 1;
    double E_0 = 1;
    double del_E = 0.5;
    double cutoff = 2.5;
    double lVstep = 4000;
    double rVstep = 4000;
    double VStepPoint_l = -1;
    double VstepPoint_r = 1;
    cout << cutoff << endl;

    gen_v(V, X, N, lVstep, rVstep, graphing_distance, VStepPoint_l, VstepPoint_r);
    gen_phi(wavefn, V, X, N, graphing_distance, phi_0, phi_1, E_0, cutoff);

    // gnuplot_two_functions("First Test", "linespoints", "X", "Potential", X, V, N, "Wavefn", X, wavefn, N, "Phi");
    gnuplot_one_function_jpg("First Test", "linespoints", "X", "wavefn", X, wavefn, N, "first_successful_test.jpg");

    for (int i = 0; i < N; i++)  {
        cout << "pot at " << X[i] << " : " << V[i] << " & wavefn : ";
        cout << wavefn[i] << "          " << i << endl;
    }


    return 0;
}
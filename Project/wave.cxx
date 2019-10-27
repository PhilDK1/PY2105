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

    if (wavefn[size - 5] < 0 && abs(wavefn[size - 5]) > cutoff) {
        last_diverge[0] = -1;
    } else {
        last_diverge[0] = 1;
    }

    return E_0, del_E, last_diverge;


}


double *gen_phi(double *wavefn,
                double *V,
                double *X,
                int size,
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
    
    // wavefn[mid] = 2*phi_1 - phi_0 - 2*del_x*del_x*(E - V[mid-1])*phi_1;
    wavefn[mid] = 2 - 1 - 2*del_x*del_x*(E - V[mid-1]);
    wavefn[mid + 1] = 2*wavefn[mid] - 1 - 2*del_x*del_x*(E - V[mid-1])*wavefn[mid];
    wavefn[mid - 1] = 2*wavefn[mid] - wavefn[mid + 1] - 2*del_x*del_x*(E - V[mid-1])*wavefn[mid];
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

    for (int i = 2; i < mid; i++) {
        phi_temp = 2*wavefn[mid - i + 1] - wavefn[mid - i + 2] - 2*del_x*del_x*(E - V[mid - i + 1])*wavefn[mid - i + 1];
        if (abs(phi_temp) < cutoff) {
            wavefn[mid - i] = phi_temp;
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



    return wavefn;
}

// Main fn
int main() {
    int N = 1001;
    double V[N], X[N], wavefn[N];
    double graphing_distance = 3;
    double E[1];
    E[0]= 1;
    double del_E[1];
    del_E[0] = 0.5;
    double cutoff = 2.5;
    double lVstep = 4000;
    double rVstep = 4000;
    double VStepPoint_l = -1;
    double VstepPoint_r = 1;
    double last_diverge[1];
    last_diverge[0] = 0;

    gen_v(V, X, N, lVstep, rVstep, graphing_distance, VStepPoint_l, VstepPoint_r);
    gen_phi(wavefn, V, X, N, graphing_distance, E[0], cutoff);

    // gnuplot_two_functions("First Test", "linespoints", "X", "Potential", X, V, N, "Wavefn", X, wavefn, N, "Phi");
    // gnuplot_one_function_jpg("First Test", "linespoints", "X", "wavefn", X, wavefn, N, "first_successful_test.jpg");

    // for (int i = 0; i < N; i++)  {
    //     cout << "pot at " << X[i] << " : " << V[i] << " & wavefn : ";
    //     cout << wavefn[i] << "          " << i << endl;
    // }
    bool cont = true;
    int count = 0;
    char prompt;
    while (cont) {
        calc(V, X, wavefn, N, last_diverge, graphing_distance, E, del_E, cutoff);
        cout << "E is: " << E[0] << endl;
        cout << "del E is: " << del_E[0] << endl;
        gen_phi(wavefn, V, X, N, graphing_distance, E[0], cutoff);
        if (count % 50 == 0) {
            gnuplot_one_function("Wave Function", "linespoints", "X", "phi", X, wavefn, N);
            cout << "Continue: ";
        
            cin >> prompt;
            switch (prompt) {
                case 'n' :
                case 'N' :
                    cont = false;
            }
        }
        count++;

    }
    char filename[30];
    char title[35];

    sprintf(filename, "plot of wavefn E=%.5lf.jpg", E[0]);
    sprintf(title, "Wave Function with Energy = %.5lf", E[0]);
    gnuplot_one_function_jpg(title, "linespoints", "X", "phi", X, wavefn, N, filename);


    return 0;
}
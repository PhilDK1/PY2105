/*
    Philip Krause
    118470776@umail.ucc.ie
    Project 9 Quantum Assignment
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
#include "Philip_Krause_Project_b.cxx"
// #include "Philip_Krause_Project_c.cxx"
#include "Philip_Krause_Project_d.cxx"



/*
const double pi = 3.14159265;
const double sq_pi = pi * pi;

int main() {
    int N = 50000;
    double V[N], X[N], wavefn[N];
    double graphing_distance = 5;
    double E[1];
    E[0]= 4;
    double del_E[1];
    del_E[0] = 0.5;
    double cutoff = 2.5;
    double Volt = 10000;
    //double rVstep = 1000;
    double L = 1;
    //double VstepPoint_r = 1;
    double last_diverge[1];
    last_diverge[0] = 0;
    double del_x = graphing_distance/N;
    double squared_function[N];
    double norm_wf[N];

    // gen_v(V, X, N, Volt, graphing_distance, L);
    mini_step(V, X, N, Volt, -2.5, 2.5, L, 0, 1);
    // gen_phi_odd(wavefn, V, X, N, 0, -del_x, graphing_distance, E[0], cutoff);
    gen_phi_odd(wavefn, V, X, N, graphing_distance, E[0], cutoff);

    bool cont = true;
    int count = 0;
    char prompt;
    while (cont) {
        calc(V, X, wavefn, N, last_diverge, graphing_distance, E, del_E, cutoff);
        cout << "E is: " << E[0] << endl;
        cout << "del E is: " << del_E[0] << endl;
        gen_phi_odd(wavefn, V, X, N, graphing_distance, E[0], cutoff);
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
    square(wavefn, squared_function, N);
    for (int i = 0; i < N; i++) {
        norm_wf[i] = wavefn[i];
    }
    double area = bound_integral(squared_function, X, -L, L, N);
    cout << "ans: " << area << endl;
    Normalise(norm_wf, area, N);
    


    char filename[35];
    char title[70];

    sprintf(filename, "plot of wavefn E=%.5lf norm.jpg", E[0]);
    sprintf(title, "Wave Function with Energy = %.5lf and normalised wavefunction", E[0]);
    gnuplot_two_functions_jpg(title, "linespoints", "X", "Wave function", X, wavefn, N, "Normalised Wavefunction", X, norm_wf, N, "Norm Wave function", filename);


    return 0;
}*/
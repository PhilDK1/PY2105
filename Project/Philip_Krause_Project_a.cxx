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
#include "Philip_Krause_Project_c.cxx"
#include "Philip_Krause_Project_d.cxx"
// #include "Philip_Krause_Project_f.cxx"


void even_squence(double energy){
    int N = 50000;
    double V[N], X[N], wavefn[N];
    double graphing_distance = 5;
    double E[1];
    E[0]= energy;
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

    gen_v(V, X, N, Volt, graphing_distance, L);
    // gen_phi_odd(wavefn, V, X, N, 0, -del_x, graphing_distance, E[0], cutoff);
    gen_phi_even(wavefn, V, X, N, graphing_distance, E[0], cutoff);

    bool cont = true;
    int count = 0;
    char prompt;
    char title[70];
    while (cont) {
        calc(V, X, wavefn, N, last_diverge, graphing_distance, E, del_E, cutoff);
        // cout << "E is: " << E[0] << endl;
        cout << "del E is: " << del_E[0] << endl;
        gen_phi_even(wavefn, V, X, N, graphing_distance, E[0], cutoff);
        if (count % 50 == 0) {
            sprintf(title, "Even Wave Function with Energy = %.5lf", E[0]);
            gnuplot_one_function(title, "linespoints", "X", "phi", X, wavefn, N);
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
    cout << "E is: " << E[0] << endl;
    cout << "\n\n\n\n" << endl;
    square(wavefn, squared_function, N);
    for (int i = 0; i < N; i++) {
        norm_wf[i] = wavefn[i];
    }
    double area = integrate(squared_function, X, L, N);
    // cout << "ans: " << area << endl;
    Normalise(norm_wf, area, N);
    


    char filename[35];
    

    sprintf(filename, "plot of Even Wave function E=%.5lf norm.jpg", E[0]);
    sprintf(title, "Even Wave Function with Energy = %.5lf and normalised wavefunction", E[0]);
    gnuplot_two_functions_jpg(title, "lines", "X", "Wave function", X, wavefn, N, "Normalised Wavefunction", X, norm_wf, N, "Norm Wave function", filename);



}


void odd_squence(double energy){
    int N = 50000;
    double V[N], X[N], wavefn[N];
    double graphing_distance = 5;
    double E[1];
    E[0]= energy;
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

    gen_v(V, X, N, Volt, graphing_distance, L);
    // gen_phi_odd(wavefn, V, X, N, 0, -del_x, graphing_distance, E[0], cutoff);
    gen_phi_odd(wavefn, V, X, N, graphing_distance, E[0], cutoff);

    bool cont = true;
    int count = 0;
    char prompt;
    char title[70];
    while (cont) {
        calc(V, X, wavefn, N, last_diverge, graphing_distance, E, del_E, cutoff);
        // cout << "E is: " << E[0] << endl;
        cout << "del E is: " << del_E[0] << endl;
        gen_phi_odd(wavefn, V, X, N, graphing_distance, E[0], cutoff);
        if (count % 50 == 0) {
            sprintf(title, "Odd Wave Function with Energy = %.5lf", E[0]);
            gnuplot_one_function(title, "linespoints", "X", "phi", X, wavefn, N);
            cout << "Continue: (Y: Yes, N:No)";
        
            cin >> prompt;
            switch (prompt) {
                case 'n' :
                case 'N' :
                    cont = false;
            }
        }
        count++;

    }
    cout << "E is: " << E[0] << endl;
    cout << "\n\n\n\n" << endl;
    square(wavefn, squared_function, N);
    for (int i = 0; i < N; i++) {
        norm_wf[i] = wavefn[i];
    }
    double area = integrate(squared_function, X, L, N);
    // cout << "ans: " << area << endl;
    Normalise(norm_wf, area, N);
    


    char filename[35];
    // char title[70];

    sprintf(filename, "plot of odd Wave function E=%.5lf norm.jpg", E[0]);
    sprintf(title, "Odd Wave Function with Energy = %.5lf and normalised wavefunction", E[0]);
    gnuplot_two_functions_jpg(title, "lines", "X", title, X, wavefn, N, "Wave function", X, norm_wf, N, "Norm Wave function", filename);
}







int main() {
    even_squence(1);
    odd_squence(4);
    even_squence(16);
    odd_squence(40);/*
    barrier_sequence(2);
    barrier_sequence(9);
    LJVsequence(-1.6, 10, 1);*/
    return 0;
}
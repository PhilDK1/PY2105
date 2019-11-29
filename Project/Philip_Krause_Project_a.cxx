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
#include "Philip_Krause_Project_f.cxx"


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



void LJVsequence(double init_E, double epsilon, double sigma){
    int N = 10000;
    double end = 5;
    // double Vstep = 4000;
    // double mVstep = 100;
    // double m_l = 0.25;
    double point = 1;
    double x_0 = 0.9;

    double V[N], X[N], wavefn[N], wavfn[N], wave_fn[N], square_fn[N], sym_V[N];
    double x_r = end - 0.005;
    double x_l = x_0 + 0.005;
    
    // double sig =1.0;
    // double ep = 10;

    LJV(V, X, x_0, end, N, sigma, epsilon);


    
    double Energy[1];
    Energy[0] = non_symadjust(V, X, wavefn, wavfn, wave_fn, N, init_E, 0.05, x_l, x_r, x_0, end);
    // cout << Energy[0] << endl;
    
    // square(wave_fn, square_fn, N);
    // double area = bound_integral(square_fn, X, x_l, x_r, N);
    // Normalise(wave_fn , area, N);
    // gnuplot_two_functions("testing", "lines", "X", "Y", X, V, N, "Lennard Jones epsilon=10, sigma =1", X, wave_fn, N, "wave function");


    char title[100];
    sprintf(title, "Wavefunction for Lennard Jones Potential E = %.5lf", Energy[0]);
    // gnuplot_one_function(title, "linespoints", "X", "phi", X, wave_fn, N);
    char file[100];
    sprintf(file, "Wavefunction_for_Lennard_Jones_Potential_E_=_%.5lf.jpg", Energy[0]);
    // gnuplot_one_function_jpg(title, "lines", "X", "Phi(x)", X, wave_fn, N, file);
    char LenJon[100];
    sprintf(LenJon, "Lennard-Jones epsilon = %.1lf, sigma = %.1lf, (scaled to fit graph)", epsilon, sigma);

    gnuplot_two_functions("Wave function for non symmetric potential", "linespoints", "X axis", "Wavefunction Amplitude", X, wave_fn, N, title , X, V, N, LenJon);
    gnuplot_one_function_jpg("Wave function for non symmetric potential", "lines", "X axis", title, X, wave_fn, N, file);
}

void barrier_sequence(double energy){
    int N = 10000;
    double end = 2;
    double Vstep = 4000;
    double mVstep = 20;
    double m_l = 0.25;
    double point = 1;
    double x_0 = -2;

    double V[N], X[N], wavefn[N], wavfn[N], wave_fn[N], square_fn[N], sym_V[N];
    double x_r = end - 0.005;
    double x_l = x_0 + 0.005;
    // cout << indexing_right(X, x_r, N);

    // double x_0 = -point - 0.1;
    // double x_f = point + 0.1;
    double sig =1.0;
    double ep = 10;

    // LJV(V, X, x_0, end, N, sig, ep);
    // mini_step(sym_V, X, S, Vstep,);
    mini_step(sym_V, X, N, Vstep, 4, 2, mVstep, 0.5);
    // cout << "yo" << endl;
    // cout << indexing_right(X, x_r, N) << endl;
    // cout << indexing_left(X, x_l, N) << endl;

    // mini_step(V, X, N, Vstep, end, point, mVstep, m_l);
    // for (int i = 0; i < N; i++) {
    //     X[i] = x_0 +i*((end -x_0)/N);
    //     V[i] = -1.61;
    // }


    // gen_phi_right(wavefn, V, X, N, -1.6, x_r, x_0, end);
    // gen_phi_left(wavfn, V, X, N, -1.6, x_l, x_0, end);
    // scale(wavefn, 0.5, N);
    // scale(wavfn, 5, N);
    double Energy[1];
    // Energy[0] = non_symadjust(V, X, wavefn, wavfn, wave_fn, N, -1.6, 0.05, x_l, x_r, x_0, end);
    // cout << Energy[0] << endl;
    Energy[0] = symadjust(sym_V, X, wavefn, wavfn, wave_fn, N, energy, 0.05, x_l, x_r, x_0, end);
    // cout << Energy[0] << endl;

    // square(wave_fn, square_fn, N);
    // double area = bound_integral(square_fn, X, x_l, x_r, N);
    // Normalise(wave_fn , area, N);


    // gnuplot_two_functions("test plot", "lines", "X", "phi", X, wavfn, N, "left start", X, wavefn, N, "right start");
    // adjust(V, X, wavefn, wavfn, N, -1.2039, 0.5, x_l, x_r, x_0, end);
    // gnuplot_two_functions("test plot", "lines", "X", "phi", X, wavfn, N, "wav", X, wavefn, N, "pot");

    // gnuplot_two_functions_jpg("Wave function for non symmetric potential", "lines", "X axis", "Wavefunction Amplitude", X, wavefn, N, "wavefunction (right)", X, wavfn, N, "wavefunction (left)", "Lennard_Jones_Potential_Wavefunction.jpg");
    char title[100];
    sprintf(title, "Wavefunction for small potential barrier E = %.5lf, potential outside = %.1lf, barrier potential = %.1lf", Energy[0], Vstep, mVstep);
    gnuplot_one_function(title, "linespoints", "X", "phi", X, wave_fn, N);
    char file[100];
    sprintf(file, "Wavefunction for small potential barrier E = %.5lf.jpg", Energy[0]);
    gnuplot_one_function_jpg(title, "lines", "X", "Phi(x)", X, wave_fn, N, file);

    // gnuplot_two_functions("sdf", "lines", "X", "phi", X, wave_fn, N,"wavefn" , X, V, N-100, "potential");
    // gnuplot_one_function("LJV", "linespoints", "x", "V", X, V,N);
}




int main() {
    even_squence(1);
    odd_squence(4);
    even_squence(16);
    odd_squence(40);
    barrier_sequence(2);
    barrier_sequence(9);
    LJVsequence(-1.6, 10, 1);
    return 0;
}
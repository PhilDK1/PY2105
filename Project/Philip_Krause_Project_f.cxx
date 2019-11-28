/*
    Philip Krause
    118470776@umail.ucc.ie
    Project 9 Quantum Assignment

    File to impliment the matching methond as per section 10.2 of computational physics 2nd edition Giordano, N.J.
    and test for lennard Jones potential

*/

// Header Files
#include <iostream>
#include <fstream>
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

#include "Philip_Krause_Project_d.cxx"
// #include "Philip_Krause_Project_c.cxx"
// #include "gnuplot.cxx"
#include "Philip_Krause_Project_b.cxx"

// int mid(N);
/*
int main() {
    int N = 10000;
    double end = 2;
    double Vstep = 4000;
    double mVstep = 100;
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

    LJV(V, X, x_0, end, N, sig, ep);
    mini_step(sym_V, X, N, 1000, 4, 2, 20, 0.5);
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

    // non_symadjust(V, X, wavefn, wavfn, wave_fn, N, -1.6, 0.05, x_l, x_r, x_0, end);
    symadjust(sym_V, X, wavefn, wavfn, wave_fn, N, -2, 0.05, x_l, x_r, x_0, end);
    square(wave_fn, square_fn, N);
    double area = bound_integral(square_fn, X, x_l, x_r, N);
    Normalise(wave_fn , area, N);


    gnuplot_two_functions("test plot", "lines", "X", "phi", X, wavfn, N, "left start", X, wavefn, N, "right start");
    // adjust(V, X, wavefn, wavfn, N, -1.2039, 0.5, x_l, x_r, x_0, end);
    // gnuplot_two_functions("test plot", "lines", "X", "phi", X, wavfn, N, "wav", X, wavefn, N, "pot");

    // gnuplot_two_functions_jpg("Wave function for non symmetric potential", "lines", "X axis", "Wavefunction Amplitude", X, wavefn, N, "wavefunction (right)", X, wavfn, N, "wavefunction (left)", "Lennard_Jones_Potential_Wavefunction.jpg");
    gnuplot_one_function("test", "lines", "X", "phi", X, wave_fn, N);


    return 0;
}*/

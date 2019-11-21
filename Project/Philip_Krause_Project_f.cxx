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
#include "Philip_Krause_Project_e.cxx"
#include "gnuplot.cxx"
#include "Philip_Krause_Project_b.cxx"

int indexing_left(double *X, double x_val, int size) {
    for (int i = 0; i < size ; i++) {
        if (x_val < X[i]) {
            return i-1;
        }    
    }
}

int indexing_right(double *X, double x_val, int size) {
    for (int i = 0; i < size ; i++) {
        if (x_val > X[size - 1 -i]) {
            return size - i;
        }    
    }
}

double *gen_phi_left(double *wavefn, double *V, double *X, int N, double E, double x_l, double start, double end){
    int ind = indexing_left(X, x_l, N);
    double del_x = (end - start)/N;
    wavefn[ind - 2] = -0.0001*del_x;
    wavefn[ind -1] = 0;

    for (int i = 0; i < N- ind; i++) {
        wavefn[ind +i] = 2*wavefn[ind +i -1] -wavefn[ind + i -2] - 2*del_x*del_x*(E- V[ind +i -1])*wavefn[ind +i -1];
        // cout << i << "      " << wavefn[ind +i] << "    " << (wavefn[ind +i] - wavefn[ind +i -1]) << "      " << V[ind +i -1] << endl;
    }


    return wavefn;
}


double *gen_phi_right(double *wavefn, double *V, double *X, int N, double E, double x_r, double start, double end) {
    int ind = indexing_right(X, x_r, N);
    double del_x = (end - start)/N;
    wavefn[ind + 2] = -0.0001*del_x;
    wavefn[ind + 1] = 0;

    for (int i = 0; i <= ind; i++) {
        wavefn[ind-i] = 2*wavefn[ind-(i -1)] - wavefn[ind-(i-2)] -2*del_x*del_x*(E-V[ind-(i-1)])*wavefn[ind-(i-1)];
    }

    return wavefn;

}




int locate_min(double *arr, int N){
    int min = 0;
    double min_val = arr[0];
    for (int i = 0; i < N; i++) {
        if (arr[i] <= min_val) {
            min = i;
            min_val = arr[i];
        }
    }
    return min;
}


double *scale(double *wavefn, double factor, int N){
    for (int i = 0; i < N; i++){
        wavefn[i] = factor*wavefn[i];
    }
    return wavefn;
}

double finite_derive(double *Y, double *X, int index) {
    double del_x = abs(X[index] -  X[index + 1]);
    double deriv = (Y[index+1]- Y[index])/del_x;
    return deriv;
}

double *adjust(double *V, double *X, double *wavefn_1, double *wavefn_2, double *wave_function, int N, double init_E, double inc_E, double x_l, double x_r, double start, double end) {
    // find minium potential of on the range
    int minimum = locate_min(V, N);
    double tolerance = ((end - start)/N)*((end - start)/N)*((end - start)/N)*((end - start)/N);
    cout << tolerance << endl;
    gen_phi_left(wavefn_1, V, X, N, init_E, x_l, start, end);
    gen_phi_right(wavefn_2, V, X, N, init_E, x_r, start, end);
    // gnuplot_two_functions("test plot", "lines", "X", "phi", X, wavefn_1, N, "left", X, wavefn_2, N, "right");
    double E = init_E;
    double ratio = (wavefn_1[minimum]/wavefn_2[minimum]);
    cout << ratio<< endl;
    scale(wavefn_2, ratio, N);
    gnuplot_two_functions("test plot", "lines", "X", "phi", X, wavefn_1, N, "left", X, wavefn_2, N, "right");
    int match = 1;
    double diff = finite_derive(wavefn_1, X, minimum) - finite_derive(wavefn_2, X, minimum); // >0 then d_phi_l/dx > d_phi_r/dx else l<r
    cout << diff << endl;

    int count = 1;
    int ans;
    if (diff >0 ) {
        match = 1;

    }else if (diff <0) {
        match = -1;
    }

    while (abs(diff) > tolerance) {
        cout << count << endl;
        if ((match == 1) && (diff >0)){
            // energy needs to increase
            // cout << "1   1  " << E << endl;
            E = E + inc_E;
            // cout << E << endl;
        } else if ((match == 1) && (diff <0)) {
            // increased last time so needs to decrease
            inc_E  = inc_E/2;
            // cout << "1  -1  " << E << endl;
            E = E - inc_E;
            // cout << E << endl;
            match = -1;
        } else if ((match == -1) && (diff >0)) {
            // decreased last time needs to increase
            inc_E  = inc_E/2;
            // cout << "-1  1  " << E << endl;
            E = E + inc_E;
            // cout << E << endl;
            match = 1;
        } else if ((match == -1) && (diff <0)) {
            // energy needs to decrease
            // cout << "-1 -1  " << E << endl;
            E = E - inc_E;
            // cout << E << endl;
        }
        gen_phi_left(wavefn_1, V, X, N, E, x_l, start, end);
        gen_phi_right(wavefn_2, V, X, N, E, x_r, start, end);
        ratio = (wavefn_1[minimum]/wavefn_2[minimum]);
        scale(wavefn_2, ratio, N);
        diff = finite_derive(wavefn_1, X, minimum) - finite_derive(wavefn_2, X, minimum); // >0 then d_phi_l/dx > d_phi_r/dx else l<r
        count++;
        // gnuplot_two_functions("test plot", "lines", "X", "phi", X, wavefn_1, N, "left", X, wavefn_2, N, "right");
        // cin >> ans;
    }
    for (int i = 0; i < N; i++) {
        if (i <= minimum) {
            wave_function[i] = wavefn_1[i];
        } else {
            wave_function[i] = wavefn_2[i];
        }
    }

}

int main() {
    int N = 1000;
    double end = 3;
    double Vstep = 4000;
    double mVstep = 100;
    double m_l = 0.25;
    double point = 1;
    double x_0 = 0.9;

    double V[N], X[N], wavefn[N], wavfn[N], wave_fn[N];
    double x_r = end - 0.05;
    double x_l = x_0 + 0.05;
    // cout << indexing_right(X, x_r, N);

    // double x_0 = -point - 0.1;
    // double x_f = point + 0.1;
    double sig =1.0;
    double ep = 250;

    LJV(V, X, x_0, end, N, sig, ep);
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

    adjust(V, X, wavefn, wavfn, wave_fn, N, -1.6, 0.05, x_l, x_r, x_0, end);


    gnuplot_two_functions("test plot", "lines", "X", "phi", X, wavfn, N, "left start", X, wavefn, N, "right start");
    // adjust(V, X, wavefn, wavfn, N, -1.2039, 0.5, x_l, x_r, x_0, end);
    // gnuplot_two_functions("test plot", "lines", "X", "phi", X, wavfn, N, "wav", X, wavefn, N, "pot");

    // gnuplot_two_functions_jpg("Wave function for non symmetric potential", "lines", "X axis", "Wavefunction Amplitude", X, wavefn, N, "wavefunction (right)", X, wavfn, N, "wavefunction (left)", "Lennard_Jones_Potential_Wavefunction.jpg");
    gnuplot_one_function("test", "lines", "X", "phi", X, wave_fn, N);


    return 0;
}
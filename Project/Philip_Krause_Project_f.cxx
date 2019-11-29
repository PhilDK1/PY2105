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
#include "Philip_Krause_Project_d.cxx"
#include "gnuplot.cxx"
#include "Philip_Krause_Project_b.cxx"

// int mid(N);

double *LJV(double *V, double *X, double starting_x_nz, double endpoint, int size, double sigma, double epsilon) {    
    // function to calculate the Lenard-Jones Potential for a range of points

    /*
    L-J potential given by the following formula
    V(x) = 4e( (sigma/x)^(12) - (sigma/x)^(6) )

    */

    // decalre del_x is the spacial step such that you have "size" number of steps between endpoint and starting_x_nz
    double del_x = (endpoint - starting_x_nz)/size; 
    // cout << "hi" << endl;

    // x is the current position i.e. x = x_0 +i*del_x kind of thing
    double x  = starting_x_nz;

    // loop through each point from starting_x_nz
    for (int i = 0; i < size; i++) {
        // assign the current calculated the position value, x = (x_0 + i*del_x)
        X[i] = x;
        // calculate the potential at the point x
        V[i] = 4*epsilon*(pow((sigma/x),12)  -pow((sigma/x), 6));
        // calculate x for the next iteration
        x += del_x;
    }
    // return the potential function
    return V;
}

double *gen_phi_left(double *wavefn, double *V, double *X, int N, double E, double x_l, double start, double end){
    int ind = indexing_left(X, x_l, N);
        double offset =abs(start - x_l);
    while (ind < 2) {
        x_l += 2*offset;
        ind = indexing_left(X, x_l, N);
    }

    double del_x = (end - start)/N;
    // cout << del_x << endl;
    // cout << abs(del_x) << endl;
    wavefn[ind - 2] = -0.0001*abs(del_x);
    wavefn[ind -1] = 0;

    for (int i = 0; i < N- ind; i++) {
        wavefn[ind +i] = 2*wavefn[ind +i -1] -wavefn[ind + i -2] - 2*del_x*del_x*(E- V[ind +i -1])*wavefn[ind +i -1];
        // cout << i << "      " << wavefn[ind +i] << "    " << (wavefn[ind +i] - wavefn[ind +i -1]) << "      " << V[ind +i -1] << endl;
    }


    return wavefn;
}


double *gen_phi_right(double *wavefn, double *V, double *X, int N, double E, double x_r, double start, double end) {
    int ind = indexing_right(X, x_r, N);
    // cout << ind << endl;
    while (ind < 2) {
        x_r = 3*x_r;
        ind = indexing_right(X, x_r, N);
    }
    double del_x = (end - start)/N;
    // cout << del_x << endl;
    // cout << abs(del_x) << endl;
    wavefn[ind + 2] = -0.0001*abs(del_x);
    wavefn[ind + 1] = 0;

    for (int i = 0; i <= ind; i++) {
        wavefn[ind-i] = 2*wavefn[ind-(i -1)] - wavefn[ind-(i-2)] -2*del_x*del_x*(E-V[ind-(i-1)])*wavefn[ind-(i-1)];
    }

    return wavefn;

}

double non_symadjust(double *V, double *X, double *wavefn_1, double *wavefn_2, double *wave_function, int N, double init_E, double inc_E, double x_l, double x_r, double start, double end) {
    // find minium potential of on the range
    int minimum = locate_min(V, N);
    double tolerance = ((end - start)/N)*((end - start)/N)*((end - start)/N)*((end - start)/N);
    double square_fn[N];
    // cout << tolerance << endl;
    gen_phi_left(wavefn_1, V, X, N, init_E, x_l, start, end);
    // cout << "left" << endl;
    gen_phi_right(wavefn_2, V, X, N, init_E, x_r, start, end);
    // gnuplot_two_functions("test plot", "lines", "X", "phi", X, wavefn_1, N, "left", X, wavefn_2, N, "right");
    double E = init_E;
    double ratio = (wavefn_1[minimum]/wavefn_2[minimum]);
    // cout << ratio<< endl;
    scale(wavefn_2, ratio, N);
    // gnuplot_two_functions("Wavefunction for Lennard-Jones Potential", "lines", "X", "phi", X, wavefn_1, N, "left", X, wavefn_2, N, "right");
    int match = 1;
    double diff = finite_derive(wavefn_1, X, minimum) - finite_derive(wavefn_2, X, minimum); // >0 then d_phi_l/dx > d_phi_r/dx else l<r
    // cout << diff << endl;

    int count = 1;
    int ans;
    if (diff >0 ) {
        match = 1;

    }else if (diff <0) {
        match = -1;
    }

    while (abs(diff) > tolerance) {
        // cout << count << endl;
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
    // cout << E<< endl;
    square(wave_function, square_fn, N);
    double area = bound_integral(square_fn, X, x_l, x_r, N);
    Normalise(wave_function , area, N);

    int maximum = locate_max(wave_function, N);
    // cout << wave_function[maximum] << endl;
    // cout <<V[2]<< endl;
    double scale_ratio = ((2*wave_function[maximum])/V[2]);
    // cout << scale_ratio << endl;
    scale(V, scale_ratio, N);
    // gnuplot_one_function("scaled V", "lines", "X", "V", X, V, N);
    return E;
}

double symadjust(double *V, double *X, double *wavefn_1, double *wavefn_2, double *wave_function, int N, double init_E, double inc_E, double x_l, double x_r, double start, double end) {
    // find minium potential of on the range
    // int minimum = locate_min(V, N);
    int mid = N/2;
    double tolerance = ((end - start)/N)*((end - start)/N)*((end - start)/N)*((end - start)/N);
    // cout << tolerance << endl;
    gen_phi_left(wavefn_1, V, X, N, init_E, x_l, start, end);
    // cout << "left" << endl;
    gen_phi_right(wavefn_2, V, X, N, init_E, x_r, start, end);
    // cout << "yo" << endl;
    // gnuplot_two_functions("test plot", "lines", "X", "phi", X, wavefn_1, N, "left", X, wavefn_2, N, "right");
    double E = init_E;
    double ratio = (wavefn_1[mid]/wavefn_2[mid]);
    // cout << ratio<< endl;
    scale(wavefn_2, ratio, N);
    // cout << "hey" << endl;
    // gnuplot_two_functions("test plot", "lines", "X", "phi", X, wavefn_1, N, "left", X, wavefn_2, N, "right");
    int match = 1;
    double diff = finite_derive(wavefn_1, X, mid) - finite_derive(wavefn_2, X, mid); // >0 then d_phi_l/dx > d_phi_r/dx else l<r
    // cout << diff << endl;

    int count = 1;
    int ans;
    if (diff >0 ) {
        match = 1;

    }else if (diff <0) {
        match = -1;
    }

    while (abs(diff) > tolerance) {
        // cout << count << endl;
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
        ratio = (wavefn_1[mid]/wavefn_2[mid]);
        scale(wavefn_2, ratio, N);
        diff = finite_derive(wavefn_1, X, mid) - finite_derive(wavefn_2, X, mid); // >0 then d_phi_l/dx > d_phi_r/dx else l<r
        count++;
        // cout << diff << endl;
        // gnuplot_two_functions("test plot", "lines", "X", "phi", X, wavefn_1, N, "left", X, wavefn_2, N, "right");
        // cin >> ans;
    }
    for (int i = 0; i < N; i++) {
        if (i <= mid) {
            wave_function[i] = wavefn_1[i];
        } else {
            wave_function[i] = wavefn_2[i];
        }
    }
    double square_fn[N];
    square(wave_function, square_fn, N);
    double area = bound_integral(square_fn, X, x_l, x_r, N);
    Normalise(wave_function , area, N);

    // int maximum = locate_max(wave_function, N);

    return E;
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





/*
int main() {
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
    Energy[0] = symadjust(sym_V, X, wavefn, wavfn, wave_fn, N, 9, 0.05, x_l, x_r, x_0, end);
    // cout << Energy[0] << endl;

    // square(wave_fn, square_fn, N);
    // double area = bound_integral(square_fn, X, x_l, x_r, N);
    // Normalise(wave_fn , area, N);


    // gnuplot_two_functions("test plot", "lines", "X", "phi", X, wavfn, N, "left start", X, wavefn, N, "right start");
    // adjust(V, X, wavefn, wavfn, N, -1.2039, 0.5, x_l, x_r, x_0, end);
    // gnuplot_two_functions("test plot", "lines", "X", "phi", X, wavfn, N, "wav", X, wavefn, N, "pot");

    // gnuplot_two_functions_jpg("Wave function for non symmetric potential", "lines", "X axis", "Wavefunction Amplitude", X, wavefn, N, "wavefunction (right)", X, wavfn, N, "wavefunction (left)", "Lennard_Jones_Potential_Wavefunction.jpg");
    char title[100];
    sprintf(title, "Wavefunction for small potential barrier E = %.5lf", Energy[0]);
    gnuplot_one_function(title, "lines", "X", "phi", X, wave_fn, N);
    char file[100];
    sprintf(file, "Wavefunction for small potential barrier E = %.5lf.jpg", Energy[0]);
    gnuplot_one_function_jpg(title, "lines", "X", "Phi(x)", X, wave_fn, N, file);

    // gnuplot_two_functions("sdf", "lines", "X", "phi", X, wave_fn, N,"wavefn" , X, V, N-100, "potential");
    // gnuplot_one_function("LJV", "linespoints", "x", "V", X, V,N);




    return 0;
}
*/


int main() {
    barrier_sequence(2);
    barrier_sequence(9);
    LJVsequence(-1.6, 10, 1);
    return 0;
}

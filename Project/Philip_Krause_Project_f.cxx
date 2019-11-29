/*
    Philip Krause
    118470776@umail.ucc.ie
    Project 9 Quantum Assignment

    File to impliment the matching methond as per section 10.2 of computational physics 2nd edition Giordano, N.J.
    and test for lennard Jones potential and for a square well with a small potential barrier

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
// #include "Philip_Krause_Project_d.cxx"
// #include "gnuplot.cxx"
// #include "Philip_Krause_Project_b.cxx"

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

// funciton to generate the wave funciton from the left hand side
double *gen_phi_left(double *wavefn, double *V, double *X, int N, double E, double x_l, double start, double end){
    // get the index of the offset on the left
    // -- code --> Philip_Krause_Project_d.cxx
    int ind = indexing_left(X, x_l, N);
    // get offset
    double offset =abs(start - x_l);
    // if the offset isn't big enough increase it an repeat above^^ until it is
    while (ind < 2) {
        // increase offset
        x_l += 2*offset;
        // get the index of the offset on the left
        // -- code --> Philip_Krause_Project_d.cxx
        ind = indexing_left(X, x_l, N);
    }

    // calculate the spacial step
    double del_x = (end - start)/N;

    // calculate initial conditions
    wavefn[ind - 2] = -0.0001*abs(del_x);
    wavefn[ind -1] = 0;

    // loop through values left to right and calculate the wave function
    for (int i = 0; i < N- ind; i++) {
        wavefn[ind +i] = 2*wavefn[ind +i -1] -wavefn[ind + i -2] - 2*del_x*del_x*(E- V[ind +i -1])*wavefn[ind +i -1];
    }

    // return the wave function
    return wavefn;
}

// funciton to generate the wave funciton from the right hand side
double *gen_phi_right(double *wavefn, double *V, double *X, int N, double E, double x_r, double start, double end) {
    // get the index of the offset on the left
    // -- code --> Philip_Krause_Project_d.cxx
    int ind = indexing_right(X, x_r, N);
    // get offset
    double offset =abs(start - x_r);
    // if the offset isn't big enough increase it an repeat above^^ until it is
    while (ind < 2) {
        x_r += 2*offset;
        // get the index of the offset on the right
        // -- code --> Philip_Krause_Project_d.cxx
        ind = indexing_right(X, x_r, N);
    }
    
    // calculate the spacial step
    double del_x = (end - start)/N;
    // calculate initial conditions
    wavefn[ind + 2] = -0.0001*abs(del_x);
    wavefn[ind + 1] = 0;

    // loop through values right to left and calculate the wave function
    for (int i = 0; i <= ind; i++) {
        wavefn[ind-i] = 2*wavefn[ind-(i -1)] - wavefn[ind-(i-2)] -2*del_x*del_x*(E-V[ind-(i-1)])*wavefn[ind-(i-1)];
    }
    // return the wavefunction
    return wavefn;

}

// routine to adjust the energy of the wave function calculated from right and left then normalise the wwave function
double non_symadjust(double *V, double *X, double *wavefn_1, double *wavefn_2, double *wave_function, int N, double init_E, double inc_E, double x_l, double x_r, double start, double end) {
    // find minium potential of on the range
    int minimum = locate_min(V, N);
    // set the tolerance as the spacial step to the 4th power
    double tolerance = ((end - start)/N)*((end - start)/N)*((end - start)/N)*((end - start)/N);
    // initialise the array to be used for the square of the wavwfunction
    double square_fn[N];
    // generate the wavefunction from the left and the right side
    gen_phi_left(wavefn_1, V, X, N, init_E, x_l, start, end);
    gen_phi_right(wavefn_2, V, X, N, init_E, x_r, start, end);
    // transfer energy to a new variable
    double E = init_E;
    // get the ratio of heights between the two wave functions when the potential is minimum
    double ratio = (wavefn_1[minimum]/wavefn_2[minimum]);
    // scale one of the wave functions so they are equal when potential is minimum
    scale(wavefn_2, ratio, N);
    // initialise variable to check difference in derivaties
    int match = 1;
    // calculate the difference in derivatives
    double diff = finite_derive(wavefn_1, X, minimum) - finite_derive(wavefn_2, X, minimum); // >0 then d_phi_l/dx > d_phi_r/dx else l<r
    // initialise count for loop
    // check if the differnce is positive or negative
    if (diff >0 ) {
        match = 1;

    }else if (diff <0) {
        match = -1;
    }

    // check if the difference in the derivatives of the functions is within the tolerance and adjust energy
    while (abs(diff) > tolerance) {
        if ((match == 1) && (diff >0)){
            // energy needs to increase
            E = E + inc_E;
            // cout << E << endl;
        } else if ((match == 1) && (diff <0)) {
            // increased last time so needs to decrease
            inc_E  = inc_E/2;
            E = E - inc_E;
            match = -1;
        } else if ((match == -1) && (diff >0)) {
            // decreased last time needs to increase
            inc_E  = inc_E/2;
            E = E + inc_E;
            match = 1;
        } else if ((match == -1) && (diff <0)) {
            // energy needs to decrease
            E = E - inc_E;
        }
        // generate the wave functions
        gen_phi_left(wavefn_1, V, X, N, E, x_l, start, end);
        gen_phi_right(wavefn_2, V, X, N, E, x_r, start, end);
        // get the ratio of heights between the two wave functions when the potential is minimum
        ratio = (wavefn_1[minimum]/wavefn_2[minimum]);
        // scale so they are equal at a point
        scale(wavefn_2, ratio, N);
        // calculate the difference in the derivatives
        diff = finite_derive(wavefn_1, X, minimum) - finite_derive(wavefn_2, X, minimum); // >0 then d_phi_l/dx > d_phi_r/dx else l<r
    }
    // loop throught the two wave functions and combine them into one continuous function
    for (int i = 0; i < N; i++) {
        // left wave function while less than the min potential point
        if (i <= minimum) {
            wave_function[i] = wavefn_1[i];
        // right wave function while greater than the min potential point
        } else {
            wave_function[i] = wavefn_2[i];
        }
    }
    // square teh wave function
    square(wave_function, square_fn, N);
    // get the area of the squared wave function
    double area = bound_integral(square_fn, X, x_l, x_r, N);
    // normalise the wave function
    Normalise(wave_function , area, N);
    // get the maximum point of the wave funciton
    int maximum = locate_max(wave_function, N);
    // get the ratio between one of the higher potential points
    double scale_ratio = ((2*wave_function[maximum])/V[2]);
    // scale the potential function so that both the wave function and the potential funciton will fit noiely on one graph
    scale(V, scale_ratio, N);
    // return energy
    return E;
}


double symadjust(double *V, double *X, double *wavefn_1, double *wavefn_2, double *wave_function, int N, double init_E, double inc_E, double x_l, double x_r, double start, double end) {
    // find minium potential of on the range
    int mid = N/2;
    double tolerance = ((end - start)/N)*((end - start)/N)*((end - start)/N)*((end - start)/N);
    gen_phi_left(wavefn_1, V, X, N, init_E, x_l, start, end);
    gen_phi_right(wavefn_2, V, X, N, init_E, x_r, start, end);
    double E = init_E;
    double ratio = (wavefn_1[mid]/wavefn_2[mid]);
    scale(wavefn_2, ratio, N);
    int match = 1;
    double diff = finite_derive(wavefn_1, X, mid) - finite_derive(wavefn_2, X, mid); // >0 then d_phi_l/dx > d_phi_r/dx else l<r

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
            E = E + inc_E;
        } else if ((match == 1) && (diff <0)) {
            // increased last time so needs to decrease
            inc_E  = inc_E/2;
            E = E - inc_E;
            match = -1;
        } else if ((match == -1) && (diff >0)) {
            // decreased last time needs to increase
            inc_E  = inc_E/2;
            E = E + inc_E;
            match = 1;
        } else if ((match == -1) && (diff <0)) {
            // energy needs to decrease
            E = E - inc_E;
        }
        gen_phi_left(wavefn_1, V, X, N, E, x_l, start, end);
        gen_phi_right(wavefn_2, V, X, N, E, x_r, start, end);
        ratio = (wavefn_1[mid]/wavefn_2[mid]);
        scale(wavefn_2, ratio, N);
        diff = finite_derive(wavefn_1, X, mid) - finite_derive(wavefn_2, X, mid); // >0 then d_phi_l/dx > d_phi_r/dx else l<r
        count++;
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


    return E;
}

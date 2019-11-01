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
#include "potential_fn.cxx"
#include "wave.cxx"

const double pi = 3.14159265;
const double sq_pi = pi * pi;



//Main files
// int main() {

//     //define variables
//     // int N = 500;
//     // int mid_point = (N-1)/2;

//     // double phi[N];

//     // double E =sq_pi/8;
//     double V[N], X[N];
//     // double L = 2;
//     // double three = 3;

//     // double phi_0 = 1.0;
//     // double phi_n1 = 1.0;

//     // double offset = 1;
//     // double V_at_L = 100000;

//     // // double graphing_distance = 2*(L+offset);
//     // double del_x = three/N;
//     // double x = 0;

// /*
//     for (int i = 0; i < N; i++) {
//         X[i] = i*del_x;
//         if (X[i] > L){
//             V[i] = V_at_L;  
//         } else {
//             V[i] = 0.0;
//         }

//     }
//     */
//     phi[0] = phi_n1;
//     phi[1] = phi_0;
//     double phi_temp;
//     int count = 2;
//     for (int i = 2; i < N; i++) {
//         phi_temp = 2*phi[i-1] - phi[i-2] - 2*del_x*del_x*(E-V[i-1])*phi[i-1];
//         if (abs(phi_temp) < 1.7) {
//             phi[i] = phi_temp;
//             count++;
//         } else {
//             break;
//         }
//     }
//     cout << count << endl;
//     double mewX[count];

//     for (int i = 0; i < count; i++ ) {
//         mewX[i] = X[i];
//     }
  

// /*
//     for (int i = 0; i < N; i++) {

//     }
// */

//      //gnuplot_one_function("Test of potential array", "linespoints", "x", "V", X, V, N);
//     gnuplot_one_function("Test of wave fn array", "linespoints", "x", "phi", mewX, phi, count);

//     return 0;
// }

// Main function 2.0


// Main fn
int main() {
    int N = 50000;
    double V[N], X[N], wavefn[N];
    double graphing_distance = 5;
    double E[1];
    E[0]= 11;
    double del_E[1];
    del_E[0] = 0.5;
    double cutoff = 2.5;
    double lVstep = 1000;
    double rVstep = 1000;
    double VStepPoint_l = -1;
    double VstepPoint_r = 1;
    double last_diverge[1];
    last_diverge[0] = 0;
    double del_x = graphing_distance/N;

    gen_v(V, X, N, lVstep, rVstep, graphing_distance, VStepPoint_l, VstepPoint_r);
    gen_phi_even(wavefn, V, X, N, 1, 1, graphing_distance, E[0], cutoff);

    bool cont = true;
    int count = 0;
    char prompt;
    while (cont) {
        calc(V, X, wavefn, N, last_diverge, graphing_distance, E, del_E, cutoff);
        cout << "E is: " << E[0] << endl;
        cout << "del E is: " << del_E[0] << endl;
        gen_phi_even(wavefn, V, X, N, 1, 1, graphing_distance, E[0], cutoff);
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
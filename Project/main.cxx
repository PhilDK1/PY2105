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

int main() {
    int N = 1001;
    double V_at_step_l = 4000;
    double V_at_step_r = 4000;
    double point = 1;
    double distance = 3;
    double V[N], X[N];

    gen_v(V, X, N, V_at_step_l, V_at_step_r, distance, -point, point);
    gnuplot_one_function("Test of potential generation", "linespoints", "x", "V", X, V, N);

    return 0;
}
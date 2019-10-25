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

const double pi = 3.14159265;
const double sq_pi = pi * pi;

//Main files
int main() {

    //define variables
    int N = 1001;
    int mid_point = (N-1)/2;

    double phi[N];

    double E =sq_pi/8;
    double V[N];
    double L = 1;

    double phi_0 = 1.0;
    double phi_n1 = 1.0;

    double offset = 1;
    double V_at_L = 100;

    double graphing_distance = 2*(L+offset);
    double del_x = graphing_distance/N;

    double x = -0.5*graphing_distance;
    double X[N];
    for (int i = 0; i < N; i++) {
        X[i] = x;
        if (x < -1*L) {
            V[i] = V_at_L;
        } else if (x > L){
            V[i] = V_at_L;
            
        } else {
            V[i] = 0.0;
        }

        x += del_x;
    }
    phi[mid_point] = phi_0;
    phi[mid_point - 1] = phi_n1;
    for (int i = mid_point; i < N; i++) {

    }
    
    // gnuplot_one_function_jpg("Test of potential array", "linespoints", "x", "V", X, V, N, "yoy.jpg" );

    return 0;
}
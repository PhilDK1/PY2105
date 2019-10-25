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

//Main files
int main() {

    //define variables
    int N = 1000;

    double phi[N];

    double E;
    double V[N];
    double L = 1;

    double phi_0;
    double phi_1;

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
    cout << abs(-x) << endl;
    cout << abs(x) << endl;
    cout << -x << endl;
    cout << x << endl;
    
    gnuplot_one_function_jpg("Test of potential array", "linespoints", "x", "V", X, V, N, "yoy.jpg" );

    return 0;
}
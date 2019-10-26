/*
    Philip Krause
    118470776@umail.ucc.ie
    Project 9 Quantum Assignment
    testing of various functions
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


double* gen_v(  double* V,
                double* X, 
                int size, 
                double lV_step, 
                double rV_step, 
                double graphing_distance, 
                double l_step_at, 
                double r_step_at) {
    
    if ((abs(l_step_at) + abs(r_step_at)) > graphing_distance) cout << "check points not all in range" << endl;

    double del_x = graphing_distance/size;
    double starting = -0.5*graphing_distance;
    
    for (int i = 0; i < size; i++) {
        X[i] = starting + (i*del_x);
        if (X[i] > r_step_at){
            V[i] = rV_step;  
        } else if (X[i] < l_step_at) { 
            V[i] = lV_step;
        } else {
            V[i] = 0.0;
        }

    }

    return V, X;

}
/*
int main() {

    int N = 1001;
    double V_at_step_l = 4000;
    double V_at_step_r = 4000;
    double point = 1;
    double distance = 3;
    double del_x = distance/N;
    double X[N], V[N], Index[N];

    gen_v(V, X, N, V_at_step_l, V_at_step_r, distance, -point, point);
/*
    for (int i = 0; i < N; i++) {
        X[i] = displacement[i];
        V[i] = potential[i];
    }



    //for (int i = 0;)

    //gnuplot_one_function("Test of potential generation", "linespoints", "x", "V", X, V, N);
}
*/

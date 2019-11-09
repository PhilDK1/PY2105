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



// double * gen_test(double *Y, int size) {
//     for (int i = 0; i < size; i++){
//         Y[i] = i*i + 2i + 1
//     }
// }

double* square(
    double *wavefn,
    double *squared_wf,
    int size) {
    

    for (int i = 0; i < size; i++ ) {
        squared_wf[i]= wavefn[i]*wavefn[i];
    }

    return squared_wf;
}

double integrate(
        double *sq_wf,
        double *X,
        double left_L,
        double right_L,
        int size,
        double del_x) {
    int count = 0;
    int upper = 0;
    // double Sum = 0;
    while (X[count] < left_L) {
        count++;
    }
    count++;

    while (X[upper] < right_L) {
        upper++;
    }
    upper--;

    double Sum = 0;
    double ans;
    Sum += sq_wf[count]+ sq_wf[upper];

    for (int i = count + 1; i < upper; i++) {
        if (i % 3 != 0) {
            Sum += 3*(sq_wf[i]);
        }

    }
    
    for (int i = count + 1; i < (upper/3); i++) {
        Sum += 2*(sq_wf[3*i]);
    }

    ans = ((3*del_x)/(8))*Sum;


}

double *Normalise(double *wavefn, double area, double size) {
    for (int i = 0; i < size; i++) {
        wavefn[i] = (1/(sqrt(area)))*wavefn[i];
    }
    return wavefn;
}


/*
int main() {
    int N = 50000;
    double V[N], X[N], wavefn[N], dupl[N];
    double graphing_distance = 5;
    double E[1];
    E[0]= 4;
    double del_E[1];
    del_E[0] = 1;
    double cutoff = 2.5;
    double lVstep = 1000;
    double rVstep = 1000;
    double VStepPoint_l = -1;
    double VstepPoint_r = 1;
    double last_diverge[1];
    last_diverge[0] = 0;
    double del_x = graphing_distance/N;
    double squared_function[N];
    double ans;

    gen_v(V, X, N, lVstep, rVstep, graphing_distance, VStepPoint_l, VstepPoint_r);
    gen_phi_odd(wavefn, V, X, N, 0, -del_x, graphing_distance, E[0], cutoff);
    gen_phi_odd(dupl, V, X, N, 0, -del_x, graphing_distance, E[0], cutoff);

    square(wavefn, squared_function, N);
    ans = integrate(squared_function, X, VStepPoint_l, VstepPoint_r, N, del_x);

    cout << ans << endl;

    Normalise(dupl, ans, N);
    square(dupl, squared_function, N);
    ans = integrate(squared_function, X, VStepPoint_l, VstepPoint_r, N, del_x);

    cout << ans << endl;
    gnuplot_two_functions("wf & sqwf", "linespoints", "X", "value", X, wavefn, N, "yurt", X, dupl, N, "norm");




    return 0;
}
*/
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


double* square(
    double *wavefn,
    double *squared_wf,
    int size) {

    for (int i = 0; i++ ; i < size) {
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
    while (count < size || X[count] < left_L) {
        count++;
    }
    count++;

    while (upper < size || X[upper] < right_L) {
        upper++;
    }

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
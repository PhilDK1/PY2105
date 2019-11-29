/*
    Philip Krause
    118470776@umail.ucc.ie
    Project 9 Quantum Assignment
    testing for wave fn

    File to contain maths functions used for various things
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




// function to calculate the square of the wavfunction
double* square(
    double *wavefn,
    double *squared_wf,
    int size) {
    
    // iterates the between 0 and size and multiples that value of the wavfn 
    // by itself and assigns it to a new value in an array
    for (int i = 0; i < size; i++ ) {
        squared_wf[i]= wavefn[i]*wavefn[i];
    }
    // returns the array of the wave function squared
    return squared_wf;
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

int locate_max(double *arr, int N){
    int max = 0;
    double max_val = arr[0];
    for (int i = 0; i < N; i++) {
        if (arr[i] >= max_val) {
            max = i;
            max_val = arr[i];
        }
    }
    return max;
}
 
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

// function designed to check between +/-L and itegrate the wavefn between those valeus and return the area
double integrate(
        double *sq_wf,
        double *X,
        double L,
        int size) {
    // intialise 2 counts, one for the left (count), and one for the right (upper)
    // they will be used to get the index of =/- L
    int count = 0;
    int upper = 0;

    double del_x = X[1]-X[0];
    
    // incriment the count value until it exceeds -L, at this point the counter is 1 short
    while (X[count] < -L) {
        count++;
    }
    // make up for the counter being one short
    count++;

    // incriment the upper value until it exceed +L, at this point the counter is 1 above the correct value
    while (X[upper] < L) {
        upper++;
    }
    // decriment the counter to get the index that gives a value less than +L
    upper--;


    // initialise a value Sum to sum the following calulations
    // below integration is using Simpson's 3/8ths rule
    double Sum = 0;
    // initialise a variable for the ans
    double ans;
    // sum first and last values who's coefficients by the 3/8th rule are 1
    Sum += sq_wf[count]+ sq_wf[upper];

    // sum up every point who's index is not divisible by 3 and multiply by the correct coefficient (3)
    for (int i = count + 1; i < upper; i++) {
        // check if i is divisible by 3 if not then 
        if (i % 3 != 0) {
            Sum += 3*(sq_wf[i]);
        }

    }
    
    // sum up every point who's index is divisible by 3 and multiply by the correct coefficient (2)
    for (int i = count + 1; i < (upper/3); i++) {
        Sum += 2*(sq_wf[3*i]);
    }
    // multply the sum by 3/8 * the width of a step
    ans = ((3*del_x)/(8))*Sum;


}

double bound_integral(double *func, 
                      double *X, 
                      double lower_bound, 
                      double upper_bound, 
                      int N) {
    int l_index = indexing_left(X, lower_bound, N);
    int r_index = indexing_left(X, upper_bound, N);
    double del_x = X[1]- X[0];

    // initialise a value Sum to sum the following calulations
    // below integration is using Simpson's 3/8ths rule
    double Sum = 0;
    // initialise a variable for the ans
    double ans;

    Sum += func[l_index] + func[r_index];


    for (int i = l_index + 1; i < r_index; i++) {
        // check if i is divisible by 3 if not then 
        if (i % 3 != 0) {
            Sum += 3*(func[i]);
        }

    }
    
    // sum up every point who's index is divisible by 3 and multiply by the correct coefficient (2)
    for (int i = l_index + 1; i < (r_index/3); i++) {
        Sum += 2*(func[3*i]);
    }
    // multply the sum by 3/8 * the width of a step
    ans = ((3*del_x)/(8))*Sum;
    return ans;
}


// function to normalise the wavefn
double *Normalise(double *wavefn, double area, double size) {
    // loop through the indexs upto size
    for (int i = 0; i < size; i++) {
        /*
        multiply each value of the wavefn by a constant factor that will satisfy that the probablilty 
        between two bounds being 1
        */
        wavefn[i] = (1/(sqrt(area)))*wavefn[i];
    }
    // return the new wavefunction
    return wavefn;
}


double finite_derive(double *Y, double *X, int index) {
    double del_x = abs(X[index] -  X[index + 1]);
    double deriv = (Y[index+1]- Y[index])/del_x;
    return deriv;
}

double *scale(double *wavefn, double factor, int N){
    for (int i = 0; i < N; i++){
        wavefn[i] = factor*wavefn[i];
    }
    return wavefn;
}
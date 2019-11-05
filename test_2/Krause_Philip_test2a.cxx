/*
    Philip Krause
    118470776

    2nd inclass test, based on numerical integration
*/

#include <iostream>
#include <string>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
using namespace std;
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>


// function to compute the integral of the fn, f(x) = 1 + b*x^2 - c*x^3, between 0 and a, by the trivial rule
double trivial(double a, double b, double c) {
    // initialize variables to be used

    // number of points to be taken
    int iterations = 100;

    // the spacial step in x direction
    double del_x = a/iterations;

    // variable to be used as the sume at the points
    double sum = 0;

    // variable for the answer
    double ans;

    // iterate from 0 to 99 and sum up f(x) at that value of x, where x(n) = n*del_x, because x(0) = 0
    for (int i = 0; i < iterations; i++) {
        sum += (1 + b*i*del_x*i*del_x - c*(i*del_x*i*del_x*i*del_x));
    }
    // multiply the sum of the y value points by the x step to get the approximate area
    ans = del_x * sum;

    // return the answer to the calling function
    return ans;
}


// function to compute the integral of the fn, f(x) = 1 + b*x^2 - c*x^3, between 0 and a, by simpsons rule
double simpson(double a, double b, double c) {
    // initialize variables to be used

    // number of points to be taken
    int N = 100;

    // the spacial step in x direction
    double del_x = a/N;

    // variable to be used as the sume at the points
    double sum = 0;

    // variable for the answer
    double ans;

    /*  iterate from 0 to 99 and sum up f(x) at that value of x, where x(n) = n*del_x, because x(0) = 0
        for the first and last values only only add the point of f(0) and f(N) without multiplying
        for the odd values of x multply the value for f(x=odd) by 4
        for the even values of x multply the value for f(x=even) by 2

        check if it's odd or even by getting the modulo of the iteration
    */
    for (int i = 0; i <= N; i++) {
        if (i == 0 || i == N){
            sum += (1 + b*i*i*del_x*del_x - c*i*i*i*del_x*del_x*del_x);
        } else if ((i % 2) == 1) {
            sum += 4*(1 + b*i*i*del_x*del_x - c*i*i*i*del_x*del_x*del_x);
        } else if ((i % 2) == 0) {
            sum += 2*(1 + b*i*i*del_x*del_x - c*i*i*i*del_x*del_x*del_x);
        }
    }

    // multiply by del_x/3  as per formula
    ans = (sum/3)*del_x;

    //return answer
    return ans;

}


// main fn
int main() {
    // declare variables

    // declare parameters
    double a, b, c;

    // declare variables to store answers
    double triv_ans, simp_ans;


    // prompt user for parameters
    cout << "Enter a value for a: ";
    cin >> a;
    cout << "Enter a value for b: ";
    cin >> b;
    cout << "Enter a value for c: ";
    cin >> c;


    // calculate trivial answer 
    triv_ans = trivial(a, b, c);


    // calculate answer as per simpsons rule
    simp_ans = simpson(a, b, c);

    // print answers
    cout << "Answer as per trivial rule is: " << triv_ans << endl;
    cout << "Answer as per Simpson rule is: " << simp_ans << endl;

    return 0;
}


/*

    Q2 analytical ans = 4
    I assume I dont have to show that 

    trivial rule has error of 0.3168
    simpsons rule mataches analytical answer exactly


*/
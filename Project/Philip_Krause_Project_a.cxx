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
#include "Philip_Krause_Project_b.cxx"
#include "Philip_Krause_Project_c.cxx"
#include "Philip_Krause_Project_d.cxx"
#include "Philip_Krause_Project_f.cxx"

// function for generating even parity solutions of the wave function in a square well
void even_squence(double energy){
    //initialise variables

    // number of points
    int N = 50000;
    // potential array,V, displacement array,X, wave function array wavefn
    double V[N], X[N], wavefn[N];
    // distance to be graphed
    double graphing_distance = 5;

    // array for energy (used arrays for pointers so it's globally modifiable)
    double E[1];
    // energy input ^^ parameter 
    E[0]= energy;
    // array for change of energy (used arrays for pointers so it's globally modifiable)
    double del_E[1];
    // initialise variable
    del_E[0] = 0.5;
    // set cut off so that divergence does not go to +/- infinity
    double cutoff = 2.5;
    //potential outside of 0 potential area
    double Volt = 10000;
    // the point outside which it is not 0 potential
    double L = 1;
    // variable to track divergence
    double last_diverge[1];
    // to make sure it's not 1 or -1
    last_diverge[0] = 0;
    // calculate spacial step
    double del_x = graphing_distance/N;
    // arrays used in normalising wave function
    double squared_function[N], norm_wf[N];

    // generate potential well
    // -- code -- > Philip_Krause_Project_b.cxx
    gen_v(V, X, N, Volt, graphing_distance, L);
    // generate wave function with initial guss for energy
    // -- code -- > Philip_Krause_Project_c.cxx
    gen_phi_even(wavefn, V, X, N, graphing_distance, E[0], cutoff);

    // start of interactive routine that checks if graph is ok else continues
    // boolean to track if it's good or not assume it's not => true, if it is => false
    bool cont = true;
    // count to keep track of iterations
    int count = 0;
    // variable for user input as a 1 letter character
    char prompt;
    // array of chars for the name used in formatting with sprinf
    char title[70];
    // loop while wavefunction is not ok
    while (cont) {
        // get the last direction of divergence and change energy accordingly, and set the new divergence
        // -- code -- > Philip_Krause_Project_c.cxx
        calc(V, X, wavefn, N, last_diverge, graphing_distance, E, del_E, cutoff);
        // print out the change of energy to see if it's small enoug
        cout << "del E is: " << del_E[0] << endl;
        // generate the wavefunction to test in next itteration
        // -- code -- > Philip_Krause_Project_c.cxx
        gen_phi_even(wavefn, V, X, N, graphing_distance, E[0], cutoff);

        // if statement to execute every 50 loops so that the program doesn't take ages of checks
        // displays values for energy, as well as graph and asks if the loop should be continued
        if (count % 50 == 0) {
            sprintf(title, "Even Wave Function with Energy = %.5lf", E[0]);
            gnuplot_one_function(title, "linespoints", "X", "phi", X, wavefn, N);
            // prompts
            cout << "Continue: (Y: Yes, N: No)";
            // gets input
            cin >> prompt;
            // switch statement rather than if else ladder to execute if value is 'N' or 'n'
            switch (prompt) {
                case 'n' :
                case 'N' :
                    cont = false;
            }
        }
        // incriment count
        count++;

    }
    // display final energy
    cout << "E is: " << E[0] << endl;
    // mulitple line breaks so it's not too close to next funciton
    cout << "\n\n\n\n" << endl;
    // square the wave function
    // -- code -- > Philip_Krause_Project_d.cxx
    square(wavefn, squared_function, N);
    // copy wavefunction into new array to be normalised elment wise
    for (int i = 0; i < N; i++) {
        norm_wf[i] = wavefn[i];
    }
    // get area under the squared function
    // -- code -- > Philip_Krause_Project_d.cxx
    double area = integrate(squared_function, X, L, N);
    // normalise the wave function so area of squared function is 1
    // -- code -- > Philip_Krause_Project_d.cxx
    Normalise(norm_wf, area, N);
    

    // array of chars for use in editting the file name of the jpg to be made
    char filename[35];
    
    // formatting of the filename and title of graphs, formats them so the energy is visable
    sprintf(filename, "plot of Even Wave function E=%.5lf norm.jpg", E[0]);
    sprintf(title, "Even Wave Function with Energy = %.5lf and normalised wavefunction", E[0]);
    // plots the wave function and the normalised wave function to a jpg
    gnuplot_two_functions_jpg(title, "lines", "X", "Wave function", X, wavefn, N, "Normalised Wavefunction", X, norm_wf, N, "Norm Wave function", filename);



}

// function for generating odd parity solutions of the wave function in a square well
void odd_squence(double energy){
    //initialise variables

    // number of points
    int N = 50000;
    // potential array,V, displacement array,X, wave function array wavefn
    double V[N], X[N], wavefn[N];
    // distance to be graphed
    double graphing_distance = 5;

    // array for energy (used arrays for pointers so it's globally modifiable)
    double E[1];
    // energy input ^^ parameter 
    E[0]= energy;
    // array for change of energy (used arrays for pointers so it's globally modifiable)
    double del_E[1];
    // initialise variable
    del_E[0] = 0.5;
    // set cut off so that divergence does not go to +/- infinity
    double cutoff = 2.5;
    //potential outside of 0 potential area
    double Volt = 10000;
    // the point outside which it is not 0 potential
    double L = 1;
    // variable to track divergence
    double last_diverge[1];
    // to make sure it's not 1 or -1
    last_diverge[0] = 0;
    // calculate spacial step
    double del_x = graphing_distance/N;
    // arrays used in normalising wave function
    double squared_function[N], norm_wf[N];



    // generate potential well
    // -- code -- > Philip_Krause_Project_b.cxx
    gen_v(V, X, N, Volt, graphing_distance, L);
    // generate wave function with initial guss for energy
    // -- code -- > Philip_Krause_Project_c.cxx
    gen_phi_odd(wavefn, V, X, N, graphing_distance, E[0], cutoff);

    // start of interactive routine that checks if graph is ok else continues
    // boolean to track if it's good or not assume it's not => true, if it is => false
    bool cont = true;
    // count to keep track of iterations
    int count = 0;
    // variable for user input as a 1 letter character
    char prompt;
    // array of chars for the name used in formatting with sprinf
    char title[70];
    // loop while wavefunction is not ok
    while (cont) {
        // get the last direction of divergence and change energy accordingly, and set the new divergence
        // -- code -- > Philip_Krause_Project_c.cxx
        calc(V, X, wavefn, N, last_diverge, graphing_distance, E, del_E, cutoff);
        // print out the change of energy to see if it's small enoug
        cout << "del E is: " << del_E[0] << endl;
        // generate the wavefunction to test in next itteration
        // -- code -- > Philip_Krause_Project_c.cxx
        gen_phi_odd(wavefn, V, X, N, graphing_distance, E[0], cutoff);
        // if statement to execute every 50 loops so that the program doesn't take ages of checks
        // displays values for energy, as well as graph and asks if the loop should be continued
        if (count % 50 == 0) {
            sprintf(title, "Odd Wave Function with Energy = %.5lf", E[0]);
            gnuplot_one_function(title, "linespoints", "X", "phi", X, wavefn, N);
            // prompts
            cout << "Continue: (Y: Yes, N: No)";
            // gets input
            cin >> prompt;
            // switch statement rather than if else ladder to execute if value is 'N' or 'n'
            switch (prompt) {
                case 'n' :
                case 'N' :
                    cont = false;
            }
        }
        // incriment count
        count++;

    }
    // display final energy
    cout << "E is: " << E[0] << endl;
    // mulitple line breaks so it's not too close to next funciton
    cout << "\n\n\n\n" << endl;
    // square the wave function
    // -- code -- > Philip_Krause_Project_d.cxx
    square(wavefn, squared_function, N);
    // copy wavefunction into new array to be normalised elment wise
    for (int i = 0; i < N; i++) {
        norm_wf[i] = wavefn[i];
    }
    // get area under the squared function
    // -- code -- > Philip_Krause_Project_d.cxx
    double area = integrate(squared_function, X, L, N);
    // normalise the wave function so area of squared function is 1
    // -- code -- > Philip_Krause_Project_d.cxx
    Normalise(norm_wf, area, N);
    

    // array of chars for use in editting the file name of the jpg to be made
    char filename[35];
    
    // formatting of the filename and title of graphs, formats them so the energy is visable
    sprintf(filename, "plot of odd Wave function E=%.5lf norm.jpg", E[0]);
    sprintf(title, "Odd Wave Function with Energy = %.5lf and normalised wavefunction", E[0]);
    // plots the wave function and the normalised wave function to a jpg
    gnuplot_two_functions_jpg(title, "lines", "X", title, X, wavefn, N, "Wave function", X, norm_wf, N, "Norm Wave function", filename);
}


// funciton for generating the wave function for the non symmetrical Lennard Jones Potential
void LJVsequence(double init_E, double epsilon, double sigma){
    // number of points to be evaluated
    int N = 10000;
    // end point for the graph
    double end = 5;
    // start point for the graph
    double x_0 = 0.9;

    // arrays for potential, displacement, as well as left, right, squared and final normalised wavefunction
    double V[N], X[N], wavefn[N], wavfn[N], wave_fn[N], square_fn[N];
    // slight offset so that phi_0, and phi_n1 (the -1th) values can be added to array on the left
    double x_r = end - 0.005;
    // as above for right side
    double x_l = x_0 + 0.005;
    
    // generate the lennard-Jones potential
    // -- code -- > Philip_Krause_Project_f.cxx
    LJV(V, X, x_0, end, N, sigma, epsilon);


    // array so that the value for energy can be accessed out side funciton
    double Energy[1];
    // function returns energy as well as goes through the routine of adjusting the values of energy and normalises it
    // -- code -- > Philip_Krause_Project_f.cxx
    Energy[0] = non_symadjust(V, X, wavefn, wavfn, wave_fn, N, init_E, 0.05, x_l, x_r, x_0, end);

    // graph file, and title name formatting 
    char title[100];
    sprintf(title, "Wavefunction for Lennard Jones Potential E = %.5lf", Energy[0]);
    char file[100];
    sprintf(file, "Wavefunction_for_Lennard_Jones_Potential_E_=_%.5lf.jpg", Energy[0]);
    char LenJon[100];
    sprintf(LenJon, "Lennard-Jones epsilon = %.1lf, sigma = %.1lf, (scaled to fit graph)", epsilon, sigma);
    //plots final results and saves it to a jpg
    gnuplot_two_functions("Wave function for non symmetric potential", "linespoints", "X axis", "Wavefunction Amplitude", X, wave_fn, N, title , X, V, N, LenJon);
    gnuplot_one_function_jpg("Wave function for non symmetric potential", "lines", "X axis", title, X, wave_fn, N, file);
}

// function for generating the wave function in a symmetrical square well with a small energy barrier also present
void barrier_sequence(double energy){
    // initialise variable
    // number of points
    int N = 10000;
    // end point of graph
    double end = 2;
    // outer potential step
    double Vstep = 4000;
    // inner potential step (barrier)
    double mVstep = 20;
    // -/+point at which the inner barrier exists
    double m_l = 0.25;
    // start point of graph
    double x_0 = -2;
    // arrays for displacement, left, right, squared and final normalised wavefunction as well as potential
    double  X[N], wavefn[N], wavfn[N], wave_fn[N], square_fn[N], sym_V[N];
    // offset for initial conditions
    double x_r = end - 0.005;
    // offset for initial conditions
    double x_l = x_0 + 0.005;

    // function for generating a square well and a small potential barrier
    // -- code --> Philip_Krause_Project_b.cxx
    mini_step(sym_V, X, N, Vstep, 4, 2, mVstep, 0.5);

    // array so that the value for energy can be accessed out side funciton
    double Energy[1];
    // function returns energy as well as goes through the routine of adjusting the values of energy and normalises it
    // -- code -- > Philip_Krause_Project_f.cxx
    Energy[0] = symadjust(sym_V, X, wavefn, wavfn, wave_fn, N, energy, 0.05, x_l, x_r, x_0, end);


    // graph file, and title name formatting
    char title[100];
    sprintf(title, "Wavefunction for small potential barrier E = %.5lf, potential outside = %.1lf, barrier potential = %.1lf", Energy[0], Vstep, mVstep);
    char file[100];
    sprintf(file, "Wavefunction for small potential barrier E = %.5lf.jpg", Energy[0]);

    //plots final results and saves it to a jpg
    gnuplot_one_function(title, "linespoints", "X", "phi", X, wave_fn, N);
    gnuplot_one_function_jpg(title, "lines", "X", "Phi(x)", X, wave_fn, N, file);

 
}



// main function
int main() {
    // lowest energy for the even solution
    even_squence(1);
    // lowest energy for the odd solution
    odd_squence(4);
    // higher energy of the even solution
    even_squence(16);
    // higher energy of the odd solution
    odd_squence(40);
    // lowest energy for the symmetrical square well with potential barrier
    barrier_sequence(2);
    // higher energy for the symmetrical square well with potential barrier
    barrier_sequence(9);
    // lowesr energy solution for the lennard jones potential
    LJVsequence(-1.6, 10, 1);
    return 0;
}
/*
    Philip Krause
    118470776@umail.ucc.ie
    Project 9 Quantum Assignment
    testing for wave fn

    File to generate wave functions for symmetrical potentials
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
// #include "gnuplot.cxx"
// #include "Philip_Krause_Project_d.cxx"
/*
//calculate fn, decides whether the function is diverging and adusts the change in energy accordingly
double* calc(double *V, 
            double *X, 
            double *wavefn,
            int size, 
            double* last_diverge,
            double graphing_distance,
            double* E_0, 
            double* del_E, 
            double cutoff) {
    
    /*
    for arguments
        double* V, is the pointer to the memory address for an array of doubles which will serve as
            the potenitial at a given point (as an array of potential for various points)
        double* X, is the pointer to the memory addess for the x axis and is used to determine where
            the potential stepts are located
        double *wavefn, is the pointer to the memory address for an array of doubles where the 
            wavefunction solutions is stored
        int size, is the number of steps
        double* last diverge is a pointer to the memory address to a double which tracks the 
            direction of the previous divergence
        double E_0, is a pointer to the location of the initial guess for energy which get changed
        double graphing distance is the total distance that is graphed using gnuplot
        del_E, is the pointer to the incriment by which the Energy is changed
        double cutoff, is the variable to ehich the function checks if the wavefunction is diverging
            beyond that and cuts off the wavefunction before it diverges to inifinity    
    */
    
    // printiting out the value of the diverging wave
    // cout << wavefn[size -5] << endl;
/*
    // if-else if ladder to check various conditions and making decisions based on that
    if (last_diverge[0] == 0 && wavefn[size - 5] > 0) {
        // first iteration, diverging up, increase energy
        E_0[0] = E_0[0] + del_E[0];
    } else if (last_diverge[0] == 0 && wavefn[size - 5] < 0) {
        // first iteration, diverging down, decrease energy
        E_0[0] = E_0[0] - del_E[0];
    } else if (last_diverge[0] == 1 && wavefn[size - 5] < 0) {
        // previously diverging up, now divering down, overshot
        del_E[0] = del_E[0] /2 ;
        E_0[0] = E_0[0] - del_E[0];
    } else if (last_diverge[0] == 1 && wavefn[size - 5] > 0) {
        // Still divering up, energy needs to increase
        E_0[0] = E_0[0] + del_E[0];
    } else if (last_diverge[0] == -1 && wavefn[size - 5] < 0) {
        // still diverging downwards, energy needs to decrease
        E_0[0] = E_0[0] - del_E[0];
    } else if (last_diverge[0] == -1 && wavefn[size - 5] > 0) {
        // previously diverging down, now diverging up
        // Energy needs to increase
        del_E[0] = del_E[0]/2;
        E_0[0] = E_0[0] + del_E[0];
    }
    // check the sign of the last divergence
    // if lessa than 0 and the absolute value is gr
    if (wavefn[size - 5] < 0 && abs(wavefn[size - 5]) > cutoff) {
        last_diverge[0] = -1;
    } else {
        last_diverge[0] = 1;
    }
    //return the calculated values
    return E_0, del_E, last_diverge;


}

// double integrate(double *wavefn,
//                  double *X,
//                  )


double *gen_phi_even(double *wavefn,
                double *V,
                double *X,
                int size,
                double graphing_distance, 
                double E,
                double cutoff) {

    /*
    for arguments
        double* V, is the pointer to the memory address for an array of doubles which will serve as
            the potenitial at a given point (as an array of potential for various points)
        double* X, is the pointer to the memory addess for the x axis and is used to determine where
            the potential stepts are located
        double *wavefn, is the pointer to the memory address for an array of doubles where the 
            wavefunction solutions is stored
        int size, is the number of steps
        double* last diverge is a pointer to the memory address to a double which tracks the 
            direction of the previous divergence
        double E, is the initial guess for the energy
        double graphing distance is the total distance that is graphed using gnuplot
        double cutoff, is the variable to ehich the function checks if the wavefunction is diverging
            beyond that and cuts off the wavefunction before it diverges to inifinity    
    */
    // for a symmetric even parity wave the initial conditions are always
    // these conditions are 1, 1, corresponging to phi(x) and d/dx(phi(x))
/*
    double phi_0 = 1;
    double phi_n1 = 1;

    // calculate the the spacial step
    double del_x = graphing_distance/size;
    
    // intialise variable to be used as temperary location for a given answer for the wave function,
    // to test if it's diverging
    double phi_temp;
    
    // variable to get the mid point of the function
    int mid;
    // check the number of points if even
    if (size % 2 == 0) {
        // the mid point is half the total number of points
        mid = (size/2);
    } else {
        // else take one from the number of points and half it to get the midpoint
        mid = size -1;
        mid = (mid/2);
    }
    
    // work out intial conditions for the wave function
    wavefn[mid] = 2*phi_0 - phi_n1 - 2*del_x*del_x*(E - V[mid-1])*phi_0;
    // work out the initial condition for one step right of 0
    wavefn[mid + 1] = 2*wavefn[mid] - phi_0 - 2*del_x*del_x*(E - V[mid-1])*wavefn[mid];
    // work out the initial condition for one step left of 0
    wavefn[mid - 1] = 2*wavefn[mid] - phi_0 - 2*del_x*del_x*(E - V[mid-1])*wavefn[mid];

    // iterate from the 2nd (greater than 0) to the Nth step and evalue at the wave function as it goes
    for (int i = 2; i < mid; i++) {
        // store the value for the wave function at that point to check for diverence
        phi_temp = 2*wavefn[mid + i - 1] - wavefn[mid + i - 2] - 2*del_x*del_x*(E - V[mid + i - 1])*wavefn[mid + i - 1];

        // check if the wave function is diverging if it's not
        if (abs(phi_temp) < cutoff) {

            // assign the value of that wave fn to the correct location
            wavefn[mid + i] = phi_temp;
        } else {
            // else it is diverging and check whether it is diverging up or down
            if (phi_temp > 0) {
                // if divering up
                for (int n = 0; n < (mid - i); n++)
                {
                    // assign every value greater than the ith as 2.5
                    wavefn[mid + i + n] = 2.5;
                }

                // break out the loop and move on to the next section
                break;
            // else it's diverging down
            } else if (phi_temp <  0) {
                // divering down
                for (int n = 0; n < (mid - i); n++)
                {
                    // assign every value less than the ith as -2.5
                    wavefn[mid + i + n] = -2.5;
                }
                // break out the loop and move on to the next section
                break;
            }
            
            
        }
        
    }


    // iterate from 2 steps less than 0 (mid-th value) to the 0th value 
    for (int i = 2; i < mid+2; i++) {
        // store the value for the wave function at that point to check for diverence
        phi_temp = 2*wavefn[mid - i + 1] - wavefn[mid - i + 2] - 2*del_x*del_x*(E - V[mid - i + 1])*wavefn[mid - i + 1];
        
        // check if the wave function is diverging if it's not
        if (abs(phi_temp) < cutoff) {
            // assign the value of that wave fn to the correct location
            wavefn[mid - i] = phi_temp;
        } else {
            // else it is diverging and check whether it is diverging up or down
            if (phi_temp > 0) {
                // if divering up
                for (int n = 0; n < (mid - i); n++)
                {
                    // assign every value less than the ith as 2.5
                    wavefn[mid - i - n] = 2.5;
                }
                break;
            // else it's diverging down
            } else if (phi_temp <  0) {
                // if divering down
                for (int n = 0; n < ( mid -i); n++)
                {
                    // assign every value less than the ith as -2.5
                    wavefn[mid - i - n] = -2.5;
                }

                // break out of the loop
                break;
            }
        }
        
    }


    // return the calculated wavefunction
    return wavefn;
}


// fn to calculate the wave fn for the odd parity ie phi(x) = -phi(-x)
double *gen_phi_odd(double *wavefn,
                double *V,
                double *X,
                int size,
                //double phi_0,
                //double phi_n1,
                double graphing_distance, 
                double E,
                double cutoff) {
    
    // calculate 
    double del_x = graphing_distance/size;
    double phi_0 = 0;
    double phi_n1 = -del_x;


    double phi_temp;

    int mid;
    if (size % 2 == 0) {
        mid = (size/2);
    } else {
        mid = size -1;
        mid = (mid/2);
    }
        
    wavefn[mid] = phi_0;
    wavefn[mid - 1] = -phi_n1;
    wavefn[mid + 1] = phi_n1;
    
    
    for (int i = 2; i < mid; i++) {
        phi_temp = 2*wavefn[mid+i-1] - wavefn[mid+i-2] - 2*del_x*del_x*(E - V[mid+i-1])*wavefn[mid+i-1];
        if (abs(phi_temp) < cutoff) {
            wavefn[mid + i] = phi_temp;
        } else {
            if (phi_temp > 0) {
                // divering up
                for (int n = 0; n < (mid - i); n++)
                {
                // stop loop and set all values to the cut off
                wavefn[mid + i + n] = 2.5;
                }
                break;
            } else if (phi_temp <  0) {
                // divering down
                for (int n = 0; n < (mid - i); n++)
                {
                // stop loop and set all values to the cut off
                wavefn[mid + i + n] = -2.5;
                }
                break;
            }
            
            
        }
        
    }

    
    for (int i = 2; i < mid+2; i++) {
        phi_temp = 2*wavefn[mid - i + 1] - wavefn[mid - i + 2] - 2*del_x*del_x*(E - V[mid - i + 1])*wavefn[mid - i + 1];
        if (abs(phi_temp) < cutoff) {
            wavefn[mid - i] = phi_temp;
        } else {
            if (phi_temp > 0) {
                // divering up
                for (int n = 0; n < (mid - i); n++)
                {
                // stop loop and set all values to the cut off
                    wavefn[mid - i - n] = 2.5;
                }
                break;
            } else if (phi_temp <  0) {
                // divering down
                for (int n = 0; n < ( mid -i); n++)
                {
                // stop loop and set all values to the cut off
                wavefn[mid - i - n] = -2.5;
                }
                break;
            }
        }
        
    }



    return wavefn;
}

double *gen_phi_left(double *wavefn, double *V, double *X, int N, double E, double x_l, double start, double end){
    int ind = indexing_left(X, x_l, N);
        double offset =abs(start - x_l);
    while (ind < 2) {
        x_l += 2*offset;
        ind = indexing_left(X, x_l, N);
    }

    double del_x = (end - start)/N;
    // cout << del_x << endl;
    // cout << abs(del_x) << endl;
    wavefn[ind - 2] = -0.0001*abs(del_x);
    wavefn[ind -1] = 0;

    for (int i = 0; i < N- ind; i++) {
        wavefn[ind +i] = 2*wavefn[ind +i -1] -wavefn[ind + i -2] - 2*del_x*del_x*(E- V[ind +i -1])*wavefn[ind +i -1];
        // cout << i << "      " << wavefn[ind +i] << "    " << (wavefn[ind +i] - wavefn[ind +i -1]) << "      " << V[ind +i -1] << endl;
    }


    return wavefn;
}


double *gen_phi_right(double *wavefn, double *V, double *X, int N, double E, double x_r, double start, double end) {
    int ind = indexing_right(X, x_r, N);
    cout << ind << endl;
    while (ind < 2) {
        x_r = 3*x_r;
        ind = indexing_right(X, x_r, N);
    }
    double del_x = (end - start)/N;
    // cout << del_x << endl;
    // cout << abs(del_x) << endl;
    wavefn[ind + 2] = -0.0001*abs(del_x);
    wavefn[ind + 1] = 0;

    for (int i = 0; i <= ind; i++) {
        wavefn[ind-i] = 2*wavefn[ind-(i -1)] - wavefn[ind-(i-2)] -2*del_x*del_x*(E-V[ind-(i-1)])*wavefn[ind-(i-1)];
    }

    return wavefn;

}



double *non_symadjust(double *V, double *X, double *wavefn_1, double *wavefn_2, double *wave_function, int N, double init_E, double inc_E, double x_l, double x_r, double start, double end) {
    // find minium potential of on the range
    int minimum = locate_min(V, N);
    double tolerance = ((end - start)/N)*((end - start)/N)*((end - start)/N)*((end - start)/N);
    // cout << tolerance << endl;
    gen_phi_left(wavefn_1, V, X, N, init_E, x_l, start, end);
    // cout << "left" << endl;
    gen_phi_right(wavefn_2, V, X, N, init_E, x_r, start, end);
    gnuplot_two_functions("test plot", "lines", "X", "phi", X, wavefn_1, N, "left", X, wavefn_2, N, "right");
    double E = init_E;
    double ratio = (wavefn_1[minimum]/wavefn_2[minimum]);
    // cout << ratio<< endl;
    scale(wavefn_2, ratio, N);
    gnuplot_two_functions("test plot", "lines", "X", "phi", X, wavefn_1, N, "left", X, wavefn_2, N, "right");
    int match = 1;
    double diff = finite_derive(wavefn_1, X, minimum) - finite_derive(wavefn_2, X, minimum); // >0 then d_phi_l/dx > d_phi_r/dx else l<r
    // cout << diff << endl;

    int count = 1;
    int ans;
    if (diff >0 ) {
        match = 1;

    }else if (diff <0) {
        match = -1;
    }

    while (abs(diff) > tolerance) {
        cout << count << endl;
        if ((match == 1) && (diff >0)){
            // energy needs to increase
            // cout << "1   1  " << E << endl;
            E = E + inc_E;
            // cout << E << endl;
        } else if ((match == 1) && (diff <0)) {
            // increased last time so needs to decrease
            inc_E  = inc_E/2;
            // cout << "1  -1  " << E << endl;
            E = E - inc_E;
            // cout << E << endl;
            match = -1;
        } else if ((match == -1) && (diff >0)) {
            // decreased last time needs to increase
            inc_E  = inc_E/2;
            // cout << "-1  1  " << E << endl;
            E = E + inc_E;
            // cout << E << endl;
            match = 1;
        } else if ((match == -1) && (diff <0)) {
            // energy needs to decrease
            // cout << "-1 -1  " << E << endl;
            E = E - inc_E;
            // cout << E << endl;
        }
        gen_phi_left(wavefn_1, V, X, N, E, x_l, start, end);
        gen_phi_right(wavefn_2, V, X, N, E, x_r, start, end);
        ratio = (wavefn_1[minimum]/wavefn_2[minimum]);
        scale(wavefn_2, ratio, N);
        diff = finite_derive(wavefn_1, X, minimum) - finite_derive(wavefn_2, X, minimum); // >0 then d_phi_l/dx > d_phi_r/dx else l<r
        count++;
        // gnuplot_two_functions("test plot", "lines", "X", "phi", X, wavefn_1, N, "left", X, wavefn_2, N, "right");
        // cin >> ans;
    }
    for (int i = 0; i < N; i++) {
        if (i <= minimum) {
            wave_function[i] = wavefn_1[i];
        } else {
            wave_function[i] = wavefn_2[i];
        }
    }

}

double *symadjust(double *V, double *X, double *wavefn_1, double *wavefn_2, double *wave_function, int N, double init_E, double inc_E, double x_l, double x_r, double start, double end) {
    // find minium potential of on the range
    // int minimum = locate_min(V, N);
    int mid = N/2;
    double tolerance = ((end - start)/N)*((end - start)/N)*((end - start)/N)*((end - start)/N);
    // cout << tolerance << endl;
    gen_phi_left(wavefn_1, V, X, N, init_E, x_l, start, end);
    // cout << "left" << endl;
    gen_phi_right(wavefn_2, V, X, N, init_E, x_r, start, end);
    gnuplot_two_functions("test plot", "lines", "X", "phi", X, wavefn_1, N, "left", X, wavefn_2, N, "right");
    double E = init_E;
    double ratio = (wavefn_1[mid]/wavefn_2[mid]);
    // cout << ratio<< endl;
    scale(wavefn_2, ratio, N);
    gnuplot_two_functions("test plot", "lines", "X", "phi", X, wavefn_1, N, "left", X, wavefn_2, N, "right");
    int match = 1;
    double diff = finite_derive(wavefn_1, X, mid) - finite_derive(wavefn_2, X, mid); // >0 then d_phi_l/dx > d_phi_r/dx else l<r
    // cout << diff << endl;

    int count = 1;
    int ans;
    if (diff >0 ) {
        match = 1;

    }else if (diff <0) {
        match = -1;
    }

    while (abs(diff) > tolerance) {
        cout << count << endl;
        if ((match == 1) && (diff >0)){
            // energy needs to increase
            // cout << "1   1  " << E << endl;
            E = E + inc_E;
            // cout << E << endl;
        } else if ((match == 1) && (diff <0)) {
            // increased last time so needs to decrease
            inc_E  = inc_E/2;
            // cout << "1  -1  " << E << endl;
            E = E - inc_E;
            // cout << E << endl;
            match = -1;
        } else if ((match == -1) && (diff >0)) {
            // decreased last time needs to increase
            inc_E  = inc_E/2;
            // cout << "-1  1  " << E << endl;
            E = E + inc_E;
            // cout << E << endl;
            match = 1;
        } else if ((match == -1) && (diff <0)) {
            // energy needs to decrease
            // cout << "-1 -1  " << E << endl;
            E = E - inc_E;
            // cout << E << endl;
        }
        gen_phi_left(wavefn_1, V, X, N, E, x_l, start, end);
        gen_phi_right(wavefn_2, V, X, N, E, x_r, start, end);
        ratio = (wavefn_1[mid]/wavefn_2[mid]);
        scale(wavefn_2, ratio, N);
        diff = finite_derive(wavefn_1, X, mid) - finite_derive(wavefn_2, X, mid); // >0 then d_phi_l/dx > d_phi_r/dx else l<r
        count++;
        // gnuplot_two_functions("test plot", "lines", "X", "phi", X, wavefn_1, N, "left", X, wavefn_2, N, "right");
        // cin >> ans;
    }
    for (int i = 0; i < N; i++) {
        if (i <= mid) {
            wave_function[i] = wavefn_1[i];
        } else {
            wave_function[i] = wavefn_2[i];
        }
    }

}


/*

int main() {
    int N = 10000;
    double wavefn[N], X[N], V[N];
    double graphing_distance = 4;
    double E[1];
    E[0] = 4;
    double cutoff = 2.5;
    gen_v(V, X, N, 10000, graphing_distance, 1);


    gen_phi_odd(wavefn, V, X, N, graphing_distance, E[0], cutoff);

    gnuplot_one_function("test", "linespoints", "x", "phi", X, wavefn, N);

    return 0;
}
*/

#ifndef __Simulation_hpp__
#define __Simulation_hpp__

#include <sstream> 
#include <string>
#include <iomanip>
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <armadillo>
#include <typeinfo>
#include <complex>
#include <utils.hpp>
#include <Box.hpp>
#include <chrono>

class Simulation{
    // member variables:
    public:

        // Constructor
        Simulation(Box BBXa, double Dta=2.5e-5, double Ta=0.008, double xca=0.25, double sigma_xa=0.05, double p_xa=200, double yca=0.5, double sigma_ya=0.05, double p_ya=0);

        double Dt;
        double T;
        double xc;
        double sigma_x;
        double p_x;
        double yc;
        double sigma_y;
        double p_y;
        int Nt;
        Box BBX;

        //  Initial states
        arma::cx_mat U0; //initial state as matrix
        arma::cx_vec u; //initial state as vector
        
        // Matrices A and B
        arma::sp_cx_mat A, B;

        // Cube to store the U matrices in time
        arma::cx_cube U;



        // Member functions
        void initialise(void);

        arma::cx_cube run_simulation(void);
};

#endif
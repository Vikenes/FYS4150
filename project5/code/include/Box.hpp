#ifndef __Box_hpp__
#define __Box_hpp__

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
// #include <Simulation.hpp>


class Box{
    private:
        friend class Simulation;
    //  member variables
    public:
        int M;
        arma::sp_mat V;
        double h;

        // constructor
        Box(double h=0.005);

        void set_up_walls(double v0, int Ns=2, double Th=0.02, double wxc=0.5, double wyc=0.5, double wall_piece_length=0.05, double aperture=0.05);


};



#endif
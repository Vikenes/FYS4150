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
#include <assert.h>


class Box{
    private:
        // Member variables

        int M;              //  number of points in each direction
        arma::sp_mat V;     //  potential 
        double h;           //  spatial step size

        friend class Simulation;
    public:
        
        // Contructor

        /**
         * Construct the 'Box' object that considers a grid of length 1 in each direction.
         * @param spatial_step_size separation between two aligned points on the grid 
        */
        Box(double spatial_step_size=0.005);

        // Member functions

        // old:
        void set_up_walls(double v0, int Ns=2, double Th=0.02, double wxc=0.5, double wyc=0.5, double wall_piece_length=0.05, double aperture=0.05);

        /**
         * Set up walls to create vertically aligned slits in box.
         * @param num_of_slits number of slits
         * @param v0 barrier potential
         * @param aperture slit opening (vertical)
         * @param wall_width width of each wall (horisontal extent)
         * @param wall_height height of each wall (vertical extent) 
         * @param horisontal_centre horisontal centre position of setup
         * @param vertical_centre vertical centre position of setup
        */
        void create_slits(int num_of_slits, double v0=1e10, double aperture=0.05, double wall_width=0.02, double wall_height=0.05, double horisontal_centre=0.5, double vertical_centre=0.5);
    
};



#endif
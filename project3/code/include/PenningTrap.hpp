#ifndef __PenningTrap_hpp__
#define __PenningTrap_hpp__

#include <cmath>
#include <vector>
#include <iostream>
#include <armadillo>
#include <assert.h>
#include <Particle.hpp>


class PenningTrap{
  private:
    // Private parameters
    double Vd2r; // ratio V/d^2
    std::vector<Particle*> particles;   //  vector to contain pointers all the particle objects
    
    /**
     * Declaration of a set of matrices that we use in the calculations.
     *  rows -> value at (x, y, z)
     *  cols -> particle number 
     */
    arma::mat Q, M;                   //  charges, masses
    arma::mat E_ext, B_ext, E_int;    //  external electric, external magnetic, internal electric fields
    arma::mat F_ext, F_int;           //  external, internal forces
    //  misc :
    arma::mat K1, K2, K3, K4, K;
    arma::mat dist;
    arma::mat norm3;

    /**
     *  rows -> (x, y, z, vx, vy, vz)
     *  cols -> particle number 
     */
    arma::mat RU; arma::mat dRU;      //  positions and velocities, change in positions and velocities


    // Member functions
    arma::mat compute_external_Efield(arma::mat R);         //  compute external E-field
    arma::mat compute_external_Bfield(arma::mat R);         //  compute external B-field
    arma::mat compute_interaction_field(arma::mat R);       //  compute force from particles

    arma::mat external_forces(arma::mat RU);    //  calculate the total forces from external fields
    arma::mat internal_forces(arma::mat RU);    //  calculate the total forces from internal fields

    arma::mat evolve_FE(double dt, arma::mat RU);   //  evolve system for one step in time using the Forward Euler scheme
    arma::mat evolve_RK4(double dt, arma::mat RU);  //  evolve system for one step in time using the Runge-Kutta 4 scheme
    arma::mat K_val(arma::mat RU);    //  helper function for RK4




  public:

    // Public parameters
    double B0;                //  magnetic field strength   [ u μs^(-1) e^(-1) ]
    double V0;                //  applied potential         [ u μm^2 μm^(-2) e^(-1) ]
    double d;                 //  characteristic dimension  [ μm ]
    int Np = 0;               //  number of particles
    bool interactions;        //  whether to include (true) interactions between particles or not (false)
    std::string filename;     //  filename of solution file
    bool isready;             //  status (are things initialised or not)
    
    /**
     * 'system' shall contain information about the position and velocity of each particle at any time
     *      slices      -> time steps
     *      cols        -> particles
     *      row(0)      -> time (t)
     *      rows(1,3)   -> position (x,y,z)
     *      rows(4,6)   -> velocity (vx,vy,vz)
     */
    arma::cube system;

    // Constructor
    PenningTrap(double B0_in, double V0_in, double d_in, bool interactions_in=true);

    // Member functions 
    void add_particle(Particle &p_in);    //  add a particle to the Penning trap
    void ready();                         //  initialise matrices etc. for a given Np
    void simulate(double T, double dt, std::string scheme="RK4", bool point=false);   //  simulate for T μs using time step dt μs using scheme
    void set_solution_filename(std::string filename);   //  define filename of soliution file
    
};

#endif
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
     * rows -> value at (x, y, z)
     * cols -> particle number 
     */
    arma::vec Q, M;
    arma::mat E_ext; arma::mat B_ext; arma::mat E_int;
    arma::mat F_ext; arma::mat F_int;
    arma::mat K1; arma::mat K2; arma::mat K3; arma::mat K4;
    arma::mat RU; arma::mat dRU;

    // Member functions
    arma::mat compute_external_Efield(arma::mat R);         //  compute external E-field
    arma::mat compute_external_Bfield(arma::mat R);         //  compute external B-field
    arma::mat compute_interaction_field(arma::mat R);       //  compute force from particles

    arma::mat external_forces(arma::mat RU);
    arma::mat internal_forces(arma::mat RU);

    arma::mat evolve_FE(double dt, arma::mat RU);
    arma::mat evolve_RK4(double dt, arma::mat RU);
    arma::mat K_val(arma::mat RU);
    


  public:

    // Public parameters
    double B0;                //  magnetic field strength   [ u μs^(-1) e^(-1) ]
    double V0;                //  applied potential         [ u μm^2 μm^(-2) e^(-1) ]
    double d;                 //  characteristic dimension  [ μm ]
    int Np = 0;               //  number of particles
    bool interactions;        //  whether to include (true) interactions between particles or not (false)
    //int dt                   //  time step length

    // Constructor
    PenningTrap(double B0_in, double V0_in, double d_in, bool interactions_in=true);

    // Member functions 
    void add_particle(Particle &p_in);    //  add a particle to the Penning trap
    void simulate(double T, double dt, std::string scheme="RK4", bool point=false);   //  simulate for T μs using time step dt μs using scheme

    
  
    
};

#endif
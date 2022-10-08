#ifndef __PenningTrap_hpp__
#define __PenningTrap_hpp__

#include <cmath>
#include <vector>
#include <iostream>
#include <armadillo>
#include <assert.h>
#include <Particle.hpp>

class PenningTrap{
    // Constructor
  PenningTrap(double B0_in, double V0_in, double d_in);

  double B0;    //  Magnetic field strength
  double V0;    //  Applied potential
  double d;     //  Characteristic dimension
  std::vector<Particle> particles;    //  Vector to contain all the particle objects. 

  // Add a particle to the trap
  void add_particle(Particle p_in);

  // External electric field at point r=(x,y,z)
  arma::vec external_E_field(arma::vec r);  

  // External magnetic field at point r=(x,y,z)
  arma::vec external_B_field(arma::vec r);  

  // Force on particle_i from particle_j
  arma::vec force_particle(int i, int j);

  // The total force on particle_i from the external fields
  arma::vec total_force_external(int i);

  // The total force on particle_i from the other particles
  arma::vec total_force_particles(int i);

  // The total force on particle_i from both external fields and other particles
  arma::vec total_force(int i);
};

#endif
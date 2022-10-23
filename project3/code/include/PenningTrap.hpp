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
    double Vd2r;                      // ratio V/d^2
    std::vector<Particle*> particles; //  vector contain pointers all the particle objects
    
    /**
     * Declaration of a set of matrices that we use in the calculations.
     *  rows -> value at (x, y, z)
     *  cols -> particle number
     */
    arma::mat Q;        //  charges
    arma::mat M;        //  masses
    arma::mat E_ext;    //  external electric field
    arma::mat B_ext;    //  external magnetic field
    arma::mat E_int;    //  internal electric field
    arma::mat F_ext;    //  external forces 
    arma::mat F_int;    //  internal forces
   
   
    arma::mat K1, K2, K3, K4, K;     //  K-values in RK4 scheme
    arma::mat Rnorm;                 //  length of position vectors
    arma::mat dist;                  //  distances from one particle to all others
    arma::mat norm3;                 //  length of vectors to the power of 3

    /**
     * Declaration of the matrices that we use to store all positions and velocities at some time step in the calculations.
     *  rows -> (x, y, z, vx, vy, vz)
     *  cols -> particle number 
     */
    arma::mat RU;       //  positions (rows 0,1,2) and velocities (rows 3,4,5) of particles
    arma::mat dRU;      //  change in positions (rows 0,1,2) and velocities (rows 3,4,5) of particles
    double t=0;         //  current time 


    // Member functions
    /* Compute external E-field at postions {R}. 
        -> result 'E_ext' */
    arma::mat compute_external_Efield(arma::mat R);
    /* Compute external B-field at postions {R}. 
        -> result 'B_ext' */
    arma::mat compute_external_Bfield(arma::mat R);
    /* Compute interaction field at postions {R}. 
        -> result 'E_int' */
    arma::mat compute_interaction_field(arma::mat R);

    /* Calculate the total forces from external fields at positions and velocities {RU}.
      -> result 'F_ext' */
    arma::mat external_forces(arma::mat RU);
    /* Calculate the total forces from internal fields at positions and velocities {RU}.
      -> result 'F_int' */
    arma::mat internal_forces(arma::mat RU); 
    /* Evolve system at positions and velocities {RU} for one step of {dt} μs in time using the forward Euler scheme. */
    arma::mat evolve_FE(double dt, arma::mat RU);
    /* Evolve system at positions and velocities {RU} for one step of {dt} μs in time using the Runge-Kutta 4 scheme. */
    arma::mat evolve_RK4(double dt, arma::mat RU);

    /* Calculate the value of Kj in the RK4-algorithm given positions and velocities {RU}.
      -> K.rows(0,2): dR = U
      -> K.rows(3,6): dV = F(R, U)/m */
    arma::mat K_val(arma::mat RU_);
    /* Find the length of the vectors inside a matrix {R} of 'Np' coloumns and 3 rows.
        -> result is a matrix of equal dimensions where the coloumns have three equal elements */
    arma::mat Pnorm(arma::mat R);

    /* Add a particle with postition {r} and velocity {v} to {list}.
        -> helps 'generate_random_identical_particles()' fill 'particles' with randomised particles */
    void generate_particle(std::vector<Particle*> &list, arma::vec r, arma::vec v); 

    // Extensions for time-dep. potential
    /* Compute external E-field at time {t} and postions {R}. 
        -> result 'E_ext' */
    arma::mat compute_external_Efield(double t, arma::mat R);
    /* Calculate the total forces from external fields at time {t} and positions and velocities {RU}.
      -> result 'F_ext' */
    arma::mat external_forces(double t, arma::mat RU);
    
    /* Evolve system at time {t} and positions and velocities {RU} for one step of {dt} μs in time using the forward Euler scheme. */
    arma::mat evolve_FE(double dt, double t, arma::mat RU);
    /* Evolve system at time {t} and positions and velocities {RU} for one step of {dt} μs in time using the Runge-Kutta 4 scheme. */
    arma::mat evolve_RK4(double dt, double t, arma::mat RU);
    /* Calculate the value of Kj in the RK4-algorithm given a time {t} and positions and velocities {RU}.
      -> K.rows(0,2): dR = U
      -> K.rows(3,6): dV = F(t, R, U)/m */
    arma::mat K_val(double t, arma::mat RU_);


  public:

    // Public parameters
    double B0;                          //  magnetic field strength [ u μs^(-1) e^(-1) ]
    double V0;                          //  applied potential [ u μm^2 μm^(-2) e^(-1) ]
    double d;                           //  characteristic dimension [ μm ]
    double f = 0;                       //  amplitude of time-dependent electric potential
    double omega_V = 0;                 //  applied angular frequency of the time-dependent electric potential [ M s^(-1) ]
    int Np = 0;                         //  number of particles to consider
    bool interactions;                  //  whether to include (true) interactions between particles or not (false)
    bool time_dep = false;              //  whether to use time-dependent potential (true) or not (false)
    bool isready = false;               //  status (are things initialised (true) or not (false))
    
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
    /* Set
      - the magnetic field strength to {B0_in} T,
      - the electric potential magnitude to {V0_in} mV, 
      - the charcteristic dimension to {d_in} μm and
      turn the Coulomb interactions on if {interactions_in}, else off. */
    PenningTrap(double B0_in, double V0_in, double d_in, bool interactions_in=true);

    // Member functions
    /* Add a particle to the Penning trap. */
    void add_particle(Particle &particle); 
    /* Switch Coulomb interaction mode from current to other.
        -> if interactions are on, they are turned off 
        -> if interactions are off, they are turned on */
    void switch_interactions();  
    /* Switch Coulomb interactions {switch_in}.
        -> turn on: switch_in = "on" or "ON"
        -> turn off: switch_in = "off" or "OFF" */
    void switch_interactions(std::string switch_in);
    /* Add a time-dependent perturbation to the applied potential, 
    i.e. set electric potential V0 -> V0*(1+f*cos(ωV*t)), where f = {amplitude} and ωV = {frequency} MHz. */
    void apply_time_dependence(double amplitude, double frequency);
    /* Initialise matrices etc. for a complete set of particles. */
    void ready();       
    /* Simulate the system for {T} μs using time step {dt} μs using {scheme} */
    void simulate(double T, double dt, std::string scheme="RK4", bool point=false);
    /* Save solution to a file named {filename}. */
    void save_solution(std::string filename);
    /* Count the particles that are located within the walls of the Penning trap. 
        -> consider |r|<d as trapped */
    int count_particles();
    /* Generate a set of {no_of_particles} identical particles with normally distributed positions and velocities. */
    void generate_random_identical_particles(double charge, double mass, int no_of_particles, int seed=69); 
    /* Print particles' positions and velocities. */
    void print_particles(); 
    
};

#endif
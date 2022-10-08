#ifndef __Particle_hpp__
#define __Particle_hpp__

#include <cmath>
#include <vector>
#include <iostream>
#include <armadillo>
#include <assert.h>


class Particle{
    private:
    //  Gives PenningTrap acces to Particle
    friend class PenningTrap;
    
    //  Define member variables
    double q_;   //  Charge
    double m_;   //  Mass
    arma::vec r_; //  Position vector
    arma::vec v_; //  Velocity vector

    public:
    //  Constructor
    Particle(int q, int m, arma::vec r, arma::vec  v);


};

#endif 
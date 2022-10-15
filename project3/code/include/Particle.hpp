#ifndef __Particle_hpp__
#define __Particle_hpp__

#include <cmath>
#include <vector>
#include <iostream>
#include <armadillo>
#include <assert.h>


class Particle{
    private:
        //  give PenningTrap access to Particle:
        friend class PenningTrap;
        
        //  define member variables:
        double q_;      //  charge
        double m_;      //  mass
        arma::vec r_;   //  position vector
        arma::vec v_;   //  velocity vector

    public:
        //  Constructor
        Particle(double q, double m, arma::vec r, arma::vec  v);

        //  Functions for returning member  variables. 
        double charge();
        double mass();
        arma::vec position();
        arma::vec velocity();
        
        //  Functions for updating to position and velocity vector
        void new_position(arma::vec r);
        void new_velocity(arma::vec v);

};


#endif 
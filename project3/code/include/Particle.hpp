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
        Particle(double q, double m, arma::vec r, arma::vec  v);

        //  Functions for returning member  variables. 
        double q();
        double m();
        arma::vec r();
        arma::vec v();

        void new_position(arma::vec new_r);

};


#endif 
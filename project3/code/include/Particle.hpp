#ifndef __Particle_hpp__
#define __Particle_hpp__

#include <cmath>
#include <vector>
#include <iostream>
#include <armadillo>
#include <assert.h>


class Particle{
    private:
        
        friend class PenningTrap;   //  give PenningTrap access to Particle:
        
        // Member variables
        double q_;      //  charge [ e ]
        double m_;      //  mass [ u ]
        arma::vec r_;   //  position vector [ μm ]
        arma::vec v_;   //  velocity vector [ μm/μs (??) ]

    public:
        // Constructor
        /* Give the particle
            - charge {q} e,
            - mass equal {m} u,
            - (initial) position {r} μm and 
            - (initial) velocity {v} m/s */
        Particle(double q, double m, arma::vec r, arma::vec  v);

        // Member functions
        /* Get the charge of the particle. */
        double charge();
        /* Get the mass of the particle. */
        double mass();
        /* Get the position of the particle. */
        arma::vec position();
        /* Get the velocity of the particle. */
        arma::vec velocity();


        /* Update the position of the particle to be {r} */
        void new_position(arma::vec r);
        /* Update the velocity of the particle to be {v} */
        void new_velocity(arma::vec v);

};


#endif 
#include <algorithm>
#include <random> 
#include "utils.hpp"
#include <omp.h>

/**
 * PROJECT 4 FYS4150
 *      By: Johan Mylius Kroken, Vetle Vikenes and Nanna Bryne
 */

inline int PBC(int idx, int L){
    /**
     * Return index with periodic boundary condition 
    */
    return (idx + L) % (L);
}

arma::mat make_Lattice(int L=20, bool ordered=false){
    /**
     * Create lattice of size LxL. 
     * Each element is +/- 1, for spin up/down.
     * All spins are aligned (+1) if ordered==true
     * If ordered==false, spins are uniformly distributed 
    */
    arma::mat Lattice;
    if(ordered==false){
        arma::arma_rng::set_seed(69);
        Lattice = arma::randi<arma::mat>(L, L, arma::distr_param(0,1));
        Lattice *= 2;
        Lattice -= 1;
    }
    else{
        // All spins point up 
        Lattice = arma::mat(L, L, arma::fill::ones);
        // Lattice *= -1; // All downwards spins 
    }
    return Lattice;
}

void compute_E_and_M(arma::mat& Lattice, double& E, double& M){
    /**
         * Compute initial energy and magnetization of Lattice  
    */
    int n_spins = Lattice.n_rows; // From Morten's lecture notes. More appropriate name should be used 
    
    for(int x=0; x<n_spins; x++){
        for(int y=0; y<n_spins; y++){
            M += (double) Lattice(x,y);
            E -= (double) Lattice(x,y) * 
                (Lattice(x, PBC(y+1,n_spins)) + Lattice(PBC(x+1,n_spins),y));
        }
    }
}

std::vector<double> DeltaE(double T){
    /**
     * Vector containing all five possible values of Delta E,
     * due to flipping a single spin in Lattice.  
    */

    double beta = 1 / T;
    std::vector<double> dE(5);
    for(int i=0; i<5; i++){
        dE[i] = exp(-beta * (4 * (i - 2)));
    }

    dE[2] = 0; // For configurations that leave energy unchanged 
    return dE;
}



void MC_cycle(arma::mat& Lattice, const std::vector<double> dE, double&E, double& M, unsigned int seed){
    /**
     * Perform one MC cycle for a given Lattice of size (LxL)
     * One cycle consists of N=L*L attempted spin flips.
     * Each spin flip attempt is performed on a random site in the Lattice  
    */
    std::mt19937_64 generator(seed);
    int L = Lattice.n_rows;

    for(int x=0; x<L; x++){
        for(int y=0; y<L; y++){
            // Pick site to consider 
            int ix = (unsigned int) generator() % L;
            int iy = (unsigned int) generator() % L;

            // Compute energy change if site is flipped. 
            int deltaE = 2*Lattice(ix,iy) * 
                        (Lattice(ix, PBC(iy+1, L)) + 
                        Lattice(PBC(ix+1, L), iy) + 
                        Lattice(ix, PBC(iy-1, L)) + 
                        Lattice(PBC(ix-1, L), iy));

            if(deltaE <= 0){
                // Accept flip, update E and M. 
                Lattice(ix,iy) *= -1;
                M += (double) 2*Lattice(ix,iy);
                E += (double) deltaE;}
    
            else{
                // Compute Boltzmann factor from energy change 
                // Accept flip if factor is larger than a random numer 0<r<1. 
                int dE_idx = deltaE/4 + 2; 
                double w = dE[dE_idx];
                std::uniform_real_distribution<double> r(0,1);
                if(r(generator) <= w){
                    Lattice(ix,iy) *= -1;
                    M += (double) 2*Lattice(ix,iy);
                    E += (double) deltaE;
                }
            }
        }
    }
}

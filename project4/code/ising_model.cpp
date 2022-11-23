#include <algorithm>
#include <random> 
#include "utils.hpp"
// #include <omp.h>

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
        arma::arma_rng::set_seed_random(); // (seed); // CHECK THIS LATER!!!!
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

double initial_E(arma::mat Lattice){
    /**
         * Compute initial energy and magnetization of Lattice  
    */
    double E = 0;
    int L = Lattice.n_rows; 
    
    for(int x=0; x<L; x++){
        for(int y=0; y<L; y++){
            E -= Lattice(x,y) * 
                (Lattice(x, PBC(y+1,L)) + Lattice(PBC(x+1,L),y));
        }
    }
    return E;
}

double initial_M(arma::mat Lattice){
    double M = 0;
    int L = Lattice.n_rows;
    for(int x=0; x < L; x++){
        for(int y=0; y < L; y++){
            M += Lattice(x,y);
        }
    }
    return M;
}


std::vector<double> Boltzmann(double T){
    /**
     * Vector containing all five possible values of Delta E,
     * due to flipping a single spin in Lattice.  
    */

    double beta = 1 / T;
    std::vector<double> dE(5);
    for(int i=0; i<5; i++){
        dE[i] = exp(-beta * (4 * (i - 2)));
    }

    return dE;
}



void MC_cycle(arma::mat& Lattice, std::vector<double> Boltz, double& E, double& M, 
            std::mt19937_64& gen){
    /**
     * Perform one MC cycle for a given Lattice of size (LxL)
     * One cycle consists of N=L*L attempted spin flips.
     * Each spin flip attempt is performed on a random site in the Lattice  
    */

    int L = Lattice.n_rows;
    std::uniform_int_distribution<int> indx(0, L-1);
    std::uniform_real_distribution<double> r(0.0, 1.0);


    for(int xy=0; xy<L*L; xy++){
        // Pick site to consider 
        int ix =  indx(gen);
        int iy =  indx(gen);

        // Compute energy change if site is flipped. 
        int deltaE = 2*Lattice(ix,iy) * 
                    (Lattice(ix, PBC(iy+1, L)) + 
                    Lattice(PBC(ix+1, L), iy) + 
                    Lattice(ix, PBC(iy-1, L)) + 
                    Lattice(PBC(ix-1, L), iy));

        int dE_idx = deltaE/4 + 2; 
        double w = Boltz[dE_idx];

        if(r(gen) <= w){
            Lattice(ix,iy) *= -1;
            M += 2*Lattice(ix,iy);
            E += deltaE;
        }
    }
}

arma::rowvec run_MC(int L, double T, int N_samples, int N_eq, 
                    unsigned int seed){
    
    std::mt19937_64 gen(seed);
    arma::rowvec averages(7, arma::fill::zeros);
    
    arma::mat Lattice = make_Lattice(L);
    double E = initial_E(Lattice);
    double M = initial_M(Lattice);

    std::vector<double> Boltzmann_ = Boltzmann(T);
    int N_tot_cycles = N_samples + N_eq; 

    double varE, varM;


    
    for(int cycle=1; cycle <= N_tot_cycles; cycle++){
        MC_cycle(Lattice, Boltzmann_, E, M, gen);
        if(cycle > N_eq){
            averages(1) += E; 
            averages(2) += E*E;
            averages(3) += abs(M);
            averages(4) += M*M; 
        }
    } 

    averages /= (double) (N_samples);
    varE = averages(2) - averages(1)*averages(1);
    varM = averages(4) - averages(3)*averages(3);

    averages(0) = T;
    averages(5) = varE;
    averages(6) = varM;

    return averages;
}

arma::mat run_MC_cumulative(int L, double T, int N_samples, int N_eq,
                    bool ordered, unsigned int seed){

    std::mt19937_64 gen(seed);
    
    arma::mat Lattice = make_Lattice(L, ordered);
    double E = initial_E(Lattice);
    double M = initial_M(Lattice);

    std::vector<double> Boltzmann_ = Boltzmann(T);
    int N_tot_cycles = N_samples + N_eq; 
    arma::mat samples(N_samples, 4, arma::fill::zeros);

    if(N_eq==0){
        // Used for estimating equilibriation time 
        // samples.insert_rows(0,1); // Add extra row on top 
        // Store initial values 
        N_tot_cycles -= 1;
    }

    for(int cycle=1; cycle <= N_tot_cycles; cycle++){
        MC_cycle(Lattice, Boltzmann_, E, M, gen);
        if(cycle > N_eq){
            samples(cycle, 0) = E; 
            samples(cycle, 1) = abs(M);
            samples(cycle, 2) = cycle;

        }
    } 
    samples.col(3) = arma::colvec (N_samples, arma::fill::value(T)); 

    samples.col(0) /= (double) (L*L);
    samples.col(1) /= (double) (L*L);


    return samples;
}
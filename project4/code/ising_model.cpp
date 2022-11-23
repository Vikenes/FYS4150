#include <algorithm>
#include <random> 
#include "utils.hpp"
// #include <omp.h>

/**
 * PROJECT 4 FYS4150
 *      By: Johan Mylius Kroken, Vetle Vikenes and Nanna Bryne
 */

/**
 * Find index along an axis of Lattice with periodic boundary conditions.
 * 
 * @param idx Index of Lattice axis.  
 * @param L length of Lattice. 
 * @return Position index of lattice in the range [0, L-1].
*/
inline int PBC(int idx, int L){
    return (idx + L) % (L);
}

/**
 * Create lattice of size (LxL) with values +1/-1, for spin up/down.
 * @param L Length of Lattice in each direction.
 * @param ordered If true: All spins are aligned (+1). 
 * If false: spins are randomly aligned.
 * @return (LxL) matrix of spins.  
*/
arma::mat make_Lattice(int L=20, bool ordered=false){
    arma::mat Lattice;
    if(ordered==false){
        // Random spin orientations
        arma::arma_rng::set_seed_random(); 
        Lattice = arma::randi<arma::mat>(L, L, arma::distr_param(0,1));
        Lattice *= 2;
        Lattice -= 1;
    }
    else{
        // All spins point up 
        Lattice = arma::mat(L, L, arma::fill::ones);
        // Lattice *= -1; // All downwards spins. Might add if desired.  
    }
    return Lattice;
}

/**
 * Compute initial energy of Lattice. 
 * @param Lattice (LxL) matrix of spins.
 * @return Total energy of Lattice    
*/
double initial_E(arma::mat Lattice){
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

/**
 * Compute initial magnetization of Lattice.
 * @param Lattice (LxL) matrix of spins. 
 * @return Total magnetization of Lattice
*/
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


/**
 * Compute the Boltzmann factors of the Lattice from a single spin flip.
 * A single spin flip yields an energy change of:
 * Delta E = {-8,-4,0,4,8}
 * 
 * @param T Temperature. 
 * 
 * @return Vector of five Boltzmann factors: exp(-Delta E / T). 
*/
std::vector<double> Boltzmann(double T){

    double beta = 1 / T;
    std::vector<double> dE(5);
    for(int i=0; i<5; i++){
        dE[i] = exp(-beta * (4 * (i - 2)));
    }

    return dE;
}



/**
 * Perform one MC cycle on a Lattice of size (LxL)
 * Attempts spin flip on a random spin in the lattice
 * One cycle consists of L*L attempted spin flips.
 * 
 * @param Lattice Adress of the (LxL) Lattice
 * @param Boltz Vector containing possible Boltzmann factor due to single spin flip
 * @param E Energy of Lattice
 * @param M Magnetization of Lattice
 * @param gen mt19937_64 generator for generating random numbers   
*/
void MC_cycle(arma::mat& Lattice, std::vector<double> Boltz, double& E, double& M, 
            std::mt19937_64& gen){

    int L = Lattice.n_rows;
    std::uniform_int_distribution<int> indx(0, L-1);    // Uniform distr of spin idices
    std::uniform_real_distribution<double> r(0.0, 1.0); // Acceptance prob. 


    for(int xy=0; xy<L*L; xy++){
        // Pick a spin at random  
        int ix =  indx(gen);
        int iy =  indx(gen);

        // Compute energy change if the spin is flipped. 
        int deltaE = 2*Lattice(ix,iy) * 
                    (Lattice(ix, PBC(iy+1, L)) + 
                    Lattice(PBC(ix+1, L), iy) + 
                    Lattice(ix, PBC(iy-1, L)) + 
                    Lattice(PBC(ix-1, L), iy));

        int dE_idx = deltaE/4 + 2; 
        double w = Boltz[dE_idx];  // Boltzmann factor of spin flip

        if(r(gen) <= w){
            // Flip spin. Update lattice
            // Update energy and magnetization
            Lattice(ix,iy) *= -1;
            M += 2*Lattice(ix,iy);
            E += deltaE;
        }
    }
}

/**
 * Compute average quantities from MC sampling 
 * by sampling after equilibriating the system
 * 
 * 
 * @param L Number of spins in each direction.
 * @param T Temperature of the system 
 * @param N_samples Number of samples over which average values are computed 
 * @param N_eq Number of steps used to equilibriate the system 
 * @param seed Random seed used for the MC cycles. 
 * 
 * @return Temperature, average of E, E^2, abs(M) and M^2. Variance of E and abs(M)
*/
arma::rowvec run_MC(int L, double T, int N_samples, int N_eq, unsigned int seed){


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


/**
 * Runs multiple MC cycles, and stores accumulated quantities from each MC cycle.
 * @param L Length of lattice
 * @param T Temperature 
 * @param N_samples Number of runs from which samples are stored 
 * @param N_eq Number of equilibriation steps to perform before storing values.  
 * @param ordered All spins aligned if true, oriented randomly if false.
 * @param seed Seed for random generator 
 * 
 * @return Energy and abs(Magnetization) per spin at each MC cycle and temperature.
 * Samples during equilibriation are not stored.  
*/
arma::mat run_MC_cumulative(int L, double T, int N_samples, int N_eq,
                    bool ordered, unsigned int seed){

    std::mt19937_64 gen(seed);
    
    arma::mat Lattice = make_Lattice(L, ordered);
    double E = initial_E(Lattice);
    double M = initial_M(Lattice);

    std::vector<double> Boltzmann_ = Boltzmann(T);
    int N_tot_cycles = N_samples + N_eq; 
    arma::mat samples(N_samples, 4, arma::fill::zeros);


    int sample_idx = 0;
    for(int cycle=1; cycle <= N_tot_cycles; cycle++){
        MC_cycle(Lattice, Boltzmann_, E, M, gen);
        if(cycle > N_eq){
            samples(sample_idx, 0) = E; 
            samples(sample_idx, 1) = abs(M);
            samples(sample_idx, 2) = cycle;
            sample_idx++;

        }
    } 
    samples.col(3) = arma::colvec (N_samples, arma::fill::value(T)); 

    samples.col(0) /= (double) (L*L);
    samples.col(1) /= (double) (L*L);


    return samples;
}
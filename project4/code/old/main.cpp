#include <algorithm>
#include <random> 
#include "utils.hpp"
#include <omp.h>

/**
 * PROJECT 4 FYS4150
 *      By: Johan Mylius Kroken, Vetle Vikenes and Nanna Bryne
 */

/**
 * Return index with periodic boundary condition 
*/
inline int PBC(int idx, int L){
    return (idx + L) % (L);
}

/**
 * Create lattice of size LxL. 
 * Each element is +/- 1, for spin up/down.
 * All spins are aligned (+1) if ordered==true
 * If ordered==false, spins are uniformly distributed 
*/
arma::mat make_Lattice(int L=20, bool ordered=false);


/**
 * Vector containing all five possible values of Delta E,
 * due to flipping a single spin in Lattice.  
*/
std::vector<double> DeltaE(double T);


/**
     * Compute initial energy and magnetization of Lattice  
*/
void compute_E_and_M(arma::mat& Lattice, double& E, double& M);



/**
 * Perform one MC cycle for a given Lattice of size (LxL)
 * One cycle consists of N=L*L attempted spin flips.
 * Each spin flip attempt is performed on a random site in the Lattice  
*/
void MC_cycle(arma::mat& Lattice, const std::vector<double> dE, double&E, double& M, unsigned int seed);



/**
 * Compute E, E^2, abs(M), M^2 as a function of MC cycles for a 2x2 Lattice. 
 * Used for comparing to analytical result. 
*/
void test_analytical(int N_MC_cycles);




/**
 * Study evolution of eps and abs(m) as a function of number of MC cycles. 
 * Used to determine number of cycles needed to equilibriate system. 
 * ordered==false: Random initial spin orientations 
 * ordered==true : All spins initially aligned
*/
void equilibriation_time(int N_MC_cycles, double T0, bool ordered);




/**
 * Sample avg eps and abs(m).
 * Equilibriate system before measuring. 
 * Makes N_samples of avg(eps) and avg(abs(m)), where.
 * average values are computed over MC_cycles_per_sample.
*/
void probability_distr(int N_samples, int equilibriation_steps, double T0, bool ordered);




int main(int argc, char* argv[]){   

    if (argc != 5){ // Expect 4 cmd-line args 
        std::string executable_name = argv[0];
        std::cerr << "Error: Wrong number of input arguments. Got " << argc << std::endl;
        std::cerr << "Usage: " << executable_name << 
                    " <Num. of MC cycles> <Num. of equilibriation steps> <Method> <filename>" << std::endl;
        
        return 1;
    }

    int N_cycles = atoi(argv[1]);
    int equi_steps = atoi(argv[2]);
    std::string method(argv[3]);
    std::string filename(argv[4]);

    auto start_time = std::chrono::high_resolution_clock::now();

    if(method=="analytical"){test_analytical()}
    // test_analytical(int(1e6));


    // equilibriation_time(100000, 1, false);
    // equilibriation_time(100000, 1, true);
    // equilibriation_time(100000, 2.4, false);
    // equilibriation_time(100000, 2.4, true);

    // probability_distr(100000, 10000, 1, false);
    // probability_distr(100000, 10000, 2.4, false);

    // testing parallel
    // #pragma omp parallel
    // {
    // int ID = omp_get_thread_num();

    // std::cout << "hello(" + std::to_string(ID) + "), ";
    // std::cout << "world(" + std::to_string(ID) + "), " << std::endl;

    // }
    

    auto stop_time = std::chrono::high_resolution_clock::now();

    double run_time = std::chrono::duration<double>(stop_time-start_time).count();

    std::cout << "Running time: " << run_time << " s (" << run_time/60.0 << " min)" << std::endl;

    return 0;
}


void MC_cycle(arma::mat& Lattice, const std::vector<double> dE, double&E, double& M, unsigned int seed){
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

void probability_distr(int N_samples, int equilibriation_steps, double T0, bool ordered){
    std::cout << "Creating N=" << N_samples << " samples at T=" << float_to_string(T0);

    std::mt19937 generator(690);
    arma::mat Lattice = make_Lattice(20, ordered);
    std::vector<double> dE = DeltaE(T0);
    std::vector<double> cycle_avg(2);
    std::vector<double> eps_samples(N_samples);
    std::vector<double> eps2_samples(N_samples);

    double E = 0;
    double M = 0;
    double N_spins = Lattice.n_rows * Lattice.n_cols; 

    compute_E_and_M(Lattice, E, M);

    for(int step=0; step <= equilibriation_steps; step++){
        // Equilibriate system. 
        MC_cycle(Lattice, dE, E, M, generator());
    }

    // Compute samples 
    for(int sample=0; sample < N_samples; sample++){
        MC_cycle(Lattice, dE, E, M, generator());

        eps_samples[sample] = E / N_spins; // avg. E per spin 
        eps2_samples[sample] = E*E / N_spins / N_spins; // avg. E^2 per spin 

        if(sample % (N_samples / 10) == 0){
            std::cout << "  N=" << sample << " samples completed" << std::endl;
        }
    }

    std::cout << "  All N=" << N_samples << " samples completed." << std::endl;

    std::string fname = "sample_eps_L20_T" + float_to_string(T0) + "_";
    if(ordered==false){fname += "un";}
    fname += "ordered";

    std::cout << "Writing output to: " << fname << ".txt" << std::endl << std::endl;

    std::string ofile_header = "eps, eps_squared";

    write_to_file_scientific(eps_samples, eps2_samples, fname, ofile_header);

}



void equilibriation_time(int N_MC_cycles, double T0, bool ordered){
    std::mt19937 generator(690);

    arma::mat Lattice = make_Lattice(20, ordered);
    std::vector<double> dE = DeltaE(T0);
    std::vector<double> average(2, 0.0);

    double E = 0;
    double M = 0;
    double N_spins = Lattice.n_rows * Lattice.n_cols; 

    compute_E_and_M(Lattice, E, M);

    // Stupid file writing method. 
    std::string path = "../output/data/";
    std::stringstream ss; 
    std::string name = "equilibriate_L20_T" + float_to_string(T0) + "_";
    if (ordered==false){name += "un";}
    name += "ordered";
    std::string fname = path + name + ".txt";

    std::ofstream avg_values;
    avg_values.open(fname.c_str());
    avg_values << "step, eps, abs_m" << std::endl;

    
    for(int cycle=1; cycle <= N_MC_cycles; cycle++){
        MC_cycle(Lattice, dE, E, M, generator());
        average[0] += E / N_spins;
        average[1] += abs(M) / N_spins; 
        avg_values << cycle << ", "
                    << scientific_format(average[0]/cycle) << ", "
                    << scientific_format(average[1]/cycle) << std::endl;
    }
    avg_values.close();
    std::cout << "Completed " << N_MC_cycles << " Monte Carlo cycles" << std::endl;
    std::cout << "  Wrote to file: " << name << std::endl;
}


void test_analytical(int N_MC_cycles){
    std::mt19937 generator(690);
    arma::mat Lattice = make_Lattice(2);
    std::vector<double> dE = DeltaE(1);

    std::vector<double> average(4, 0.0);


    double E = 0;
    double M = 0;

    compute_E_and_M(Lattice, E, M);

    // Stupid method for writing to file. Will be improved later. 
    std::string path = "../output/data/";
    std::string name = "test_analytical_" + std::to_string(N_MC_cycles) + "_cycles";
    std::string fname = path + name + ".txt";

    std::ofstream avg_values;  
    avg_values.open(fname.c_str()); // , std::ios_base::app | std::ios_base::in);
    avg_values << "step, E, E2, absM, M2" << std::endl;

    for(int cycle=1; cycle <= N_MC_cycles; cycle++){
        MC_cycle(Lattice, dE, E, M, generator());
        average[0] += E; average[1] += E*E;
        average[2] += abs(M); average[3] += M*M;

        avg_values << cycle << ", ";
        avg_values << scientific_format(average[0]) << ", "
                    << scientific_format(average[1]) << ", "
                    << scientific_format(average[2]) << ", "
                    << scientific_format(average[3]) << std::endl;
    }

    avg_values.close();
    
    std::cout << "Completed " << N_MC_cycles << " Monte Carlo cycles" << std::endl;
}


arma::mat make_Lattice(int L, bool ordered){
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

std::vector<double> DeltaE(double T){

    double beta = 1 / T;
    std::vector<double> dE(5);
    for(int i=0; i<5; i++){
        dE[i] = exp(-beta * (4 * (i - 2)));
    }

    dE[2] = 0; // For configurations that leave energy unchanged 
    return dE;
}

void compute_E_and_M(arma::mat& Lattice, double& E, double& M){
    int n_spins = Lattice.n_rows; // From Morten's lecture notes. More appropriate name should be used 
    
    for(int x=0; x<n_spins; x++){
        for(int y=0; y<n_spins; y++){
            M += (double) Lattice(x,y);
            E -= (double) Lattice(x,y) * 
                (Lattice(x, PBC(y+1,n_spins)) + Lattice(PBC(x+1,n_spins),y));
        }
    }
}
#ifndef __Simulation_hpp__
#define __Simulation_hpp__

#include <utils.hpp>
#include <Box.hpp>


class Simulation{

    private:

        // Member variables

        double dt_;
        double T_;
        int Nt_;
        Box box_;

        double xc_, yc_;            //  gaussian centre position
        double sigmax_, sigmay_;    //  gaussian extent (stdv.)
        double px_, py_;            //  gaussian momentum (p-vec)

        // initial states:
        arma::cx_mat U0_;    //  initial state as matrix
        arma::cx_vec u_;     //  initial state as vector
        
        // matrices A and B:
        arma::sp_cx_mat A_, B_;   // matrix A, matrix B

        // system:
        arma::cx_cube U_;    // cube to store the U matrices in time 


        // Member functions
        /** 
         * @brief Set up an initial 2D Gaussian wave packet
         * @param extent sigma-vector in Gaussian
         * @param centre_position xc-vector in Gaussian
         * @param momentum p-vector in Gaussian
        */
        void gaussian_wavepacket(std::tuple<double,double> extent, std::tuple<double,double> centre_position, std::tuple<double,double> momentum);
    
    public:

        // Member variables
        std::tuple<double, double, double, double, double, double> gaussian_params_;    //  parameters in the Gaussian 

        // Constructor
        /**
         * @brief Construct simulation object.
         * @param duration simulation duration
         * @param time_step_size temporal step size
         * @param gaussian_extent sigma-vector in Gaussian
         * @param gaussian_centre_position xc-vector in Gaussian
         * @param gaussian_momentum p-vector in Gaussian
        */
        Simulation(Box box, double duration=0.008, double time_step_size=2.5e-5, std::tuple<double,double> gaussian_extent=std::make_tuple(0.05,0.05), std::tuple<double,double> gaussian_centre_position=std::make_tuple(0.25,0.50), std::tuple<double,double> gaussian_momentum=std::make_tuple(200.0,0.0));
        

        // Member functions
        /**
         * @brief Initialise the simulation.
        */
        void initialise();
        /**
         * @brief Run the simulation. 
         * @note Solves system using Armadillo's spsolve.
         * @returns arma::cx_cube 'U' for which a slice n contains the solutions u_{i,j} for time point n
        */
        arma::cx_cube run_simulation();
        /**
         * @brief Easily extend the initial Gaussian in x- and/or y-direction.
         * @param vertical_extent standard deviation of Gaussian in y-direction
         * @param horisontal_extent standard deviation of Gaussian in x-direction
        */
        void extend_wavepacket(double vertical_extent=0.05, double horisontal_extent=0.05);
};

#endif
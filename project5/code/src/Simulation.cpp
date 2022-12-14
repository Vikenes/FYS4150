
#include "Simulation.hpp"


Simulation::Simulation(Box box, double duration, double time_step_size, std::tuple<double,double> gaussian_extent, std::tuple<double,double> gaussian_centre_position, std::tuple<double,double> gaussian_momentum){
    dt_ = time_step_size;
    T_ = duration;
    Nt_ = int(T_/dt_)+1;
    box_ = box;
    
    //  set parameters in Gaussian:
    gaussian_wavepacket(gaussian_extent, gaussian_centre_position, gaussian_momentum);

    //  set up matices and cube:
    U0_ = arma::cx_mat(box_.M, box_.M);
    A_ = arma::sp_cx_mat((box_.M-2)*(box_.M-2), (box_.M-2)*(box_.M-2));
    B_ = arma::sp_cx_mat((box_.M-2)*(box_.M-2), (box_.M-2)*(box_.M-2));
    U_ = arma::cx_cube(box_.M, box_.M, Nt_);
}

void Simulation::gaussian_wavepacket(std::tuple<double,double> extent, std::tuple<double,double> centre_position, std::tuple<double,double> momentum){
    std::tie(sigmax_, sigmay_) = extent;
    std::tie(xc_, yc_) = centre_position;
    std::tie(px_, py_) = momentum;
    //  update all member vairables:
    gaussian_params_ = std::make_tuple(xc_, yc_, sigmax_, sigmay_, px_, py_);
}

void Simulation::extend_wavepacket(double vertical_extent, double horisontal_extent){
    //  update parameters:
    gaussian_wavepacket(std::make_tuple(horisontal_extent, vertical_extent), std::make_tuple(xc_, yc_), std::make_tuple(px_, py_));
}

void Simulation::initialise(){
    
    //  initial state:
    std::cout << " M = " << box_.M << std::endl;
    std::cout << "Nt = " << Nt_ << std::endl;
    std::cout << "Initialising state ..." << std::endl;
    
    initialise_state(U0_, box_.M, box_.h, xc_, yc_, sigmax_, sigmay_, px_, py_); //  from 'utils.cpp'
    u_ = make_column_vector(U0_, box_.M); //  from 'utils.cpp'

    std::cout << "-> State initialised" << std::endl;
    
    //  set up A and B matrices:
    fill_AB_matrix(box_.M, box_.h, dt_, box_.V, A_, B_);  //  from 'utils.cpp'
    U_.slice(0) = U0_;
    std::cout << "-> U-cube set up" << std::endl;
}


arma::cx_cube Simulation::run_simulation(){

    std::cout << "Running simulation:" << std::endl;

    auto t1 = std::chrono::high_resolution_clock::now();
    
    //  solve for each time step:
    for(int n=0; n<Nt_-1; n++){
        std::cout << "-> " << n+1 << " / " << Nt_-1 << std::endl;
        arma::cx_vec bvec = B_*u_;
        u_ = arma::spsolve(A_, bvec);
        U_.slice(n+1) = make_matrix(u_, box_.M);   //  from 'utils.cpp'
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    double duration_seconds = std::chrono::duration<double>(t2-t1).count();
    
    //  update runner:
    std::cout << "Simulation finished" << std::endl;
    std::cout << "Run time: " << duration_seconds << "s (" << duration_seconds / 60.0 << " mins)\n\n" << std::endl<< std::endl;

    return U_;
}
// #include "Box.hpp"
#include "Simulation.hpp"
// #include "utils.hpp"

Simulation::Simulation(Box BBXa, double Dta, double Ta, double xca, double sigma_xa, double p_xa, double yca, double sigma_ya, double p_ya){
    Dt = Dta;
    T = Ta;
    xc = xca;
    sigma_x = sigma_xa;
    p_x = p_xa;
    yc = yca;
    sigma_y = sigma_ya;
    p_y = p_ya;
    Nt = int(Ta/Dta)+1;
    BBX = BBXa;

    U0 = arma::cx_mat(BBX.M, BBX.M);
    A = arma::sp_cx_mat((BBX.M-2)*(BBX.M-2), (BBX.M-2)*(BBX.M-2));
    B = arma::sp_cx_mat((BBX.M-2)*(BBX.M-2), (BBX.M-2)*(BBX.M-2));
    U = arma::cx_cube(BBX.M, BBX.M, Nt);
}

void Simulation::initialise(void){
    // Initial state:
    std::cout<<"M: "<< BBX.M <<std::endl;
    std::cout<<"Nt: "<< Nt <<std::endl;
    std::cout<<"Initialising state..."<<std::endl;
    initialise_state(U0, BBX.M, BBX.h, xc, yc, sigma_x, sigma_y, p_x, p_y);
    u = make_column_vector(U0, BBX.M);
    std::cout<<"-> State initialised"<<std::endl;
    //  Set up A and B matrices:
    fill_AB_matrix(BBX.M, BBX.h, Dt, BBX.V, A, B);
    std::cout<<"-> A and B set up"<<std::endl;
    U.slice(0) = U0;
    std::cout<<"-> U-cube set up"<<std::endl;
}

arma::cx_cube Simulation::run_simulation(void){
    std::cout<<"Running simulation:"<<std::endl;
    auto t1 = std::chrono::high_resolution_clock::now();
    for(int n=0; n<Nt-1; n++){
        std::cout<<"-> "<<n+1<<" / "<<Nt-1<<std::endl;
        arma::cx_vec bvec = B*u;
        u = arma::spsolve(A,bvec);
        U.slice(n+1) = make_matrix(u, BBX.M);
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    double duration_seconds = std::chrono::duration<double>(t2-t1).count();
    std::cout<<"Simulation successfully finished"<<std::endl;
    std::cout<<"Run time: "<<duration_seconds<< "s, or " << duration_seconds / 60.0<<" minutes\n\n"<<std::endl;
    return U;
}
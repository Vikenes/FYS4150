// #include "Box.hpp"
#include "Simulation.hpp"
// #include "utils.hpp"

Simulation::Simulation(Box BBX, double Dt, double T, double xc, double sigma_x, double p_x, double yc, double sigma_y, double p_y){
    Dt = Dt;
    T = T;
    xc = xc;
    sigma_x = sigma_x;
    p_x = p_x;
    yc = yc;
    sigma_y = sigma_y;
    p_y = p_y;
    Nt = int(T/Dt)+1;
    BBX = BBX;

    U0 = arma::cx_mat(BBX.M, BBX.M);
    A = arma::sp_cx_mat((BBX.M-2)*(BBX.M-2), (BBX.M-2)*(BBX.M-2));
    B = arma::sp_cx_mat((BBX.M-2)*(BBX.M-2), (BBX.M-2)*(BBX.M-2));
    U = arma::cx_cube(BBX.M, BBX.M, Nt);
}

void Simulation::initialise(void){
    // Initial state:
    initialise_state(U0, BBX.M, BBX.h, xc, yc, sigma_x, sigma_y, p_x, p_y);
    u = make_column_vector(U0, BBX.M);

    //  Set up A and B matrices:
    fill_AB_matrix(BBX.M, BBX.h, Dt, BBX.V, A, B);

    U.slice(0) = U0;
}

arma::cx_cube Simulation::run_simulation(void){
    for(int n=0; n<Nt-1; n++){
        std::cout<<"Time step: "<<n+1<<" of "<<Nt-1<<std::endl;

        arma::cx_vec bvec = B*u;
        u = arma::spsolve(A,bvec);
        U.slice(n+1) = make_matrix(u, BBX.M);
    }
    return U;
}
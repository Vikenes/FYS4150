// Include stuff
#include <iostream>
#include "utils.hpp"
#include "Box.hpp"
#include "Simulation.hpp"


/**
 * @brief Main (testing) program for FYS4150 project5 code.
 * Vetle A. Vikenes, Nanna Bryne, Johan Mylius Kroken
 * 
 */

void no_slit(void){
    // double h = 0.005;
    // double Dt = 2.5e-5;
    // double T = 0.008;
    // double xc = 0.25;
    // double sigma_x = 0.05;
    // double p_x = 200;
    // double yc = 0.5;
    // double sigma_y = 0.05;
    // double p_y = 0;
    double v0 = 0;
    Box b1 = Box(0.005);
    b1.set_up_walls(v0);
    // Simulation s1 = Simulation(b1);
    // s1.initialise();
    
    // std::cout<<"Running no slit experiment"<<std::endl;
    // arma::cx_cube U = s1.run_simulation();
    // U.save("../output/binfiles/no_slit_arma_cube.bin");
}

// void double_slit_broad_sigma_y(void){
//     double h = 0.005;
//     double Dt = 2.5e-5;
//     double T = 0.008;
//     double xc = 0.25;
//     double sigma_x = 0.05;
//     double p_x = 200;
//     double yc = 0.5;
//     double sigma_y = 0.10;  //  Broader sigma
//     double p_y = 0;
//     double v0 = 1e10;   // High potential to set up the slits

//     std::cout<<"Running double slit experiment"<<std::endl;
//     arma::cx_cube U = simulation(h, Dt, T, xc, sigma_x, p_x, yc, sigma_y, p_y, v0);
//     U.save("../output/binfiles/double_slit_arma_cube.bin");
// }



int main(){

    // //  Preliminary testing
    // std::complex<double> z = 5;
    // std::complex<double> r = std::complex<double>(1,1);
    // int M = 5;
    // arma::cx_vec cvec = arma::cx_vec((M-2)*(M-2)).fill(z);
    // std::cout << get_AB_matrix(M, cvec, r) << std::endl;
    // int matrix_size = (M-2)*(M-2);
    // arma::sp_cx_mat AB = arma::sp_cx_mat(matrix_size, matrix_size);
    // std::cout << AB << std::endl;
    // get_AB_matrix(AB, M, cvec, r);
    // std::cout << AB << std::endl;
    // // std::cout << size(AB)(0) << std::endl;

    // double z = 5;
    // double r = 1;
    // int M = 6;
    // arma::vec vec = arma::vec((M-2)*(M-2)).fill(z);
    // arma::mat A = get_AB_matrix(M, vec, -r);
    // arma::mat B = get_AB_matrix(M, vec, r);
    // std::cout << "A:" << std::endl;
    // std::cout << A << std::endl;

    // std::cout << "B:" << std::endl;
    // std::cout << B << std::endl;

    // Testing complex number properties

    // std::complex<double> z1 = std::complex<double>(3.0, 5.0);
    // std::complex<double> z2 = std::complex<double>(7.0, 1.0);

    // std::cout << z1 << "  " << z2 << std::endl;
    // std::cout << sqrt(real(conj(z1)*z1)) << std::endl;
    // std::cout << sqrt(norm(z1)) << std::endl;

    // int M = 100;
    // double h = 1.0/(M-1);
    // double v0 = 1.0;
    // std::cout<<h<<std::endl;
    // arma::sp_mat V = arma::sp_mat(M,M);
    // std::cout << V << std::endl;
    // set_up_walls(V, v0, M, h);
    // std::cout << V << std::endl;

    no_slit();
    // double_slit_broad_sigma_y();


    // Testing index functions
    // int M = 7;
    // int i = 5;
    // int j = 2;
    // int k = idx_k(i,j,M);
    // int ii, jj;
    // std::tie(ii,jj) = idx_ij(k,M);

    // std::cout << k << std::endl;
    // std::cout << ii << std::endl;
    // std::cout << jj << std::endl;
}

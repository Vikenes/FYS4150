// Include stuff
#include <iostream>
#include <utils.hpp>


/**
 * @brief Main (testing) program for FYS4150 project5 code.
 * Vetle A. Vikenes, Nanna Bryne, Johan Mylius Kroken
 * 
 */

// void simulation()



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

    int M = 100;
    double h = 1.0/(M-1);
    double v0 = 1.0;
    std::cout<<h<<std::endl;
    arma::sp_mat V = arma::sp_mat(M,M);
    std::cout << V << std::endl;
    set_up_walls(V, v0, M, h);
    std::cout << V << std::endl;
}

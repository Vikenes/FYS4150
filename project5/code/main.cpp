// Include stuff
#include <iostream>
#include <utils.hpp>


/**
 * @brief Main (testing) program for FYS4150 project5 code.
 * Vetle A. Vikenes, Nanna Bryne, Johan Mylius Kroken
 * 
 */


int main(){

    //  Preliminary testing
    std::complex<double> z = 5;
    std::complex<double> r = std::complex<double>(1,1);
    int M = 5;
    arma::cx_vec cvec = arma::cx_vec((M-2)*(M-2)).fill(z);
    std::cout << get_AB_matrix(M, cvec, r) << std::endl;
    int matrix_size = (M-2)*(M-2);
    arma::sp_cx_mat AB = arma::sp_cx_mat(matrix_size, matrix_size);
    std::cout << AB << std::endl;
    get_AB_matrix(AB, M, cvec, r);
    std::cout << AB << std::endl;
    // std::cout << size(AB)(0) << std::endl;


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

}

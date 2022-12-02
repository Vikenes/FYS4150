#ifndef __utils_hpp__
#define __utils_hpp__

#include <sstream> 
#include <string>
#include <iomanip>
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <armadillo>
#include <typeinfo>
#include <complex>

// std::string scientific_format(const double d, const int width=15, const int prec=10);

// /**
//      * Convert simple float to string with a given precision (default precision=2)
//      * Used to simplify file naming for different parameter choices.
// */
// std::string float_to_string(const double d, const int prec=2);

// int write_to_file_scientific(std::vector<double> col1, std::vector<double> col2, std::string filename, int width=15, int prec=10);
// int write_to_file_scientific(std::vector<double> col1, std::vector<double> col2, std::string filename, std::string header, int width=15, int prec=10);

int idx_k(int i, int j, int M);

// void submatrix_diag(arma::sp_cx_mat &AB, int idx, int M, arma::cx_vec vec, std::complex<double> r);

arma::sp_cx_mat get_AB_matrix(int M, arma::cx_vec cvec, std::complex<double> r);

// arma::mat get_AB_matrix(int M, arma::vec vec, double r);

void get_AB_matrix(arma::sp_cx_mat &AB, int M, arma::cx_vec cvec, std::complex<double> r);

void fill_AB_matrix(int M, double h, double Dt, const arma::mat &V, arma::sp_cx_mat &A, arma::sp_cx_mat &B);

void initialise_state(arma::cx_mat &u0, int M, double h, double xc, double yc, double sigma_x, double sigma_y, double p_x, double p_y);

std::complex<double> unnormalised_gaussian(double x, double y, double xc, double yc, double sigma_x, double sigma_y, double p_x, double p_y);

void set_up_walls(arma::sp_mat &V, double v0, int M, double h, int Ns=2, double T=0.02, double xc=0.5, double yc=0.5, double Sw=0.05, double Sa=0.05);


#endif
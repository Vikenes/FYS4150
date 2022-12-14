#ifndef __utils_hpp__
#define __utils_hpp__

#include <sstream> 
#include <string.h>
#include <iomanip>
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <armadillo>
#include <typeinfo>
#include <complex>
// #include<stdio.h>

// std::string scientific_format(const double d, const int width=15, const int prec=10);

// /**
//      * Convert simple float to string with a given precision (default precision=2)
//      * Used to simplify file naming for different parameter choices.
// */
// std::string float_to_string(const double d, const int prec=2);

// int write_to_file_scientific(std::vector<double> col1, std::vector<double> col2, std::string filename, int width=15, int prec=10);
// int write_to_file_scientific(std::vector<double> col1, std::vector<double> col2, std::string filename, std::string header, int width=15, int prec=10);

int idx_k(int i, int j, int M);

std::tuple<int, int> idx_ij(int k, int M);

// void submatrix_diag(arma::sp_cx_mat &AB, int idx, int M, arma::cx_vec vec, std::complex<double> r);

arma::sp_cx_mat get_AB_matrix(int M, arma::cx_vec cvec, std::complex<double> r);

// arma::mat get_AB_matrix(int M, arma::vec vec, double r);

void get_AB_matrix(arma::sp_cx_mat &AB, int M, arma::cx_vec cvec, std::complex<double> r);

void fill_AB_matrix(int M, double h, double Dt, const arma::sp_mat &V, arma::sp_cx_mat &A, arma::sp_cx_mat &B);

arma::cx_vec make_column_vector(arma::cx_mat U, int M);

arma::cx_mat make_matrix(arma::cx_vec column_vector, int M);

void initialise_state(arma::cx_mat &u0, int M, double h, double xc, double yc, double sigma_x, double sigma_y, double p_x, double p_y);

std::complex<double> unnormalised_gaussian(double x, double y, double xc, double yc, double sigma_x, double sigma_y, double p_x, double p_y);

// void set_up_walls_utils(arma::sp_mat &V, double v0, int M, double h, int Ns=2, double T=0.02, double xc=0.5, double yc=0.5, double Sw=0.05, double Sa=0.05);

// arma::cx_cube simulation(double h, double Dt, double T, double xc, double sigma_x, double p_x, double yc, double sigma_y, double p_y, double v0);

#endif
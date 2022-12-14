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
#include <chrono>
#include <tuple>

// std::string scientific_format(const double d, const int width=15, const int prec=10);

// /**
//      * Convert simple float to string with a given precision (default precision=2)
//      * Used to simplify file naming for different parameter choices.
// */
// std::string float_to_string(const double d, const int prec=2);

// int write_to_file_scientific(std::vector<double> col1, std::vector<double> col2, std::string filename, int width=15, int prec=10);
// int write_to_file_scientific(std::vector<double> col1, std::vector<double> col2, std::string filename, std::string header, int width=15, int prec=10);



// Functions to help Simulation

/*
A, B, a, b, r, etc.: Following the notation in the report.
*/

/** 
 * @brief Get index k for flattened vector u
 * @param i horisontal index
 * @param j vertical index
 * @param M number of points in each direction
 * @returns index k
*/
int idx_k(int i, int j, int M);

/**
 * @brief Get indices i, j for matrix U
 * @param k flattened index
 * @param M number of points in each direction
 * @returns tuple of indices (i, j)
*/
std::tuple<int, int> idx_ij(int k, int M);

/**
 * @brief Set up A/B matrix of size ((M-2)**2, (M-2)**2)
 * @param M number of points in each direction
 * @param cvec complex vector a/b
 * @param r ratio r
 * @returns complex matrix A/B
*/
arma::sp_cx_mat get_AB_matrix(int M, arma::cx_vec cvec, std::complex<double> r);

/**
 * @brief Set up the A and B matrix of size ((M-2)**2, (M-2)**2)
 * @param AB reference to A/B matrix to fill
 * @param M number of points in each direction
 * @param cvec complex vector a/b
 * @param r ratio r
 */
void get_AB_matrix(arma::sp_cx_mat &AB, int M, arma::cx_vec cvec, std::complex<double> r);

/**
 * @brief Fill A- and B-matrices.
 * @param M number of points in each direction
 * @param h spatial step size
 * @param dt temporal step size
 * @param V reference to a sparse matrix containing the grid potential
 * @param A reference to a sparse, complex matrix A
 * @param B reference to a sparse, complex matrix B
*/
void fill_AB_matrix(int M, double h, double dt, const arma::sp_mat &V, arma::sp_cx_mat &A, arma::sp_cx_mat &B);

/**
 * @brief Make a flatted column vector u out of matrix U
 * @param U complex matrix of u_{i,j} values
 * @param M number of points in each direction
 * @returns complex column vector u of length (M-2)*(M-2)
*/
arma::cx_vec make_column_vector(arma::cx_mat U, int M);

/**
 * @brief Set up a matrix of size M*M from the flattened column vector
 * @note Imposes periodic boundary conditions.
 * @param column_vector flattened column vector u 
 * @param M number of points in each direction
*/
arma::cx_mat make_matrix(arma::cx_vec column_vector, int M);

/**
 * @brief Compute the unnormalised 2D Gaussian.
 * @returns Gaussian evaluated at (x,y), a complex number
*/
std::complex<double> unnormalised_gaussian(double x, double y, double xc, double yc, double sigma_x, double sigma_y, double p_x, double p_y);

/**
 * @brief Initialise the particle state on the lattice with a 2D Gaussian.
 * @param x horisontal position for which the Gaussian is evaluated
 * @param y vertical position for which the Gaussian is evaluated
*/
void initialise_state(arma::cx_mat &u0, int M, double h, double xc, double yc, double sigma_x, double sigma_y, double p_x, double p_y);




#endif
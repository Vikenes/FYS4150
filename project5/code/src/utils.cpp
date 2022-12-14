#include "utils.hpp"
#include "Box.hpp"
#include "Simulation.hpp"

// std::string scientific_format(const double d, const int width, const int prec){
//     std::stringstream ss;
//     ss << std::setw(width) << std::setprecision(prec) << std::scientific << d; 
//     return ss.str();
// }

// std::string float_to_string(const double d, const int prec){
//     std::stringstream ss;
//     ss << std::setprecision(prec) << d; 
//     return ss.str();
// }


// int write_to_file_scientific(std::vector<double> col1, std::vector<double> col2, std::string filename, int width, int prec){
                        
//     std::string path = "../output/data/"; // path for .txt files
//     std::string file = path + filename + ".txt";
//     std::ofstream ofile;

//     ofile.open(file.c_str());

//     int n = col1.size();
//     for(int i=0; i<n; i++){
//         ofile << scientific_format(col1[i], width, prec) << ", "
//               << scientific_format(col2[i], width, prec)
//               << std::endl;
//     }
//     std::cout << "Finished writing to file: " << file << std::endl;
//     ofile.close();

//     return 0;
// }



// int write_to_file_scientific(std::vector<double> col1, std::vector<double> col2, std::string filename, std::string header, int width, int prec){
                        
//     std::string path = "../output/data/"; // path for .txt files
//     std::string file = path + filename + ".txt";
//     std::ofstream ofile;

//     ofile.open(file.c_str());

//     ofile << header << std::endl;

//     int n = col1.size();
//     for(int i=0; i<n; i++){
//         ofile << scientific_format(col1[i], width, prec) << ", "
//               << scientific_format(col2[i], width, prec)
//               << std::endl;
//     }

//     ofile.close();

//     return 0;
// }


// int write_arma_to_file_scientific(arma::cube R, std::string filename, int width, int prec){
//     // Write armadillo cube to file                         
//     std::string path = "../output/data/"; // path for .txt files
//     std::ofstream ofile;

//     // std::cout << "test:" << std::endl;
//     for(int p=0; p<R.n_cols; p++){
        
//         std::string file = path + filename + "_p" + std::to_string(p+1) + ".txt";
//         ofile.open(file.c_str());
//         std::string header = "t, x, y, z, vx, vy, vz \n";
//         ofile << header;
//         // std::cout << "sl" << std::endl;
        
//         for(int i=0; i<R.n_slices; i++){
//             ofile << scientific_format(R(0,p,i));
//             for(int j=1; j<R.n_rows; j++){
//                 //x = R.slice(i).col(p).row(j);
//                 ofile << ", " << scientific_format(R(j,p,i)); // position, velocity
//             }
//             ofile << std::endl;
//         }
        
//         ofile.close();
//     }


//     return 0;
// }



int idx_k(int i, int j, int M){
    return (j-1)*(M-2) + (i-1);
}

std::tuple<int, int> idx_ij(int k, int M){
    return std::make_tuple(k%(M-2)+1, k/(M-2)+1);
}

arma::sp_cx_mat get_AB_matrix(int M, arma::cx_vec cvec, std::complex<double> r){

    //  declare the matrix:
    int matrix_size = (M-2)*(M-2);
    arma::sp_cx_mat AB = arma::sp_cx_mat(matrix_size, matrix_size);
    AB.diag(M-2).fill(r);
    AB.diag(-(M-2)).fill(r);
    AB.diag() = cvec;
    //  iterate and fill:
    for(int q=0; q<M-2;q++){
        for(int p=0; p<M-2-1;p++){
            int idx = q*(M-2)+p;
            AB.diag(-1)(idx) = r;
            AB.diag(1)(idx) = r;
        }
    }
    return AB;
}

void get_AB_matrix(arma::sp_cx_mat &AB, int M, arma::cx_vec cvec, std::complex<double> r){
    
    AB.diag(M-2).fill(r);
    AB.diag(-(M-2)).fill(r);
    AB.diag() = cvec;

    //  iterate and fill:
    for(int q=0; q<M-2;q++){
        for(int b=0; b<M-2-1;b++){
            int idx = q*(M-2)+b;
            AB.diag(-1)(idx) = r;
            AB.diag(1)(idx) = r;
        }
    }
}

void fill_AB_matrix(int M, double h, double dt, const arma::sp_mat &V, arma::sp_cx_mat &A, arma::sp_cx_mat &B){
    
    //  define r, a and b:
    std::complex<double> r = std::complex<double>(0, dt/(2*h*h));
    arma::cx_vec a((M-2)*(M-2));
    arma::cx_vec b((M-2)*(M-2));
    
    //  iterate through grid and fill a and b:
    for(int i=1; i<=M-2; i++){
        for(int j=1; j<=M-2; j++){
            int k = idx_k(i,j,M);
            double abval = dt/2*V(i,j);
            a(k) = std::complex<double>(1,0) + std::complex<double>(4,0)*r + std::complex<double>(0,abval);
            b(k) = std::complex<double>(1,0) - std::complex<double>(4,0)*r - std::complex<double>(0,abval);
        }
    }

    //  set up the matrices:
    get_AB_matrix(A, M, a, -r);
    get_AB_matrix(B, M, b, r);
}

arma::cx_vec make_column_vector(arma::cx_mat matrix, int M){

    arma::cx_vec column_vector = arma::cx_vec((M-2)*(M-2));
    
    //  iterate through the lattice and fill the column vector with matrix elements;
    for(int i=1;i<=M-2;i++){
        for(int j=1;j<=M-2;j++){
            int k = idx_k(i,j, M);
            column_vector(k) = matrix(i,j);
        }
    }

    return column_vector;
}

arma::cx_mat make_matrix(arma::cx_vec column_vector, int M){

    arma::cx_mat matrix = arma::cx_mat(M,M);
    //  iterate through the flattened grid and fill matrix with vector elements:
    for(int k=0;k<(M-2)*(M-2);k++){
        int i,j;
        std::tie(i,j) = idx_ij(k, M);
        matrix(i,j) = column_vector(k);
    }

    //  impose boundary conditions:
    matrix.col(0).fill(0);
    matrix.col(M-1).fill(0);
    matrix.row(0).fill(0);
    matrix.row(M-1).fill(0);

    return matrix;
}

std::complex<double> unnormalised_gaussian(double x, double y, double xc, double yc, double sigma_x, double sigma_y, double p_x, double p_y){
    //  compute real and imaginary part of exponent:
    double real_exp_part = -((x-xc)*(x-xc)/(2*sigma_x*sigma_x))-((y-yc)*(y-yc)/(2*sigma_y*sigma_y));
    double imag_exp_part = p_x*(x-xc) + p_y*(y-yc);
    
    //  compute the exponent:
    std::complex<double> exponent = std::complex<double>(real_exp_part, imag_exp_part);
    
    return exp(exponent);
}

void initialise_state(arma::cx_mat &u0, int M, double h, double xc, double yc, double sigma_x, double sigma_y, double p_x, double p_y){

    double pnorm = 0.0;
    //  iterate through grid and fill the initial state with the Gaussian exaluated at that point:
    for(int i=1; i<=M-2; i++){
        for(int j=1; j<=M-2;j++){
            double x = i*h;
            double y = j*h;
            u0(i,j) = unnormalised_gaussian(x, y, xc, yc, sigma_x, sigma_y, p_x, p_y);
            pnorm += std::norm(u0(i,j));
        }   
    }

    //  boundary conditions:
    u0.col(0).fill(0);
    u0.col(M-1).fill(0);
    u0.row(0).fill(0);
    u0.row(M-1).fill(0);

    //  normalise:
    u0 /= sqrt(pnorm);

}

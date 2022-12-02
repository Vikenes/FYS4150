#include "utils.hpp"

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
    // Valid for i,j in [1,M-2]
    // (j-1)->j for i,j in [0, M-2)
    int k = (j-1)*(M-2) +i-1;
    return k;
}

// void submatrix_diag(arma::sp_cx_mat &AB, int idx, int M, arma::cx_vec vec, std::complex<double> r){
//     // Better way of doing this by slicing incoming vector.
//     // Could also take a large matrix as reference and fill it with a routine. 
//     int idx0 = idx * (M-2);
//     int idx1 = idx0 + (M-3);
//     arma::cx_mat ret_mat = arma::cx_mat(M-2, M-2);
//     ret_mat.diag() = vec.subvec(idx0, idx1);
//     ret_mat.diag(1).fill(r);
//     ret_mat.diag(-1).fill(r);
//     int first = idx * (M-2);
//     AB.submat(first, first+(M-2), first, first+(M-2)) = ret_mat;
// }

arma::sp_cx_mat get_AB_matrix(int M, arma::cx_vec cvec, std::complex<double> r){
    /**
     * @brief Sets up and return the A and B matrix of size ((M-2)**2, (M-2)**2)
     * @param M
     * @param cvec 
     * @param r
     * 
     * @return AB
     * 
     */
    int matrix_size = (M-2)*(M-2);
    arma::sp_cx_mat AB = arma::sp_cx_mat(matrix_size, matrix_size);
    AB.diag(M-2).fill(r);
    AB.diag(-(M-2)).fill(r);
    AB.diag() = cvec;
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
    /**
     * @brief Sets up and return the A and B matrix of size ((M-2)**2, (M-2)**2)
     * @param AB reference to AB matrix to fill
     * @param M
     * @param cvec 
     * @param r
     * 
     */
    // int M = std::pow(size(AB)(0),1/2) + 2;
    // std::cout<<M<<std::endl;
    AB.diag(M-2).fill(r);
    AB.diag(-(M-2)).fill(r);
    AB.diag() = cvec;
    for(int q=0; q<M-2;q++){
        for(int b=0; b<M-2-1;b++){
            int idx = q*(M-2)+b;
            AB.diag(-1)(idx) = r;
            AB.diag(1)(idx) = r;
        }
    }
}

// arma::mat get_AB_matrix(int M, arma::vec vec, double r){
//     int matrix_size = (M-2)*(M-2);
//     arma::mat AB = arma::mat(matrix_size, matrix_size);
//     // Check the below diagonal indices
//     AB.diag(M-2).fill(r);
//     AB.diag(-(M-2)).fill(r);
//     AB.diag() = vec;
//     for(int q=0; q<M-2;q++){
//         for(int b=0; b<M-2-1;b++){
//             int idx = q*(M-2)+b;
//             AB.diag(-1)(idx) = r;
//             AB.diag(1)(idx) = r;
//         }
//     }
//     return AB;
// }

void fill_AB_matrix(int M, double h, double Dt, const arma::mat &V, arma::sp_cx_mat &A, arma::sp_cx_mat &B){
    // Define r, a and b
    std::complex<double> r = std::complex<double>(0,Dt/(2*h*h));
    arma::cx_vec a((M-2)*(M-2));
    arma::cx_vec b((M-2)*(M-2));
    for(int i=1; i<=M-2;i++){
        for(int j=1; j<=M-2; j++){
            int k = idx_k(i,j,M);
            double abval = Dt/2*V(i,j);
            a(k) = std::complex<double>(1,0) + std::complex<double>(4,0)*r + std::complex<double>(0,abval);
            b(k) = std::complex<double>(1,0) - std::complex<double>(4,0)*r - std::complex<double>(0,abval);
        }
    }
    get_AB_matrix(A, M, a, -r);
    get_AB_matrix(B, M, b, r);
}

std::complex<double> unnormalised_gaussian(double x, double y, double xc, double yc, double sigma_x, double sigma_y, double p_x, double p_y){
    double real_exp_part = -((x-xc)*(x-xc)/(2*sigma_x*sigma_x))-((y-yc)*(y-yc)/(2*sigma_y*sigma_y));
    double imag_exp_part = p_x*(x-xc) + p_y*(y-yc);
    std::complex<double> exponent = std::complex<double>(real_exp_part, imag_exp_part);
    return exp(exponent);
}

void initialise_state(arma::cx_mat &u0, int M, double h, double xc, double yc, double sigma_x, double sigma_y, double p_x, double p_y){
    double sqrt_norm = 0.0;
    for(int i=1; i<=M-2; i++){
        for(int j=1; i<=M-2;j++){
            double x = i*h;
            double y = j*h;
            u0(i,j) = unnormalised_gaussian(x, y, xc, yc, sigma_x, sigma_y, p_x, p_y);
            sqrt_norm += sqrt(norm(u0(i,j)));
        }   
    }
    // Boundary conditions 
    // Not sure if this will work
    u0.col(0).fill(0);
    u0.col(M-1).fill(0);
    u0.row(0).fill(0);
    u0.row(M-1).fill(0);

    // Normalise
    u0 = u0/sqrt_norm;
}


void set_up_walls(arma::sp_mat &V, double v0, int M, double h, int Ns, double T, double xc, double yc, double Sw, double Sa){
    /**
     * @brief Sets up the wall used to generate the slit apertures.
     * @param V reference to potential matrix V.
     * @param v0 constant potential value of wall.
     * @param M number of points in x and y direction.
     * @param h step length between points.
     * @param Ns number of slits.
     * @param T thickness of slits.
     * @param xc slit wall x-centre.
     * @param yc slit wall y-centre.
     * @param Sw slit wall segment length in y-direction
     * @param Sa slit aperture width in y-direction
     * 
     */
    // std::cout<<M<<std::endl;
    for(int i=1; i<=M-2; i++){
        for(int j=1; j<=M-2; j++){
            double x = i*h;
            double y = j*h;
            // Ignore number of slits for now, just assume that it is two
            
            // First set x-criterion
            if(x<=xc+T/2 && x>xc-T/2){
                // Then set the various y-criteria
                // Will be a loop here if you include varying number of slits (must differ between odd and even)
                // First slit wall in centre
                if(y<=yc+Sw/2 && y>yc-Sw/2){
                    V(i,j) = v0;
                }
                // Wall on side of larger y-vals
                else if(y<=yc+3/2*Sw+Sa && y>yc+1/2*Sw+Sa){
                    V(i,j) = v0;
                }
                // Wall on side of smaller y-vals
                else if(y<=yc-1/2*Sw-Sa && y>yc-3/2*Sw-Sa){
                    V(i,j) = v0;
                }
            }
        }
    }
    // Not sure if this needs to be imposed here
    V.col(0).fill(0);
    V.col(M-1).fill(0);
    V.row(0).fill(0);
    V.row(M-1).fill(0);
}


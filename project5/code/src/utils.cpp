#include "utils.hpp"

std::string scientific_format(const double d, const int width, const int prec){
    std::stringstream ss;
    ss << std::setw(width) << std::setprecision(prec) << std::scientific << d; 
    return ss.str();
}

std::string float_to_string(const double d, const int prec){
    std::stringstream ss;
    ss << std::setprecision(prec) << d; 
    return ss.str();
}


int write_to_file_scientific(std::vector<double> col1, std::vector<double> col2, std::string filename, int width, int prec){
                        
    std::string path = "../output/data/"; // path for .txt files
    std::string file = path + filename + ".txt";
    std::ofstream ofile;

    ofile.open(file.c_str());

    int n = col1.size();
    for(int i=0; i<n; i++){
        ofile << scientific_format(col1[i], width, prec) << ", "
              << scientific_format(col2[i], width, prec)
              << std::endl;
    }
    std::cout << "Finished writing to file: " << file << std::endl;
    ofile.close();

    return 0;
}



int write_to_file_scientific(std::vector<double> col1, std::vector<double> col2, std::string filename, std::string header, int width, int prec){
                        
    std::string path = "../output/data/"; // path for .txt files
    std::string file = path + filename + ".txt";
    std::ofstream ofile;

    ofile.open(file.c_str());

    ofile << header << std::endl;

    int n = col1.size();
    for(int i=0; i<n; i++){
        ofile << scientific_format(col1[i], width, prec) << ", "
              << scientific_format(col2[i], width, prec)
              << std::endl;
    }

    ofile.close();

    return 0;
}


int write_arma_to_file_scientific(arma::cube R, std::string filename, int width, int prec){
    // Write armadillo cube to file                         
    std::string path = "../output/data/"; // path for .txt files
    std::ofstream ofile;

    // std::cout << "test:" << std::endl;
    for(int p=0; p<R.n_cols; p++){
        
        std::string file = path + filename + "_p" + std::to_string(p+1) + ".txt";
        ofile.open(file.c_str());
        std::string header = "t, x, y, z, vx, vy, vz \n";
        ofile << header;
        // std::cout << "sl" << std::endl;
        
        for(int i=0; i<R.n_slices; i++){
            ofile << scientific_format(R(0,p,i));
            for(int j=1; j<R.n_rows; j++){
                //x = R.slice(i).col(p).row(j);
                ofile << ", " << scientific_format(R(j,p,i)); // position, velocity
            }
            ofile << std::endl;
        }
        
        ofile.close();
    }


    return 0;
}

int idx_k(int i, int j, int M){
    int k = j*(M-2) +1;
    return k;
}

arma::cx_mat submatrix_diag(int idx, int M, arma::cx_vec vec, std::complex<double> r){
    // Better way of doing this by slicing incoming vector.
    // Could also take a large matrix as reference and fill it with a routine. 
    int idx0 = idx * (M-2);
    int idx1 = idx0 + (M-3);
    arma::cx_mat ret_mat = arma::cx_mat(M-2, M-2);
    ret_mat.diag() = vec.subvec(idx0, idx1);
    ret_mat.diag(1).fill(r);
    ret_mat.diag(-1).fill(r);
    
}
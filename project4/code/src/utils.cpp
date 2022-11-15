#include "utils.hpp"

std::string scientific_format(const double d, const int width, const int prec){
    std::stringstream ss;
    ss << std::setw(width) << std::setprecision(prec) << std::scientific << d; 
    return ss.str();
}

// std::string scientific_format(double d, const int width, const int prec){
    // std::stringstream ss;
    // ss << std::setw(width) << std::setprecision(prec) << std::scientific << d; 
    // return ss.str();
// }

int write_to_file_scientific(std::vector<double> x, std::vector<double> v, std::string filename, int width, int prec){
                        
    std::string path = "../output/data/"; // path for .txt files
    std::string file = path + filename + ".txt";
    std::ofstream ofile;

    ofile.open(file.c_str());

    int n = x.size();
    for(int i=0; i<n; i++){
        ofile << scientific_format(x[i], width, prec) << ", "
              << scientific_format(v[i], width, prec)
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




int write_to_file(std::vector<double> a, std::vector<double> b, std::string filename, int width){
                        
    std::string path = "../output/data/"; // path for .txt files
    std::string file = path + filename + ".txt";
    std::ofstream ofile;

    ofile.open(file.c_str());

    int n = a.size();
    for(int i=0; i<n; i++){
        ofile << std::setw(width) << a[i] << ", "
              << std::setw(width) << b[i]
              << std::endl;
    }

    ofile.close();

    return 0;
}

int write_to_file(std::vector<double> a, std::vector<int> b, std::string filename, int width){
                        
    std::string path = "../output/data/"; // path for .txt files
    std::string file = path + filename + ".txt";
    std::ofstream ofile;

    ofile.open(file.c_str());

    int n = a.size();
    for(int i=0; i<n; i++){
        ofile << std::setw(width) << a[i] << ", "
              << std::setw(width) << b[i]
              << std::endl;
    }

    ofile.close();

    return 0;
}



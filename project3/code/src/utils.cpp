#include "utils.hpp"

std::string scientific_format(const double d, const int width, const int prec){
    std::stringstream ss;
    ss << std::setw(width) << std::setprecision(prec) << std::scientific << d; 
    return ss.str();
}


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

//  Is it allowed to write it like this:?

// double k_e = 1.38935333 * 10 ** 5;
// double T = 9.64852558 * 10 ** 1;
// double V = 9.64852558* 10 ** 7;
// double B0 = 1 * T;
// double V0 = 10 * V;
// double d = 10 ** 4;
// double Vdr = 9.65;
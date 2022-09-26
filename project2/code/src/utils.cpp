#include "utils.hpp"

std::string scientific_format(const double d, const int width, const int prec){
    std::stringstream ss;
    ss << std::setw(width) << std::setprecision(prec) << std::scientific << d; 
    return ss.str();
}

int write_to_file(std::vector<double> x, std::vector<double> v, std::string filename, int width, int prec){
                        
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
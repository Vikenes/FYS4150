#include <iostream>
#include <string>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>
#include <chrono>

double f(double x){
    return 100 * exp(-10*x);
}

std::vector<double> Thomas(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> g){
    int m = b.size();
    std::vector<double> btilde(m);
    std::vector<double> gtilde(m);

    btilde[0] = b[0];
    gtilde[0] = g[0];
     
    // Step 1 - Forward

    for(int i=1; i<m; i++){
        double K = a[i]/btilde[i-1];
        btilde[i] = b[i] - K*c[i-1];
        gtilde[i] = g[i] - K*gtilde[i-1];
    }

    // Step 2 - Backward

    std::vector<double> vstar(m);
    vstar[-1] = gtilde[-1]/btilde[-1];

    for(int i=m-2; i>=0; i--){
        vstar[i] = gtilde[i] - c[i]*vstar[i+1]/btilde[i];
    }

    return vstar;

}

std::vector<double> Thomas_special(std::vector<double> g){
    int m = g.size();

    std::vector<double> gtilde(m);

    gtilde[0] = g[0];
    
    // Step 1 - Forward

    for(int i=1; i<m; i++){
        double ii = i-1;
        gtilde[i] = g[i] + (ii+1)/(ii+2)*gtilde[i-1];
    }

    // Step 2 - Backward

    std::vector<double> vstar(m);
    double mm = m-1;
    vstar[-1] = gtilde[-1]* (mm+1)/(mm+2);

    for(int i=m-2; i>=0; i--){
        vstar[i] = gtilde[i] + vstar[i+1]*(i+1.0)/(i+2.0);
    }

    return vstar;

}

int main(){

    auto start_time = std::chrono::high_resolution_clock::now();

    double x_min = 0;
    double x_max = 1;
    int n_steps = 100;
    int n_points = n_steps + 1;
    double h = (x_max - x_min) / n_steps;

    int m = n_points - 2;
    std::vector<double> x(n_points);
    x[0] = x_min;
    x[n_steps] = x_max;
    std::vector<double> g(m);
   
    for(int i=1; i<=m; i++){
        int j = i-1;
        x[i] = x_min + i*h;
        g[j] = std::pow(h,2) * f(x[i]); 
    }

    std::vector<double> a(m, -1);
    std::vector<double> b(m, 2);
    std::vector<double> c(m, -1);
    
    // std::vector<double> vstar = Thomas(a, b, c, g);
    std::vector<double> vstar = Thomas_special(g);
    std::vector<double> v(n_points);

    // Include endpoints

    v[0] = 0;
    v[n_steps] = 0;
    for(int i=1; i<=m; i++){
        int j = i-1;
        v[i] = vstar[j];
    }

    // Write to file

    std::string path="folder/";

    std::string filename="num_sol_" + std::to_string(n_steps) + "steps.txt";
    std::ofstream ofile;
    std::string file = path + filename;
    ofile.open(file.c_str());

    int width=15;
    int prec=5;

    for(int i=0; i<=n_steps; i++){
        ofile << std::setw(width) << std::setprecision(prec) << std::scientific << x[i] 
              << std::setw(width) << std::setprecision(prec) << std::scientific << v[i]
              << std::endl;
    }

    ofile.close();

    auto stop_time = std::chrono::high_resolution_clock::now();

    double run_time = std::chrono::duration<double>(stop_time-start_time).count();

    std::cout << "Running time: " << run_time << " s" << std::endl;



}
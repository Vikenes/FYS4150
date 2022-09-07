#include <iostream>
#include <string>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>
#include <chrono>



double u(double x){
    // Analytical solution
    return 1 - (1 - exp(-10))*x - exp(-10*x);
}

double f(double x){
    // Right hand side of diff. eq.
    return 100 * exp(-10*x);
}

std::vector<double> generalThomas(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> g){
    // Thomas algorithm (general)
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
    vstar[m-1] = gtilde[m-1]/btilde[m-1];

    for(int i=m-2; i>=0; i--){
        vstar[i] = (gtilde[i] - c[i]*vstar[i+1]) / btilde[i];
    }

    return vstar;
}

std::vector<double> specialThomas(std::vector<double> g){
    // Thomas algorithm (general)
    int m = g.size();
    std::vector<double> btilde(m);
    std::vector<double> gtilde(m);

    btilde[0] = 2;
    gtilde[0] = g[0];
     
    // Step 1 - Forward

    for(int i=1; i<m; i++){
        btilde[i] = (i+2.0)/(i+1);
        gtilde[i] = g[i] + gtilde[i-1]/btilde[i-1];
    }

    // Step 2 - Backward

    std::vector<double> vstar(m);
    vstar[m-1] = gtilde[m-1]/btilde[m-1];

    for(int i=m-2; i>=0; i--){
        vstar[i] = (gtilde[i] + vstar[i+1]) / btilde[i];
    }

    return vstar;

}




int writeto_file(std::vector<double> x, std::vector<double> v, std::string filename){

    int width=15;
    int prec=5;
    std::string path = "../output/data/"; // path for .txt files
    std::string file = path + filename + ".txt";
    std::ofstream ofile;
    ofile.open(file.c_str());

    int n = x.size();
    for(int i=0; i<n; i++){
        ofile << std::setw(width) << std::setprecision(prec) << std::scientific << x[i] 
              << std::setw(width) << std::setprecision(prec) << std::scientific << v[i]
              << std::endl;
    }

    ofile.close();

    return 0;
}


// Defining global variables (don't know if this is OK...)
const double x_min = 0; 
const double x_max = 1;
const double v0 = 0;
const double v1 = 0;


// declare functions
int problem2(int n_steps);
int problem7(int n_steps);
int problem9(int n_steps);
int problem8(int n_steps, bool absolute_error, bool maximum_error);

int main(){

    auto start_time = std::chrono::high_resolution_clock::now();
 
    problem2(100);

    problem7(10);
    problem7(100);
    problem7(1000);

    problem9(10);
    problem9(100);
    problem9(1000);

    problem8(10, true, false);
    problem8(100, true, false);
    problem8(1000, true, false);

    problem8(10, false, false);
    problem8(100, false, false);
    problem8(1000, false, false);

    auto stop_time = std::chrono::high_resolution_clock::now();

    double run_time = std::chrono::duration<double>(stop_time-start_time).count();

    std::cout << "Running time: " << run_time << " s" << std::endl;

    return 0;

}

int problem2(int n_steps){
    
    // Problem 2

    int n_points = n_steps+1;
    double h = (x_max - x_min) / n_steps;

    std::vector<double> x(n_points);
    std::vector<double> v(n_points);

    for(int i=0; i<=n_steps; i++){
        x[i] = x_min + i*h;
        v[i] = u(x[i]);
    }
    writeto_file(x, v, "x_u");
    return 0;
}

int problem7(int n_steps){


    // Problem 7

    int n_points = n_steps+1;
    double h = (x_max - x_min) / n_steps;
    int m = n_points - 2;
    std::vector<double> x(n_points);
    x[0] = x_min;
    x[n_points-1] = x_max;
    std::vector<double> g(m);
   
    for(int i=1; i<=m; i++){
        x[i] = x_min + i*h;
        g[i-1] = std::pow(h,2) * f(x[i]); 
    }
    g[0] += v0;
    g[n_points-1] += v1;

    std::vector<double> a(m, -1);
    std::vector<double> b(m, 2);
    std::vector<double> c(m, -1);
    
    std::vector<double> vstar = generalThomas(a, b, c, g);
    
    // Include endpoints

    std::vector<double> v(n_points);
    v[0] = v0;
    v[n_points-1] = v1;
    for(int i=1; i<=m; i++){
        v[i] = vstar[i-1];
    }

    writeto_file(x, v, "num_sol_" + std::to_string(n_steps) + "steps");

    return 0;
}

int problem9(int n_steps){
    
    // Problem 9

    int n_points = n_steps+1;
    double h = (x_max - x_min) / n_steps;
    int m = n_points - 2;
    std::vector<double> x(n_points);
    x[0] = x_min;
    x[n_points-1] = x_max;
    std::vector<double> g(m);
   
    for(int i=1; i<=m; i++){
        x[i] = x_min + i*h;
        g[i-1] = std::pow(h,2) * f(x[i]); 
    }
    g[0] += v0;
    g[n_points-1] += v1;
    
    std::vector<double> vstar = specialThomas(g);
    
    // Include endpoints

    std::vector<double> v(n_points);
    v[0] = v0;
    v[n_points-1] = v1;
    for(int i=1; i<=m; i++){
        v[i] = vstar[i-1];
    }

    writeto_file(x, v, "special_num_sol_" + std::to_string(n_steps) + "steps");

    return 0;
}

int problem8(int n_steps, bool absolute_error, bool maximum_error){
    int n_points = n_steps+1;
    double h = (x_max - x_min) / n_steps;
    int m = n_points - 2;
    std::vector<double> x(n_points);
    x[0] = x_min;
    x[n_points-1] = x_max;
    std::vector<double> g(m);
   
    for(int i=1; i<=m; i++){
        x[i] = x_min + i*h;
        g[i-1] = std::pow(h,2) * f(x[i]); 
    }
    g[0] += v0;
    g[n_points-1] += v1;

    std::vector<double> a(m, -1);
    std::vector<double> b(m, 2);
    std::vector<double> c(m, -1);
    
    std::vector<double> vstar = generalThomas(a, b, c, g);
    
    // Include endpoints

    std::vector<double> v(n_points);
    v[0] = v0;
    v[n_points-1] = v1;
    for(int i=1; i<=m; i++){
        v[i] = vstar[i-1];
    }
    
    std::vector<double> u_exact(n_points);
    for(int i=0; i<=n_steps; i++){
        u_exact[i] = u(x[i]);
    }

    if (absolute_error==true && maximum_error==false){
        std::vector<double> absolute_error(n_points);
        absolute_error[0] = 0;
        absolute_error[n_points-1] = 0;
        for(int i=1; i<n_steps; i++){
            double abs_err = std::abs(u(x[i])-v[i]);
            absolute_error[i] = std::log10(abs_err);
        }
        writeto_file(x, absolute_error, "absolute_error" + std::to_string(n_steps) + "steps");
    }
    else if (absolute_error==false && maximum_error==false){
        std::vector<double> relative_error(n_points);
        relative_error[0] = 0;
        relative_error[n_points-1] = 0;
        for(int i=1; i<n_steps; i++){
            double rel_err = std::abs((u(x[i])-v[i])/u(x[i]));
            relative_error[i] = std::log10(rel_err);
        }
        writeto_file(x, relative_error, "relative_error" + std::to_string(n_steps) + "steps");
    }
    // MAXIMUM RELATIVE ERROR FOR n_steps up to 10^7.
    return 0;
}
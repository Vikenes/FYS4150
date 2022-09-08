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
double problem8(int n_steps, bool absolute_error, bool maximum_error);
int problem8c(int number_of_n_steps);

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

    problem8c(7);

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
    g[m-1] += v1;

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
    g[m-1] += v1;
    
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

double problem8(int n_steps, bool absolute_error, bool maximum_error){
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
    g[m-1] += v1;

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
    
    //std::vector<double> u_exact(n_points);
    //for(int i=0; i<=n_steps; i++){
    //    u_exact[i] = u(x[i]);
    //}

    

    double return_val;
    if (absolute_error==true){
        std::vector<double> absolute_error(m);
        for(int i=0; i<m; i++){
            double abs_err = std::abs(u(x[i+1])-v[i+1]);
            absolute_error[i] = std::log10(abs_err);
        }
        // WRONG! x includes end-points
        writeto_file(x, absolute_error, "absolute_error" + std::to_string(n_steps) + "steps");
        return_val = 0;
    }
    else if (absolute_error==false){
        std::vector<double> relative_error(m);
        std::vector<double> log_relative_error(m);
        for(int i=0; i<m; i++){
            double rel_err = std::abs((u(x[i+1])-v[i+1])/u(x[i+1]));
            relative_error[i] = rel_err;
            log_relative_error[i] = std::log10(rel_err);
        }
        if(maximum_error==true){
            double max_val = 0;
            int index;
            for(int i=0; i<m; i++){
                double current = relative_error[i];
                //double prior = relative_error[i-1];
                if(current>max_val){
                    index = i;
                    max_val = current;
                }
            }
            //std::cout << index << ' ' << max_val << std::endl;
            return_val = max_val;
        }
        else{
            writeto_file(x, log_relative_error, "relative_error" + std::to_string(n_steps) + "steps");
            return_val = 0;
        }
    }
    //std::cout << return_val << std::endl;

    return return_val;
}

int problem8c(int number_of_n_steps){
    //  problem 8c
    // int n_step_max = 1e7;
    // int n_step_min = 1e1;
    // int step_difference = (n_step_max-n_step_min)/(number_of_n_steps-1);
    std::vector<double> list_of_n_steps(number_of_n_steps); //should be int, but that messes with the writeto_file() function
    std::vector<double> max_error_per_step_number(number_of_n_steps);
    for(int i=1; i<=number_of_n_steps; i++){
        list_of_n_steps[i-1] = std::pow(10,i);
        int n_steps = std::pow(10,i);
        max_error_per_step_number[i-1] = problem8(n_steps, false, true);
        //std::cout << n_steps << "  " << max_error_per_step_number[i-1] << std::endl;
    }
    writeto_file(list_of_n_steps, max_error_per_step_number, "max_relative_error");
    return 0;
}
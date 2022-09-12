#include <algorithm>
#include <chrono>

#include "utils.hpp"
#include "algorithms.hpp"


const double x_min = 0;
const double x_max = 1; 
const double v_min = 0;
const double v_max = 0;


void write_u(int n_steps, std::string filename);

double step_size(int n_steps){
    double h=(x_max - x_min)/n_steps;
    return h;
}

std::vector<double> x_array(int n_steps){
    // Create x_array with n_points=n_steps+1 between x_min and x_max
    int n_points = n_steps + 1;
    double h = step_size(n_steps);

    std::vector<double> x(n_points);
    x[0] = x_min;
    x[n_steps] = x_max;

    for(int i=1; i<n_steps; i++){x[i] = x_min + i*h;}

    return x;
}

std::vector<double> g_array(int n_steps){
    int m = n_steps-1;
    std::vector<double> g(m);
    double h=step_size(n_steps);
    std::vector<double> x=x_array(n_steps);
    for (int i=1; i<=m; i++){
        g[i-1] = std::pow(h,2) * f(x[i]);
    }
    g[0] += v_min;
    g[m-1] += v_max;

    return g;
}

std::vector<double> compute_generalThomas(int n_steps, bool write=false, int timing=0){
    int m = n_steps - 1;
    std::vector<double> a(m, -1);
    std::vector<double> b(m, 2);
    std::vector<double> c(m, -1);

    std::vector<double> x=x_array(n_steps);
    std::vector<double> g=g_array(n_steps);

    if (timing!=0){
        std::vector<double> time_mean_stddev(2);
        std::vector<double> times(timing);
        double mean = 0;
        double variance = 0;
        double stddev = 0;

        for(int j=0; j<timing; j++){
            auto t1 = std::chrono::high_resolution_clock::now();
            std::vector<double> v_time = generalThomas(a, b, c, g);
            auto t2 = std::chrono::high_resolution_clock::now();
            times[j] = std::chrono::duration<double>(t2 - t1).count();
            mean += times[j];
        }
        mean /= timing; // finding the actual mean, not the sum.
        
        for(int j=0; j<timing; j++){
            variance += std::pow((times[j] - mean), 2);
        }

        variance /= timing;

        stddev = std::pow(variance, .5);

        time_mean_stddev[0] = mean;
        time_mean_stddev[1] = stddev;

        return time_mean_stddev;
    }

    else{
        std::vector<double> v = generalThomas(a, b, c, g);
        if (write==true){
            write_to_file(x, v, "generalThomas_" + std::to_string(n_steps) + "steps");
        }
        return v;
    }

} 

std::vector<double> compute_specialThomas(int n_steps, bool write=false, int timing=0){
    
    std::vector<double> g=g_array(n_steps);

    if (timing!=0){
        std::vector<double> time_mean_stddev(2);
        std::vector<double> times(timing);
        double mean = 0;
        double variance = 0;
        double stddev = 0;

        for(int j=0; j<timing; j++){
            auto t1 = std::chrono::high_resolution_clock::now();
            std::vector<double> v_time = specialThomas(g);
            auto t2 = std::chrono::high_resolution_clock::now();
            times[j] = std::chrono::duration<double>(t2 - t1).count();
            mean += times[j];
        }
        mean /= timing; // finding the actual mean, not the sum.
        
        for(int j=0; j<timing; j++){
            variance += std::pow((times[j] - mean), 2);
        }

        variance /= timing;

        stddev = std::pow(variance, .5);

        time_mean_stddev[0] = mean;
        time_mean_stddev[1] = stddev;

        return time_mean_stddev;
    }

    else{   
    std::vector<double> v_star = specialThomas(g);
    if (write==true){
        std::vector<double> x=x_array(n_steps);
        write_to_file(x, v_star, "ST_" + std::to_string(n_steps) + "steps");
    }
    return v_star;
    
    }

}

void absolute_error(int n_steps){
    // Make either utils or algorithms script for this computation?
    int m = n_steps - 1;

    std::vector<double> x=x_array(n_steps);
    std::vector<double> v=compute_generalThomas(n_steps);


    std::vector<double> abs_error(m);
    std::vector<double> x_red = std::vector<double>(x.begin()+1, x.end()-1);
    
    for(int i=0; i<m; i++){
        double abs_err = std::abs(u(x_red[i])-v[i+1]);
        abs_error[i] = std::log10(abs_err);
        
    }
    write_to_file(x_red, abs_error, "abs_error_" + std::to_string(n_steps) + "steps");
}

std::vector<double> relative_error(int n_steps, bool write=false){
    int m = n_steps - 1;
    double return_val; 

    std::vector<double> x=x_array(n_steps);
    std::vector<double> v=compute_generalThomas(n_steps);

    std::vector<double> rel_error(m);
    std::vector<double> log_rel_error(m);
    std::vector<double> x_red = std::vector<double>(x.begin()+1, x.end()-1);

    for(int i=0; i<m; i++){
        // x_reduced[i] = x[i+1];
        double rel_err = std::abs((u(x_red[i])-v[i+1])/u(x_red[i]));
        rel_error[i] = rel_err;
        log_rel_error[i] = std::log10(rel_err);
    }

    if(write==true){
        write_to_file(x_red, log_rel_error, "rel_error_" + std::to_string(n_steps) + "steps", 20,15);
    }

    return rel_error;
}

int max_rel_err(int number_of_n_steps){
    std::vector<double> n_steps_array(number_of_n_steps); //should be int, but that messes with the writeto_file() function
    std::vector<double> max_error_n_step(number_of_n_steps);
    for(int i=1; i<=number_of_n_steps; i++){
        n_steps_array[i-1] = std::pow(10,i);
        int n_steps = std::pow(10,i);
        std::vector<double> rel_err = relative_error(n_steps);
        max_error_n_step[i-1] = *std::max_element(rel_err.begin(), rel_err.end());
    }
    write_to_file(n_steps_array, max_error_n_step, "max_rel_error");
    return 0;
}

void timing(){
    int n_step_sizes = 6;
    int n_simulations = 500;
    std::vector<double> list_of_n_steps(n_step_sizes); // double for writing to file 
    std::vector<double> general_Thomas_timed(n_step_sizes);
    std::vector<double> special_Thomas_timed(n_step_sizes);
    std::vector<double> general_Thomas_stddev(n_step_sizes);
    std::vector<double> special_Thomas_stddev(n_step_sizes);
    for(int i=0; i<n_step_sizes; i++){
        list_of_n_steps[i] = std::pow(10,i+1);
    }

    std::string filename = "thomas_timed";
    int width=15;
    int prec=5;
    std::string path = "../output/data/"; // path for .txt files
    std::string file = path + filename + ".txt";
    std::ofstream ofile;
    ofile.open(file.c_str());

    for(int i=0; i<n_step_sizes; i++){
        int n_steps = list_of_n_steps[i];

        std::vector<double> general_mean_stddev = compute_generalThomas(n_steps, false, n_simulations);
        std::vector<double> special_mean_stddev = compute_specialThomas(n_steps, false, n_simulations);

        general_Thomas_timed[i] = general_mean_stddev[0];
        special_Thomas_timed[i] = special_mean_stddev[0];
        general_Thomas_stddev[i] = general_mean_stddev[1];
        special_Thomas_stddev[i] = special_mean_stddev[1];

        ofile << std::setw(width) << std::setprecision(prec) << std::scientific << n_steps << ","
            << std::setw(width) << std::setprecision(prec) << std::scientific << general_Thomas_timed[i] << ","
            << std::setw(width) << std::setprecision(prec) << std::scientific << special_Thomas_timed[i] << ","
            << std::setw(width) << std::setprecision(prec) << std::scientific << general_Thomas_stddev[i] << ","
            << std::setw(width) << std::setprecision(prec) << std::scientific << special_Thomas_stddev[i]
            << std::endl;
    }
    
    // write_to_file(list_of_n_steps, general_Thomas_timed, "general_thomas_timed");
    // write_to_file(list_of_n_steps, special_Thomas_timed, "special_thomas_timed");
    // write_to_file(list_of_n_steps, general_Thomas_stddev, "general_thomas_stddev");
    // write_to_file(list_of_n_steps, special_Thomas_stddev, "special_thomas_stddev");

    ofile.close();

    return;
}

int main(){

    auto start_time = std::chrono::high_resolution_clock::now();
    
    // problem 2
    write_u(1000, "analytical_x_u");

    // problem 7, 8a) and 8b)
    int max_order = 3;
    for(int i=1; i<=max_order; i++){
        int n = std::pow(10, i);
        compute_generalThomas(n, true);
        absolute_error(n);
        relative_error(n, true);
    }

    //  problem 8c)
    max_rel_err(7);

    //  problem10
    timing();

    auto stop_time = std::chrono::high_resolution_clock::now();

    double run_time = std::chrono::duration<double>(stop_time-start_time).count();

    std::cout << "Running time: " << run_time << " s" << std::endl;

    return 0;

}

void write_u(int n_steps, std::string filename){
    // Computes analytical u(x) for plotting  
    int n_points = n_steps + 1; 

    std::vector<double> x=x_array(n_steps);
    std::vector<double> u_arr(n_points);

    for(int i=0; i<n_points; i++){u_arr[i] = u(x[i]);}

    write_to_file(x, u_arr, filename);
}
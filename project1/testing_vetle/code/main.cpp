#include <algorithm>

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



std::vector<double> compute_generalThomas(int n_steps, bool write=false){
    int m = n_steps - 1;
    std::vector<double> a(m, -1);
    std::vector<double> b(m, 2);
    std::vector<double> c(m, -1);

    std::vector<double> x=x_array(n_steps);
    std::vector<double> g=g_array(n_steps);

    std::vector<double> v = generalThomas(a, b, c, g);

    if (write==true){
        write_to_file(x, v, "GT_" + std::to_string(n_steps) + "steps");
    }
    return v;
} 


void compute_specialThomas(int n_steps, bool write=false){
    int m = n_steps - 1;
    
    std::vector<double> x=x_array(n_steps);
    std::vector<double> g=g_array(n_steps);

    std::vector<double> v_star = specialThomas(g);

    if (write==true){
        write_to_file(x, v_star, "ST_" + std::to_string(n_steps) + "steps");
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
    write_to_file(n_steps_array, max_error_n_step, "test_max_relative_error");
    return 0;
}

int main(){
    
    // write_u(50, "analytical_xu_50");

    // std::vector<double> a()
    // compute_generalThomas(50);
    // compute_specialThomas(50);
    // absolute_error(10);
    // absolute_error(100);
    // absolute_error(1000);
    
    // relative_error(100);
    // relative_error(1000);

    max_rel_err(7);

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
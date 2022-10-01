#include <algorithm>
#include "utils.hpp"
#include "algorithms.hpp"

/**
 * PROJECT 3 FYS4150
 *      By: Johan Mylius Kroken, Vetle Vikenes and Nanna Bryne
 */

// define boundary conditions:
const double xhat_0 = 0;
const double xhat_n = 1;
const double v_0 = 0;
const double v_n = 0;

// parameters for Jacobi algorithm:
const double eps = 1e-8;
const int M_max = int(1e6);


double step_size(int n){
    /**
     * Return step size for a given number of discretisation steps.
     */
    double h=(xhat_n - xhat_0)/n;
    return h;
}

int main(){
    
    return 0;
}
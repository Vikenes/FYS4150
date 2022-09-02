#include <iostream>
#include <string>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>

double f(double x){
    return 100 * exp(-10*x);
}

int main(){

    double x_min = 0;
    double x_max = 1;
    int n_steps = 100;
    int n_points = n_steps + 1;
    double h = (x_max - x_min) / n_steps;

    int m = n_points - 2;
    std::vector<double> x(m);
    std::vector<double> g(m);
   
    for(int i=1; i<=m; i++){
        int j = i-1;
        x[j] = x_min + i*h;
        g[j] = std::pow(h,2) * f(x[j]); 
    }

    /// Special case:

    double a = -1;
    double b = 2;
    double c = -1;

    std::vector<double> btilde(m);
    std::vector<double> gtilde(m);

    btilde[0] = b;
    gtilde[0] = g[0];

    for(int i=2; i<=m; i++){
        int j = i-1;
        double K = a/btilde[j-1];
        btilde[j] = b - K*c;
        gtilde[j] = g[j] - K*gtilde[j-1];
    }

    std::vector<double> vstar(m);
    vstar[m-1] = gtilde[m-1]/btilde[m-1];

    for(int i=m-1; i>0; i--){
        int j = i - 1;
        vstar[j] = gtilde[j] - c*vstar[j+1]/btilde[j];
    }

    std::vector<double> v(n_points);
    
    v[0] = 0;
    v[n_points-1] = 0;
    for(int i=1; i<=m; i++){
        int j = i-1;
        v[i] = vstar[j];
        std::cout << v[i] << std::endl;
    }

}
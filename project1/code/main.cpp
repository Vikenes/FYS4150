#include <iostream>
#include <string>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>



double u(double x){
    return 1 - (1 - exp(-10))*x - exp(-10*x);
}



int main(){

    // task 2
    std::string filename="x_u.txt";
    std::ofstream ofile;
    ofile.open(filename);

    int width=15;
    int prec=5;



    double x_min = 0; 
    double x_max = 1;
    int n_steps = 100;
    int n_points = n_steps+1;
    double h = (x_max - x_min) / n_steps;

    std::vector<double> x(n_points);
    std::vector<double> v(n_points);

    x[0] = x_min;
    v[0] = u(x[0]);

    for(int i=0; i<n_steps; i++){
        x[i+1] = x[i] + h;
        v[i+1] = u(x[i+1]);
    }

    for(int i=0; i<=n_steps; i++){
        ofile << std::setw(width) << std::setprecision(prec) << std::scientific << x[i] 
              << std::setw(width) << std::setprecision(prec) << std::scientific << v[i]
              << std::endl;
    }


    ofile.close();


    // task 7 THIS IS PROBABLY WRONG
    std::vector<double> a(n_points, 2);
    std::vector<double> b(n_points-1, -1);
    std::vector<double> c(n_points-1, -1);
    std::vector<double> g(n_points);

    for(int i=1; i<n_points-1; i++){
        g[i] = a[i-1]*v[i-1]+ b[i]*v[i] + c[i]*v[i+1];
    }
    
    //boundary points
    g[0]= 

    return 0;

}
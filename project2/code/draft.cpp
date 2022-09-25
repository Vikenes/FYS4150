#include <sstream> 
#include <string>
#include <iomanip>
#include <vector>
#include <fstream>
#include <iostream>
#include <armadillo>
#include <cmath>



const double pi = 3.14159265359;




// Create tridiagonal matrix from vectors.
// - lower diagonal: vector a, lenght N-1
// - main diagonal:  vector d, lenght N
// - upper diagonal: vector e, lenght N-1
arma::mat create_tridiagonal(const arma::vec& a, const arma::vec& d, const arma::vec& e){
    // Start from identity matrix
    int N = d.size();
    arma::mat A = arma::mat(N, N, arma::fill::eye);
    
    // Fill first row (row index 0)
    A(0,0) = d[0];
    A(0,1) = e[0];
    
    // Loop that fills rows 2 to N-1 (row indices 1 to N-2)
    for(int i=1; i<N-1; i++){
        A(i, i) = d[i];
        A(i, i-1) = a[i];
        A(i, i+1) = e[i];
    }

    // Fill last row (row index N-1)
    A(N-1,N-1) = d[N-1];
    A(N-1,N-2) = a[N-2];
    
    return A;
}


// Create a tridiagonal matrix tridiag(a,d,e) of size n*n
// from scalar input a, d and e
arma::mat create_tridiagonal(int N, double a, double d, double e)
{
  arma::vec a_vec = arma::vec(N-1, arma::fill::ones) * a;
  arma::vec d_vec = arma::vec(N, arma::fill::ones) * d;
  arma::vec e_vec = arma::vec(N-1, arma::fill::ones) * e;

  // Call the vector version of this function and return the result
  return create_tridiagonal(a_vec, d_vec, e_vec);
}


// Create a symmetric tridiagonal matrix tridiag(a,d,a) of size n*n
// from scalar input a and d.
arma::mat create_symmetric_tridiagonal(int N, double a, double d)
{
  // Call create_tridiagonal and return the result
  return create_tridiagonal(N, a, d, a);
}


int main(){

    // PROBLEM 2
    int N = 6;
    int n = N + 1;
    double xhat_max = 1;
    double xhat_min = 0;
    double h = (xhat_max - xhat_min)/n;
    double h2 = std::pow(h, 2);
    double a = -1/h2;
    double d = 2/h2;


    // print values a, d
    std::cout << a << ' ' << d << std::endl;

    // create A = triadiag(a,d,a)
    arma::mat A = create_symmetric_tridiagonal(N, a, d);


    // print to check
    std::cout << A << std::endl;

    // find eigenvals, eigenvecs
    arma::vec eigval; 
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, A);

    
    std::cout << eigval << std::endl;
    std::cout << eigvec << std::endl;

    // analytical solutions

    arma::vec lambda = arma::vec(N);
    arma::mat v = arma::mat(N, N);
    


    for(int i=1; i<=N; i++){
        lambda(i-1) = d + 2*a*std::cos(i*pi/(N+1));
        for(int j=1; j<=N; j++){
            v(j-1,i-1) = std::sin(j * i*pi/(N+1));
        }
    }
    arma::mat vnorm =arma::normalise(v);

    std::cout << lambda << std::endl;
    std::cout << vnorm << std::endl;


    return 0;

}
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



// A function that finds the max off-diag element of a symmetric matrix A.
// - The matrix indices of the max element are returned by writing to the  
//   int references k and l (row and column, respectively)
// - The value of the max element A(k,l) is returned as the function
//   return value
double max_offdiag_symmetric(const arma::mat& A, int& k, int& l){
    // Get size of the matrix A. Use e.g. A.n_rows, see the Armadillo documentation
    int N = A.n_rows;
    
    // Possible consistency checks:
    // Check that A is square and larger than 1x1. Here you can for instance use A.is_square(), 
    // see the Armadillo documentation.
    
    //  
    
    // The standard function 'assert' from <assert.h> can be useful for quick checks like this
    // during the code development phase. Use it like this: assert(some condition),
    // e.g assert(a==b). If the condition evaluates to false, the program is killed with 
    // an assertion error. More info: https://www.cplusplus.com/reference/cassert/assert/

    // Initialize references k and l to the first off-diagonal element of A
    k = 0;
    l = 1;

    // Initialize a double variable 'maxval' to A(k,l). We'll use this variable
    // to keep track of the largest off-diag element.
    double maxval = A(k, l);

    // Loop through all elements in the upper triangle of A (not including the diagonal)
    // When encountering a matrix element with larger absolute value than the current value of maxval
    // update k, l and max accordingly.
    for(int i=0; i<N-1; i++){
        for(int j=i+1; j<N; j++){
            double element = std::abs(A(i, j));
            if(element > maxval){
                maxval = element;
                k = i;
                l = j;
            }
        }
    }

    // Return maxval 
    return maxval;
}



// Performs a single Jacobi rotation, to "rotate away"
// the off-diagonal element at A(k,l).
// - Assumes symmetric matrix, so we only consider k < l
// - Modifies the input matrices A and R
void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l){
    double tau = (A(l,l)-A(k,k))/(2*A(k,l));
    double t;
    double c;
    double s;
    double rootoneplustausquared = std::pow(1+std::pow(tau,2),1/2);
    if(tau>0){
        t = 1/(tau + rootoneplustausquared);
    }
    else if(tau<0){
        t = -1/(-tau + rootoneplustausquared);
    }
    //  Could assert tau=0 or something here. 
    c = 1 / (std::pow(1+std::pow(t,2),1/2));
    s = c*t;

    //  Transform A-matrix
    A(k,k) = A(k,k) * std::pow(c,2) - 2*A(k,l)*c*s + A(l,l) * std::pow(s,2);
    A(l,l) = A(l,l) * std::pow(c,2) + 2*A(k,l)*c*s + A(k,k) * std::pow(s,2);
    A(k,l) = 0;
    A(l,k) = 0;

    for(int i=0; i<A.n_rows; i++){
        if(i!=k && i!= l){
            double A_ik = A(i,k);
            double A_il = A(i,l);
            A(i,k) = A_ik * c - A_il*s;
            A(k,i) = A(i,k);
            A(i,l) = A_il * c + A_ik * s;
            A(l,i) = A(i,l);
        }
        //  Transform R-matrix
        double R_ik = R(i,k);
        double R_il = R(i,l);
        R(i,k) = R_ik * c - R_il * s;
        R(i,l) = R_il * c - R_ik * s;
    }
}

// Jacobi method eigensolver:
// - Runs jacobo_rotate until max off-diagonal element < eps
// - Writes the eigenvalues as entries in the vector "eigenvalues"
// - Writes the eigenvectors as columns in the matrix "eigenvectors"
//   (The returned eigenvalues and eigenvectors are sorted using arma::sort_index)
// - Stops if it the number of iterations reaches "maxiter"
// - Writes the number of iterations to the integer "iterations"
// - Sets the bool reference "converged" to true if convergence was reached before hitting maxiter
void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, const int maxiter, int& iterations, bool& converged);





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
    //std::cout << a << ' ' << d << std::endl;

    // create A = triadiag(a,d,a)
    arma::mat A = create_symmetric_tridiagonal(N, a, d);


    // print to check
    //std::cout << A << std::endl;

    // find eigenvals, eigenvecs
    arma::vec eigval; 
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, A);
    arma::mat eigvecnorm = arma::normalise(eigvec);
    
    //std::cout << eigval << std::endl;
    //std::cout << eigvecnorm << std::endl;

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

    //std::cout << lambda << std::endl;
    //std::cout << vnorm << std::endl;

    // PROBLEM 3


    arma::mat A3 = arma::mat(4, 4, arma::fill::eye);
    A3(0,3) = 0.5;
    A3(1,2) = -0.7;
    A3(2,1) = -0.7;
    A3(3,0) = 0.5;
    int k = 0;
    int l = 1;
    double maxval = max_offdiag_symmetric(A3, k, l);

    std::cout << maxval << ' ' << k << ' ' << l << std::endl;

    return 0;

}
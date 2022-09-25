#include <iostream>
#include <armadillo>
#include <cmath>
#include <assert.h>


const double pi = 3.14159265359;

// Boundary conditions
const double xhat_min = 0;
const double xhat_max = 1;
const double v_min = 0; // dont know if these names make sense
const double v_max = 0; // maybe v_xhat_max


// Notation:
//      n: number of steps in our discretization
//      N: size of quadratic matrix A
//  N = n - 1
// (no. of points: n+1)



double step_size(int n){
    // Return step size for a given number of steps 
    double h=(xhat_max - xhat_min)/n;
    return h;
}



// Create tridiagonal matrix from vectors.
//  - subdiagonal:      vector a, lenght N-1
//  - main diagonal:    vector d, lenght N
//  - superdiagonal:    vector e, lenght N-1
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


// Create a tridiagonal matrix tridiag(a,d,e) of size N*N from scalar input a, d and e
arma::mat create_tridiagonal(int N, double a, double d, double e){
    arma::vec a_vec = arma::vec(N-1, arma::fill::ones) * a;
    arma::vec d_vec = arma::vec(N,   arma::fill::ones) * d;
    arma::vec e_vec = arma::vec(N-1, arma::fill::ones) * e;
    
    // Call the vector version of this function and return the result
    return create_tridiagonal(a_vec, d_vec, e_vec);
}

// Create a symmetric tridiagonal matrix tridiag(a,d,a) of size N*N from scalar input a and d.
arma::mat create_symmetric_tridiagonal(int N, double a, double d){
    // Call create_tridiagonal and return the result
    return create_tridiagonal(N, a, d, a);
}


double max_offdiag_symmetric(const arma::mat& A, int& k, int& l){
    // Get size of the matrix A.
    int N = A.n_rows;

    // Consistency checks:
    assert(A.is_square());   
    assert(A.is_symmetric()); 

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


// Test funtion for max_offdiag_symmetric()
int test_max_offdiag_symmetric(){
    // Creating test matrix
    arma::mat A = arma::mat(4, 4, arma::fill::eye);
    A(0,3) = 0.5;
    A(1,2) = -0.7;
    A(2,1) = -0.7;
    A(3,0) = 0.5;
    
    // Calling the function to test
    int k;
    int l;
    double maxval = max_offdiag_symmetric(A, k, l);

    // Asserting expected vs. computet values
    assert(maxval == 0.7);
    assert(k == 1);
    assert(l == 2);

    return 0;
}


// Find eigenvalues and -vectors for a tridiagonal symmetric matrix A with signature (a,d,a)
int analytical_eigenproblem(arma::mat A, arma::vec& eigval, arma::mat& eigvec){
    int N = A.n_rows;
    // Consistency checks:
    assert(A.is_square());   
    assert(A.is_symmetric()); 

    // Get signature
    double d = A(0,0);
    double a = A(0,1);

    for(int i=1; i<=N; i++){
        eigval(i-1) = d + 2*a*std::cos(i*pi/(N+1));
        arma::vec v = arma::vec(N);
        for(int j=1; j<=N; j++){
            v(j-1) = std::sin(j * i*pi/(N+1));
        eigvec.col(i-1) = arma::normalise(v);
        }
    }

    return 0;
}

// Check results for the 6x6 tridiagonal symmetric matrix A with signature (a,d,a)
//      (Comparing computed eigenvalues and -vectors using Armadillo- and/or Jacobi-solver.)
int check_for_babycase(std::string which="arma"){

    // Matrix A for N = 6
    int N = 6;
    double h = step_size(N+1);
    double h2 = std::pow(h, 2);
    double a = -1/h2;
    double d = 2/h2;
    arma::mat A = create_symmetric_tridiagonal(N, a, d);

    // Eigenvectors, eigenvalues with analytical expressions
    arma::vec eigval = arma::vec(N);
    arma::mat eigvec = arma::mat(N, N);
    analytical_eigenproblem(A, eigval, eigvec);


    arma::vec eigval_test; 
    arma::mat eigvec_test;
    // Check with Armadillo
    if(which == "arma"){
        arma::eig_sym(eigval_test, eigvec_test, A);   
        std::cout << "Checking if we use arma::eig_sym correctly." << std::endl;
    }
    // Check with Jacobi algorithm
    else if(which == "jacobi"){
        int p; // do something else
        std::cout << "Checking if we implement Jacobi rotation method correctly." << std::endl;
    }
    else if(which == "both"){
        int p; // do both
        std::cout << "Checking if we use arma::eig_sym and implement Jacobi rotation method correctly." << std::endl;
    }
    else{
        std::cout << "Provide valid argument" << std::endl;
    }

    // Check if they are equal
    arma::vec vals = eigval_test/eigval;
    arma::mat vecs = eigvec_test/eigvec;
    vecs = arma::abs(vecs);

    double tol = 0.0000001;
    bool is_vecs = arma::approx_equal(vecs, arma::mat(N,N).fill(1.), "absdiff",  tol);
    bool is_vals = arma::approx_equal(vals, arma::vec(N).fill(1.),   "absdiff", tol);

    assert(is_vecs);
    assert(is_vals);

    return 0;

}




int main(){

    // PROBLEM 2
    check_for_babycase();

    // PROBLEM 3
    test_max_offdiag_symmetric();

    

    return 0;
}
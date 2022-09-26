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
    A(0,0) = d(0);
    A(0,1) = e(0);
    
    // Loop that fills rows 2 to N-1 (row indices 1 to N-2)
    for(int i=1; i<N-1; i++){
        A(i, i) = d(i);
        A(i, i-1) = a(i);
        A(i, i+1) = e(i);
    }

    // Fill last row (row index N-1)
    A(N-1,N-1) = d(N-1);
    A(N-1,N-2) = a(N-2);
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
    // std::cout << "val:" << maxval << ", k:" << k << ", l:" << l << std::endl;
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
int analytical_eigenproblem(const arma::mat A, arma::vec& eigval, arma::mat& eigvec){
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

// Performs a single Jacobi rotation, to "rotate away"
// the off-diagonal element at A(k,l).
// - Assumes symmetric matrix, so we only consider k < l
// - Modifies the input matrices A and R
void jacobi_rotate(arma::mat &A, arma::mat &R, int k, int l){
    double tau;
    double t;
    double c;
    double s;
    // std::cout << "DEBUG: " << __FILE__ << ":" << __LINE__ << std::endl;
    // std::cout << k << " " << l << std::endl;
    if (A(k, l) == 0){
        // std::cout << "if statement" << std::endl;
        c = 1;
        s = 0;
        t = 0;
    }
    else{
        // std::cout << "else" << std::endl;
        tau = (A(l, l) - A(k, k)) / (2 * A(k, l));
        if (tau > 0){
            // t = -tau + sqrt(1+tau*tau);
            t = 1 / (tau + sqrt(1 + tau*tau));
        }
        else{
            // t = -tau - sqrt(1+tau*tau);
            t = -tau - sqrt(1+tau*tau);
        }
        c = 1 / sqrt(1 + t*t);
        s = c * t;
    }
    // std::cout << "tau: " << tau << std::endl;
    // std::cout << "t: " << t << std::endl;
    // std::cout << "c: " << c << std::endl;
    // std::cout << "s: " << s << std::endl;

    //  Transform A-matrix
    double A_kk = A(k,k);
    double A_ll = A(l,l);
    A(k, k) = A_kk * c*c - 2 * A(k, l) * c * s + A_ll * s * s;
    A(l, l) = A_ll * c*c + 2 * A(k, l) * c * s + A_kk * s * s;
    A(k, l) = 0;
    A(l, k) = 0;

    for (int i = 0; i < A.n_rows; i++){
        if (i != k && i != l){
            // std::cout << i << " " << l << " " << k << std::endl;
            double A_ik = A(i, k);
            double A_il = A(i, l);
            A(i, k) = A_ik * c - A_il * s;
            A(k, i) = A(i, k);
            A(i, l) = A_il * c + A_ik * s;
            A(l, i) = A(i, l);
        }
        //  Transform R-matrix
        double R_ik = R(i, k);
        double R_il = R(i, l);
        R(i, k) = R_ik * c - R_il * s;
        R(i, l) = R_il * c + R_ik * s;
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
void jacobi_eigensolver(arma::mat A, double eps, arma::vec &eigenvalues, arma::mat &eigenvectors, const int maxiter, int &iterations, bool &converged){
    int k;
    int l;
    int N = A.n_rows;
    arma::mat R = arma::mat(N, N, arma::fill::eye);
    int iter = 0;
    // double calc_maxval = max_offdiag_symmetric(A, k, l);
    // loop
    // std::cout << "DEBUG: " << __FILE__ << ":" << __LINE__ << std::endl;
    while (max_offdiag_symmetric(A, k, l) > eps){
        if (iter >= maxiter){
            break;
        }
        else{
            // std::cout << "DEBUG: " << __FILE__ << ":" << __LINE__ << std::endl;
            jacobi_rotate(A, R, k, l);
            // std::cout << "Iterations: " << iter+1 << std::endl;
            // std::cout << A << std::endl;
            // std::cout << R << std:: endl;
        }
        // std::cout << iter << std::endl;
        iter++;
        // calc_maxval = max_offdiag_symmetric(A, k, l);
        // std::cout << "Maxval: " << calc_maxval << std::endl;
    }
    // std::cout << "out of loop" << std::endl;
    iterations = iter;
    // Check for convergens or iteration stop.
    if (iterations < maxiter){
        converged = true;
    }
    else{
        converged = false;
    }
    // std::cout << "after bool" << std::endl;
    // std::cout << "DEBUG: " << __FILE__ << ":" << __LINE__ << std::endl;
    for (int i = 0; i < N; i++){
        eigenvalues(i) = A(i, i);
        // std::cout << i << std::endl;
        eigenvectors.col(i) = R.col(i);
    }
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
    int iterations;
    bool converged;
    // std::cout << A << std::endl;

    // Eigenvectors, eigenvalues with analytical expressions
    arma::vec eigval = arma::vec(N);
    arma::mat eigvec = arma::mat(N, N);
    analytical_eigenproblem(A, eigval, eigvec);
    // arma::eig_sym(eigval, eigvec, A);

    // std::cout << "DEBUG: " << __FILE__ << ":" << __LINE__ << std::endl;
    arma::vec eigval_test = arma::vec(N); 
    arma::mat eigvec_test = arma::mat(N, N);
    // Check with Armadillo
    if(which == "arma"){
        arma::eig_sym(eigval_test, eigvec_test, A);   
        std::cout << "Checking if we use arma::eig_sym correctly." << std::endl;
    }
    // Check with Jacobi algorithm
    else if(which == "jacobi"){
        // std::cout << "DEBUG: " << __FILE__ << ":" << __LINE__ << std::endl;
        jacobi_eigensolver(A, 1e-9, eigval_test, eigvec_test, 10000, iterations, converged);
        // std::cout << A << std::endl;
        std::cout << "Checking if we implement Jacobi rotation method correctly." << std::endl;
    }
    else if(which == "both"){
        int p; // do both
        std::cout << "Checking if we use arma::eig_sym and implement Jacobi rotation method correctly." << std::endl;
    }
    else{
        std::cout << "Provide valid argument" << std::endl;
    }
    std::cout << "Using jacobi with iter: " << iterations << std::endl;
    // std::cout << A << std::endl;
    std::cout << eigval_test << std::endl;
    std::cout << eigvec_test << std::endl;

    std::cout << "Analytical:" << std::endl;
    std::cout << eigval << std::endl;
    std::cout << eigvec << std::endl;

    // Check if they are equal
    arma::vec vals = eigval_test/eigval;
    arma::mat vecs = eigvec_test/eigvec;
    vecs = arma::abs(vecs);

    double tol = 0.0000001;
    bool is_vecs = arma::approx_equal(vecs, arma::mat(N,N).fill(1.), "absdiff",  tol);
    bool is_vals = arma::approx_equal(vals, arma::vec(N).fill(1.),   "absdiff", tol);

    // assert(is_vecs);
    assert(is_vals);

    return 0;

}




int main(){

    // PROBLEM 2
    // std::cout << "DEBUG: " << __FILE__ << ":" << __LINE__ << std::endl;
    check_for_babycase("jacobi");

    // PROBLEM 3
    // test_max_offdiag_symmetric();

    

    return 0;
}
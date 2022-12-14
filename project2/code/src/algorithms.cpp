#include "algorithms.hpp"

const double pi = 3.14159265359;


arma::mat create_tridiagonal(const arma::vec& a, const arma::vec& d, const arma::vec& e){
    /**
     * Create tridiagonal matrix from vectors.
     *      subdiagonal:      vector a, lenght N-1
     *      main diagonal:    vector d, lenght N
     *      superdiagonal:    vector e, lenght N-1
     */
    
    //  start from identity matrix:
    int N = d.size();
    arma::mat A = arma::mat(N, N, arma::fill::eye);
    
    //  fill first row:
    A(0,0) = d(0);
    A(0,1) = e(0);
    
    //  loop that fills rows 2 to N-1 (row indices 1 to N-2):
    for(int i=1; i<N-1; i++){
        A(i, i) = d(i);
        A(i, i-1) = a(i);
        A(i, i+1) = e(i);
    }

    //  fill last row (row index N-1):
    A(N-1,N-1) = d(N-1);
    A(N-1,N-2) = a(N-2);
    return A;
}


arma::mat create_tridiagonal(int N, double a, double d, double e){
    /**
     * Create a tridiagonal matrix tridiag(a,d,e) of size N*N from scalar input a, d and e
     */
    arma::vec a_vec = arma::vec(N-1, arma::fill::ones) * a;
    arma::vec d_vec = arma::vec(N,   arma::fill::ones) * d;
    arma::vec e_vec = arma::vec(N-1, arma::fill::ones) * e;
    return create_tridiagonal(a_vec, d_vec, e_vec);
}


arma::mat create_symmetric_tridiagonal(int N, double a, double d){
    /**
     * Create a symmetric tridiagonal matrix tridiag(a,d,a) of size N*N from scalar input a and d.
     */
    return create_tridiagonal(N, a, d, a);
}


double max_offdiag_symmetric(const arma::mat& A, int& k, int& l){
    /**
     * Find maximal value, not counting the diagonal, of a symmetric matrix.
     */

    //  get size of the matrix A:
    int N = A.n_rows;

    // Consistency checks
    assert(A.is_square());   
    assert(A.is_symmetric()); 

    // Find max value
    //  initialize references k and l to the first off-diagonal element of A:
    k = 0;
    l = 1;

    //  initialize a double variable 'maxval' to A(k,l):
    double maxval = std::abs(A(k, l));
   
    /**
     * Loop through all elements in the upper triangle of A (not including the diagonal)
     *  When encountering a matrix element with larger absolute value than the current value of maxval
     *  update k, l and max accordingly.
     */
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
    return maxval;
}


int analytical_eigenproblem(const arma::mat A, arma::vec& eigval, arma::mat& eigvec){
    /**
     * Find eigenvalues and -vectors for a tridiagonal symmetric matrix A with signature (a,d,a)
     */

    // Consistency checks
    assert(A.is_square());   
    assert(A.is_symmetric()); 

    // Analytical formulas from project description
    int N = A.n_rows;

    //  get signature:
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


void jacobi_rotate(arma::mat &A, arma::mat &R, int k, int l){
    /**
     * Performs a single Jacobi rotation, to "rotate away" the off-diagonal element at A(k,l).
     *      assumes symmetric matrix, so we only consider k < l
     *      modifies the input matrices A and R
     */
    double tau;
    double t;
    double c;
    double s;
    if (A(k, l) == 0){
        c = 1;
        s = 0;
        t = 0;
    }
    else{
        tau = (A(l, l) - A(k, k)) / (2 * A(k, l));
        if (tau > 0){
            t = 1 / (tau + sqrt(1 + tau*tau));
        }
        else{
            t = -1 /(-tau + sqrt(1 + tau*tau));
        }
        c = 1 / sqrt(1 + t*t);
        s = c * t;
    }

    //  transform A-matrix:
    double A_kk = A(k,k);
    double A_ll = A(l,l);
    A(k, k) = A_kk * c*c - 2 * A(k, l) * c * s + A_ll * s * s;
    A(l, l) = A_ll * c*c + 2 * A(k, l) * c * s + A_kk * s * s;
    A(k, l) = 0;
    A(l, k) = 0;

    for (int i = 0; i < A.n_rows; i++){
        if (i != k && i != l){
            double A_ik = A(i, k);
            double A_il = A(i, l);
            A(i, k) = A_ik * c - A_il * s;
            A(k, i) = A(i, k);
            A(i, l) = A_il * c + A_ik * s;
            A(l, i) = A(i, l);
        }
        //  transform R-matrix:
        double R_ik = R(i, k);
        double R_il = R(i, l);
        R(i, k) = R_ik * c - R_il * s;
        R(i, l) = R_il * c + R_ik * s;
    }
}



void jacobi_eigensolver(arma::mat A, double eps, arma::vec &eigenvalues, arma::mat &eigenvectors, const int maxiter, int &iterations, bool &converged){
    /**
     * Jacobi method eigensolver:
     *  - Runs 'jacobo_rotate' until max off-diagonal element < 'eps'
     *  - Writes the eigenvalues as entries in the vector 'eigenvalues'
     *  - Writes the eigenvectors as columns in the matrix 'eigenvectors'
     *  - (The returned eigenvalues and eigenvectors are sorted using arma::sort_index)
     *  - Stops if it the number of iterations reaches 'maxiter'
     *  - Writes the number of iterations to the integer 'iterations'
     *  - Sets the bool reference 'converged' to true if convergence was reached before hitting maxiter
     */
    int k;
    int l;
    int N = A.n_rows;
    arma::mat R = arma::mat(N, N, arma::fill::eye);
    int iter = 0;
    while (max_offdiag_symmetric(A, k, l) > eps){
        if (iter >= maxiter){
            break;
        }
        else{
            jacobi_rotate(A, R, k, l);
            iter++;
        }
    }
    iterations = iter;
    
    //  check for convergens or iteration stop:
    if (iterations < maxiter){
        converged = true;
    }
    else{
        converged = false;
    }
    for (int i = 0; i < N; i++){
        eigenvalues(i) = A(i, i);
    }
    arma::uvec indecies = sort_index(eigenvalues);
    eigenvalues = sort(eigenvalues);
    for (int i =0; i < N; i++){
        eigenvectors.col(i) = R.col(indecies(i));
    }

}
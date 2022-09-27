#include <algorithm>
#include "utils.hpp"
#include "algorithms.hpp"

/**
 * PROJECT 2 FYS4150
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


arma::mat problem_matrix_A(int N){
    /**
     * Create problem specific, given N.
     *      A is NxN, tridiagonal, symmetric with signature (a,d,a) where
     *          a = -1/h^2
     *          d = 2/h^2
     *      and h is the step size.
     *  
     * @note N = n-1, where n is the number of discretisation steps
     */
    arma::mat A = arma::mat(N,N);
    double h = step_size(N+1);
    double h2 = std::pow(h, 2);
    double a = -1/h2;
    double d = 2/h2;
    A = create_symmetric_tridiagonal(N, a, d);
    return A;
}


int test_max_offdiag_symmetric(){
    /**
     * Test functrion for 'max_offdiag_symmetric()'.
     * 
     * Runs silently when successful.
     */

    //  create test matrix:
    arma::mat A = arma::mat(4, 4, arma::fill::eye);
    A(0,3) = 0.5;
    A(1,2) = -0.7;
    A(2,1) = -0.7;
    A(3,0) = 0.5;
    
    //  calling the function:
    int k;
    int l;
    double maxval = max_offdiag_symmetric(A, k, l);

    //  asserting expected vs. computet values
    assert(maxval == 0.7);
    assert(k == 1);
    assert(l == 2);

    return 0;
}


int check_for_babycase(std::string which="arma"){
    /**
     * Check results for problem specific A when A is 6x6
     * 
     * Comparing computed eigenvalues and -vectors using Armadillo- or Jacobi-solver.
     */

    //  declare matrix A for N = 6:
    int N = 6;
    arma::mat A = problem_matrix_A(N); 


    // declare and compute eigenvectors, eigenvalues with analytical expressions:
    arma::vec eigval = arma::vec(N);
    arma::mat eigvec = arma::mat(N, N);
    analytical_eigenproblem(A, eigval, eigvec);

    // declare eigenvectors, eigenvalues
    arma::vec eigval_test = arma::vec(N); 
    arma::mat eigvec_test = arma::mat(N, N);

    if(which == "arma"){
        //  compute eigenvectors, eigenvalues with Armadillo:
        arma::eig_sym(eigval_test, eigvec_test, A);   
        std::cout << "\nChecking if we use arma::eig_sym correctly." << std::endl;
    }
    else if(which == "jacobi"){
        int M; 
        bool converged; 
        //  compute eigenvectors, eigenvalues with Armadillo:
        jacobi_eigensolver(A, 1e-8, eigval_test, eigvec_test, 10000, M, converged);
        std::cout << "\nChecking if we implement Jacobi rotation method correctly." << std::endl;
    }
    else{
        std::cout << "Provide valid argument" << std::endl;
    }

    /**
     * Check to see if they are equal.
     *      ! Account for oppositely directed eigenvectors.
     */
    arma::mat r = eigvec_test/eigvec;
    for(int i=0; i<N; i++){
        bool opposite = arma::all(r.col(i)<0);
        if(opposite){
            //  change sign:
            eigvec_test.col(i) = - eigvec_test.col(i);
        }
    }

    //  assert:
    bool is_vecs = arma::approx_equal(eigvec_test, eigvec, "absdiff", eps);
    bool is_vals = arma::approx_equal(eigval_test, eigval, "absdiff", eps);
    assert(is_vecs);
    assert(is_vals);

    std::cout << "Check OK." << std::endl;

    return 0;
}


int run_jacobi_algorithms(int N_max, bool tridiag=true, int N_min=2){
    /**
     * Runs Jacobi algorithm for choices of N betweeen N_min and N_max, 
     * and saves corresponding M.
     * 
     * Option to run with a random symmetric matrix if 'tridiag' is false.
     */

    std::string type;
    if(tridiag==true){
        type = "tridiag";   
    }
    else{
        type = "dense"; 
    }

    std::string filename = "transformations_per_" + type + "_N_matrix";
    std::cout << "\nRunning Jacobi algorithm for N = " << N_min << ", ..., " << N_max << " (A is " << type << ")." << std::endl;

    //  declare empty vectors: 
    int L = N_max - N_min + 1;
    std::vector<double> Ns(L); // choices of N (N_min, N_min+1, ..., N_max-1, N_max)
    std::vector<double> Ms(L); // number of transformations needed
    
    for(int i=0; i<L; i++){
        int N = i + N_min;
        //  prepare solver:
        arma::mat A = arma::mat(N,N);
        if(tridiag==true){
            A = problem_matrix_A(N);
        }
        else{
            A = arma::mat(N,N).randn();
            A = arma::symmatu(A);
        }
        int M;
        bool converged;
        arma::vec eigval_j = arma::vec(N);
        arma::mat eigvec_j = arma::mat(N,N);
        //  solve and save relevant parameters:
        jacobi_eigensolver(A, eps, eigval_j, eigvec_j, M_max, M, converged);
        Ns[i] = N;
        Ms[i] = M;
        assert(converged); // to make sure the algorithm actually converged
    }

    write_to_file(Ns, Ms, filename);
    std::cout << "\nWrote result to '" << filename << ".txt'.\n\n"  << std::endl;

    return 0;

}


int write_solution_to_file(std::string method, arma::vec x, arma::mat V){
    /**
     * Write solution(s) of eigenvalue problem to file.
     */
    int npoints = V.n_rows;
    int n = npoints - 1;
    int save_first = V.n_cols;

    std::string filename = method + "_solution_" + std::to_string(n) + "steps";
    std::string path = "../output/data/"; // path for .txt files
    std::string file = path + filename + ".txt";
    std::ofstream ofile;

    int w = 20;
    int p = 15;

    ofile.open(file.c_str());

    std::string deli = ", ";
    for(int i=0; i<npoints; i++){
        ofile << scientific_format( x(i), w, p) << deli;
        for(int j=0; j<save_first; j++){
            ofile << scientific_format( V(i, j), w, p);
            if(j!=save_first-1){
                ofile << deli;
            }
            else{
                ofile << std::endl;
            }
        }
    }

    ofile.close();

    std::cout << "\nWrote solution to '" << filename << ".txt'.\n"<< std::endl;

    return 0;
}


int compute_solution(int n, int save_first=3){
    /**
     * Computes solution of the problem specific eigenvalue problem,
     * both using Jacobi algorithm and analytical formula.
     * 
     * Saves the first three eigenvectors to files.
     */

    //  initialise with boundary points
    arma::vec xhat = arma::vec(n+1);
    xhat(0) = xhat_0;
    xhat(n) = xhat_n;
    double h = step_size(n);
    for(int i=1; i<n; i++){
        xhat(i) = i*h + xhat_0;
    }

    //  create matrix A:
    int N = n-1;
    arma::mat A = problem_matrix_A(N);
    int M;
    bool converged;

    // Analytical solution
    arma::vec eigval_a = arma::vec(N);
    arma::mat eigvec_a = arma::mat(N, N);
    analytical_eigenproblem(A, eigval_a, eigvec_a);


    // Jacobi solution
    arma::vec eigval_J = arma::vec(N);
    arma::mat eigvec_J = arma::mat(N,N);
    jacobi_eigensolver(A, eps, eigval_J, eigvec_J, M_max, M, converged);

    /**
     * Find first (3) eigenvectors and save them 
     * (already sorted with increasing eigenvalues)
     */

    //  initialise:
    arma::mat V_a = arma::mat(n+1, save_first);
    arma::mat V_J = arma::mat(n+1, save_first);

    for(int i=0; i<save_first; i++){
        //  implement boundary conditions:
        V_a(0, i) = v_0;
        V_J(0, i) = v_0;
        V_a(n, i) = v_n;
        V_J(n, i) = v_n;
        //  fill with computed eigenvectors:
        V_a(arma::span(1, N), i) = eigvec_a.col(i);
        V_J(arma::span(1, N), i) = eigvec_J.col(i);
        
        //  account for oppositely directed vectors:
        arma::vec r = V_a.col(i)/V_J.col(i);
        r.shed_row(n);
        r.shed_row(0);
        bool opposite = arma::all(r<0);
        if(opposite){
            V_J.col(i) = - V_J.col(i);
        }  
    }
    //  save solutions to file:
    write_solution_to_file("analytical", xhat, V_a);
    write_solution_to_file("Jacobi", xhat, V_J);
    
    return 0;
}


int main(){

    // PROBLEM 2
    check_for_babycase("arma");

    // PROBLEM 3
    test_max_offdiag_symmetric();

    // PROBLEM 4
    check_for_babycase("jacobi");

    // PROBLEM 5
    //  a)
    run_jacobi_algorithms(100);
    //  b)
    run_jacobi_algorithms(100, false); 

    // PROBLEM 6
    //  a)
    compute_solution(10);
    //  b)
    compute_solution(100);


    return 0;
}
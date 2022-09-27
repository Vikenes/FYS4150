#include <algorithm>
#include "utils.hpp"
#include "algorithms.hpp"


//  boundary conditions
const double xhat_min = 0;
const double xhat_max = 1;
const double v_min = 0;
const double v_max = 0;

double step_size(int n){
    // Return step size for a given number of steps 
    double h=(xhat_max - xhat_min)/n;
    return h;
}

// Create problem specific tridiagonal symmetric matrix, given N
arma::mat create_tridiag(int N){
    arma::mat A = arma::mat(N,N);
    double h = step_size(N+1);
    double h2 = std::pow(h, 2);
    double a = -1/h2;
    double d = 2/h2;
    A = create_symmetric_tridiagonal(N, a, d);
    return A;
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


// Check results for the 6x6 tridiagonal symmetric matrix A with signature (a,d,a)
//      (Comparing computed eigenvalues and -vectors using Armadillo- and/or Jacobi-solver.)
int check_for_babycase(std::string which="arma"){

    // Matrix A for N = 6
    int N = 6;
    arma::mat A = create_tridiag(N); 
    int iterations;
    bool converged;


    // Eigenvectors, eigenvalues with analytical expressions
    arma::vec eigval = arma::vec(N);
    arma::mat eigvec = arma::mat(N, N);
    analytical_eigenproblem(A, eigval, eigvec);


    arma::vec eigval_test = arma::vec(N); 
    arma::mat eigvec_test = arma::mat(N, N);
    // Check with Armadillo
    if(which == "arma"){
        arma::eig_sym(eigval_test, eigvec_test, A);   
        std::cout << "\nChecking if we use arma::eig_sym correctly." << std::endl;
    }
    // Check with Jacobi algorithm
    else if(which == "jacobi"){
        jacobi_eigensolver(A, 1e-8, eigval_test, eigvec_test, 10000, iterations, converged);
        std::cout << "\nChecking if we implement Jacobi rotation method correctly." << std::endl;
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

    std::cout << "Check OK.\n" << std::endl;

    return 0;
}


int run_jacobi_algorithms(int N_max, bool tridiag=true, int N_min=2, double eps=(1e-8), int maxiter=int(1e6)){

    std::string type;

    if(tridiag==true){
        type = "tridiag";   
    }
    else{
        type = "dense"; 
    }

    std::string filename = "transformations_per_" + type + "_N_matrix";
    std::cout << "\nRunning Jacobi algorithm for N = " << N_min << ", ..., " << N_max << " (A is " << type << ")." << std::endl;

    int L = N_max - N_min + 1;
    std::vector<double> Ns(L); // choices of N (N_min, N_min+1, ..., N_max-1, N_max)
    std::vector<double> n_transf(L);  // number of transformations needed
    
    
    for(int i=0; i<L; i++){
        int N = i + N_min;
        arma::mat A = arma::mat(N,N);
        
        if(tridiag==true){
            A = create_tridiag(N);
        }
        else{
            A = arma::mat(N,N).randn();
            A = arma::symmatu(A);
        }
        int iterations;
        bool converged;

        arma::vec eigval_j = arma::vec(N);
        arma::mat eigvec_j = arma::mat(N,N);
        jacobi_eigensolver(A, eps, eigval_j, eigvec_j, maxiter, iterations, converged);
        Ns[i] = N;
        n_transf[i] = iterations;
        //assert(converged) //?
    }

    write_to_file(Ns, n_transf, filename);
    std::cout << "\nWrote result to '" << filename << ".txt'.\n\n"  << std::endl;

    return 0;

}



int write_solution_to_file(std::string method, arma::vec x, arma::mat V){
    int npoints = V.n_rows;
    int n = npoints - 1;
    int save_first = V.n_cols;

    std::string filename = method + "_solution_" + std::to_string(n) + "steps";

    std::string path = "../output/data/"; // path for .txt files
    std::string file = path + filename + ".txt";
    std::ofstream ofile;

    int w = 15;
    int p = 10;

    ofile.open(file.c_str());

    std::stringstream ss;
    ss << std::setw(w) << std::setprecision(p) << std::scientific;

    for(int i=0; i<npoints; i++){
        ofile << scientific_format( x(i), w, p) << ", ";
        for(int j=0; j<save_first; j++){
            ofile << scientific_format( V(i, j), w, p) << ", ";
        }
        ofile << std::endl;
    }

    ofile.close();
    std::cout << "\nWrote solution to '" << filename << ".txt'.\n"<< std::endl;

    return 0;
}


int compute_solution(int n=10, int save_first=3, double eps=(1e-8), int maxiter=int(1e6)){
    // Initialise with boundary points
    arma::vec xhat = arma::vec(n+1);
    xhat(0) = xhat_min;
    xhat(n) = xhat_max;
    double h = step_size(n);
    for(int i=1; i<n; i++){
        xhat(i) = i*h + xhat_min;
    }

    int N = n-1;
    arma::mat A = create_tridiag(N);
    int iterations;
    bool converged;


    // Jacobi solution
    arma::vec eigval_J = arma::vec(N);
    arma::mat eigvec_J = arma::mat(N,N);
    jacobi_eigensolver(A, eps, eigval_J, eigvec_J, maxiter, iterations, converged);

    // Analytical solution
    arma::vec eigval_a = arma::vec(N);
    arma::mat eigvec_a = arma::mat(N, N);
    analytical_eigenproblem(A, eigval_a, eigvec_a);

    // Already sorted with increasing eigenvalues




    arma::mat V_J = arma::mat(n+1, save_first);
    arma::mat V_a = arma::mat(n+1, save_first);
    
    for(int i=0; i<save_first; i++){
        V_J( arma::span(1, N), i ) = eigvec_J.col(i);
        V_a( arma::span(1, N), i ) = eigvec_a.col(i);
    }

    write_solution_to_file("Jacobi", xhat, V_J);
    write_solution_to_file("analytical", xhat, V_a);
    
    return 0;
}

int main(){

    // PROBLEM 2
    //check_for_babycase("arma");

    // PROBLEM 3
    //test_max_offdiag_symmetric();

    // PROBLEM 4
    //check_for_babycase("jacobi");

    // PROBLEM 5
    //run_jacobi_algorithms(100);
    //run_jacobi_algorithms(100, false);

    // PROBLEM 6
    compute_solution();

    return 0;
}
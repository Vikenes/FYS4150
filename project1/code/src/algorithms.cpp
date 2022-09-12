#include "algorithms.hpp"

double u(const double x){
    return 1 - (1 - exp(-10))*x - exp(-10*x);
}

double f(const double x){
    return 100 * exp(-10*x);
}



std::vector<double> generalThomas(std::vector<double> a, 
                                std::vector<double> b, 
                                std::vector<double> c, 
                                std::vector<double> g,
                                double v_min, double v_max){
    // Thomas algorithm (general)
    int m = b.size();
    std::vector<double> btilde(m);
    std::vector<double> gtilde(m);

    btilde[0] = b[0];
    gtilde[0] = g[0];
     
    // Step 1 - Forward

    for(int i=1; i<m; i++){
        double K = a[i]/btilde[i-1];
        btilde[i] = b[i] - K*c[i-1];
        gtilde[i] = g[i] - K*gtilde[i-1];
    }

    // Step 2 - Backward

    std::vector<double> vstar(m+2);
    // std::vector<double> vstar(m);
    vstar[0] = v_min; 
    vstar[m+1] = v_max; 

    vstar[m] = gtilde[m-1]/btilde[m-1];

    for(int i=m-2; i>=0; i--){
        vstar[i+1] = (gtilde[i] - c[i]*vstar[i+2]) / btilde[i];
    }

    return vstar;
}

std::vector<double> specialThomas(std::vector<double> g, 
                                    double v_min,
                                    double v_max, 
                                    int b_tilde0){
    // Thomas algorithm (general)
    int m = g.size();
    std::vector<double> btilde(m);
    std::vector<double> gtilde(m);

    btilde[0] = b_tilde0;
    gtilde[0] = g[0];
     
    // Step 1 - Forward

    for(int i=1; i<m; i++){
        btilde[i] = (i+2.0)/(i+1);
        gtilde[i] = g[i] + gtilde[i-1]/btilde[i-1];
    }

    // Step 2 - Backward

    std::vector<double> vstar(m+2);

    vstar[0] = v_min;
    vstar[m+1] = v_max;

    vstar[m] = gtilde[m-1]/btilde[m-1];

    for(int i=m-2; i>=0; i--){
        vstar[i+1] = (gtilde[i] + vstar[i+2]) / btilde[i];
    }

    return vstar;

}

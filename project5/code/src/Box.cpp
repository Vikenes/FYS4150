#include "Box.hpp"


Box::Box(double spatial_step_size){
    h = spatial_step_size;
    M = int(1/h) + 1; 
    V = arma::sp_mat(M,M);
}


void Box::create_slits(int num_of_slits, double v0, double aperture, double wall_width, double wall_height, double horisontal_centre, double vertical_centre){

    assert(num_of_slits > 0);

    int num_of_walls = num_of_slits + 1;
    std::cout << num_of_walls << std::endl;

    double sep = aperture + wall_height;
    double x_c = horisontal_centre;
    /*
    Stupid way of solving this problem:
    (There has to be a smoother way??)c
    */
    std::vector<double> y_c(num_of_walls); 
    std::vector<double> y_ll(num_of_walls); 
    std::vector<double> y_ur(num_of_walls); 

    y_c[0] = vertical_centre - sep*num_of_slits/2;

    for(int k=0; k<num_of_walls; k++){
        y_c[k] = y_c[0] + k*sep;
        y_ll[k] = y_c[k] - wall_height/2;   // (y) lower left corner
        y_ur[k] = y_c[k] + wall_height/2;   // (y) upper right corner
    }

    double x_ll = x_c-wall_width/2;     // (x) lower left corner
    double x_ur = x_c+wall_width/2;     // (x) upper right corner

    arma::vec x = arma::linspace(0, 1, M);
    arma::vec y = arma::linspace(0, 1, M);

    arma::uvec x_indices = arma::find(x >= x_ll && x <= x_ur); // <= gives right dimensions!
    // check actual width:
    arma::vec xwall = x(x_indices);
    arma::vec ywall;
    double width = xwall(x_indices.n_elem-1) - xwall(0);
    double height;
    

    // Woops, correct misunderstanding:
    y_ll[0] = 0;
    y_ur[num_of_slits] = 1;

    /* 
    If somebody could help me find a way to run through these uvec-s, that would be sprutnice :p
    */
    for(int w=0; w<num_of_walls; w++){
        arma::uvec y_indices = arma::find(y >= y_ll[w] && y <= y_ur[w]);
        ywall = y(y_indices);
        height = ywall(y_indices.index_max()) - ywall(y_indices.index_min());
        for(int i=0; i<x_indices.n_elem; i++){
            for(int j=0; j<y_indices.n_elem; j++){
                V(x_indices[i],y_indices[j]) = v0;
            }
        }
    }

    
    std::cout << "Sat up " << num_of_slits << " slits using " << num_of_walls << " walls." << std::endl;
    std::cout << "  -> width  (x-dir): " << width  << std::endl;
    std::cout << "  -> height (y-dir): " << height << std::endl;


}
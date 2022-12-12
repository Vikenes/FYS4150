#include "Box.hpp"
// #include "Simulation.hpp"
// #include "utils.hpp"

Box::Box(double ha){
    h = ha;
    M = int(1/ha) +1;
    V = arma::sp_mat(M,M);
}


// void Box::set_up_walls(void){
//     set_up_walls_utils(V, v0, M, h);
// }

void Box::set_up_walls(double v0, int Ns, double Th, double wxc, double wyc, double ws_length, double aperture){
    /**
     * @brief Sets up the wall used to generate the slit apertures.
     * @param V reference to potential matrix V.
     * @param v0 constant potential value of wall.
     * @param M number of points in x and y direction.
     * @param h step length between points.
     * @param Ns number of slits.
     * @param Th thickness of slits.
     * @param wxc slit wall x-centre.
     * @param wyc slit wall y-centre.
     * @param ws_length slit wall segment length in y-direction
     * @param aperture slit aperture width in y-direction
     * 
     */
    std::cout<<"Setting up wall with "<< Ns << " slits, and v0: "<<v0<<std::endl;



    for(int i=1; i<=M-2; i++){
        for(int j=1; j<=M-2; j++){
            double x = i*h;
            double y = j*h;
            // Ignore number of slits for now, just assume that it is two
            
            // First set x-criterion
            if(x<=wxc+Th/2 && x>wxc-Th/2){
                // Then set the various y-criteria
                // Will be a loop here if you include varying number of slits (must differ between odd and even)
                // First slit wall in centre
                if(y<=wyc+ws_length/2 && y>wyc-ws_length/2){
                    V(i,j) = v0;
                }
                // Wall on side of larger y-vals
                else if(y<=wyc+3/2*ws_length+aperture && y>wyc+1/2*ws_length+aperture){
                    V(i,j) = v0;
                }
                // Wall on side of smaller y-vals
                else if(y<=wyc-1/2*ws_length-aperture && y>wyc-3/2*ws_length-aperture){
                    V(i,j) = v0;
                }
            }
        }
    }
    // Not sure if this needs to be imposed here
    V.col(0).fill(0);
    V.col(M-1).fill(0);
    V.row(0).fill(0);
    V.row(M-1).fill(0);
    std::cout<<"-> potential wall set up"<<std::endl;
}


void Box::create_slits(int num_of_slits, double v0, double aperture, double wall_width, double wall_height, double horisontal_centre, double vertical_centre){

    assert(num_of_slits > 0);

    int num_of_walls = num_of_slits + 1;
    std::cout << num_of_walls << std::endl;

    double sep = aperture + wall_height;
    double x_c = horisontal_centre;
    /*
    Stupid way of solving this problem:
    (There has to be a smoother way??)
    */
    std::vector<double> y_c(num_of_walls); 
    std::vector<double> y_ll(num_of_walls); 
    std::vector<double> y_ur(num_of_walls); 

    y_c[0] = vertical_centre - sep*num_of_slits/2;
    // std::cout << y_c[0] << std::endl;
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

    // arma::uvec y_indices = arma::find(y >= y_ll[0] && y <= y_ur[0]);

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
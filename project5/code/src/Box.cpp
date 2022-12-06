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

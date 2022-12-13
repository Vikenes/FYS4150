// Include stuff
#include <iostream>
#include "utils.hpp"
#include "Box.hpp"
#include "Simulation.hpp"


/**
 * @brief Main (testing) program for FYS4150 project5 code.
 * Vetle A. Vikenes, Nanna Bryne, Johan Mylius Kroken
 * 
 */

void no_slit(void){
    // double h = 0.005;
    // double Dt = 2.5e-5;
    // double T = 0.008;
    // double xc = 0.25;
    // double sigma_x = 0.05;
    // double p_x = 200;
    // double yc = 0.5;
    // double sigma_y = 0.05;
    // double p_y = 0;
    double v0 = 0;
    Box b1 = Box(0.005);
    b1.set_up_walls(v0);
    Simulation s1 = Simulation(b1);
    s1.initialise();
    
    std::cout<<"\nExperiment: NO SLIT\n"<<std::endl;
    arma::cx_cube U1 = s1.run_simulation();
    U1.save("../output/binfiles/NS_arma_cube.bin");
}

void double_slit_broad_sigma_y(void){
    double h = 0.005;
    double Dt = 2.5e-5;
    double T = 0.008;
    double xc = 0.25;
    double sigma_x = 0.05;
    double p_x = 200;
    double yc = 0.5;
    double sigma_y = 0.10;  //  Broader sigma
    double p_y = 0;
    double v0 = 1e10;   // High potential to set up the slits
    // double v0 = 1e10;
    Box b2 = Box(h);
    b2.set_up_walls(v0);
    Simulation s2 = Simulation(b2, Dt, T, xc, sigma_x, p_x, yc, sigma_y, p_y);
    s2.initialise();

    std::cout<<"\nExperiment: DOUBLE-SLIT\n"<<std::endl;
    arma::cx_cube U2 = s2.run_simulation();
    U2.save("../output/binfiles/DS1_arma_cube.bin");
}

void double_slit_broader_sigma_y_short_time(void){
    double h = 0.005;
    double Dt = 2.5e-5;
    double T = 0.002;
    double xc = 0.25;
    double sigma_x = 0.05;
    double p_x = 200;
    double yc = 0.5;
    double sigma_y = 0.20;  //  Broader sigma
    double p_y = 0;
    double v0 = 1e10;   // High potential to set up the slits
    // double v0 = 1e10;
    Box b3 = Box(h);
    b3.set_up_walls(v0);
    Simulation s3 = Simulation(b3, Dt, T, xc, sigma_x, p_x, yc, sigma_y, p_y);
    s3.initialise();

    std::cout<<"\nExperiment: DOUBLE-SLIT with broader sigma_y and shorter time\n"<<std::endl;
    arma::cx_cube U3 = s3.run_simulation();
    U3.save("../output/binfiles/DS2_arma_cube.bin");
}




// Må "rydde opp", er noe som ikke helt stemmer:


void no_slits(){
    Box b1 = Box();
    Simulation s1 = Simulation(b1, 2.5e-5, 0.008, 0.25, 0.05, 200.0, 0.5, 0.05, 0.0);
    s1.initialise();

    std::cout<<"\n1st experiment: NO SLITS\n"<<std::endl;
    arma::cx_cube U1 = s1.run_simulation();
    U1.save("../output/binfiles/NS_arma_cube.bin");
}


void double_slit_first(){
    Box b2 = Box();
    b2.create_slits(2);
    Simulation s2 = Simulation(b2, 2.5e-5, 0.008, 0.25, 0.05, 200.0, 0.5, 0.10, 0.0);
    s2.initialise();

    std::cout<<"\n2nd experiment: DOUBLE-SLIT (1)\n"<<std::endl;
    arma::cx_cube U2 = s2.run_simulation();
    U2.save("../output/binfiles/DS1_arma_cube.bin");
}

void double_slit_second(){
    Box b3 = Box();
    b3.create_slits(2);
    Simulation s3 = Simulation(b3, 2.5e-5, 0.002, 0.25, 0.05, 200.0, 0.5, 0.20, 0.0);
    s3.initialise();

    std::cout<<"\n3rd experiment: DOUBLE-SLIT (2)\n"<<std::endl;
    arma::cx_cube U3 = s3.run_simulation();
    U3.save("../output/binfiles/DS2_arma_cube.bin");
}



void single_slit(){

    Box b4 = Box();
    b4.create_slits(1);
    Simulation s4 = Simulation(b4, 2.5e-5, 0.004, 0.25, 0.05, 200.0, 0.5, 0.20, 0.0);
    s4.initialise();

    std::cout<<"\n4th experiment: SINGLE-SLIT \n"<<std::endl;
    arma::cx_cube U4 = s4.run_simulation();
    U4.save("../output/binfiles/SS_arma_cube.bin");

}

void triple_slit(){

    Box b5 = Box();
    b5.create_slits(3);
    Simulation s5 = Simulation(b5, 2.5e-5, 0.004, 0.25, 0.05, 200.0, 0.5, 0.20, 0.0);
    s5.initialise();

    std::cout<<"\n5th experiment: TRIPLE-SLIT \n"<<std::endl;
    arma::cx_cube U5 = s5.run_simulation();
    U5.save("../output/binfiles/TS_arma_cube.bin");

}




int main(){
    // no_slit();
    // double_slit_broad_sigma_y();
    // double_slit_broader_sigma_y_short_time();

    no_slits();
    double_slit_first();
    double_slit_second();
    single_slit();
    triple_slit();
}

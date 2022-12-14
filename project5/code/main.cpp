//  Include stuff

#include "utils.hpp"
#include "Box.hpp"
#include "Simulation.hpp"


/**
 * @brief Main program for FYS4150 project 5 code.
 * @author Nanna Bryne
 * @author Johan M. Kroken
 * @author Vetle A. Vikenes
 */



//  Functions for running specific simulations

void no_slits(){

    Box b = Box();
    Simulation S = Simulation(b);
    S.extend_wavepacket(0.05);
    S.initialise();
    
    std::cout << "\n1st experiment: NO SLITS\n" << std::endl;
    arma::cx_cube U = S.run_simulation();
    U.save("../output/binfiles/NS_arma_cube.bin");
}

void double_slit_first(){
    Box b = Box();
    b.create_slits(2);
    Simulation S = Simulation(b);
    S.extend_wavepacket(0.10);
    S.initialise();

    std::cout << "\n2nd experiment: DOUBLE-SLIT (1)\n" << std::endl;
    arma::cx_cube U = S.run_simulation();
    U.save("../output/binfiles/DS1_arma_cube.bin");
}

void double_slit_second(){
    Box b = Box();
    b.create_slits(2);
    Simulation S = Simulation(b, 0.002);
    S.extend_wavepacket(0.20);
    S.initialise();

    std::cout << "\n2nd experiment: DOUBLE-SLIT (2)\n" << std::endl;
    arma::cx_cube U = S.run_simulation();
    U.save("../output/binfiles/DS2_arma_cube.bin");
}


void single_slit(){
    Box b = Box();
    b.create_slits(1);
    Simulation S = Simulation(b, 0.004);
    S.extend_wavepacket(0.20);
    S.initialise();

    std::cout << "\n2nd experiment: SINGLE-SLIT\n" << std::endl;
    arma::cx_cube U = S.run_simulation();
    U.save("../output/binfiles/SS_arma_cube.bin");
}


void triple_slit(){
    Box b = Box();
    b.create_slits(3);
    Simulation S = Simulation(b, 0.004);
    S.extend_wavepacket(0.20);
    S.initialise();

    std::cout << "\n2nd experiment: TRIPLE-SLIT\n" << std::endl;
    arma::cx_cube U = S.run_simulation();
    U.save("../output/binfiles/TS_arma_cube.bin");
}






/**
 * main that runs all five simulations unless command line argument is given.
 * @param argv what simulation to run by code, can be {NS, DS1, DS2, SS, TS, ALL}
*/
int main(int argc, char **argv){


    // if argument is not given, run all:
    std::string arg1;
    if(argc==1){arg1 = "ALL";} 
    else{arg1.assign(argv[1]); }

    // run commanded simulation(s):
    if(arg1=="NS"){no_slits();}
    else if(arg1=="DS1"){double_slit_first();}
    else if(arg1=="DS2"){double_slit_second();}
    else if(arg1=="SS"){single_slit();}
    else if(arg1=="TS"){triple_slit();}
    else if(arg1=="ALL"){ 
        no_slits();
        double_slit_first(); 
        double_slit_second();    
        single_slit();   
        triple_slit();   
    }
    else{
        std::cout << "Please provide a valid command line argument." << std::endl;
        std::cout << "It can be one of the following: {NS, DS1, DS2, SS, TS, ALL}" << std::endl; 
    }

}

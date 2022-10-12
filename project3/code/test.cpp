// g++ test.cpp -o test.o -larmadillo && ./test.o

#include <iostream>
#include<armadillo>
#include<vector>

using namespace std;

class Particle {
    public:
    int charge_;
    double mass_;
    arma::vec position_;
    arma::vec velocity_;

    Particle(double charge_in, double mass_in, arma::vec position_in, arma::vec velocity_in);

    void new_position(arma::vec new_position);
    void print();
};


Particle::Particle(double charge_in, double mass_in, arma::vec position_in, arma::vec velocity_in){
    charge_ = charge_in;
    mass_ = mass_in;
    position_ = position_in;
    velocity_ = velocity_in;
}

void Particle::new_position(arma::vec new_position){
    position_ = new_position;
}

void Particle::print() {
    cout << "r " << endl <<position_ << endl << "v " << endl  << velocity_ << endl;
}

int main()  {
    arma::vec r0(3, arma::fill::ones);
    arma::vec v0(3, arma::fill::zeros);
    arma::vec r_new(3, arma::fill::randu);


    cout << "Metod 1" <<endl << "Før:" << endl;
    Particle p1(1, 0, r0, v0);
    p1.print();
    vector<Particle> p1_list;
    p1_list.push_back(p1); // UPSI! Her legger vi til en kopi av p1, ikke faktisk p1.
    cout << "p1: " << &p1 << " != p1_list.at(i): " << &p1_list.at(0) << endl;
    Particle& p1_ref = p1_list.at(0);
    p1_ref.new_position(r_new); 
    cout << "Etter" <<endl;
    p1_list.at(0).print(); // Endret p1_list.at(0) men IKKE p1

    cout << endl << "Metode 2" << endl << "Før:" <<endl;
    Particle p2(1, 0, r0, v0);

    p2.print();
    vector<Particle*> p2_list;
    p2_list.push_back(&p2); // Her legger vi til en peker som peker på p2
    cout << "p2: " << &p2 << " == p2_list.at(i): " << p2_list.at(0) << endl;
    p2_list.at(0)->new_position(r_new); // p2.new_position(r_new); vil gjøre det samme
    p2.new_position(r_new);
    cout << "Etter" << endl;
    p2.print(); // p2_list.at(0)->print(); vil gjøre det samme 

    return 0;
}
#include <iostream>
#include "tensor_type.h"
#include <math.h>  //fabs(double x) =|x|
#include <algorithm>
#include <stdlib.h>  //for div(q,n).rem(quot),rand
#include <time.h>
#include <fstream>
#include <limits>
#include <iomanip>
#include <stdio.h>

using namespace std;


class Hartree_Fock_Engine{

public:

    //Declarations must go here----------------------------------//

    int Length_x, Length_y, Length_z, Size, Hamil_Size; //No of sites in x,y,z dir
    int N_orb;    //No of orbitals
    int N_e_PerSite, N_e_Total;    //No of  e's per site and total
    Mat_2_doub Hamiltonian;
    Mat_2_doub Hopping_connections;


    double U, U_p, J_H;   //U,U_p=Onsite coloumbic int, J_H=Onsite Hunds coupling(See Notes)
    Mat_1_doub onsite_potential; //Crystal_field_splitting,disorder etc goes here

    //----------------------------------------------------------------//


    //declaration of Engine parts go below--------------------------------//
    void read_INPUT();
    void Initialize_parameters();
    void Initialize_Tensors();

    void Add_Kinetic_E_part();
    void Add_Spin_orbit_part();
    void Add_Interaction_part();
    //-------------------------------------------------------------------//

};

void Hartree_Fock_Engine::read_INPUT(){

    //later must be read and done whatever by reading input funtion---//
    Length_x = 10;
    Length_y = 10;
    Length_z = 1;
    N_orb = 3;
    N_e_PerSite = 4;
    onsite_potential.resize(2*N_orb*Length_x*Length_y*Length_z);
    for(int i=0;i<2*N_orb*Length_x*Length_y*Length_z;i++){
        onsite_potential[i]=0;
    }

    //---------------------------------------------------------------//
}

void Hartree_Fock_Engine::Initialize_parameters(){

    Size=Length_x*Length_y*Length_z;
    N_e_Total= N_e_PerSite*Size;
    Hamil_Size=2*N_orb*Size;

}

void Hartree_Fock_Engine::Initialize_Tensors(){

    Hamiltonian.resize(Hamil_Size);
    for(int i=0;i<Hamil_Size;i++){
        Hamiltonian[i].resize(Hamil_Size);
    }

    Hopping_connections.resize(Hamil_Size);
    for(int i=0;i<Hamil_Size;i++){
        Hopping_connections[i].resize(Hamil_Size);
    }

}


void Hartree_Fock_Engine::Add_Kinetic_E_part(){

    /* basis is in following order, where notation is c_{site}{orbital}{spin}
   [c_00up c_01up c_02up......c_00dn c01dn...... ]
    */


}

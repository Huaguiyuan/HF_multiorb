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
    /*At present moment this Engine only works for
 3 orbital model and square geometry,
 should make general in future.
*/
public:

    //Declarations must go here----------------------------------//

    int Length_x, Length_y, Length_z, Size, Hamil_Size, Hamil_Sizeby2; //No of sites in x,y,z dir
    int N_orb;    //No of orbitals
    int N_e_PerSite, N_e_Total;    //No of  e's per site and total
    Mat_2_doub Hamiltonian;
    Mat_2_doub Hopping_connections;


    double U, U_p, J_H;   //U,U_p=Onsite coloumbic int, J_H=Onsite Hunds coupling(See Notes)
    Mat_1_doub onsite_potential; //Crystal_field_splitting,disorder etc goes here

    double Convergence_error_n_up, Convergence_error_n_dn; //Convergence_errors for order parameters
    double diff_n_up, diff_n_dn;
    double Alpha_n_up, Alpha_n_dn;


    bool Simple_Mixing_bool, Broyden_bool;

    Mat_1_doub n_up_in, n_dn_in, s_plus_in, s_minus_in;
    Mat_5_doub CdagC;
    string n_up_in_file_str, n_dn_in_file_str; //files from where input order params are being read
    Mat_1_doub n_up_out, n_dn_out, s_plus_out, s_minus_out;

    //Classical energies
    double H1_classical, H2_classical, H3_classical, H4_classical;

    //----------------------------------------------------------------//


    //declaration of Engine parts go below--------------------------------//
    void read_INPUT();
    void Initialize_parameters();
    void Initialize_Tensors();

    void Construct_Hamiltonian();
    void Diagonalize_Hamiltonian();
    void Calculate_order_parameters_out();
    void Add_Kinetic_E_part();
    void Add_Spin_orbit_part();
    void Add_Interaction_part();
    void Add_H1_part();
    void Add_H2_part();
    void Add_H3_part();
    void Add_H4_part();
    void Classical_energy();
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
        onsite_potential[i ]= 0;
    }

    Simple_Mixing_bool = true;
    Broyden_bool = false;
    Convergence_error_n_up=10e-4;
    Convergence_error_n_dn=10e-4;

    //read n_up_in_file_str, n_dn_in_file_str
    // and then n_up_in, n_dn_in of size N_orb*Length_x*Length_y*Length_z//
    //---------------------------------------------------------------//
}

void Hartree_Fock_Engine::Initialize_parameters(){

    Size=Length_x*Length_y*Length_z;
    N_e_Total= N_e_PerSite*Size;
    Hamil_Size=2*N_orb*Size;
    Hamil_Sizeby2=N_orb*Size;

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

void Hartree_Fock_Engine::Construct_Hamiltonian(){


    for(int i=0;i<Hamil_Size;i++){
        for(int j=0;j<Hamil_Size;j++){
            Hamiltonian[i][j]=0;
        }
    }

    Add_Kinetic_E_part();
    Add_Interaction_part();

}

void Hartree_Fock_Engine::Diagonalize_Hamiltonian(){

}

void Hartree_Fock_Engine::Calculate_order_parameters_out(){

}

void Hartree_Fock_Engine::Add_Kinetic_E_part(){

    /* basis is in following order, where notation is c_{site_x}{site_y}{orbital}{spin}
    [c_{0,0,0,up} c_{1,0,0,up}....c_{Lx-1,0,0,up}...c_{x_i,y_i,orb_i,spin_i}...c_xy0dn cxy1dn...... ]
    in future we can add z dir also.
    */

    /*
    spin=up <=> spin_i=0, spin=dn <=> spin=1
    orb_i \belongs {0,1,2}
    x_i \belongs {0,1..,Length_x-1}
    y_i \belongs {0,1..,Length_y-1}

    Total elements in above sequence above are = 2*N_orb*Length_x*Length_y*Length_z i.e Hamil_Size
    jth term in above sequence is related to spin_i,orb_i,x_i,y_i as below:
    j=x_i + Length_x*y_i + Length_x*Length_y*orb_i +  Length_x*Length_y*N_orb*spin_i
    */

    /*
      Hopping_connections[j][jp] will give connection b/w j and jp( i.e j')
     */


    for(int j=0;j<Hamil_Size;j++){
        for(int jp=0;jp<Hamil_Size;jp++){
            Hamiltonian[j][jp] = Hamiltonian[j][jp] + Hopping_connections[j][jp];
        }
    }


}


void Hartree_Fock_Engine::Add_H1_part(){

    //Hartree term
    for(int j=0;j<Hamil_Size;j++){
        if(j<Hamil_Sizeby2){
            Hamiltonian[j][j] = Hamiltonian[j][j] + U*n_up_in[j];
        }
        else{
            Hamiltonian[j][j] = Hamiltonian[j][j] + U*n_up_in[j-Hamil_Sizeby2];
        }

    }


    //Fock term
    for(int j=0;j<Hamil_Sizeby2;j++){
        Hamiltonian[j][j+Hamil_Sizeby2] = Hamiltonian[j][j+Hamil_Sizeby2] - U*s_minus_in[j];
        Hamiltonian[j+Hamil_Sizeby2][j] = Hamiltonian[j+Hamil_Sizeby2][j] - U*s_plus_in[j];

    }


}


void Hartree_Fock_Engine::Add_H2_part(){

    //H_2 part 1(see notes)

    int spin_alpha, spin_beta;
    for(int x_i=0;x_i<Length_x;x_i++){
        for(int y_i=0;y_i<Length_y;y_i++){
            for(int orb_beta=0;orb_beta<N_orb;orb_beta++){
                for(int orb_alpha=0;orb_alpha<orb_beta;orb_alpha++){

                  //n_{i,alpha,up}<n_{i,beta,up}>
                  spin_alpha=0;//up
                  spin_beta=0;//up
                  j=x_i + Length_x*y_i + Length_x*Length_y*orb_alpha +
                          Length_x*Length_y*N_orb*spin_alpha;
                  jp=x_i + Length_x*y_i + Length_x*Length_y*orb_beta +
                          Length_x*Length_y*N_orb*spin_beta;

                  Hamiltonian[j][j] = Hamiltonian[j][j] +
                                      (U_p - J_H*0.5)*n_up_in[jp];

                  //n_{i,alpha,up}<n_{i,beta,dn}>
                  spin_alpha=0;//up
                  spin_beta=1;//dn
                  j=x_i + Length_x*y_i + Length_x*Length_y*orb_alpha +
                          Length_x*Length_y*N_orb*spin_alpha;
                  jp=x_i + Length_x*y_i + Length_x*Length_y*orb_beta +
                          Length_x*Length_y*N_orb*spin_beta;

                  Hamiltonian[j][j] = Hamiltonian[j][j] +
                                      (U_p - J_H*0.5)*n_dn_in[jp-Hamil_Sizeby2];

                  //n_{i,alpha,dn}<n_{i,beta,up}>
                  spin_alpha=1;//dn
                  spin_beta=0;//up
                  j=x_i + Length_x*y_i + Length_x*Length_y*orb_alpha +
                          Length_x*Length_y*N_orb*spin_alpha;
                  jp=x_i + Length_x*y_i + Length_x*Length_y*orb_beta +
                          Length_x*Length_y*N_orb*spin_beta;

                  Hamiltonian[j][j] = Hamiltonian[j][j] +
                                      (U_p - J_H*0.5)*n_up_in[jp];

                  //n_{i,alpha,dn}<n_{i,beta,dn}>
                  spin_alpha=1;//dn
                  spin_beta=1;//dn
                  j=x_i + Length_x*y_i + Length_x*Length_y*orb_alpha +
                          Length_x*Length_y*N_orb*spin_alpha;
                  jp=x_i + Length_x*y_i + Length_x*Length_y*orb_beta +
                          Length_x*Length_y*N_orb*spin_beta;

                  Hamiltonian[j][j] = Hamiltonian[j][j] +
                                      (U_p - J_H*0.5)*n_dn_in[jp-Hamil_Sizeby2];


                  //n_{i,beta,up}<n_{i,alpha,up}>
                  spin_alpha=0;//up
                  spin_beta=0;//up
                  j=x_i + Length_x*y_i + Length_x*Length_y*orb_alpha +
                          Length_x*Length_y*N_orb*spin_alpha;
                  jp=x_i + Length_x*y_i + Length_x*Length_y*orb_beta +
                          Length_x*Length_y*N_orb*spin_beta;

                  Hamiltonian[jp][jp] = Hamiltonian[jp][jp] +
                                      (U_p - J_H*0.5)*n_up_in[j];


                  //n_{i,beta,up}<n_{i,alpha,dn}>
                  spin_alpha=1;//dn
                  spin_beta=0;//up
                  j=x_i + Length_x*y_i + Length_x*Length_y*orb_alpha +
                          Length_x*Length_y*N_orb*spin_alpha;
                  jp=x_i + Length_x*y_i + Length_x*Length_y*orb_beta +
                          Length_x*Length_y*N_orb*spin_beta;

                  Hamiltonian[jp][jp] = Hamiltonian[jp][jp] +
                                      (U_p - J_H*0.5)*n_dn_in[j-Hamil_Sizeby2];

                  //n_{i,beta,dn}<n_{i,alpha,up}>
                  spin_alpha=0;//up
                  spin_beta=1;//dn
                  j=x_i + Length_x*y_i + Length_x*Length_y*orb_alpha +
                          Length_x*Length_y*N_orb*spin_alpha;
                  jp=x_i + Length_x*y_i + Length_x*Length_y*orb_beta +
                          Length_x*Length_y*N_orb*spin_beta;

                  Hamiltonian[jp][jp] = Hamiltonian[jp][jp] +
                                      (U_p - J_H*0.5)*n_up_in[j];

                  //n_{i,beta,dn}<n_{i,alpha,dn}>
                  spin_alpha=1;//dn
                  spin_beta=1;//dn
                  j=x_i + Length_x*y_i + Length_x*Length_y*orb_alpha +
                          Length_x*Length_y*N_orb*spin_alpha;
                  jp=x_i + Length_x*y_i + Length_x*Length_y*orb_beta +
                          Length_x*Length_y*N_orb*spin_beta;

                  Hamiltonian[jp][jp] = Hamiltonian[jp][jp] +
                                      (U_p - J_H*0.5)*n_dn_in[j-Hamil_Sizeby2];


                }

            }


        }


    }


    //H2 part2


}

void Hartree_Fock_Engine::Classical_energy(){
    H1_classical=0;

    //Hartree term classical energy
    for(int j=0;j<Hamil_Sizeby2;j++){
        H1_classical = H1_classical - U*n_up_out[j]*n_dn_out[j];
    }


}

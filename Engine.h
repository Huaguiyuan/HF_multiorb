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

    //Models-----------------------------------------------
    string _MODEL;
    //-----------------------------------------------------


    //Declarations must go here----------------------------------//

    int Length_x, Length_y, Length_z, Size, Hamil_Size, Hamil_Sizeby2; //No of sites in x,y,z dir
    int N_orb;    //No of orbitals
    int N_e_PerSite, N_e_Total;    //No of  e's per site and total
    Mat_2_doub Hamiltonian;
    Mat_2_doub Eigenvectors;
    Mat_1_doub Evals;
    Mat_2_doub Hopping_connections;
    int _SPIN_UP;
    int _SPIN_DN;


    double U, U_p, J_H;   //U,U_p=Onsite coloumbic int, J_H=Onsite Hunds coupling(See Notes)
    Mat_1_doub onsite_potential; //Crystal_field_splitting,disorder etc goes here
    Mat_2_doub T_matrix;    //Nearest Neighbour hoppin only.

    double Convergence_error_n_up, Convergence_error_n_dn; //Convergence_errors for order parameters
    double diff_n_up, diff_n_dn;
    double Alpha_n_up, Alpha_n_dn;


    bool Simple_Mixing_bool, Broyden_bool;


    Mat_6_doub CdagC_in;//[xi][yi][orb1][orb2][spin1][spin2]
    Mat_6_doub CdagC_out;//[xi][yi][orb1][orb2][spin1][spin2]
    string n_up_in_file_str, n_dn_in_file_str; //files from where input order params are being read
    Mat_1_doub n_up_out, n_dn_out, s_plus_out, s_minus_out;

    //Classical energies
    double H1_classical, H2_classical, H3_classical, H4_classical;


    //----------------------------------------------------------------//


    //declaration of Engine parts go below--------------------------------//
    void read_INPUT();
    void Initialize_parameters();
    void Initialize_Tensors();
    void Create_Hopping_connections();
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
    Length_x = 20;
    Length_y = 20;
    Length_z = 1;
    N_orb = 3;
    N_e_PerSite = 4;
    onsite_potential.resize(2*N_orb*Length_x*Length_y*Length_z);
    for(int i=0;i<2*N_orb*Length_x*Length_y*Length_z;i++){
        onsite_potential[i]= 0;
    }

    Simple_Mixing_bool = true;
    Broyden_bool = false;
    Convergence_error_n_up=10e-4;
    Convergence_error_n_dn=10e-4;
    _MODEL="Pnictides_3orb";

    //if Pnictides_3orb read T_matrix

    //T_matrix values
    T_matrix.resize(3);
    for(int i=0;i<3;i++){
        T_matrix[i].resize(3);
    }
    T_matrix[0][0]=0;T_matrix[0][1]=0;T_matrix[0][2]=0;
    T_matrix[1][0]=0;T_matrix[1][1]=0;T_matrix[1][2]=0;
    T_matrix[2][0]=0;T_matrix[2][1]=0;T_matrix[2][2]=-1;



    //read n_up_in_file_str, n_dn_in_file_str
    // and then n_up_in, n_dn_in of size N_orb*Length_x*Length_y*Length_z//
    //---------------------------------------------------------------//
}

void Hartree_Fock_Engine::Initialize_parameters(){

    _SPIN_UP=0;
    _SPIN_DN=1;
    Size=Length_x*Length_y*Length_z;
    N_e_Total= N_e_PerSite*Size;
    Hamil_Size=2*N_orb*Size;
    Hamil_Sizeby2=N_orb*Size;

}

void Hartree_Fock_Engine::Initialize_Tensors(){

    Hamiltonian.resize(Hamil_Size);
    Eigenvectors.resize(Hamil_Size);
    Evals.resize(Hamil_Size);

    for(int i=0;i<Hamil_Size;i++){
        Hamiltonian[i].resize(Hamil_Size);
        Eigenvectors[i].resize(Hamil_Size);
    }

    Hopping_connections.resize(Hamil_Size);
    for(int i=0;i<Hamil_Size;i++){
        Hopping_connections[i].resize(Hamil_Size);
    }

}

void Hartree_Fock_Engine::Create_Hopping_connections(){

    //Hopping_connections
    //Right now working only for 3 orbital Hubbard model (Nearest neighbour hopping)
    //j=x_i + Length_x*y_i + Length_x*Length_y*orb_i +  Length_x*Length_y*N_orb*spin_i

    for(int i =0;i<Hamil_Size;i++){
        for(int j =0;j<Hamil_Size;j++){
            Hopping_connections[i][j]=0;
        }
    }


    int N_neigh_xy=4;

    double neigh_xi[4]={1,0,-1,0};
    double neigh_yi[4]={0,1,0,-1};

    int xi_neigh, yi_neigh;

    int N_neigh_orb=3;
    double orb_neigh[3]={0,1,2};


    int j_neigh,j;


    for(int x_i=0;x_i<Length_x;x_i++){
        for(int y_i=0;y_i<Length_y;y_i++){
            for(int orb_i=0;orb_i<N_orb;orb_i++){
                for (int spin_i=0;spin_i<2;spin_i++){
                    j=x_i + Length_x*y_i +
                            Length_x*Length_y*orb_i +
                            Length_x*Length_y*N_orb*spin_i;

                    for(int neigh_xy_i=0;neigh_xy_i<N_neigh_xy;neigh_xy_i++){
                        xi_neigh=x_i+neigh_xi[neigh_xy_i];
                        yi_neigh=y_i+neigh_yi[neigh_xy_i];

                        if(xi_neigh==-1){
                            xi_neigh=Length_x-1;
                        }
                        if(xi_neigh==Length_x){
                            xi_neigh=0;
                        }

                        if(yi_neigh==-1){
                            yi_neigh=Length_y-1;
                        }
                        if(yi_neigh==Length_y){
                            yi_neigh=0;
                        }



                        for(int neigh_orb_i=0;neigh_orb_i<N_neigh_orb;neigh_orb_i++){
                            j_neigh=xi_neigh+ Length_x*yi_neigh +
                                    Length_x*Length_y*orb_neigh[neigh_orb_i] +
                                    Length_x*Length_y*N_orb*spin_i;
                            Hopping_connections[j][j_neigh]=T_matrix[orb_i][orb_neigh[neigh_orb_i]];
                        }
                    }
                }
            }

        }


    }


}

void Hartree_Fock_Engine::Construct_Hamiltonian(){


    for(int i=0;i<Hamil_Size;i++){
        for(int j=0;j<Hamil_Size;j++){
            Hamiltonian[i][j]=0;
        }
    }

    Add_Kinetic_E_part();
    //Add_Interaction_part();

}

void Hartree_Fock_Engine::Diagonalize_Hamiltonian(){


    int LDA=Hamil_Size, info;
    /* Local arrays */


    double* eval = new double[Hamil_Size];
    double* mat = new double[Hamil_Size*Hamil_Size];

    for(int i=0;i<Hamil_Size;i++){
        for(int j=0;j<=i;j++){

            mat[i*Hamil_Size+j]=Hamiltonian[i][j];

        }
    }

    info=LAPACKE_dsyev(LAPACK_ROW_MAJOR,  'V', 'L',Hamil_Size , mat , LDA, eval );


    /* Check for convergence */
    if( info > 0 ) {

        cout<< "The LAPACKE_dsyev failed to diagonalize."<<endl;

    }

    //Eigenvectors
    for(int i=0;i<Hamil_Size;i++){
        for(int j=0;j<Hamil_Size;j++){

            Eigenvectors[i][j]=mat[j*Hamil_Size+i];

        }
    }

    //Evals
    for(int j=0;j<Hamil_Size;j++){
        Evals[j]=eval[j];
    }



    free(eval);
    free(mat);

}

void Hartree_Fock_Engine::Calculate_order_parameters_out(){

}

void Hartree_Fock_Engine::Add_Kinetic_E_part(){

    /* basis is in following order, where notation is c_{site_x}{site_y}{orbital}{spin}
    [c_{0,0,0,up} c_{1,0,0,up}....c_{Lx-1,0,0,up}...c_{x_i,y_i,orb_i,spin_i}...c_xy0dn cxy1dn...... ]
    in future we can add z dir also.
    */

    /*
    spin=up <=> spin_i=_SPIN_UP, spin=dn <=> spin=_SPIN_DN
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

    int j;
    //Hartree term
    for(int x_i=0;x_i<Length_x;x_i++){
        for(int y_i=0;y_i<Length_y;y_i++){
            for(int orb=0;orb<N_orb;orb++){
                for(int spin=0;spin<2;spin++){
                    j=x_i + Length_x*y_i + Length_x*Length_y*orb +
                            Length_x*Length_y*N_orb*spin;
                    Hamiltonian[j][j] = Hamiltonian[j][j] + U*CdagC_in[x_i][y_i][orb][orb][spin][spin];
                }
            }
        }
    }



    //Fock term
    for(int x_i=0;x_i<Length_x;x_i++){
        for(int y_i=0;y_i<Length_y;y_i++){
            for(int orb=0;orb<N_orb;orb++){
                j=x_i + Length_x*y_i + Length_x*Length_y*orb;

                Hamiltonian[j][j+Hamil_Sizeby2] =
                        Hamiltonian[j][j+Hamil_Sizeby2] +
                        U*CdagC_in[x_i][y_i][orb][orb][_SPIN_DN][_SPIN_UP];
                Hamiltonian[j+Hamil_Sizeby2][j] =
                        Hamiltonian[j+Hamil_Sizeby2][j] +
                        U*CdagC_in[x_i][y_i][orb][orb][_SPIN_UP][_SPIN_DN];



            }
        }
    }


}


void Hartree_Fock_Engine::Add_H2_part(){

    int j,jp;
    double _factor;
    _factor = (U_p - J_H*0.5);
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
                            _factor*
                            CdagC_in[x_i][y_i][orb_beta][orb_beta][_SPIN_UP][_SPIN_UP];


                    //n_{i,alpha,up}<n_{i,beta,dn}>
                    spin_alpha=0;//up
                    spin_beta=1;//dn
                    j=x_i + Length_x*y_i + Length_x*Length_y*orb_alpha +
                            Length_x*Length_y*N_orb*spin_alpha;
                    jp=x_i + Length_x*y_i + Length_x*Length_y*orb_beta +
                            Length_x*Length_y*N_orb*spin_beta;

                    Hamiltonian[j][j] = Hamiltonian[j][j] +
                            _factor*
                            CdagC_in[x_i][y_i][orb_beta][orb_beta][_SPIN_DN][_SPIN_DN];

                    //n_{i,alpha,dn}<n_{i,beta,up}>
                    spin_alpha=1;//dn
                    spin_beta=0;//up
                    j=x_i + Length_x*y_i + Length_x*Length_y*orb_alpha +
                            Length_x*Length_y*N_orb*spin_alpha;
                    jp=x_i + Length_x*y_i + Length_x*Length_y*orb_beta +
                            Length_x*Length_y*N_orb*spin_beta;

                    Hamiltonian[j][j] = Hamiltonian[j][j] +
                            _factor*
                            CdagC_in[x_i][y_i][orb_beta][orb_beta][_SPIN_UP][_SPIN_UP];

                    //n_{i,alpha,dn}<n_{i,beta,dn}>
                    spin_alpha=1;//dn
                    spin_beta=1;//dn
                    j=x_i + Length_x*y_i + Length_x*Length_y*orb_alpha +
                            Length_x*Length_y*N_orb*spin_alpha;
                    jp=x_i + Length_x*y_i + Length_x*Length_y*orb_beta +
                            Length_x*Length_y*N_orb*spin_beta;

                    Hamiltonian[j][j] = Hamiltonian[j][j] +
                            _factor*
                            CdagC_in[x_i][y_i][orb_beta][orb_beta][_SPIN_DN][_SPIN_DN];


                    //n_{i,beta,up}<n_{i,alpha,up}>
                    spin_alpha=0;//up
                    spin_beta=0;//up
                    j=x_i + Length_x*y_i + Length_x*Length_y*orb_alpha +
                            Length_x*Length_y*N_orb*spin_alpha;
                    jp=x_i + Length_x*y_i + Length_x*Length_y*orb_beta +
                            Length_x*Length_y*N_orb*spin_beta;

                    Hamiltonian[jp][jp] = Hamiltonian[jp][jp] +
                            _factor*
                            CdagC_in[x_i][y_i][orb_alpha][orb_alpha][_SPIN_UP][_SPIN_UP];


                    //n_{i,beta,up}<n_{i,alpha,dn}>
                    spin_alpha=1;//dn
                    spin_beta=0;//up
                    j=x_i + Length_x*y_i + Length_x*Length_y*orb_alpha +
                            Length_x*Length_y*N_orb*spin_alpha;
                    jp=x_i + Length_x*y_i + Length_x*Length_y*orb_beta +
                            Length_x*Length_y*N_orb*spin_beta;

                    Hamiltonian[jp][jp] = Hamiltonian[jp][jp] +
                            _factor*
                            CdagC_in[x_i][y_i][orb_alpha][orb_alpha][_SPIN_DN][_SPIN_DN];

                    //n_{i,beta,dn}<n_{i,alpha,up}>
                    spin_alpha=0;//up
                    spin_beta=1;//dn
                    j=x_i + Length_x*y_i + Length_x*Length_y*orb_alpha +
                            Length_x*Length_y*N_orb*spin_alpha;
                    jp=x_i + Length_x*y_i + Length_x*Length_y*orb_beta +
                            Length_x*Length_y*N_orb*spin_beta;

                    Hamiltonian[jp][jp] = Hamiltonian[jp][jp] +
                            _factor*
                            CdagC_in[x_i][y_i][orb_alpha][orb_alpha][_SPIN_UP][_SPIN_UP];

                    //n_{i,beta,dn}<n_{i,alpha,dn}>
                    spin_alpha=1;//dn
                    spin_beta=1;//dn
                    j=x_i + Length_x*y_i + Length_x*Length_y*orb_alpha +
                            Length_x*Length_y*N_orb*spin_alpha;
                    jp=x_i + Length_x*y_i + Length_x*Length_y*orb_beta +
                            Length_x*Length_y*N_orb*spin_beta;

                    Hamiltonian[jp][jp] = Hamiltonian[jp][jp] +
                            _factor*
                            CdagC_in[x_i][y_i][orb_alpha][orb_alpha][_SPIN_DN][_SPIN_DN];


                }

            }


        }


    }


    _factor= - (U_p - J_H*0.5);
    //H2 part2


    //Remember orb1_OP,spin1_OP is for Cdag_{orb1_OP,spin1_OP}  inside <CdadC>
    //Remember orb2_OP,spin2_OP is for C_{orb2_OP,spin2_OP}  inside <CdadC>
    //Remember orb1,spin1 is for Cdag_{orb1,spin1}  inside CdadC
    //Remember orb2,spin2 is for C_{orb2,spin2}  inside CdadC


    for(int x_i=0;x_i<Length_x;x_i++){
        for(int y_i=0;y_i<Length_y;y_i++){
            for(int orb_beta=0;orb_beta<N_orb;orb_beta++){
                for(int orb_alpha=0;orb_alpha<orb_beta;orb_alpha++){


                    int orb1_OP[8]={orb_alpha,orb_beta,orb_alpha,orb_beta,orb_alpha,orb_beta,orb_alpha,orb_beta};
                    int orb2_OP[8]={orb_beta,orb_alpha,orb_beta,orb_alpha,orb_beta,orb_alpha,orb_beta,orb_alpha};

                    int spin1_OP[8]={_SPIN_UP,_SPIN_UP,_SPIN_UP,_SPIN_DN,_SPIN_DN,_SPIN_UP,_SPIN_DN,_SPIN_DN};
                    int spin2_OP[8]={_SPIN_UP,_SPIN_UP,_SPIN_DN,_SPIN_UP,_SPIN_UP,_SPIN_DN,_SPIN_DN,_SPIN_DN};



                    for(int index=0;index<8;index++){

                        j=x_i + Length_x*y_i + Length_x*Length_y*orb2_OP[index] +
                                Length_x*Length_y*N_orb*spin2_OP[index];
                        jp=x_i + Length_x*y_i + Length_x*Length_y*orb1_OP[index] +
                                Length_x*Length_y*N_orb*spin1_OP[index];

                        Hamiltonian[j][jp] = Hamiltonian[j][jp] +
                                _factor*
                                CdagC_in[x_i][y_i][orb1_OP[index]][orb2_OP[index]][spin1_OP[index]][spin2_OP[index]];


                    }

                }
            }
        }
    }




}


void Hartree_Fock_Engine::Add_H3_part(){

    int j,jp;
    double _factor;


    //Remember orb1_OP,spin1_OP is for Cdag_{orb1_OP,spin1_OP}  inside <CdadC>
    //Remember orb2_OP,spin2_OP is for C_{orb2_OP,spin2_OP}  inside <CdadC>
    //Remember orb1,spin1 is for Cdag_{orb1,spin1}  inside CdadC
    //Remember orb2,spin2 is for C_{orb2,spin2}  inside CdadC


    int orb2[12];



    //H3 part1
    _factor = -J_H/2;
    for(int x_i=0;x_i<Length_x;x_i++){
        for(int y_i=0;y_i<Length_y;y_i++){
            for(int orb_beta=0;orb_beta<N_orb;orb_beta++){
                for(int orb_alpha=0;orb_alpha<orb_beta;orb_alpha++){


                    int orb1_OP[12]={orb_alpha,orb_beta,orb_alpha,orb_beta,orb_alpha,orb_beta,
                             orb_alpha,orb_beta,orb_alpha,orb_beta,orb_alpha,orb_beta};
                    int orb1[12]={orb_beta,orb_alpha,orb_beta,orb_alpha,orb_beta,orb_alpha,
                          orb_beta,orb_alpha,orb_beta,orb_alpha,orb_beta,orb_alpha};


                    int spin1_OP[12]={_SPIN_UP,_SPIN_DN,_SPIN_DN,_SPIN_UP,_SPIN_UP,_SPIN_UP,
                              _SPIN_DN,_SPIN_UP,_SPIN_UP,_SPIN_DN,_SPIN_DN,_SPIN_DN};
                    int spin2_OP[12]={_SPIN_DN,_SPIN_UP,_SPIN_UP,_SPIN_DN,_SPIN_UP,_SPIN_UP,
                              _SPIN_DN,_SPIN_UP,_SPIN_UP,_SPIN_DN,_SPIN_DN,_SPIN_DN};

                    int spin1[12]={_SPIN_DN,_SPIN_UP,_SPIN_UP,_SPIN_DN,_SPIN_UP,_SPIN_UP,
                           _SPIN_UP,_SPIN_DN,_SPIN_DN,_SPIN_UP,_SPIN_DN,_SPIN_DN};
                    int spin2[12]={_SPIN_UP,_SPIN_DN,_SPIN_DN,_SPIN_UP,_SPIN_UP,_SPIN_UP,
                           _SPIN_UP,_SPIN_DN,_SPIN_DN,_SPIN_UP,_SPIN_DN,_SPIN_DN};




                    for(int index=0;index<12;index++){

                        if(index<4){
                            _factor = 2*_factor;
                        }

                        if(index>5 && index<10){
                            _factor = -1.0*_factor;
                        }

                        j=x_i + Length_x*y_i + Length_x*Length_y*orb1[index] +
                                Length_x*Length_y*N_orb*spin1[index];
                        jp=x_i + Length_x*y_i + Length_x*Length_y*orb1[index] +
                                Length_x*Length_y*N_orb*spin2[index];

                        Hamiltonian[j][jp] = Hamiltonian[j][jp] +
                                _factor*
                                CdagC_in[x_i][y_i][orb1_OP[index]][orb1_OP[index]][spin1_OP[index]][spin2_OP[index]];


                    }

                }
            }
        }
    }


    //H3 part2
    _factor = J_H/2;
    for(int x_i=0;x_i<Length_x;x_i++){
        for(int y_i=0;y_i<Length_y;y_i++){
            for(int orb_beta=0;orb_beta<N_orb;orb_beta++){
                for(int orb_alpha=0;orb_alpha<orb_beta;orb_alpha++){


                    int orb1_OP[12]={orb_alpha,orb_beta,orb_alpha,orb_beta,orb_alpha,orb_beta,
                             orb_alpha,orb_beta,orb_alpha,orb_beta,orb_alpha,orb_beta};
                    int orb2_OP[12]={orb_beta,orb_alpha,orb_beta,orb_alpha,orb_beta,orb_alpha,
                             orb_beta,orb_alpha,orb_beta,orb_alpha,orb_beta,orb_alpha};


                    int spin1_OP[12]={_SPIN_UP,_SPIN_DN,_SPIN_DN,_SPIN_UP,_SPIN_UP,_SPIN_UP,
                              _SPIN_DN,_SPIN_UP,_SPIN_UP,_SPIN_DN,_SPIN_DN,_SPIN_DN};
                    int spin2_OP[12]={_SPIN_UP,_SPIN_DN,_SPIN_DN,_SPIN_UP,_SPIN_UP,_SPIN_UP,
                              _SPIN_UP,_SPIN_DN,_SPIN_DN,_SPIN_UP,_SPIN_DN,_SPIN_DN};


                    for(int index=0;index<12;index++){

                        if(index<4){
                            _factor = 2*_factor;
                        }

                        if(index>5 && index<10){
                            _factor = -1.0*_factor;
                        }

                        j=x_i + Length_x*y_i + Length_x*Length_y*orb2_OP[index] +
                                Length_x*Length_y*N_orb*spin2_OP[index];
                        jp=x_i + Length_x*y_i + Length_x*Length_y*orb1_OP[index] +
                                Length_x*Length_y*N_orb*spin1_OP[index];

                        Hamiltonian[j][jp] = Hamiltonian[j][jp] +
                                _factor*
                                CdagC_in[x_i][y_i][orb1_OP[index]][orb2_OP[index]][spin1_OP[index]][spin2_OP[index]];


                    }

                }
            }
        }
    }



}


void Hartree_Fock_Engine::Add_H4_part(){

    int j,jp;
    double _factor;

    //Remember orb1_OP,spin1_OP is for Cdag_{orb1_OP,spin1_OP}  inside <CdadC>
    //Remember orb2_OP,spin2_OP is for C_{orb2_OP,spin2_OP}  inside <CdadC>
    //Remember orb1,spin1 is for Cdag_{orb1,spin1}  inside CdadC
    //Remember orb2,spin2 is for C_{orb2,spin2}  inside CdadC







    int orb1[4];
    int orb2[4];
    int spin1[4];
    int spin2[4];

    //H4 part1
    _factor = J_H;
    for(int x_i=0;x_i<Length_x;x_i++){
        for(int y_i=0;y_i<Length_y;y_i++){
            for(int orb_beta=0;orb_beta<N_orb;orb_beta++){
                for(int orb_alpha=0;orb_alpha<orb_beta;orb_alpha++){


                    int orb1_OP[4]={orb_alpha,orb_alpha,orb_beta,orb_beta};
                    int orb2_OP[4]={orb_beta,orb_beta,orb_alpha,orb_alpha};


                    int spin1_OP[4]={_SPIN_UP,_SPIN_DN,_SPIN_UP,_SPIN_DN};
                    int spin2_OP[4]={_SPIN_DN,_SPIN_UP,_SPIN_DN,_SPIN_UP};



                    for(int index=0;index<4;index++){

                        _factor = -1.0*_factor;


                        j=x_i + Length_x*y_i + Length_x*Length_y*orb1_OP[index] +
                                Length_x*Length_y*N_orb*spin2_OP[index];
                        jp=x_i + Length_x*y_i + Length_x*Length_y*orb2_OP[index] +
                                Length_x*Length_y*N_orb*spin1_OP[index];

                        Hamiltonian[j][jp] = Hamiltonian[j][jp] +
                                _factor*
                                CdagC_in[x_i][y_i][orb1_OP[index]][orb2_OP[index]][spin1_OP[index]][spin2_OP[index]];


                    }

                }
            }
        }
    }


    //H4 part1
    _factor = J_H;
    for(int x_i=0;x_i<Length_x;x_i++){
        for(int y_i=0;y_i<Length_y;y_i++){
            for(int orb_beta=0;orb_beta<N_orb;orb_beta++){
                for(int orb_alpha=0;orb_alpha<orb_beta;orb_alpha++){


                    int orb1_OP[4]={orb_alpha,orb_alpha,orb_beta,orb_beta};
                    int orb2_OP[4]={orb_beta,orb_beta,orb_alpha,orb_alpha};


                    int spin1_OP[4]={_SPIN_UP,_SPIN_DN,_SPIN_DN,_SPIN_UP};
                    int spin1[4]={_SPIN_DN,_SPIN_UP,_SPIN_UP,_SPIN_DN};



                    for(int index=0;index<4;index++){


                        j=x_i + Length_x*y_i + Length_x*Length_y*orb1_OP[index] +
                                Length_x*Length_y*N_orb*spin1[index];
                        jp=x_i + Length_x*y_i + Length_x*Length_y*orb2_OP[index] +
                                Length_x*Length_y*N_orb*spin1[index];

                        Hamiltonian[j][jp] = Hamiltonian[j][jp] +
                                _factor*
                                CdagC_in[x_i][y_i][orb1_OP[index]][orb2_OP[index]][spin1_OP[index]][spin1_OP[index]];


                    }

                }
            }
        }
    }





}

void Hartree_Fock_Engine::Classical_energy(){
    H1_classical=0;

    //Hartree term classical energy
    for(int j=0;j<Hamil_Sizeby2;j++){
        H1_classical = H1_classical - U*n_up_out[j]*n_dn_out[j];
    }


}

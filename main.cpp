#include <iostream>  //for cin and cout
#include <math.h>  // for pow
#include <stdlib.h>  //for div(q,n).rem(quot),abs(int n)
#include <time.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <complex>
#include <assert.h>
#include "Engine.h"
#include "iconprint.h"

using namespace std; 


int main(int argc, char** argv )
{
    iconprint();

    cout<<"Using Hamiltonian parameters from input file = "<<argv[1]<<endl;
    cout<<"Using Initial guess for order parameters from input file = "<<argv[2]<<endl;



    Hartree_Fock_Engine HF_Engine;


    HF_Engine.FILE_HAMIL_PARAM = argv[1];
    HF_Engine.FILE_ORDER_PARAM = argv[2];
    HF_Engine.FILE_ORDER_PARAM_OUT = argv[3];



    HF_Engine.read_INPUT();

    HF_Engine.Initialize_parameters();



    HF_Engine.Initialize_Tensors();
    HF_Engine.Create_Hopping_connections();




    bool Convergence_gained;
    Convergence_gained=false;
    //Start doing Hartree Fock calculation--------------------------------------//
    int iteration=0;
    HF_Engine.diff_energy=1.0;
    while(Convergence_gained==false)
    {
        HF_Engine.Construct_Hamiltonian();

        HF_Engine.Diagonalize_Hamiltonian();

        HF_Engine.Calculate_order_parameters_out();

        HF_Engine.E_old=HF_Engine.E_new;

        if(HF_Engine._CONVERGENCE_GLOBAL==false){
    //checking if Engine should continue to run or stop----------------------------//
        if( (HF_Engine.diff_n_up >= HF_Engine.Convergence_error_n_up)
                || (HF_Engine.diff_n_up >= HF_Engine.Convergence_error_n_dn) )
        {
            Convergence_gained=false;
            HF_Engine.Guess_new_input();
        }
        else
        {
            Convergence_gained=true;
        }

     //---------------------------------------------------------------------------//
        }
        else{

            if(
               (HF_Engine.diff_OP_global < HF_Engine.Convergence_error_global)
                    &&
               (HF_Engine.diff_energy < HF_Engine.Convergence_error_energy)
              )
            {
                Convergence_gained=true;

            }
            else
            {
                Convergence_gained=false;
                HF_Engine.Guess_new_input();
            }
        }

    cout.precision(20);
    cout<<iteration<<"    "<<HF_Engine.diff_OP_global<<"    ";
    cout<<HF_Engine.Total_energy<<"    "<<HF_Engine.diff_energy<<"     "<<HF_Engine.Total_N_exp<<endl;
    iteration++;
    }

    //Hartree Fock convergence achieved-----------------------------------------//
    HF_Engine.Write_Order_Params();

    //Observables calulation starts----------------------------------------------//


    //---------------------------------------------------------------------------//

    return 0;
}

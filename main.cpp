#include <iostream>  //for cin and cout
#include <math.h>  // for pow
#include <stdlib.h>  //for div(q,n).rem(quot),abs(int n)
#include <time.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <complex>
#include "Engine.h"

using namespace std; 


int main(int argc, char** argv)
{


    Hartree_Fock_Engine HF_Engine;

    HF_Engine.read_INPUT();
    HF_Engine.Initialize_parameters();
    HF_Engine.Initialize_Tensors();

    bool Convergence_gained;
    Convergence_gained=false;
    //Start doing Hartree Fock calculation--------------------------------------//

    while(Convergence_gained==false)
    {
        HF_Engine.Construct_Hamiltonian();

        HF_Engine.Diagonalize_Hamiltonian();

        HF_Engine.Calculate_order_parameters_out();





    //checking if Engine should continue to run or stop----------------------------//
        if( (HF_Engine.diff_n_up > HF_Engine.Convergence_error_n_up)
                || (HF_Engine.diff_n_up > HF_Engine.Convergence_error_n_dn) )
        {
            Convergence_gained=false;
        }
        else
        {
            Convergence_gained=true;
        }

     //---------------------------------------------------------------------------//

    }


    //Hartree Fock calculation done----------------------------------------------//


    return 0;
}

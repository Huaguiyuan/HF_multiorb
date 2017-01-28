#include <iostream>  //for cin and cout
#include <math.h>  // for pow
#include <stdlib.h>  //for div(q,n).rem(quot),abs(int n)
#include <time.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <complex>
using namespace std; 

typedef vector< complex<double> >  Mat_1_Complex_doub;
typedef vector<Mat_1_Complex_doub> Mat_2_Complex_doub;

int main(int argc, char** argv){

//Declarations must go here----------------------------------

    int Length_x, Length_y, Length_z, Size, Hamil_Size; //No of sites in x,y,z dir
    int N_orb;    //No of orbitals
    int N_e_PerSite, N_e_Total;     //No of  e's per site and total
    Mat_2_Complex_doub Hamiltonian;

    double U, U_p, J ;   //U,U_p=Onsite coloumbic int, J=Hunds coupling(See Notes)

//----------------------------------------------------------------


    //later from input---
     dim = 2;
     Length_x = 10;
     Length_y = 10;
     Length_z = 1;
     N_orb = 3;
     N_e_PerSite = 4;
    //later from input---

    Size=Length_x*Length_y*Length_x;
    N_e_Total= N_e_PerSite*Size;
    Hamil_Size=2*N_orb*Size;

    Hamiltonian.resize(Hamil_size);
    for(int i=0;i<Hamil_size;i++){
    Hamiltonian[i].resize(Hamil_Size);
    }


    Construct_KE_part(Hamiltonian);



    return 0;
}

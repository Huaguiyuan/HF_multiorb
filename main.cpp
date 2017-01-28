#include <iostream>  //for cin and cout
#include <math.h>  // for pow
#include <stdlib.h>  //for div(q,n).rem(quot),abs(int n)
#include <time.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <complex>

using namespace std; 


int main(int argc, char** argv){






    Hamiltonian.resize(Hamil_size);
    for(int i=0;i<Hamil_size;i++){
    Hamiltonian[i].resize(Hamil_Size);
    }


    Construct_KE_part(Hamiltonian);



    return 0;
}

#include <iostream>
#include <fstream>
#include <string> 
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include "tensor_type.h"
#include <assert.h>
using namespace std;



void read_order_parameters(Mat_6_doub & CdagC_in, string FILE_ORDER_PARAM,
                           int Length_x,int Length_y, int N_orb,
                           bool RANDOM_INITAL_GUESS,int N_e_Total)
{



    ifstream inputfile( FILE_ORDER_PARAM.c_str());

    /*******************************resizing***********************/
    CdagC_in.resize(Length_x);
    for(int xi=0;xi<Length_x;xi++){
        CdagC_in[xi].resize(Length_y);
         for(int yi=0;yi<Length_y;yi++){
             CdagC_in[xi][yi].resize(N_orb);
             for(int o1=0;o1<N_orb;o1++){
                 CdagC_in[xi][yi][o1].resize(N_orb);
                 for(int o2=0;o2<N_orb;o2++){
                     CdagC_in[xi][yi][o1][o2].resize(2);
                     for(int s1=0;s1<2;s1++){
                         CdagC_in[xi][yi][o1][o2][s1].resize(2);
                     }
                 }
             }
         }
    }
    /****************************************************************************/

    if(RANDOM_INITAL_GUESS==false){
    int x,y, o1, o2, s1, s2;
    int line_no=0;
    double val;
    for(int x_i=0;x_i<Length_x;x_i++){
        for(int y_i=0;y_i<Length_y;y_i++){
            for(int orb1_i=0;orb1_i<N_orb;orb1_i++){
                for(int orb2_i=0;orb2_i<N_orb;orb2_i++){
                    for (int spin1_i=0;spin1_i<2;spin1_i++){
                        for (int spin2_i=0;spin2_i<2;spin2_i++){

                            inputfile>>x>>y>>o1>>o2>>s1>>s2>>val;
                            CdagC_in[x][y][o1][o2][s1][s2]=val;
                            line_no++;
                        }
                    }
                }
            }
        }
    }


    if(line_no!=Length_x*Length_y*N_orb*N_orb*4){
        cout<<"no. of lines in "<<FILE_ORDER_PARAM<<" must be "<<
              Length_x*Length_y*N_orb*N_orb*4<<endl;
        assert (false);
    }
    cout<<"Order parameters were read"<<endl;
       }

    else{
        double den_temp=0;
        double val;
        srand (time(NULL));
        for(int xi=0;xi<Length_x;xi++){
            for(int yi=0;yi<Length_y;yi++){
                for(int o1=0;o1<N_orb;o1++){
                    for(int o2=0;o2<N_orb;o2++){
                        for(int s1=0;s1<2;s1++){
                            for(int s2=0;s2<2;s2++){
                                val=rand();
                                val=val/RAND_MAX;
                                CdagC_in[xi][yi][o1][o2][s1][s2]=val;
                                if(o1==02 && s1==s2){
                                    den_temp=den_temp+val;
                                }

                            }
                        }
                    }
                }
            }
        }

        for(int xi=0;xi<Length_x;xi++){
            for(int yi=0;yi<Length_y;yi++){
                for(int o1=0;o1<N_orb;o1++){
                    for(int o2=0;o2<N_orb;o2++){
                        for(int s1=0;s1<2;s1++){
                            for(int s2=0;s2<2;s2++){

                                if(o1==02 && s1==s2){

                                    CdagC_in[xi][yi][o1][o2][s1][s2]=
                                       N_e_Total*(CdagC_in[xi][yi][o1][o2][s1][s2]/den_temp);
                                }

                            }
                        }
                    }
                }
            }
        }


       cout<<"Random guess for order parameters"<<endl;
    }

}

#include "Continuum.h"

void Runge_Kutta_Df_fixed_Cont(vec2x& bgs, double& EF1, double& EF2, vec1C& ArrayCont,int ContMode)
{
    complexd pulse1, pulse2;
    for(int i=0;i<ArrayCont.size();i++){      
        int mmax = ArrayCont[i].Mmax; // As I call it several times, it will be faster to have it in a local variable
        int lm; // As I call it several times, it will be faster to have it in a local variable
        if(ContMode==0){
            for (int eps=0; eps<ArrayCont[i].NE; eps++)
            {
                for (int L=0; L<=ArrayCont[i].Lmax; L++)
                {
                    int M1=( L<=mmax ? L : mmax);
                    for (int M=-M1; M<=M1; M++)
                    {
                        lm = spH(L,M,mmax);
                        fill(ArrayCont[i].Vte[eps][lm].begin(), ArrayCont[i].Vte[eps][lm].end(), 0.); // Reset summatory to 0
                        for(auto j:ArrayCont[i].UniqueStates){ // Iterating throw different couplings
                            if(ArrayCont[i].Allow[0][j] == 1){
                                pulse1 = (ArrayCont[i].DIPpump)[ArrayCont[i].Indexes[0][j]][eps][lm][0]*EF1;
                            }
                            else pulse1 = (0,0);

                            if(ArrayCont[i].Allow[1][j] == 1){
                                pulse2 = (ArrayCont[i].DIPprobe)[ArrayCont[i].Indexes[1][j]][eps][lm][0]*EF2;
                            }
                            else pulse2 = (0,0);

                            ArrayCont[i].Vte[eps][lm][0]+=c1*(pulse1 + pulse2)*bgs[j][0]; // Now it has only one state (1x0 matrix)
                        }
                        ArrayCont[i].bev[eps][lm][0] = (-c1*(ArrayCont[i].PES[0]+ArrayCont[i].E[eps])- 0.5*ArrayCont[i].Gamma)*ArrayCont[i].be[eps][lm][0] - ArrayCont[i].Vte[eps][lm][0];
                    }
                }
            }
        }
        else{
            for (int eps=0; eps<ArrayCont[i].NE; eps++)
            {
                for (int L=0; L<=ArrayCont[i].Lmax; L++)
                {
                    int M1=( L<=mmax ? L : mmax);
                    for (int M=-M1; M<=M1; M++)
                    {
                        lm = spH(L,M,mmax);
                        fill(ArrayCont[i].Vte[eps][lm].begin(), ArrayCont[i].Vte[eps][lm].end(), 0.); // Reset summatory to 0
                        for(auto j:ArrayCont[i].UniqueStates){ // Iterating throw different couplings
                            if(ArrayCont[i].Allow[0][j] == 1){
                                pulse1 = (ArrayCont[i].DIPpump)[ArrayCont[i].Indexes[0][j]][eps][lm][0]*EF1;
                            }
                            else pulse1 = (0,0);

                            if(ArrayCont[i].Allow[1][j] == 1){
                                pulse2 = (ArrayCont[i].DIPprobe)[ArrayCont[i].Indexes[1][j]][eps][lm][0]*EF2;
                            }
                            else pulse2 = (0,0);

                            ArrayCont[i].Vte[eps][lm][0]+=c1*(pulse1 + pulse2)*bgs[j][0]; // Now it has only one state (1x0 matrix)
                        }
                        ArrayCont[i].bev[eps][lm][0] = (-c1*(ArrayCont[i].PES[0]+ArrayCont[i].E[eps])- 0.5*ArrayCont[i].Gamma)*ArrayCont[i].be2[eps][lm][0] - ArrayCont[i].Vte[eps][lm][0];
                    }
                }
            }
        }
    }
}

void Runge_Kutta_Ac_Cont(double& dt, vec1C& ArrayCont)
{
    int NR=1;
    if (NR>1) {    
        for(int i=0;i<ArrayCont.size();i++){          
            int mmax = ArrayCont[i].Mmax;
            int lm;
            for (int iR=0; iR<NR; iR++)
            {
                for (int eps=0; eps<ArrayCont[i].NE; eps++)
                {
                    for (int L=0; L<=ArrayCont[i].Lmax; L++)
                    {
                        int M1=( L<=mmax ? L : mmax);
                        for (int M=-M1; M<=M1; M++)
                        {
                            lm = spH(L,M,mmax);
                            ArrayCont[i].be1[eps][lm][iR] += ArrayCont[i].bev[eps][lm][iR]*dt;
                        }
                    }
                }
            }
        }
    }
    else{    
        for(int i=0;i<ArrayCont.size();i++){          
            int mmax = ArrayCont[i].Mmax;
            int lm;
            
            for (int eps=0; eps<ArrayCont[i].NE; eps++)
            {
                for (int L=0; L<=ArrayCont[i].Lmax; L++)
                {
                    int M1=( L<=mmax ? L : mmax);
                    for (int M=-M1; M<=M1; M++)
                    {
                        lm = spH(L,M,mmax);
                        ArrayCont[i].be1[eps][lm][0] += ArrayCont[i].bev[eps][lm][0]*dt;
                    }
                }
            }
        }
    }
}

void Runge_Kutta_Ad_Cont(double& dt, vec1C& ArrayCont, int ContMode)
{
    int NR=1;
    
    if (NR>1) {
        for(int i=0;i<ArrayCont.size();i++){          
            int mmax = ArrayCont[i].Mmax;
            int lm;
            for (int iR=0; iR<NR; iR++) // CHANGE CONTMODE BEFORE THAT LOOP FOR INCREASING VELOCITY
            {
                if(ContMode==1){ // If here much faster than inside the other loops
                    for (int eps=0; eps<ArrayCont[i].NE; eps++){
                        for (int L=0; L<=ArrayCont[i].Lmax; L++)
                        {
                            int M1=( L<=mmax ? L : mmax);
                            for (int M=-M1; M<=M1; M++)
                            {
                                lm = spH(L,M,mmax);
                                ArrayCont[i].be1[eps][lm][iR] = ArrayCont[i].be[eps][lm][iR] + ArrayCont[i].bev[eps][lm][iR]*dt;
                            }
                        }
                    }
                }
                else if(ContMode==2){
                    for (int eps=0; eps<ArrayCont[i].NE; eps++){
                        for (int L=0; L<=ArrayCont[i].Lmax; L++)
                        {
                            int M1=( L<=mmax ? L : mmax);
                            for (int M=-M1; M<=M1; M++)
                            {
                                lm = spH(L,M,mmax);
                                ArrayCont[i].be2[eps][lm][iR] = ArrayCont[i].be[eps][lm][iR] + ArrayCont[i].bev[eps][lm][iR]*dt; 
                            }
                        }
                    }
                }
                else{
                    for (int eps=0; eps<ArrayCont[i].NE; eps++){
                        for (int L=0; L<=ArrayCont[i].Lmax; L++)
                        {
                            int M1=( L<=mmax ? L : mmax);
                            for (int M=-M1; M<=M1; M++)
                            {
                                lm = spH(L,M,mmax);
                                ArrayCont[i].be[eps][lm][iR] = ArrayCont[i].be1[eps][lm][iR] + ArrayCont[i].bev[eps][lm][iR]*dt;
                            }
                        }
                    }
                }
            }
        }
    }
    else{    
        for(int i=0;i<ArrayCont.size();i++){          
            int mmax = ArrayCont[i].Mmax;
            int lm;
            if(ContMode==1){ // If here much faster than inside the other loops
                for (int eps=0; eps<ArrayCont[i].NE; eps++){
                    for (int L=0; L<=ArrayCont[i].Lmax; L++)
                    {
                        int M1=( L<=mmax ? L : mmax);
                        for (int M=-M1; M<=M1; M++)
                        {
                            lm = spH(L,M,mmax);
                            ArrayCont[i].be1[eps][lm][0] = ArrayCont[i].be[eps][lm][0] + ArrayCont[i].bev[eps][lm][0]*dt;
                        }
                    }
                }
            }
            else if(ContMode==2){
                for (int eps=0; eps<ArrayCont[i].NE; eps++){
                    for (int L=0; L<=ArrayCont[i].Lmax; L++)
                    {
                        int M1=( L<=mmax ? L : mmax);
                        for (int M=-M1; M<=M1; M++)
                        {
                            lm = spH(L,M,mmax);
                            ArrayCont[i].be2[eps][lm][0] = ArrayCont[i].be[eps][lm][0] + ArrayCont[i].bev[eps][lm][0]*dt;
                        }
                    }
                }
            }
            else{
                for (int eps=0; eps<ArrayCont[i].NE; eps++){
                    for (int L=0; L<=ArrayCont[i].Lmax; L++)
                    {
                        int M1=( L<=mmax ? L : mmax);
                        for (int M=-M1; M<=M1; M++)
                        {
                            lm = spH(L,M,mmax);
                            ArrayCont[i].be[eps][lm][0] = ArrayCont[i].be1[eps][lm][0] + ArrayCont[i].bev[eps][lm][0]*dt;
                        }
                    }
                }
            }
        }
    }
}
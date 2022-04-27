/*
 *  RK.h
 *  
 *  Created by APA on 13/02/2019
 *
 *
 */
#include "Continuum.h"

class Kinetic
{
public:
    Kinetic(string name,double dR);
    double _dR;
    double Energy;
};
Kinetic::Kinetic(string name,double dR): _dR(dR){
    if(name=="CO"){
        double MassC=12.00000000;
        double MassO=15.99491464;
        //double MassN=14.0030740052;
        double massAU=MassC*MassO*1822.88839/(MassC+MassO); //unit "u" or Dalton and a.u.
        Energy=1./(2.*massAU*_dR*_dR);
    }
    else if(name=="H"){
        double MassH=1.0078250322;
        double massAU=MassH*1822.88839;
        Energy=1./(2.*massAU*_dR*_dR);
    }
    else if(name=="N2O"){
        double MassN=14.0030740052;
        double MassO=15.99491464;
        double massAU=MassN*MassN*MassO*1822.88839/(2.*MassN+MassO);
        Energy=1./(2.*massAU*_dR*_dR);
    }
}
/*
double kintf(int iR) //NOTE: we need to include dR in the input file and transfer that variable to this function
{
    double MassC=12.00000000;
    double MassO=15.99491464;
    //double MassN=14.0030740052;
    double massAU=MassC*MassO*1822.88839/(MassC+MassO); //unit "u" or Dalton and a.u.
    double dR = 0.01;
    return 1./(2.*massAU*dR*dR);
}
 */

void Runge_Kutta_Df_fixed(vec1x& Vt, vec2x& bgs, vec2x& bgsv, vec2d& PES, vec3d& Dip1,vec3d& Dip2, vec1d& Gamma, double& EF1, double& EF2, vec1C& ArrayCont,int ContMode)

{
    int NEi = PES.size();
    
//#pragma omp parallel for schedule(dynamic) default(shared) -- it causes a problem due to Vt array, we should use Vt[Ei][iR]
    for (int Ei=0; Ei<NEi; Ei++)
    {
        fill(Vt.begin(), Vt.end(), 0.);
        
        for (int Ej=0; Ej<NEi; Ej++)
        {
            Vt[0]+=c1*(Dip1[Ei][Ej][0]*EF1+Dip2[Ei][Ej][0]*EF2)*bgs[Ej][0];
        }
        bgsv[Ei][0] = (-c1*PES[Ei][0]- 0.5*Gamma[Ei])*bgs[Ei][0] - Vt[0];
        // cout << "Ei: " << Ei << " " << PES[Ei][0] << endl;
        // cout << Gamma[Ei] << endl;
    }

    if(ArrayCont.size()>=0){
        complexd pulse1, pulse2;
        for(int i=0;i<ArrayCont.size();i++){      
            // for(auto a:ArrayCont[i].Allow[0]) cout << a << " ";
            // cout << endl;
            // for(auto a:ArrayCont[i].Allow[1]) cout << a << " ";
            // cout << endl;
            // for(auto a:ArrayCont[i].Indexes[0]) cout << a << " ";
            // cout << endl;
            // for(auto a:ArrayCont[i].Indexes[1]) cout << a << " ";
            // cout << endl;
            // exit(1);  
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
                            //for(int Ej=0; Ej<=ArrayCont[i].maxBoundState; Ej++){
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
                                // cout << "j: " << j << " Correct index: " << ArrayCont[i].Indexes[0][j] <<  " allowance: " << ArrayCont[i].Allow[0][j] << endl;
                                // cout << "j: " << j << " Correct index: " << ArrayCont[i].Indexes[1][j] <<  " allowance: " << ArrayCont[i].Allow[1][j] << endl;
                            }
                            // cout << ArrayCont[i].Vte[eps][lm][0] << endl;

                            ArrayCont[i].bev[eps][lm][0] = (-c1*(ArrayCont[i].PES[0]+ArrayCont[i].E[eps])- 0.5*ArrayCont[i].Gamma)*ArrayCont[i].be[eps][lm][0] - ArrayCont[i].Vte[eps][lm][0];
                            // cout << ArrayCont[i].Gamma << endl;
                            // cout << "0: " << ArrayCont[i].bev[eps][lm][0] << endl;
                            // cout << "PES: " << ArrayCont[i].PES[0] << endl;
                            // cout << "EPS: "<< ArrayCont[i].E[eps] << endl;
                            // cout << "PES+EPS: "<< ArrayCont[i].PES[0]+ArrayCont[i].E[eps] << endl;
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
                            //for(int Ej=0; Ej<=ArrayCont[i].maxBoundState; Ej++){
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
                            // cout << "1: " << ArrayCont[i].bev[eps][lm][0] << endl;
                        }
                    }
                }
            }
        }
    }
}

void Runge_Kutta_Df(Kinetic NuclearKE,vec1x& Vt, vec2x& bgs, vec2x& bgsv, vec2d& PES, vec3d& Dip1,vec3d& Dip2, vec1d& Gamma, double& EF1, double& EF2, vec1C& ArrayCont, int ContMode)
{
    
    int NEi = PES.size();
    int NR = PES[0].size();
    
//#pragma omp parallel for schedule(dynamic) default(shared) -- it causes a problem due to Vt array, we should use Vt[Ei][iR]
    for (int Ei=0; Ei<NEi; Ei++)
    {
        fill(Vt.begin(), Vt.end(), 0.);
        
        for (int Ej=0; Ej<NEi; Ej++)
        {
            #pragma omp simd
            for (int iR=1; iR<NR-1; iR++)
            {
                Vt[iR]+=c1*(Dip1[Ei][Ej][iR]*EF1+Dip2[Ei][Ej][iR]*EF2)*bgs[Ej][iR];
            }
        }
        #pragma omp simd
        for (int iR=1; iR<NR-1; iR++)
        {
            bgsv[Ei][iR] = -c1*(-bgs[Ei][iR+1]+2.*bgs[Ei][iR]-bgs[Ei][iR-1])*NuclearKE.Energy +(-c1*PES[Ei][iR]- 0.5*Gamma[Ei])*bgs[Ei][iR] - Vt[iR];
        }
    }
    if(ArrayCont.size()>=0){
        complexd pulse1, pulse2;
        for(int i=0;i<ArrayCont.size();i++){          
            int mmax = ArrayCont[i].Mmax; // As I call it several times, it will be faster to have it in a local variable
            int lm; // As I call it several times, it will be faster to have it in a local variable
            bool allow; // Counting couplings
            if(ContMode==0){
                for (int iR=1; iR<NR-1; iR++)
                {
                    for (int eps=0; eps<ArrayCont[i].NE; eps++)
                    {
                        for (int L=0; L<=ArrayCont[i].Lmax; L++)
                        {
                            int M1=( L<=mmax ? L : mmax);
                            for (int M=-M1; M<=M1; M++)
                            {
                                lm = spH(L,M,mmax);
                                fill(ArrayCont[i].Vte[eps][lm].begin(), ArrayCont[i].Vte[eps][lm].end(), 0.); // Reset summatory to 0
                                //for(int Ej=0; Ej<=ArrayCont[i].maxBoundState; Ej++){
                                for(auto j:ArrayCont[i].UniqueStates){ // Iterating throw different couplings
                                    allow=0;
                                    if(ArrayCont[i].Allow[0][j] == 1){
                                        pulse1 = (ArrayCont[i].DIPpump)[ArrayCont[i].Indexes[0][j]][eps][lm][iR]*EF1;
                                        allow = 1;
                                    }
                                    else pulse1 = (0,0);

                                    if(ArrayCont[i].Allow[1][j] == 1){
                                        pulse2 = (ArrayCont[i].DIPprobe)[ArrayCont[i].Indexes[1][j]][eps][lm][iR]*EF2;
                                    }
                                    else pulse2 = (0,0);

                                    ArrayCont[i].Vte[eps][lm][iR]+=c1*(pulse1 + pulse2)*bgs[j][iR]; // Now it has only one state (1x0 matrix)
                                }
                                ArrayCont[i].bev[eps][lm][iR] = (-c1*(ArrayCont[i].PES[iR]+ArrayCont[i].E[eps])- 0.5*ArrayCont[i].Gamma)*ArrayCont[i].be[eps][lm][iR] - ArrayCont[i].Vte[eps][lm][0];
                            }
                        }
                    }
                }
            }
            else{
                for (int iR=1; iR<NR-1; iR++)
                {
                    for (int eps=0; eps<ArrayCont[i].NE; eps++)
                    {
                        for (int L=0; L<=ArrayCont[i].Lmax; L++)
                        {
                            int M1=( L<=mmax ? L : mmax);
                            for (int M=-M1; M<=M1; M++)
                            {
                                lm = spH(L,M,mmax);
                                fill(ArrayCont[i].Vte[eps][lm].begin(), ArrayCont[i].Vte[eps][lm].end(), 0.); // Reset summatory to 0
                                //for(int Ej=0; Ej<=ArrayCont[i].maxBoundState; Ej++){
                                for(auto j:ArrayCont[i].UniqueStates){ // Iterating throw different couplings
                                    allow=0;
                                    if(ArrayCont[i].Allow[0][j] == 1){
                                        pulse1 = (ArrayCont[i].DIPpump)[ArrayCont[i].Indexes[0][j]][eps][lm][iR]*EF1;
                                        allow = 1;
                                    }
                                    else pulse1 = (0,0);

                                    if(ArrayCont[i].Allow[1][j] == 1){
                                        pulse2 = (ArrayCont[i].DIPprobe)[ArrayCont[i].Indexes[1][j]][eps][lm][iR]*EF2;
                                    }
                                    else pulse2 = (0,0);

                                    ArrayCont[i].Vte[eps][lm][iR]+=c1*(pulse1 + pulse2)*bgs[j][iR]; // Now it has only one state (1x0 matrix)
                                }
                                ArrayCont[i].bev[eps][lm][iR] = (-c1*(ArrayCont[i].PES[iR]+ArrayCont[i].E[eps])- 0.5*ArrayCont[i].Gamma)*ArrayCont[i].be2[eps][lm][iR] - ArrayCont[i].Vte[eps][lm][0];
                            }
                        }
                    }
                }
            }
        }
    }
}

void Runge_Kutta_Ac(vec2x& bgs1, vec2x& bgsv, double& dt, vec1C& ArrayCont)
{
    int NEi=bgs1.size();
    int NR=bgs1[0].size();
    
    if (NR>1) {
        for (int Ei=0; Ei<NEi; Ei++)
        {
            #pragma omp simd
            for (int iR=0; iR<NR; iR++)
            {
                bgs1[Ei][iR]+=bgsv[Ei][iR]*dt;
            }
        }
        if(ArrayCont.size()>=0){     
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
    }
    else{
        for (int Ei=0; Ei<NEi; Ei++) bgs1[Ei][0]+=bgsv[Ei][0]*dt;

        if(ArrayCont.size()>=0){     
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
}


void Runge_Kutta_Ad(vec2x& bgs0, vec2x& bgs1, vec2x& bgsv, double& dt, vec1C& ArrayCont, int ContMode)
{
    int NEi=bgs1.size();
    int NR=bgs1[0].size();
    
    if (NR>1) {
        for (int Ei=0; Ei<NEi; Ei++)
        {
            #pragma omp simd
            for (int iR=0; iR<NR; iR++)
            {
                bgs1[Ei][iR]=bgs0[Ei][iR] + bgsv[Ei][iR]*dt;
            }
        }
        if(ArrayCont.size()>=0){     
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
    }
    else{
        for (int Ei=0; Ei<NEi; Ei++){
            bgs1[Ei][0]=bgs0[Ei][0] + bgsv[Ei][0]*dt;
        }
        if(ArrayCont.size()>=0){     
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
                                // cout << ArrayCont[i].be[eps][lm][0] << endl;
                            }
                        }
                    }
                }
            }
        }
    }
}

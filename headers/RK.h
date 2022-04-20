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

//void Runge_Kutta_Df_fixed(vec1x& Vt, vec2x& bgs, vec2x& bgsv, vec2d& PES, vec3d& Dip1,vec3d& Dip2, vec1d& Gamma, double& EF1, double& EF2)
void Runge_Kutta_Df_fixed(vec1x& Vt, vec2x& bgs, vec2x& bgsv, vec2d& PES, vec3d& Dip1,vec3d& Dip2, vec1d& Gamma, double& EF1, double& EF2, vec1C ArrayCont, vec1x& Vti )

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
    }
    // ,vec1C& ArrayCont,vec1x& Vti, vec4x& biv
    // ADD NContStat, ArrayCont and Runge-Kutta vectors to function input
    if(ArrayCont.size()>=0){
        complexd pulse1, pulse2;
        for(int i=0;i<ArrayCont.size();i++){
            // ArrayCont[i].positionEj(NEi); 

            int Ncoupl = (ArrayCont[i]).StatCouplPump.size() + (ArrayCont[i]).StatCouplProbe.size();
            int maxcpl= max((ArrayCont[i]).StatCouplPump.size(), (ArrayCont[i]).StatCouplProbe.size());
            int pumpsize = (ArrayCont[i].StatCouplPump).size();
            int mmax = ArrayCont[i].Mmax; // As I call it several times, it will be faster to have it in a local variable
            // Can a state couple with both lasers?? YES
            fill(Vti.begin(), Vti.end(), 0.);
            for (int eps=0; eps<ArrayCont[i].NE; eps++)
            {
                for (int L=0; L<=ArrayCont[i].Lmax; L++)
                {
                    int M1=( L<=mmax ? L : mmax);
                    for (int M=-M1; M<=M1; M++)
                    {
                        int lm = spH(L,M,mmax); // As I call it several times, it will be faster to have it in a local variable
                        // cout << "here it works" << endl;
                        // cout << "Allow is a matrixi " << ArrayCont[i].Allow.size() << " x " << ArrayCont[i].Allow[0].size() << endl;
                        // for(auto a:ArrayCont[i].Allow[0]) cout << a << " ";
                        // cout << endl;
                        // for(auto a:ArrayCont[i].Allow[1]) cout << a << " ";
                        // cout << endl;
                        // cout << "Max bound state " << ArrayCont[i].maxBoundState << endl;
                        for(int Ej=0; Ej<=ArrayCont[i].maxBoundState; Ej++){
                            if(ArrayCont[i].Allow[0][Ej] == 1){
                                pulse1 = (ArrayCont[i].DIPpump)[ArrayCont[i].Indexes[0][Ej]][eps][lm][0]*EF1;
                                cout << "Ej for pump: "<< ArrayCont[i].Indexes[0][Ej] << endl;
                            }
                            else pulse1 = (0,0);

                            // if(ArrayCont[i].Allow[1][Ej] == 1){
                            //     //pulse2 = (ArrayCont[i].DIPpump)[ArrayCont[i].Indexes[1][Ej]][eps][lm][0]*EF2;
                            //     cout << "Ej for probe: "<< ArrayCont[i].Indexes[1][Ej] << endl;
                            // }
                            // else pulse2 = (0,0);
                            // cout << "size of bgs: " << bgs.size() << endl;
                            // if(ArrayCont[i].Allow[0][Ej] == 1 || ArrayCont[i].Allow[1][Ej] == 1){
                            //     Vti[0]+=c1*(pulse1 + pulse2)*bgs[Ej][0];
                            // }
                        }
                        // biv[i][eps][lm][0] = (-c1*ArrayCont[i].PES[i][0]- 0.5*ArrayCont[i].Gamma[i])*bgs[i][0] - Vti[0];


                        // if(Ej<bar)        (ArrayCont[i].DIPpump)[Epsj][eps][SpH(L,M,mmax)][0]      = x*u1[0]+y*u1[1]+z*u1[2];
                        // else if(bar == 0) (ArrayCont[i].DIPprobe)[Epsj][eps][SpH(L,M,mmax)][0]     = x*u2[0]+y*u2[1]+z*u2[2];
                        // else              (ArrayCont[i].DIPprobe)[bar-Epsj][eps][SpH(L,M,mmax)][0] = x*u2[0]+y*u2[1]+z*u2[2];
                        // PES path with double core-hole energies
                       
                    }
                }
            }
        }
    }
    // cout << "bg " << bgsv[0][0] << " " << bgsv[20][0] << endl;
}

void Runge_Kutta_Df(Kinetic NuclearKE,vec1x& Vt, vec2x& bgs, vec2x& bgsv, vec2d& PES, vec3d& Dip1,vec3d& Dip2, vec1d& Gamma, double& EF1, double& EF2)
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
}

void Runge_Kutta_Ac(vec2x& bgs1, vec2x& bgsv, double& dt)
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
    }
    else for (int Ei=0; Ei<NEi; Ei++) bgs1[Ei][0]+=bgsv[Ei][0]*dt;
}


void Runge_Kutta_Ad(vec2x& bgs0, vec2x& bgs1, vec2x& bgsv, double& dt)
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
    }
    else for (int Ei=0; Ei<NEi; Ei++) bgs1[Ei][0]=bgs0[Ei][0] + bgsv[Ei][0]*dt;
}

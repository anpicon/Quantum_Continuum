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

void Runge_Kutta_Df_fixed(vec1x& Vt, vec2x& bgs, vec2x& bgsv, vec2d& PES, vec3d& Dip1,vec3d& Dip2, vec1d& Gamma, double& EF1, double& EF2)

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
    else{
        for (int Ei=0; Ei<NEi; Ei++) bgs1[Ei][0]+=bgsv[Ei][0]*dt;
    }
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
    else{
        for (int Ei=0; Ei<NEi; Ei++){
            bgs1[Ei][0]=bgs0[Ei][0] + bgsv[Ei][0]*dt;
        }
    }
}

/*
 *  Observables.h
 *  
 *  Created by APA on 13/02/2019
 *
 *
 */



bool check_conservation(ofstream& fp_population,vec2x& bgs, double time)
{
    
    int NEi=bgs.size();
    int NR=bgs[0].size();
    double sum=0.;
    vec1d sumgs(NEi,0.);
    
    for (int Ei=0; Ei<NEi; Ei++)
    {
        for (int iR=0; iR<NR; iR++)
        {
            sumgs[Ei]+=norm(bgs[Ei][iR]);
        }
        sum+=sumgs[Ei];
    }
    
    
    //cout << "Norm " << sum << endl;
    fp_population << time*time_au_fs;
    for (int Ei=0; Ei<NEi; Ei++)
    {
        fp_population << " " << sumgs[Ei] ;
    }
    fp_population << endl;
    
    if (abs(sum - 1.0) < 0.01)
    return true;
    else
    return false;
}

void PrintEF(ofstream& fp_E, double& EF1, double& EF2,double& time)
{
    fp_E << time*time_au_fs << " " << EF1 << " " << EF2 << endl;
} //End PrintEF


void PrintWFD(ofstream& fp_E, vec1d& R, vec2x& b0)
{
    int NEi=b0.size();
    int NR=b0[0].size();
    
    //Printing wavefunctions
    for(int iR=0; iR<NR; iR++){
        fp_E << R[iR] << " ";
        for (int Ei=0; Ei<NEi; Ei++){
           fp_E << norm(b0[Ei][iR]) << " ";
        }
        fp_E << endl;
    }
} //End PrintWFD

void PrintAmpl(ofstream& fp_E, vec2d& PES, vec2x& b0)
{
    int NEi = PES.size();
    int NR = PES[0].size();
    
    for (int Ei=0; Ei<NEi; Ei++){
        for(int iR=0; iR<NR; iR++){
            fp_E << PES[Ei][iR] << " " << real(b0[Ei][iR]) << " " << imag(b0[Ei][iR]) << " ";
        }
        fp_E << endl;
    }
} //End PrintAmpl

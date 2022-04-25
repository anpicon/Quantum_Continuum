/*
 *  Observables.h
 *  
 *  Created by APA on 13/02/2019
 *
 *
 */



bool check_conservation(ofstream& fp_population,vec2x& bgs, double time, vec1C& ArrayCont,ofstream& fp_Contpopulation)
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
        // cout << " aaa " << endl;
        // cout << "size: " << ArrayCont.size() << endl;
    }
    fp_population << endl;

    if(ArrayCont.size()>=0){     
        for(int i=0;i<ArrayCont.size();i++){ 
            double ContPopulation =0.0;     
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
                            ContPopulation += norm(ArrayCont[i].be[eps][lm][iR]);
                        }
                    }
                }
            }
            fp_Contpopulation << time*time_au_fs;
            fp_Contpopulation << " " << ContPopulation;
        }
        fp_Contpopulation << endl;
    }

    if (abs(sum - 1.0) < 0.01) return true;
    else return false;
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

#include "Continuum.h"

void Runge_Kutta_Df_fixed_Cont(vec2x& bgs, double& EF1, double& EF2, vec1C& ArrayCont,int ContMode)
{
    for(int i=0;i<ArrayCont.size();i++){      
        if(ContMode==0){
            #pragma omp parallel for
            for (int eps=0; eps<ArrayCont[i].NE; eps++)
            {
                complexd pulse1, pulse2;
                for(auto lm: ArrayCont[i].lm)
                {
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
        else{
            #pragma omp parallel for
            for (int eps=0; eps<ArrayCont[i].NE; eps++)
            {
                complexd pulse1, pulse2;
                for(auto lm: ArrayCont[i].lm)
                {
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
void Runge_Kutta_Df_Cont(vec2x& bgs, double& EF1, double& EF2, vec1C& ArrayCont, int ContMode)
{
    int NR = ArrayCont[0].Vte[0][0].size();
    
    //#pragma omp parallel for schedule(dynamic) default(shared) -- it causes a problem due to Vt array, we should use Vt[Ei][iR]
    complexd pulse1, pulse2;
    for(int i=0;i<ArrayCont.size();i++){
        if(ContMode==0){
            for (int iR=1; iR<NR-1; iR++)
            {
                for (int eps=0; eps<ArrayCont[i].NE; eps++)
                {
                    for(auto lm: ArrayCont[i].lm)
                    {
                        fill(ArrayCont[i].Vte[eps][lm].begin(), ArrayCont[i].Vte[eps][lm].end(), 0.); // Reset summatory to 0
                        for(auto j:ArrayCont[i].UniqueStates){ // Iterating throw different couplings
                            if(ArrayCont[i].Allow[0][j] == 1){
                                pulse1 = (ArrayCont[i].DIPpump)[ArrayCont[i].Indexes[0][j]][eps][lm][iR]*EF1;
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
        else{
            for (int iR=1; iR<NR-1; iR++)
            {
                for (int eps=0; eps<ArrayCont[i].NE; eps++)
                {
                    for(auto lm: ArrayCont[i].lm)
                    {
                        fill(ArrayCont[i].Vte[eps][lm].begin(), ArrayCont[i].Vte[eps][lm].end(), 0.); // Reset summatory to 0
                        for(auto j:ArrayCont[i].UniqueStates){ // Iterating throw different couplings
                            if(ArrayCont[i].Allow[0][j] == 1){
                                pulse1 = (ArrayCont[i].DIPpump)[ArrayCont[i].Indexes[0][j]][eps][lm][iR]*EF1;
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

void Runge_Kutta_Ac_Cont(double& dt, vec1C& ArrayCont)
{
    int NR = ArrayCont[0].Vte[0][0].size();
    if (NR>1) {    
        for(int i=0;i<ArrayCont.size();i++){
            for (int iR=0; iR<NR; iR++)
            {
                for (int eps=0; eps<ArrayCont[i].NE; eps++)
                {
                    for(auto lm: ArrayCont[i].lm)
                    {
                        ArrayCont[i].be1[eps][lm][iR] += ArrayCont[i].bev[eps][lm][iR]*dt;
                    }
                }
            }
        }
    }
    else{    
        for(int i=0;i<ArrayCont.size();i++){          
            #pragma omp parallel for
            for (int eps=0; eps<ArrayCont[i].NE; eps++)
            {
                for(auto lm: ArrayCont[i].lm)
                {
                    ArrayCont[i].be1[eps][lm][0] += ArrayCont[i].bev[eps][lm][0]*dt;
                }
            }
        }
    }
}

void Runge_Kutta_Ad_Cont(double& dt, vec1C& ArrayCont, int ContMode)
{
    int NR = ArrayCont[0].Vte[0][0].size();
    
    if (NR>1) {
        for(int i=0;i<ArrayCont.size();i++){
            for (int iR=0; iR<NR; iR++) // CHANGE CONTMODE BEFORE THAT LOOP FOR INCREASING VELOCITY
            {
                if(ContMode==1){ // If here much faster than inside the other loops
                    for (int eps=0; eps<ArrayCont[i].NE; eps++){
                        for(auto lm: ArrayCont[i].lm)
                        {
                            ArrayCont[i].be1[eps][lm][iR] = ArrayCont[i].be[eps][lm][iR] + ArrayCont[i].bev[eps][lm][iR]*dt;
                        }
                    }
                }
                else if(ContMode==2){
                    for (int eps=0; eps<ArrayCont[i].NE; eps++){
                        for(auto lm: ArrayCont[i].lm)
                        {
                            ArrayCont[i].be2[eps][lm][iR] = ArrayCont[i].be[eps][lm][iR] + ArrayCont[i].bev[eps][lm][iR]*dt;
                        }
                    }
                }
                else{
                    for (int eps=0; eps<ArrayCont[i].NE; eps++){
                        for(auto lm: ArrayCont[i].lm)
                        {
                            ArrayCont[i].be[eps][lm][iR] = ArrayCont[i].be1[eps][lm][iR] + ArrayCont[i].bev[eps][lm][iR]*dt;
                        }
                    }
                }
            }
        }
    }
    else{    
        for(int i=0;i<ArrayCont.size();i++){          
            if(ContMode==1){ // If here much faster than inside the other loops
                // #pragma omp parallel for
                for (int eps=0; eps<ArrayCont[i].NE; eps++)
                {
                    for(auto lm: ArrayCont[i].lm)
                    {
                        ArrayCont[i].be1[eps][lm][0] = ArrayCont[i].be[eps][lm][0] + ArrayCont[i].bev[eps][lm][0]*dt;
                    }
                }
            }
            else if(ContMode==2){
                // #pragma omp parallel for
                for (int eps=0; eps<ArrayCont[i].NE; eps++){
                    for(auto lm: ArrayCont[i].lm)
                    {
                        ArrayCont[i].be2[eps][lm][0] = ArrayCont[i].be[eps][lm][0] + ArrayCont[i].bev[eps][lm][0]*dt;
                    }
                }
            }
            else{
                // #pragma omp parallel for
                for (int eps=0; eps<ArrayCont[i].NE; eps++){
                    for(auto lm: ArrayCont[i].lm)
                    {
                        ArrayCont[i].be[eps][lm][0] = ArrayCont[i].be1[eps][lm][0] + ArrayCont[i].bev[eps][lm][0]*dt;
                    }
                }
            }
        }
    }
}

void Calculate_Cont_Ampl(double& EF1, double& EF2,vec1C& ArrayCont){
    int NR = ArrayCont[0].Vte[0][0].size();
    for(int i=0;i<ArrayCont.size();i++){      
        ArrayCont[i].load_dip_BS(NR); // Reset summatory
        int count =0;
        for (int iR=0; iR<NR; iR++) // CHANGE CONTMODE BEFORE THAT LOOP FOR INCREASING VELOCITY
        {
            for (int eps=0; eps<ArrayCont[i].NE; eps++)
            {
                complexd pulse1, pulse2;
                for(auto lm: ArrayCont[i].lm)
                {
                    count=0;
                    for(auto j:ArrayCont[i].UniqueStates){ // Iterating throw different couplings
                        if(ArrayCont[i].Allow[0][j] == 1){
                            pulse1 = (ArrayCont[i].DIPpump)[ArrayCont[i].Indexes[0][j]][eps][lm][0]*EF1;
                            cout << "pulse1 " << pulse1 << endl;
                        }
                        else pulse1 = (0,0);

                        if(ArrayCont[i].Allow[1][j] == 1){
                            pulse2 = (ArrayCont[i].DIPprobe)[ArrayCont[i].Indexes[1][j]][eps][lm][0]*EF2;
                        }
                        else pulse2 = (0,0);

                        ArrayCont[i].dip_BS[count][0] += c1*(pulse1 + pulse2)*ArrayCont[i].be[eps][lm][0]*ArrayCont[i].dE;
                        // ArrayCont[i].dip_BS[count][0] += c1*(pulse1 + pulse2)*(ArrayCont[i].be[eps][lm][0] / complexd(ArrayCont[i].NE));
                        count++;
                    }
                }
            }
        }
    }
}
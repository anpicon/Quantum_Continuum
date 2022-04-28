/*
 *  loading_parameters.h
 *  
 *  Created by APA on 13/02/2019
 *
 *
 */

#include "tk_spline.h"
#include "Continuum.h"

//Notation for saving spherical harmonic basis
int SpH(int l, int m, int mmax){
    int A=0;
    if (l<=mmax) {
        for (int j=0; j<l; j++) A+=(2*j+1);
        A--;
        if (m==0) A+=l+1;
        else A+=l+1+m;
    }
    else {
        for (int j=0; j<=mmax; j++) A+=(2*j+1);
        for (int j=mmax+1; j<l; j++) A+=(2*mmax+1);
        A--;
        if (m==0) A+=mmax+1;
        else A+=mmax+1+m;
    }
    return A;
}

void Read_PES(ifstream& fp_input,vec1d& R, vec2d& PESi)
{
    if (!fp_input.is_open())
    {
        printf("Error opening PES file ...\n");
        exit(1);
    }
    printf("Reading PES.txt file ...\n");
    
    //double dummy=0.;
    string dummystring;
    
    int NEi; fp_input >> dummystring >> NEi; //read number electronics states printed in file
    int NR; fp_input >> dummystring >> NR; //read number of nuclear grid printed in file
    
    printf("#NEi= %d\n",NEi);
    
    vec2d    PES(NEi,vec1d(NR,0.));
    vec1d    Rgrid(NR,0.);
    
    if (NEi != int(PESi.size()))
    {
        cout << "Number of electronic states in input is different than in the printed file" << endl;
        exit(1);
    }
    
    fp_input >> dummystring;
    for (int Ei=0; Ei<NEi; Ei++) {
        fp_input >> dummystring;
    }
    
    if (NR>1) {
        for (int iR=0; iR<NR; iR++)
        {
            fp_input >> Rgrid[iR];
            for (int Ei=0; Ei<NEi; Ei++) fp_input >> PES[Ei][iR];
        }
        
        printf("Initial grid R=(%3.3f, %3.3f), #R=%d , dR= %1.3f\n", Rgrid[0],Rgrid[NR-1],NR,abs(Rgrid[1]-Rgrid[0]));
        
        int NRi=PESi[0].size();
        if (NR > NRi)
        {
            cout << "Number of interpolated points for nuclear grid is lower than current data" << endl;
            exit(1);
        }
        //We find the minimum energy in the ground state and we shift all energies
        double Eshift=10000000.;
        for (int iR=0; iR<NR; iR++) if(PES[0][iR]<Eshift) Eshift=PES[0][iR];
        for (int Ei=0; Ei<NEi; Ei++){
            for (int iR=0; iR<NR; iR++){
                PES[Ei][iR]-=Eshift;
            }
        }
        printf("Energies have been shifted by the minimum energy of the ground state: Eshift=%3.5f a.u.\n", Eshift);
        //Starting the interpolation for a finer nuclear grid
        tk::spline s;
        for (int Ei=0; Ei<NEi; Ei++)
        {
            s.set_points(Rgrid,PES[Ei]);
            for (int iR=0; iR<NRi; iR++)
            {
                double x=Rgrid[0]+iR*(Rgrid[NR-1]-Rgrid[0])/(NRi-1);
                PESi[Ei][iR]=s(x);
                R[iR]=x;
            }
        }
        printf("Interpolated grid R=(%3.3f, %3.3f), #R=%d , dR= %1.3f\n", R[0],R[NRi-1],NRi,abs(R[1]-R[0]));
    }
    else if (NR==1) {
        for (int iR=0; iR<NR; iR++)
        {
            fp_input >> Rgrid[iR];
            for (int Ei=0; Ei<NEi; Ei++) fp_input >> PESi[Ei][iR];
        }
        //We find the minimum energy in the ground state and we shift all energies
        double Eshift=10000000.;
        for (int iR=0; iR<NR; iR++) if(PES[0][iR]<Eshift) Eshift=PESi[0][iR];
        for (int Ei=0; Ei<NEi; Ei++){
            for (int iR=0; iR<NR; iR++){
                PESi[Ei][iR]-=Eshift;
            }
        }
        printf("Energies have been shifted by the minimum energy of the ground state: Eshift=%3.5f a.u.\n", Eshift);
        
        printf("Initial grid R=(%3.3f, %3.3f), #R=%d , dR= %1.3f\n", Rgrid[0],Rgrid[NR-1],NR,0.);
    }
    
}//End Read_PES

void Calculate_GS(Kinetic NuclearKE,vec1x& b0, vec1d& PES)
{
    int NR=PES.size();
    arma::mat H(NR,NR);
    arma::vec lambda(NR);
    arma::mat U(NR,NR);
    for (int iR=0; iR<NR; iR++){
        H(iR,iR)=(2.*NuclearKE.Energy)+PES[iR];
        if(iR!=NR-1) H(iR,iR+1)=-1.*NuclearKE.Energy;
    }
    H=arma::symmatu(H);
    arma::eig_sym(lambda,U,H);
    bool printEigenStates=true; //Printing first 20 vibrational states
    if (printEigenStates==true) {
        ofstream fp_output; fp_output.open("check_eigenvectors.txt");
        for (int iR=0; iR<NR; iR++){
            fp_output << iR << " " << PES[iR] << " ";
            for (int Ei=0; Ei<20; Ei++){
                fp_output << U(iR,Ei) << " ";
            }
            fp_output << endl;
        }
        fp_output.close();
    }
    if (printEigenStates==true) {
        ofstream fp_output; fp_output.open("check_hamiltonian.txt");
        for (int iR=0; iR<20; iR++){
            fp_output << iR << " " << PES[iR] << " ";
            for (int Ei=0; Ei<20; Ei++){
                fp_output << H(iR,Ei) << " ";
            }
            fp_output << endl;
        }
        fp_output.close();
    }
    
    printf("Energies of the first vibrational states:\n");
    for (int iR=0; iR<10; iR++) printf("v%d   %3.6f a.u.\n", iR, lambda[iR]);
    
    
    for (int iR=0; iR<NR; iR++) b0[iR]=U(iR,0); //we are passing the first eigenvector
}//End Calculate_GS



// Dipoles for bound states
void Read_Dipoles(ifstream& fp_input,vec3d& Dip1i,vec3d& Dip2i,vec1d& u1,vec1d& u2)
{
    int NRi=Dip1i[0][0].size();
    
    if (!fp_input.is_open())
    {
        cout << "Error opening Dipoles file" << endl;
        exit(1);
    }
    
    string dummystring;
    
    printf("Reading Dipoles.txt file ...\n");
    int NEi; fp_input >> dummystring >> NEi; //read number electronics states printed in file
    int NR; fp_input >> dummystring >> NR; //read number of nuclear grid printed in file
    if (NEi != int(Dip1i.size()))
    {
        cout << "Number of electronic states in input is different than in the printed file" << endl;
        exit(1);
    }
    
    vec3d    Dip1(NEi,vec2d(NEi,vec1d(NR,0.))); //Dipole transitions pump: Dip[Ei][Ej][R]
    vec3d    Dip2(NEi,vec2d(NEi,vec1d(NR,0.))); //Dipole transitions Stokes: Dip[Ei][Ej][R]
    vec1d    Rgrid(NR,0.);
    
    fp_input >> dummystring;
    for (int Ei=0; Ei<NEi; Ei++) for (int Ej=Ei+1; Ej<NEi; Ej++) { fp_input >> dummystring >> dummystring >> dummystring;}
    
    if (NR>1) {
        for (int iR=0; iR<NR; iR++)
        {
            fp_input >> Rgrid[iR];
            for (int Ei=0; Ei<NEi; Ei++)
            {
                for (int Ej=Ei+1; Ej<NEi; Ej++)
                {
                    double x,y,z; fp_input >> x >> y >> z;
                    Dip1[Ei][Ej][iR]=x*u1[0]+y*u1[1]+z*u1[2];
                    Dip2[Ei][Ej][iR]=x*u2[0]+y*u2[1]+z*u2[2];
                }
            }
        }

        //Starting the interpolation for a finer nuclear grid
        tk::spline s1;
        tk::spline s2;
        for (int Ei=0; Ei<NEi; Ei++)
        {
            for (int Ej=Ei+1; Ej<NEi; Ej++)
            {
                s1.set_points(Rgrid,Dip1[Ei][Ej]);
                s2.set_points(Rgrid,Dip2[Ei][Ej]);
                for (int iR=0; iR<NRi; iR++)
                {
                    double x=Rgrid[0]+iR*(Rgrid[NR-1]-Rgrid[0])/(NRi-1);
                    Dip1i[Ei][Ej][iR]=s1(x);
                    Dip2i[Ei][Ej][iR]=s2(x);
                }
            }
        }
    }
    else if (NR==1){
        for (int iR=0; iR<NR; iR++)
        {
            fp_input >> Rgrid[iR];
            for (int Ei=0; Ei<NEi; Ei++)
            {
                for (int Ej=Ei+1; Ej<NEi; Ej++)
                {
                    double x,y,z; fp_input >> x >> y >> z;
                    Dip1i[Ei][Ej][iR]=x*u1[0]+y*u1[1]+z*u1[2];
                    Dip2i[Ei][Ej][iR]=x*u2[0]+y*u2[1]+z*u2[2];
                }
            }
        }
    }
    
    //we assume dipole transition to be symmetric
    for (int Ei=0; Ei<NEi; Ei++)
    {
        for (int Ej=Ei+1; Ej<NEi; Ej++)
        {
            for (int iR=0; iR<NRi; iR++){
                Dip1i[Ej][Ei][iR]=Dip1i[Ei][Ej][iR];
                Dip2i[Ej][Ei][iR]=Dip2i[Ei][Ej][iR];
            }
        }
    }
    
}//End Read_Dipoles

void Read_Cont_PES(ifstream& fp_input,vec1C& ArrayCont,int& ind, vec2d& PES)
{
    if (!fp_input.is_open())
    {
        printf("Error opening Continuum PES file ...\n");
        exit(1);
    }
    printf("Reading Continuum_PES.txt file ...\n");
    
    string dummystring;
    
    int NR; fp_input >> dummystring >> NR; //read number of nuclear grid printed in file
    
    ArrayCont[ind].load_PES(NR); // Give dimensions to PES
    vec1d Rgrid(NR,0.);
    
    fp_input >> dummystring >> dummystring; // Reading headers
    if (NR==1) {
        fp_input >> Rgrid[0];
        fp_input >> ArrayCont[ind].PES[0];
        printf("Initial grid R=(%3.3f, %3.3f), #R=%d , dR= %1.3f\n", Rgrid[0],Rgrid[NR-1],NR,0.);
    }
    
    // if (NR>1) {
    //     for (int iR=0; iR<NR; iR++)
    //     {
    //         fp_input >> Rgrid[iR];
    //         for (int Ei=0; Ei<NEi; Ei++) fp_input >> PES[Ei][iR];
    //     }
        
    //     printf("Initial grid R=(%3.3f, %3.3f), #R=%d , dR= %1.3f\n", Rgrid[0],Rgrid[NR-1],NR,abs(Rgrid[1]-Rgrid[0]));
        
    //     int NRi=ArrayCont[ind].PES.size();
    //     if (NR > NRi)
    //     {
    //         cout << "Number of interpolated points for nuclear grid is lower than current data" << endl;
    //         exit(1);
    //     }
    //     //We find the minimum energy in the ground state and we shift all energies
    //     double Eshift=10000000.;
    //     for (int iR=0; iR<NR; iR++) if(PES[0][iR]<Eshift) Eshift=PES[0][iR];
    //     for (int Ei=0; Ei<NEi; Ei++){
    //         for (int iR=0; iR<NR; iR++){
    //             PES[Ei][iR]-=Eshift;
    //         }
    //     }
    //     printf("Energies have been shifted by the minimum energy of the ground state: Eshift=%3.5f a.u.\n", Eshift);
    //     //Starting the interpolation for a finer nuclear grid
    //     tk::spline s;
    //     for (int Ei=0; Ei<NEi; Ei++)
    //     {
    //         s.set_points(Rgrid,PES[Ei]);
    //         for (int iR=0; iR<NRi; iR++)
    //         {
    //             double x=Rgrid[0]+iR*(Rgrid[NR-1]-Rgrid[0])/(NRi-1);
    //             PESi[Ei][iR]=s(x);
    //             R[iR]=x;
    //         }
    //     }
    //     printf("Interpolated grid R=(%3.3f, %3.3f), #R=%d , dR= %1.3f\n", R[0],R[NRi-1],NRi,abs(R[1]-R[0]));
    // }        
}//End Read_Cont_PES

 //Read Dipoles for continuum states
void Read_Cont_Dipoles(ifstream& fp_input,vec1C& contstate,int& ind,int& Ej,vec1d& u1,vec1d& u2)
{
    //int NRi=(contstate[ind].DIP).size();
    if (!fp_input.is_open())
    {
        cout << "Error opening Contimuum Dipoles " << ind << "file" << endl;
        exit(1);
    }
 
    string dummystring;
    printf("Reading Continuum_Dipoles.txt file from Continuum State %i ...\n",ind);
    int NE; fp_input >> dummystring >> NE; //read number of Continuum energies in file
    int NR; fp_input >> dummystring >> NR; //read number of nuclear grid printed in file
    int lmax; fp_input >> dummystring >> lmax; //read number of the maximum l quantum number
    int mmax; fp_input >> dummystring >> mmax; //read number of the maximum m quantum number
    
    vec1d Rgrid(NR,0.);

    contstate[ind].NE=int(((contstate[ind].Emax-contstate[ind].Emin)/contstate[ind].dE)); //Change to read parameters
    cout << "   Number of Continuum energies : " << contstate[ind].NE << endl;  // DEBUG
    int foo = SpH(contstate[ind].Lmax,contstate[ind].Mmax,contstate[ind].Mmax);
    // cout << "DIMENSION OF HEADER: " << (foo+1)*3+1 << endl; /// DEBUG

    if (NE != contstate[ind].NE)
    {
        cout << "Number of Continuum Energies in input is different than in the printed file" << endl;
        exit(1);
    }
    for(int i=0;i<(foo+1)*3+1;i++) fp_input >> dummystring; // Read headers
    if(Ej==0){
        int sizepump = (contstate[ind].StatCouplPump).size();
        int sizeprobe = (contstate[ind].StatCouplProbe).size();
        int lm= SpH(lmax,mmax,mmax)+1;
        contstate[ind].load_DIP(sizepump, contstate[ind].NE, NR, lm, contstate[ind].DIPpump); // Give dimension DIP Pump
        contstate[ind].load_DIP(sizeprobe, contstate[ind].NE, NR, lm, contstate[ind].DIPprobe); // Give dimension DIP Probe
        contstate[ind].load_E(); // Load Continuum photoelectron energies
        contstate[ind].load_RKvariables(contstate[ind].NE, NR, lm); // Initialize Runge-Kutta vectors
    }

    int bar = (contstate[ind].StatCouplPump).size();

    if (NR==1){   
        for (int iR=0; iR<NR; iR++)
        {
            fp_input >> Rgrid[iR];
            for (int eps=0; eps<contstate[ind].NE; eps++)
            {
                for (int L=0; L<=lmax; L++)
                {
                    int M1=( L<=mmax ? L : mmax);
                    for (int M=-M1; M<=M1; M++)
                    {
                        complexd x,y,z; fp_input >> x >> y >> z;
                        if(Ej<bar)        (contstate[ind].DIPpump)[Ej][eps][SpH(L,M,mmax)][0]      = x*u1[0]+y*u1[1]+z*u1[2];
                        else if(bar == 0) (contstate[ind].DIPprobe)[Ej][eps][SpH(L,M,mmax)][0]     = x*u2[0]+y*u2[1]+z*u2[2];
                        else              (contstate[ind].DIPprobe)[bar-Ej][eps][SpH(L,M,mmax)][0] = x*u2[0]+y*u2[1]+z*u2[2];
                    }
                }
            }
        }
    }
//     if (NR>1){   
//         for (int iR=0; iR<NR; iR++)
//         {
//             fp_input >> Rgrid[iR];
//             for (int eps=0; eps<contstate[ind].NE; eps++)
//             {
//                 complexd dip = (0,0);
//                 for (int L=0; L<=lmax; L++)
//                 {
//                     int M1=( L<=mmax ? L : mmax);
//                     for (int M=-M1; M<=M1; M++)
//                     {
//                         cout << "iR: " << iR << endl;              // DEBUG
//                         cout << "eps: " << eps << endl;            // DEBUG
//                         cout << "l: " << L << " m: " << M << endl; // DEBUG
//                         complexd x,y,z; fp_input >> x >> y >> z;
//                         cout << x << " " << y << " " << z << endl; // DEBUG
//                         if(Ej<bar)        (contstate[ind].DIPpump)[Ej][eps][SpH(L,M,mmax)][iR]      = x*u1[0]+y*u1[1]+z*u1[2];
//                         else if(bar == 0) (contstate[ind].DIPprobe)[Ej][eps][SpH(L,M,mmax)][iR]     = x*u2[0]+y*u2[1]+z*u2[2];
//                         else              (contstate[ind].DIPprobe)[bar-Ej][eps][SpH(L,M,mmax)][iR] = x*u2[0]+y*u2[1]+z*u2[2];
//                     }
//                 }
//             }
//         }
//     }
//         // //Starting the interpolation for a finer nuclear grid
} // End Read Continuum Dipoles

// void Read_Cont_Decays(ifstream& fp_input,vec1C& ArrayCont, int& ind)
// {
//     if (!fp_input.is_open())
//     {
//         cout << "Error opening Decays file" << endl;
//         exit(1);
//     }
    
//     printf("Reading Decays.txt file ...\n");
    
//     double dummy;
//     int size = (ArrayCont[ind].StatCouplPump).size() + (ArrayCont[ind].StatCouplProbe).size();
//     for (int i=0; i<size; i++)
//     {
//         fp_input >> dummy >> ArrayCont[ind].Gamma[i];
//         //printf("For state %d the decay width is: %3.4f eV -> Lifetime: %3.4f fs\n",Ei,Gamma[Ei],0.659/Gamma[Ei]);
//         ArrayCont[ind].Gamma[i]*=energy_eV_au; //we convert it from eV to a.u.
//     }
// }//End Read_Cont_Decays

void Read_Decays(ifstream& fp_input,vec1d& Gamma, int NEi)
{
    if (!fp_input.is_open())
    {
        cout << "Error opening Decays file" << endl;
        exit(1);
    }
    
    printf("Reading Decays.txt file ...\n");   
    double dummy;
    for (int Ei=0; Ei<NEi; Ei++)
    {
        fp_input >> dummy >> Gamma[Ei];
        printf("For state %d the decay width is: %3.4f eV -> Lifetime: %3.4f fs\n",Ei,Gamma[Ei],0.659/Gamma[Ei]);
        Gamma[Ei]*=energy_eV_au; //we convert it from eV to a.u.
    }
}//End Read_Decays



void Read_Initial_State(ifstream& fp_input,vec2x& bgs)
{
    if (!fp_input.is_open())
    {
        cout << "Error opening Initial_State file" << endl;
        exit(1);
    }
    
    printf("Reading Initial_State.txt file ...\n");
    
    int NEi=bgs.size();
    int NR=bgs[0].size();
    double dummy;
    
    for (int iR=0; iR<NR; iR++)
    {
        fp_input >> dummy;
        for (int Ei=0; Ei<NEi; Ei++)
        {
            //complexd c;
            fp_input >> bgs[Ei][iR];
        }
    }
    
    //check normalization
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
    
    
    cout << "Norm initial state: ";
    for (int Ei=0; Ei<NEi; Ei++)
    {
        cout << Ei << ") " << sumgs[Ei] << ", ";
    }
    cout << endl;
    cout << "Total norm: " << sum << endl;
    
}//End Read_Initial_State

//%%%%%%%%%%%%%%%%%%% HERE WE READ ELECTRONIC STRUCTURE AND PES

fp_input.open(path[0].c_str()); //"INPUTS/PES.txt"
Read_PES(fp_input,R,PES);
fp_input.close();

if (NR!=1)
{
    printf("Nuclear propagation activated.\n");
    dR=R[1]-R[0];
}
else printf("Nuclear propagation deactivated.\n");

Kinetic NuclearKE(molecule,dR);

if (iInitialState=="Read") {
    fp_input.open(path[3].c_str());//"INPUTS/Initial_State.txt"
    Read_Initial_State(fp_input,b0);
    fp_input.close();
}
else if (iInitialState=="Calculate" && NR>1) {
    printf("Calculating ground state ...\n");
    Calculate_GS(NuclearKE,b0[0],PES[0]);
}
else if (NR==1) {
    b0[0][0]=1.;
}
else{
    printf("We can only calculate initial state for NR>1.\n");
    exit(1);
}

// Bound states part
fp_input.open(path[2].c_str()); //path[2]="INPUTS/Dipoles.txt";
Read_Dipoles(fp_input,Dip1,Dip2,pulse1.pol,pulse2.pol);
fp_input.close();

fp_input.open(path[1].c_str()); //"INPUTS/Decays.txt"
Read_Decays(fp_input,Gamma,NEi);
fp_input.close();

// Continuum states part
if(NContStat>=0){
    for(int i=0;i<=NContStat;i++){
        fp_input.open((ArrayCont[i].PESpath).c_str()); // Read double core-hole states
        Read_Cont_PES(fp_input, ArrayCont, i, PES);
        cout << "Read PES of Continuum State " << i << endl;
        fp_input.close();
        
        int Ej = (ArrayCont[i]).StatCouplPump.size() + (ArrayCont[i]).StatCouplProbe.size();

        for(int j=0;j<Ej;j++){
            cout << ((ArrayCont[i].DIPpath)[j]).c_str() << endl;
            fp_input.open(((ArrayCont[i].DIPpath)[j]).c_str());
            Read_Cont_Dipoles(fp_input,ArrayCont,i,j,pulse1.pol,pulse2.pol);  // Read continuum dipoles
            // ArrayCont[i].statesINcommon(); // Array management for RK
            // ArrayCont[i].positionEj(NEi);
            cout << "Read dipoles of Continuum State " << i << ", coupling number: " << j+1 << endl;
            fp_input.close();
        }
        ArrayCont[i].positionEj(); // Prepare coupled arrays for Runge-Kutta

        fp_input.open((ArrayCont[i].DECAYSpath).c_str()); // Read continuum decays
        ArrayCont[i].load_Gamma(Ej);
        Read_Cont_Decays(fp_input, ArrayCont, i);
        cout << "Read Continuum Decays " << i << endl;
        fp_input.close();
    }
}
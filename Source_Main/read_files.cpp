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

//path[2]="INPUTS/Dipoles.txt";
fp_input.open(path[2].c_str());
Read_Dipoles(fp_input,Dip1,Dip2,pulse1.pol,pulse2.pol);
fp_input.close();

fp_input.open(path[1].c_str());//"INPUTS/Decays.txt"
Read_Decays(fp_input,Gamma,NEi);
fp_input.close();

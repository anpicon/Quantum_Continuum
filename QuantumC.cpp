/* Created by APA and Solene Oberli 13/02/2019
 Time dynamics program in real space for quantum control simulations
 */

#include "Source_Main/include_headers.cpp"
using namespace std;

int main (int argc, char* argv[])
{
    //%%%%%%%%%%%%%%%%%%% SIMULATION VARIABLES
    //clock_t clock1,clock2; clock1=clock();
    //double dt = 0.01;
    //string iWFD="Print"; //"Print"
    //string iInitialState="Calculate"; //"Read" or "Calculate"
    //printf("Time step RK: %1.3f a.u.\n", dt);
    //double time, dummy; // Dummy variable used to read unused data
    
    ifstream fp_input;
    if ( argc != 2 )
    {
        cout<<"usage: "<< argv[0] <<" <filename>" << endl;
        exit(1);
    }
    else
    {
        fp_input.open(argv[1]);
        if (!fp_input.is_open())
        {
            cout << "error opening file " << argv[1] << endl;
            exit(1);
        }
    }
    
    time_t time1 = time(0);
    tm* now = localtime(&time1);
    int ihour = now->tm_hour;
    int iday = now->tm_mday;
    int imin = now-> tm_min;
    int isec = now-> tm_sec;

    printf("\n\n\nCalculation started day %4d/%02d/%02d at %02d.%02d.%02d \n", (now->tm_year + 1900), (now->tm_mon + 1), (now->tm_mday), (now->tm_hour), (now->tm_min), (now->tm_sec));
    
    //**********************************************************
    //*  Definition of variables                               *
    //**********************************************************
    #include "Source_Main/variables.cpp"
    
    //**********************************************************
    //*  Read INPUT                                            *
    //**********************************************************
    Read_Input(fp_input,dt,iInitialState,NR,tf,NEi,path,pulse1,sigma1,gaussian1,pulse2,sigma2,gaussian2,DELAY,EnergyShift,molecule,iObservables,ArrayCont,NContStat);
    fp_input.close();
    
    printf("molecule:          %s\n", molecule.c_str());
    printf("# states:          %d\n", NEi);
    printf("# nuclear points:  %d\n", NR);
    
    //**********************************************************
    //*  Print laser parameters                                *
    //**********************************************************
    #include "Source_Main/print_laser.cpp"
    
    //**********************************************************
    //*  Allocating memory for large arrays                    *
    //**********************************************************
    #include "Source_Main/allocate_arr_big.cpp"
    
    //**********************************************************
    //*  Reading files with electronic structure               *
    //**********************************************************
    #include "Source_Main/read_files.cpp"
    
    //**********************************************************
    //*  TIME PROPAGATION: RUNGE KUTTA                         *
    //**********************************************************
    EF1=0.; EF2=0.;
    int ncyc=4;
    if(tf<1.e-8){
        if(gaussian1) tf = 2.*ncyc*sigma1;
        else          tf = sigma1*pulse1.get_Period();
    }
    int nstep = int(tf/dt) ;
    printf("Number of steps: %3d \nFinal time (fs): %5.2f\n", nstep, tf*time_au_fs);
    //double contador=sigma1/25.;
    double contador=tf/100.;
    
    //for(int it=0;it<=1;it++)
    for(int it=0;it<=nstep;it++)
    {
        if ( it %  (nstep/10) == 0)         printf("Step # %5d             ---------->  %3.2f %% \n", it, it*100./nstep);
        double time=it*dt;
        
        //----- Print Electric Fields and Amplitudes
        // sum of the probability, shoud be one
        if(time>=icont*contador)
        {
            if (!check_conservation(fp_population,b0,time,ArrayCont,fp_Contpopulation))
            {
                // printf("probability not conserved\n");
                // exit(1);
            }
            if (iObservables[0])
            {
                string name_file; stringstream sname;
                sname.seekp(0,ios::beg); sname << icont;
                name_file= "Output/wfd_" + sname.str() + ".txt";
                ofstream fp_wf; fp_wf.open(name_file.c_str());
                PrintWFD(fp_wf,R,b0);
                fp_wf.close();
            }
            if (iObservables[1])
            {
                string name_file; stringstream sname;
                sname.seekp(0,ios::beg); sname << icont;
                name_file= "Output/ampl_" + sname.str() + ".txt";
                ofstream fp_wf; fp_wf.open(name_file.c_str());
                PrintAmpl(fp_wf,PES,b0);
                fp_wf.close();
            }
            if (iObservables[2])
            {
                string name_file; stringstream sname;
                sname.seekp(0,ios::beg); sname << icont;
                name_file= "Output/XPS/xps_" + sname.str() + ".txt";
                ofstream fp_wf; fp_wf.open(name_file.c_str());
                PrintXPS(fp_wf, time, ArrayCont);
                fp_wf.close();
                // cout << time << " " << time*time_au_fs << endl;
            }
            PrintEF(fp_EF,EF1,EF2,time);
            icont++;
        }
        // Include if(ArrayCont.size()>0)
        if (NR>1) {
            #include "Source_Main/RK_nuclear.cpp"
        }
        else {
            #include "Source_Main/RK_fixed.cpp"
            #include "Source_Main/RK_fixed_Cont.cpp"
        }
        
    }//--------- END TIME LOOP
    
    //**********************************************************
    //*  Ending program                                        *
    //**********************************************************
    fp_population.close(); fp_EF.close(); fp_Contpopulation.close();
    
    /*
    if (iObservables[1]){
        fp_ampl.open("Output/FinalAmpl.txt");
        PrintAmpl(fp_ampl,PES,b0);
        fp_ampl.close();
    }*/
    
    time_t time_ = time(0);
    tm* now_ = localtime(&time_);

    printf("Calculation completed day %4d/%02d/%02d at %02d.%02d.%02d \n", (now_->tm_year + 1900), (now_->tm_mon + 1), (now_->tm_mday), (now_->tm_hour), (now_->tm_min), (now_->tm_sec));

    int deltad =  (now_->tm_mday) - iday;
    int deltah =  (now_->tm_hour) - ihour;

    int deltam =  (now_->tm_min)  - imin;
    int deltas =  (now_->tm_sec)  - isec;


    if (deltas < 0)
    {
        deltas += 60;
        deltam -= 1;
    }
    if (deltam < 0)
    {
        deltam += 60;
        deltah -= 1;
    }
    if (deltah < 0)
    {
        deltah += 24;
        deltad -= 1;
    }

    printf("Elapsed time:  %02d days, %02d hours, %02d minutes and %02d seconds.\n", deltad, deltah, deltam, deltas);
    
    clock2=clock();
    printf("Time of computation: %8.4f s\n", (double) (clock2-clock1)/CLOCKS_PER_SEC);
    
return 0;
}

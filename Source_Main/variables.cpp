
//**********************************************************
//*  General variables                                     *
//**********************************************************
clock_t clock1,clock2; clock1=clock();
string molecule;
double dummy;

string iInitialState="Calculate"; //"Read" or "Calculate"

//**********************************************************
//*  Variables TDSE                                        *
//**********************************************************
double dt = 0.01;  // time step
int NR=1;          // # of nuclear points
double dR=1.;      // resolution of the internuclear grid
vec1d R(1,0.);     // grid of internuclear distance
double tf=0.;      // final time

double dt2=(dt/2.0); // for RK
double dt3=(dt/3.0); // for RK
double dt6=(dt/6.0); // for RK

int icont=0;       // for RK, it controlls when observables are printed

//**********************************************************
//*  Electronic states                                     *
//**********************************************************
int NEi;           //Number of electronic states
vec1s path;        //Array of addresses

double EnergyShift=0.; // Energy shift in order to make 0 in the min point of PES ground state

path.resize(4);
//path[0] Array of addresses pointing to PES
//path[1] Array of addresses pointing to decays for each state
//path[2] Array of addresses pointing to dipole moments
//path[3] Array of addresses pointing to initial state population

//**********************************************************
//*  Continuum States                                      *
//**********************************************************
int NPump;       // Number of states coupled with Pump pulse
int NProbe;      // Number of states coupled with Probe pulse
int NContStat;   // (Npump+Nprobe - 1, as we start counting from 0) It will index the ArrayCont vector
vec1C ArrayCont; // Array of Continuum states

//**********************************************************
//*  Observables                                           *
//**********************************************************
vector<bool> iObservables(3);                // Vector to keep track of the observables to print or compute
iObservables[0]=false;                       // Print wave functions
iObservables[1]=false;                       // Print amplitudes and energies after pulse
iObservables[2]=false;                       // Print XPS spectrum


//**********************************************************
//*  Variables of the laser pulses                         *
//**********************************************************
Laser pulse1,pulse2;
//vec1d u1(3,0.);                                  //polarization pulse 1
//vec1d u2(3,0.);                                  //polarization pulse 2
//double wx1,wx2;                                  //frequencies pulse 1 and 2
double EF1,EF2;                                  //electric field amplitude for pulse 1 and 2
double sigma1,sigma2;                            //define width of time pulse profile
bool gaussian1=true;                             //define sin2 (false) or gaussian (true,default) envelope
bool gaussian2=true;                             //define sin2 (false) or gaussian (true,default) envelope
double radiusf=1.e-6; radiusf*=space_m_au;       //radius focus
double areafocus=pi*radiusf*radiusf;             //area focus
double GAMMA1x=0.0;                              //Photoionization Cross Sections: Gamma_x = (EF1*EF1/2.) * (sum_i |Dip_x|^2)
double DELAY=0.;                                 //Delay between pump and probe pulse

pulse1.set_CEP(0.);
pulse2.set_CEP(0.);

//**********************************************************
//*  Output files                                          *
//**********************************************************
system("mkdir -p Output");
ofstream fp_population; fp_population.open("Output/Population.txt");
ofstream fp_EF;  fp_EF.open("Output/EF.txt");
ofstream fp_ampl;
system("mkdir -p Output/XPS");
ofstream fp_Contpopulation; fp_Contpopulation.open("Output/ContPopulation.txt");

//ofstream fp_PES;  fp_PES.open("Output/PES.txt");
//ofstream fp_Initial;  fp_Initial.open("Output/Initial.txt");





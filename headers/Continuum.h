#include <vector>
#pragma once
using namespace std;

typedef complex<double>                         complexd;
typedef vector<complexd>                        vec1x;
typedef vector<vec1x>                           vec2x;
typedef vector<vec2x>                           vec3x;
typedef vector<vec3x>                           vec4x;

//Notation for saving spherical harmonic basis
int spH(int l, int m, int mmax){
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

class Continuum{
public:
	Continuum(){}
	int State;                      //Number of the continuum state (starts from 0)
	Continuum(int var){             //Initialize the object
		State=var;
	}
	vector<string> DIPpath;         // Continuum Dipole Moments Path (First the ones coupled with Pump, then Probe)
	string PESpath;                 // Double core-hole energies Path
	string DECAYSpath;
	vector<int> StatCouplPump;      // States that couples with Pump
	vector<int> StatCouplProbe;     // States that couples with Probe
	// vector<int> CommonStates;       // States that couples with both Pump and Probe
	vector<vector<int>> Indexes;	// Indexes of couplings
	vector<vector<double>> Allow;	    // Indexes of couplings
	vector<double> PES;	            // Double-core hole energies
	vector<double> Gamma;	// Decays for each continuum state
	double dE;
	vector<double> Emin;
	vector<double> Emax;
	vector<double> BE;                      // Binding Energy   
	int    Lmax;
	int    Mmax;
	int    NE;                      // Number of energies: we define it after having Emax,Emin and dE values NE=int((Emax-Emin)/dE)
	vec4x DIPpump;					// Dipoles that couple with pump laser
	vec4x DIPprobe;                 // Dipoles that couple with probe pulse
	void load_DIP(int& NEj, int& NEps, int& NR, int& lm, vec4x& DIP){  // Gives dimensionality to DIP array   
		vec4x v(NEj,vec3x(NEps,vec2x(lm,vec1x(NR,complexd(0,0))))); // Dimension[NEj][NEps][l/m][R] (R has to be the last to interpolate it better)
		DIP = v;
	}
	void load_PES(int& NR){  // Gives dimensionality to PES array   
		PES.resize(NR);
	}
	void load_Gamma(int& Ej){  // Gives dimensionality to Gamma array   
		Gamma.resize(Ej);
	}


	// void statesINcommon(){ // Use it after reading dipoles
	// 	vector<int> v(StatCouplPump.size()+StatCouplProbe.size());
	// 	vector<int>::iterator it, common;
	// 	it = set_intersection(StatCouplPump.begin(), StatCouplPump.end(), StatCouplProbe.begin(), StatCouplProbe.end(), v.begin());
	// 	for(common = v.begin(); common != it; ++common) CommonStates.push_back(*common);
	// 	for(int i=0;i<CommonStates.size();i++){
	// 		for(int j=0;j<StatCouplPump.size();j++){
	// 			if(CommonStates[i]==StatCouplPump[j]) StatCouplPump.erase(StatCouplPump.begin()+j);
	// 		}
	// 		for(int k=0;k<StatCouplProbe.size();k++){
	// 			if(CommonStates[i]==StatCouplProbe[k]) StatCouplProbe.erase(StatCouplProbe.begin()+k);
	// 		}
	// 	}
	// }

	void positionEj(int& NEi){
		int count=-1;
		for(int i=0;i<2;i++){ //Looping in pump and probe arrays
			for(int j=0;j<NEi;j++){
				if(i==0 && j!=StatCouplPump[j]){
					Indexes[i].push_back(0);
					Allow[i].push_back(0.0);
				}
				else{
					Indexes[i].push_back(++count);
					Allow[i].push_back(1.0);
				}
			}
		}
	}

	void check_DIP(int& Ej, vector<int>& StatCoupl, double *allow, int *ind){
		for(int i=0;i<StatCoupl.size();i++){
			if(StatCoupl[i]==Ej){
				*ind = i;
				*allow = 1.0;
				break;
			}
			else{
				*ind = 0;
				*allow = 0.0;
			}
		}
	}

};
/*
If we want to use circularly polarized laser, we shoul increase dimension to DIP vectors by separating x, y and z components
*/
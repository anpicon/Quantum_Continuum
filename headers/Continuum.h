#include <vector>
#pragma once
using namespace std;

typedef complex<double>                         complexd;
typedef vector<complexd>                        vec1x;
typedef vector<vec1x>                           vec2x;
typedef vector<vec2x>                           vec3x;
typedef vector<vec3x>                           vec4x;
typedef vector<double>                          vec1d;
typedef vector<vec1d>                           vec2d;
typedef vector<vec2d>                           vec3d;

class Continuum{
private:
	void load_UniqueStates(){ // Create an array with non repeated elements of indexes of coupling
		UniqueStates.insert(UniqueStates.begin(), StatCouplPump.begin(), StatCouplPump.end());
    	UniqueStates.insert(UniqueStates.end(), StatCouplProbe.begin(), StatCouplProbe.end());
		sort(UniqueStates.begin(),UniqueStates.end());
		vector<int>::iterator it;
    	it=unique(UniqueStates.begin(),UniqueStates.end());
    	UniqueStates.resize(distance(UniqueStates.begin(),it));
	}
public:
	Continuum(){}
	int State;                      //Number of the continuum state (starts from 0)
	Continuum(int var){             //Initialize the object
		State=var;
	}
	vector<string> DIPpath;         // Continuum Dipole Moments Path (First the ones coupled with Pump, then Probe)
	string PESpath;                 // Double core-hole energies Path
	vector<int> StatCouplPump;      // Bound states that couples with Pump from input
	vector<int> StatCouplProbe;     // Bound states that couples with Probe from input
	vector<int> UniqueStates;		// Indexes of couplings non repeated
	vector<vector<int>> Indexes;	// Indexes of couplings
	vector<vector<bool>> Allow;	    // Allowance of couplings
	vector<double> PES;	            // Double-core hole energies
	double Gamma;	                // Double-core hole state decay
	vector<double> E;				// Array containing the
	double Emax;
	double Emin;
	double dE;
	int maxBoundState;              // Index of the maximum bound state coupled
	int    Lmax;
	int    Mmax;
	vector<int> lm;
	int    NE;                      // Number of energies: we define it after having Emax,Emin and dE values NE=int((Emax-Emin)/dE)
	vec4x DIPpump;					// Dipoles that couple with pump laser
	vec4x DIPprobe;                 // Dipoles that couple with probe pulse
	// Runge-Kutta variables
	vec3x Vte;						// Dipole term
	vec3x bev;						// Derivative term
	vec3x be;						// Continuum amplitude
	vec3x be1;						// Aux1 variable
	vec3x be2;						// Aux2 variable
	vec2x dip_BS;					// Dipole term that removes population from Bound States
	void load_E(){ // It calculates the whole range of Continuum energies
		vector<double> v(NE);
		E = v;
		for(int i=0;i<NE;i++){
			E[i] = (Emin + double(i)*dE);
		}
	}
	void load_RKvariables(int& NEps, int& NR, int& lm){// Initializing RK variables
		load_UniqueStates();
		vec3x v(NEps,vec2x(lm,vec1x(NR,complexd(0,0))));
		vec2x w(UniqueStates.size(),vec1x(NR,complexd(0,0)));
		Vte = v;
		bev = v;
		be  = v;
		be1 = v;
		be2 = v;
		dip_BS = w;
	}
	void load_dip_BS(int& NR){
		load_UniqueStates();
		vec2x w(UniqueStates.size(),vec1x(NR,complexd(0,0)));
		dip_BS = w;
	}
	void load_DIP(int& NEj, int& NEps, int& NR, int& lm, vec4x& DIP){  // Initializing DIP array   
		vec4x v(NEj,vec3x(NEps,vec2x(lm,vec1x(NR,complexd(0,0))))); // Dimension[NEj][NEps][l/m][R] (R has to be the last to interpolate it better)
		DIP = v;
	}
	void load_PES(int& NR){  // Initializing PES array   
		PES.resize(NR);
	}
	void positionEj(){ // It organizes couplings indexes and makes a boolean matrix to indicate allowance of coupligs in each pulse
		int max1 = (StatCouplPump.size() > 0) ? *max_element(StatCouplPump.begin(), StatCouplPump.end()) : 0;
		int max2 = (StatCouplProbe.size() > 0) ? *max_element(StatCouplProbe.begin(), StatCouplProbe.end()) : 0;
		maxBoundState = max(max1,max2); // max bound state to construct a matrix with size: 2 x max1  
		// For Pump pulse
		bool check;
		Indexes.resize(2);
		Allow.resize(2);
		for(int j=0;j<=maxBoundState;j++){
			check=0;
			for(int Ei=0;Ei<StatCouplPump.size();Ei++){
				if(j==StatCouplPump[Ei]){
					Indexes[0].push_back(Ei);
					Allow[0].push_back(1);
					check=1;
				}
			}
			if(check==0){
					Indexes[0].push_back(0);
					Allow[0].push_back(0);
				}
		}
		// For Probe pulse
		for(int j=0;j<=maxBoundState;j++){
			check=0;
			for(int Ei=0;Ei<StatCouplProbe.size();Ei++){
				if(j==StatCouplProbe[Ei]){
					Indexes[1].push_back(Ei);
					cout << Ei << endl;
					Allow[1].push_back(1);
					check=1;
				}
			}
			if(check==0){
				Indexes[1].push_back(0);
				Allow[1].push_back(0);
			}
		}
	}
	void load_sphVec(){
		int count = 0;
		if(Lmax==Mmax) for(int i=0; i<((Lmax+1)*(Lmax+1)); i++) lm.push_back(count++);
		else{
			for(int i=0;i<Lmax+1;i++){
				if(i<=Mmax) for(int j=-i;j<=i;j++) lm.push_back(count++);
				else for(int j=-Mmax;j<=Mmax;j++) lm.push_back(count++);
    		}
		}
		cout << "VECTORRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR" << endl;
		for(auto a:lm) cout << a << " ";
		cout << endl;
	}
};

/*
If we want to use circularly polarized laser, we should increase dimension to DIP vectors by separating x, y and z components
*/
/*
 *  Laser.h
 *  
 *
 *  Created by Antonio Picon on 4/15/11.
 *  Copyright 2011 JILA, University of Colorado. All rights reserved.
 *
 */

#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

//const double pi=4.*atan(1.); Define it in Constants.h
//const double c=137.036; Define it in Constants.h

class Laser{
protected:
	double lambdanm,wl,Period,Intensity,A0,E0,CEP;
public:
    vec1d  pol;
	void set_freq (double a) {wl=a; lambdanm=2.*pi*c*0.0529/wl; Period=2*pi/wl;}
	void set_lambdanm (double a){lambdanm=a; wl=2.*pi*c*0.0529/lambdanm; Period=2*pi/wl;}
	void set_IntensityWcm2 (double a){Intensity=a; A0=sqrt(Intensity/(3.51e+16))/wl; E0=A0*wl;}
	void set_A0 (double a){A0=a; Intensity=A0*A0*(3.51e+16)*wl*wl; E0=A0*wl;}
	void set_E0 (double a){E0=a; A0=E0/wl; Intensity=A0*A0*(3.51e+16)*wl*wl;}
    void set_CEP (double CEP_){CEP=CEP_;}
    void set_pol (vec1d pol_){
        pol.resize(3);
        pol[0]=pol_[0]; pol[1]=pol_[1]; pol[2]=pol_[2];
    }
	void print_par(string&);
	double get_Period(void){return Period;}
	double get_Intensity(void){return Intensity;}
	double get_E0(void){return E0;}
	double get_wl(void){return wl;}
	double PlaneWave (double,double,double);
	double PlaneWave_A (double,double,double,double);
	double PseudoPW (double,double,double,double,double);
	double PseudoPW_A (double,double,double,double,double,double);
	double PseudoPW_Envelope (double,double,double,double,double);
	double Sin2 (double,double,double);
	double Sin2_A (double,double,double,double);
	double Sin2_Envelope (double,double,double);
	double Gaussian(double,double,double);
	double Gaussian_A(double,double,double,double);
	double Gaussian_Envelope(double,double,double);
	double Gaussian_Chirp(double,double,double,double);
	friend class Laser_FT;
};

void Laser::print_par(string& title)
{
    cout << title << "\n";
    printf("*  Energy                                %11.5f au   %11.5f eV                    *\n", wl, wl*energy_au_eV);
    printf("*  Wavelength                            %11.5f au   %11.5f nm                    *\n", 10.*lambdanm*space_A_au,lambdanm);
    printf("*  Period                                %11.5f au   %11.5f fs                    *\n", Period, Period*time_au_fs);
    printf("*  A0                                    %11.5e au   %11.5e V s/m                 *\n", A0, A0);
    printf("*  E0                                    %11.5e au   %11.5e V/ang                 *\n", E0, E0/(Efield_au_Vm*1.e10));
    printf("*  Intensity                             %11.5e au   %11.5e W/cm^2                *\n", Intensity*intensity_Wcm2_au,Intensity);
    printf("*  CEP                                   %11.5e rad                               *\n", CEP);
    printf("*  polarization vectors ->              (%2.4f,%2.4f,%2.4f)                       *\n", pol[0],  pol[1],  pol[2]);
    //printf("*  polarization vectors ->         (%2.4f,%2.4f,%2.4f) (%2.4f,%2.4f,%2.4f)            *\n", pol1.cart[0],  pol1.cart[1],  pol1.cart[2], pol2.cart[0],  pol2.cart[1], pol2.cart[2]);
}

double Laser::PlaneWave(double tinitial, double tfinal,double time){
	if(time>=tinitial && time<=tfinal)return (E0*sin(wl*(time-tinitial)+CEP));
	else return 0.;
}
double Laser::PlaneWave_A(double tinitial, double tfinal,double time, double phase){
	if(time>=tinitial && time<=tfinal)return (A0*sin(wl*(time-tinitial)+phase));
	else return 0.;
}
double Laser::PseudoPW(double tinitial,double t1,double t2,double tfinal,double time){
	if(time>=tinitial && time<t1)return (E0*(time-tinitial)*sin(wl*(time-tinitial)+CEP)/(t1-tinitial));
	if(time>=t1 && time<t2)return (E0*sin(wl*(time-tinitial)+CEP));
	if(time>=t2 && time<=tfinal)return (E0*(time-tfinal)*sin(wl*(time-tinitial)+CEP)/(t2-tfinal));
	else return 0.;
}
double Laser::PseudoPW_A(double tinitial,double t1,double t2,double tfinal,double time,double phase){
	if(time>=tinitial && time<t1)return (A0*(time-tinitial)*sin(wl*(time-tinitial)+phase)/(t1-tinitial));
	if(time>=t1 && time<t2)return (A0*sin(wl*(time-tinitial)+phase));
	if(time>=t2 && time<=tfinal)return (A0*(time-tfinal)*sin(wl*(time-tinitial)+phase)/(t2-tfinal));
	else return 0.;
}
double Laser::PseudoPW_Envelope(double tinitial,double t1,double t2,double tfinal,double time){
	if(time>=tinitial && time<t1)return (E0*(time-tinitial)/(t1-tinitial));
	if(time>=t1 && time<t2)return (E0);
	if(time>=t2 && time<=tfinal)return (E0*(time-tfinal)/(t2-tfinal));
	else return 0.;
}
double Laser::Sin2(double tinitial, double tfinal,double time){
	if(time>=tinitial && time<=tfinal)return E0*sin(wl*(time-tinitial)+CEP)*pow(sin(pi*(time-tinitial)/(tfinal-tinitial)),2);
	else return 0.;
}
double Laser::Sin2_A(double tinitial, double tfinal,double time, double phase){
	if(time>=tinitial && time<=tfinal)return A0*sin(wl*(time-tinitial)+phase)*pow(sin(pi*(time-tinitial)/(tfinal-tinitial)),2);
	else return 0.;
}
double Laser::Sin2_Envelope(double tinitial, double tfinal,double time){
	if(time>=tinitial && time<=tfinal)return E0*pow(sin(pi*(time-tinitial)/(tfinal-tinitial)),2);
	else return 0.;
}
double Laser::Gaussian(double b, double sigma,double time){
	if(pow(time-b,2)<40*sigma*sigma){
		// cout << "WWWWWWWWWWWWWWWWWWWWLLLLLLLLLLLLLLLLLLLLL"<< wl << endl;
		return E0*sin(wl*(time-b)+CEP)*exp(-pow(time-b,2)/(2*sigma*sigma));
	}
	else return 0.;
}
double Laser::Gaussian_A(double b, double sigma,double time, double phase){
	if(pow(time-b,2)<40*sigma*sigma)return A0*cos(wl*(time-b)+phase)*exp(-pow(time-b,2)/(2*sigma*sigma));
	else return 0.;
}
double Laser::Gaussian_Envelope(double b, double sigma,double time){
	if(pow(time-b,2)<40*sigma*sigma)return E0*exp(-pow(time-b,2)/(2*sigma*sigma));
	else return 0.;
}
double Laser::Gaussian_Chirp(double b, double sigma,double time, double chirp){
	if(pow(time-b,2)<16.*sigma*sigma)return E0*sin((wl+chirp)*(time-b))*exp(-pow(time-b,2)/(2*sigma*sigma));
	else return 0.;
}


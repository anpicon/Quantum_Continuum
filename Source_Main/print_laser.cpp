/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*
*  Parameters of the lasers, Pump                          *
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
string title1="\nX rays parameters - Pump:";
pulse1.print_par(title1);

printf("Gaussian profile: Sigma(fs):          %8.4f\n", sigma1*time_au_fs                       );
printf("                  FWHM(fs):           %8.4f\n", 2.0*sqrt(2.0*log(2.0))*sigma1*time_au_fs);
printf("                  FWHM Intensity(fs): %8.4f\n", 2.0*sqrt(log(2.0))*sigma1*time_au_fs    );

printf("Assuming a radius of R = %8.4e(m)\n", radiusf*space_au_m);

dummy=0.0;
for(int it=0;(it*dt)<=(8.0*sigma1);it++)
{
    double time=it*dt;
    EF1=pulse1.Gaussian(4.0*sigma1,sigma1,time+dt/2);         // Gaussian in time, phase is not important if frequency is high
    dummy+=EF1*EF1*dt;
}
printf("Energy required per pulse: %8.4e(mJ)\n", dummy*areafocus*energy_au_mJ);

//printf("Polarization in cartesian coordinates: (ux,uy,uz) -> (%1.2f, %1.2f, %1.2f)\n",ue1[0],ue1[1],ue1[2]);printf("Polarization in cartesian coordinates: (ux,uy,uz) -> (%1.2f, %1.2f, %1.2f)\n",ue1[0],ue1[1],ue1[2]);


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*
*  Variables of the x rays, Stokes                        *
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
string title2="X rays parameters - Stokes:";
pulse2.print_par(title2);

printf("Gaussian profile: Sigma(fs):          %8.4f\n", sigma2*time_au_fs                       );
printf("                  FWHM(fs):           %8.4f\n", 2.0*sqrt(2.0*log(2.0))*sigma2*time_au_fs);
printf("                  FWHM Intensity(fs): %8.4f\n", 2.0*sqrt(log(2.0))*sigma2*time_au_fs    );

printf("Assuming a radius of R = %8.4e(m)\n", radiusf*space_au_m);


dummy=0.0;
for(int it=0;(it*dt)<=(8.0*sigma2);it++)
{
    double time=it*dt;
    EF2=pulse2.Gaussian(4.0*sigma2,sigma2,time+dt/2);         // Gaussian in time, phase is not important if frequency is high
    dummy+=EF2*EF2*dt;
}
printf("Energy required per pulse: %8.4e(mJ)\n", dummy*areafocus*energy_au_mJ);


printf("Delay between pump and Stokes: %8.4e(fs)\n", DELAY*time_au_fs);

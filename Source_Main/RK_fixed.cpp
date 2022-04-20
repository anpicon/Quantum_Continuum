
        //------ 1st step Runge-Kutta
        EF1=pulse1.Gaussian(ncyc*sigma1,sigma1,time);
        EF2=pulse2.Gaussian(ncyc*sigma1+DELAY,sigma2,time);
        // Runge_Kutta_Df_fixed(Vt,b0,bv,PES,Dip1,Dip2,Gamma,EF1,EF2);
        Runge_Kutta_Df_fixed(Vt,b0,bv,PES,Dip1,Dip2,Gamma,EF1,EF2, ArrayCont, Vti);
        Runge_Kutta_Ad(b0,b1,bv,dt6);
        Runge_Kutta_Ad(b0,b2,bv,dt2);
        
        //------ 2nd step Runge-Kutta
        time=it*dt+dt2;
        EF1=pulse1.Gaussian(ncyc*sigma1,sigma1,time);
        EF2=pulse2.Gaussian(ncyc*sigma1+DELAY,sigma2,time);
       
        // Runge_Kutta_Df_fixed(Vt,b2,bv,PES,Dip1,Dip2,Gamma,EF1,EF2);
        Runge_Kutta_Df_fixed(Vt,b2,bv,PES,Dip1,Dip2,Gamma,EF1,EF2, ArrayCont, Vti);
        Runge_Kutta_Ac(b1,bv,dt3);
        Runge_Kutta_Ad(b0,b2,bv,dt2);

        //------ 3rd step Runge-Kutta
        // Runge_Kutta_Df_fixed(Vt,b2,bv,PES,Dip1,Dip2,Gamma,EF1,EF2);
        Runge_Kutta_Df_fixed(Vt,b2,bv,PES,Dip1,Dip2,Gamma,EF1,EF2, ArrayCont, Vti);
        Runge_Kutta_Ac(b1,bv,dt3);
        Runge_Kutta_Ad(b0,b2,bv,dt);

        //------ 4th step Runge-Kutta
        time=(it+1)*dt;
        EF1=pulse1.Gaussian(ncyc*sigma1,sigma1,time);
        EF2=pulse2.Gaussian(ncyc*sigma1+DELAY,sigma2,time);
        
        // Runge_Kutta_Df_fixed(Vt,b2,bv,PES,Dip1,Dip2,Gamma,EF1,EF2);
        Runge_Kutta_Df_fixed(Vt,b2,bv,PES,Dip1,Dip2,Gamma,EF1,EF2, ArrayCont, Vti);
        Runge_Kutta_Ad(b1,b0,bv,dt6);

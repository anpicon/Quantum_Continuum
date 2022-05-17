        // dt2=(dt/2.0);
        // dt3=(dt/3.0);
        // dt6=(dt/6.0);
        //------ 1st step Runge-Kutta
        EF1=pulse1.Gaussian(ncyc*sigma1,sigma1,time);
        EF2=pulse2.Gaussian(ncyc*sigma1+DELAY,sigma2,time);
        Runge_Kutta_Df_Cont(b0,EF1,EF2,ArrayCont,0);
        Runge_Kutta_Ad_Cont(dt6,ArrayCont,1);
        Runge_Kutta_Ad_Cont(dt2,ArrayCont,2);
        
        //------ 2nd step Runge-Kutta
        time = it*dt+(dt2);
        EF1=pulse1.Gaussian(ncyc*sigma1,sigma1,time);
        EF2=pulse2.Gaussian(ncyc*sigma1+DELAY,sigma2,time);
       
        Runge_Kutta_Df_Cont(b2,EF1,EF2,ArrayCont,1);
        Runge_Kutta_Ac_Cont(dt3,ArrayCont);
        Runge_Kutta_Ad_Cont(dt2,ArrayCont,2);

        //------ 3rd step Runge-Kutta
        Runge_Kutta_Df_fixed_Cont(b2,EF1,EF2,ArrayCont,1);
        Runge_Kutta_Ac_Cont(dt3, ArrayCont);
        Runge_Kutta_Ad_Cont(dt,ArrayCont,2);

        //------ 4th step Runge-Kutta
        time=(it+1)*dt; // CHANGE
        EF1=pulse1.Gaussian(ncyc*sigma1,sigma1,time);
        EF2=pulse2.Gaussian(ncyc*sigma1+DELAY,sigma2,time);
        
        Runge_Kutta_Df_Cont(b2,EF1,EF2,ArrayCont,1);
        Runge_Kutta_Ad_Cont(dt6,ArrayCont,3);

        // ----- Calculates cont amplituds at the end of each step
        Calculate_Cont_Ampl(EF1, EF2, ArrayCont);
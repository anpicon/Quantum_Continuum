
void Separate_string(string& s, vector<string>& str, size_t& pos)
{
        s = s.substr(0,pos); //cout << s << endl;
        str.push_back("");
        size_t count = 0;
        size_t w=0;  
        while(w != pos)            
        {
            if( s[w] != ' ')
            {
                str[count] += s[w];
            }
            else
            {
              if( 0 != strcasecmp(str[count].c_str(), "") )
              {
                count++;
                str.push_back("");
              }
            }
            w++;
        }
        if(str.back()=="")
          str.pop_back();                   //     for(int i=0; i<str.size(); i++) cout << str[i].c_str() << endl;
        //if(str.size() > 2)
        //{ //cout << str.size();
        //   cout << str[str.size() -1] << endl;
        //  printf("ERROR: string \" %s \": I can't recognize all the parameters\n", s.c_str());
        //  exit(1);
        //}
    
}



bool blank(string& s)
{
  size_t pos = s.length();
  size_t w = 0;
  bool blank = true;
  while(w!=pos) 
  {
      if(s[w]!= ' ')
      {
           blank = false;
           break;
      }
      w++;
  }
  return blank;
}

void End_Section(ifstream& fp_input, string str)
{
  string s;
  while((getline(fp_input,s)))
  {
   // cout << s << endl;

    if(0 == strcasecmp(s.c_str(), "")){}
    else
    {
          
          if(!blank(s))
          {
              size_t pos = s.find("}");
              if(pos == -1)
              { 
                  printf("Cannot find end of %s. Correct your input.\n", str.c_str());
                  exit(1);
              }
              else
              {
                 // printf("%s section read.\n", str.c_str());
                  return;
              } 
          }
    } 
  }
}

void flaser(ifstream& fp_input, Laser& pulse1, double& sigma1, bool& gaussian1,double& DELAY,string& s,vector<string>& str, size_t pos)
{
    s="";
    while(s==""){getline(fp_input, s);}
    
    while(s.find("}") == -1 || 0 == strcasecmp( s.c_str(), ""))
    {
        str.clear();
        pos = s.length();
        Separate_string(s, str, pos);
        if( 0 ==strcasecmp(str[0].c_str(), "frequency") )
        {
            if(str.size()==2)
            {
                pulse1.set_freq(atof(str[1].c_str()));
            }
            else
            {
                if( 0 ==strcasecmp(str[2].c_str(), "ev") )
                {
                    pulse1.set_freq(energy_eV_au*atof(str[1].c_str()));
                }
                else
                {
                    printf("Cannot understand %s as units.", str[2].c_str());
                    exit(1);
                }
            }
        }
        
        else if( 0 ==strcasecmp(str[0].c_str(), "delay") )
        {
            if(str.size()==2)
            {
                DELAY=atof(str[1].c_str());
            }
            else
            {
                if( 0 ==strcasecmp(str[2].c_str(), "fs") )
                {
                    DELAY=time_fs_au*atof(str[1].c_str());
                }
                else if( 0 ==strcasecmp(str[2].c_str(), "au") )
                {
                    DELAY=atof(str[1].c_str());
                }
                else
                {
                    printf("Cannot understand %s as units.", str[2].c_str());
                    exit(1);
                }
            }
        }
        
        else if ( 0 == strcasecmp( str[0].c_str(),  "envelope") )
        {
            if(str.size()==2)
            {
                if( 0 == strcasecmp(str[1].c_str(), "gaussian") )
                  gaussian1 = true;
                else if(str[1] == "sin2")
                  gaussian1 = false;
                else
                {
                   printf("Allowed types for laser envelope: gaussian, sin2. Your type \" %s \" is not available.\n", str[1].c_str());
                   exit(1);
                }
            }
            else
            {
                printf("Problem reading envelope.");
                exit(1);
            }
        }
                     
        else if ( 0 == strcasecmp( str[0].c_str(),  "wavelength") )
        {
            if(str.size()==2)
            {
                pulse1.set_lambdanm(atof(str[1].c_str()));
            }
            else if( 0 == strcasecmp( str[2].c_str(), "nm") )
                pulse1.set_lambdanm(atof(str[1].c_str()));
            else if( 0 == strcasecmp( str[2].c_str(), "au") || 0 == strcasecmp( str[2].c_str(), "atomicunit") )
                pulse1.set_lambdanm(atof(str[1].c_str())*10.*space_au_A);
            else
            {
                printf("Cannot understand %s as units.", str[2].c_str());
                exit(1);
            }
        }
        else if(0 == strcasecmp( str[0].c_str(),  "intensity") )
        {
            if(str.size()==2)
            {
                  pulse1.set_IntensityWcm2(atof(str[1].c_str()));
            }
            else if(0 == strcasecmp( str[2].c_str(), "wcm2") )
                  pulse1.set_IntensityWcm2(atof(str[1].c_str()));
            else if(0 == strcasecmp( str[2].c_str(), "au") || 0 == strcasecmp( str[2].c_str(),  "atomicunit") )
                  pulse1.set_IntensityWcm2(atof(str[1].c_str())*intensity_au_Wcm2);
            else
            {
                  printf("Cannot understand %s as units.", str[2].c_str());
                  exit(1);
            }
        }
        else if(0 == strcasecmp( str[0].c_str(), "polarization") )
        {
            double value1 = atof(str[1].c_str());
            double value2 = atof(str[2].c_str()); //cout << value2 << endl;
            double value3 = atof(str[3].c_str()); //cout << value3 << endl;
            
            double value=sqrt(value1*value1 + value2*value2 + value3*value3);
            value1/=value; value2/=value; value3/=value;
            vec1d temp(3,0.); temp[0] = value1; temp[1] = value2; temp[2] = value3;
            //u1[0]=value1; u1[1]=value2; u1[2]=value3;
            pulse1.set_pol(temp);
        }
        else if(0 == strcasecmp( str[0].c_str(), "sigma") )
        {
            if(gaussian1 == false)
            {
                printf("You can't define sigma for sin2 wave. Change your input.\n");
                exit(1);
            }
            else
            {
                if(str.size()<3)
                {
                    sigma1 = atof(str[1].c_str());
                }
                else if(0 == strcasecmp( str[2].c_str(), "fs"))
                {//cout << "si " << endl;
                    sigma1 = atof(str[1].c_str());
                    sigma1 *= time_fs_au;
                }
                else if(0 == strcasecmp( str[2].c_str(), "au") || 0 == strcasecmp( str[2].c_str(), "atomicunit"))
                {
                    sigma1 = atof(str[1].c_str());
                }
            }
        }
        else if(0 == strcasecmp( str[0].c_str(), "cycles") )
        {
            if(gaussian1 == true)
            {
                printf("You can't define cycles for gaussian wave. Change your input.\n");
                exit(1);
            }
            sigma1 = atof(str[1].c_str());
        }
        else if(0 == strcasecmp( str[0].c_str(), "CEP") )
        {
            pulse1.set_CEP(atof(str[1].c_str()));
        }
        else if(!blank(str[0]) || (0 != strcasecmp(s.c_str(), "")))
        {
            printf("cannot recognize key %s", s.c_str());
            exit(1);
        }
        if(!(getline(fp_input,s)))
        {
            printf("Can't find the end of laser section.\n");
            exit(1);
        }
        

    }//end of while
}//end of laser block

void ftdse(ifstream& fp_input,double& dt, string& iInitialState, vec1s& path, int& NR,double& tf,string& s,vector<string>& str, size_t pos)
{
    s="";
    while(s==""){getline(fp_input, s);}
    
    while(s.find("}") == -1 || 0 == strcasecmp( s.c_str(), ""))
    {
        str.clear();
        pos =s.length();
        Separate_string(s, str, pos);
        if( 0 ==strcasecmp(str[0].c_str(), "dt") )
        {
            if(str.size() < 2)
            {
                dt = atof(str[1].c_str());
            }
            else if(0 ==strcasecmp(str[2].c_str(), "fs"))
            {
                dt = atof(str[1].c_str())*time_fs_au;
            }
            else if(0 ==strcasecmp(str[2].c_str(), "au") || (0 ==strcasecmp(str[0].c_str(), "atomicunit")))
            {
                dt = atof(str[1].c_str());
            }
            else
            {
                printf("Units %s aren't allowed for dt.\n", str[2].c_str());
                exit(1);
            }
        }//end dt
        else if( 0 ==strcasecmp(str[0].c_str(), "InitialState") )
        {
            if( 0 ==strcasecmp(str[1].c_str(), "calculate") )
            {
                iInitialState="Calculate";
            }
            else if( 0 ==strcasecmp(str[1].c_str(), "read") )
            {
                iInitialState="Read";
                if (str.size() == 3) {
                    path[3]=str[2].c_str();
                }
                else
                {
                    printf("Specify path to initial state.\n");
                    exit(1);
                }
            }
            else
            {
                printf("%s is not correct, only calculate or read is allowed.\n", str[1].c_str());
                exit(1);
            }
        } //end InitialState
        else if( 0 ==strcasecmp(str[0].c_str(), "NR") )
        {
            NR = atof(str[1].c_str());
        }
        if( 0 ==strcasecmp(str[0].c_str(), "tf") )
        {
            if(str.size() < 2)
            {
                tf = atof(str[1].c_str());
            }
            else if(0 ==strcasecmp(str[2].c_str(), "fs"))
            {
                tf = atof(str[1].c_str())*time_fs_au;
            }
            else if(0 ==strcasecmp(str[2].c_str(), "au") || (0 ==strcasecmp(str[0].c_str(), "atomicunit")))
            {
                tf = atof(str[1].c_str());
            }
            else
            {
                printf("Units %s aren't allowed for dt.\n", str[2].c_str());
                exit(1);
            }
        }//end dt
        
        if(!(getline(fp_input, s)))
        {
             printf("Can't find the end of tdse section.\n");
             exit(1);
        }
    }//end reading tdse
    
} //end of tdse block


void felectronicstates(ifstream& fp_input,int& NEi,vec1s& path,double& EnergyShift,string& s,vector<string>& str, size_t pos)
{
    s="";
    while(s==""){getline(fp_input, s);}
    
    while(s.find("}") == -1 || 0 == strcasecmp( s.c_str(), ""))
    {
        str.clear();
        pos =s.length();
        Separate_string(s, str, pos);
        if( 0 ==strcasecmp(str[0].c_str(), "NEi") )
        {
            NEi=atof(str[1].c_str());
        }
        else if( 0 ==strcasecmp(str[0].c_str(), "PESpath") )
        {
            path[0]=str[1];
        }
        else if( 0 ==strcasecmp(str[0].c_str(), "DECAYSpath") )
        {
            path[1]=str[1];
        }
        else if( 0 ==strcasecmp(str[0].c_str(), "DIPpath") )
        {
            path[2]=str[1];
        }
        else if( 0 ==strcasecmp(str[0].c_str(), "energyshift") )
        {
            if(str.size() < 3)
            {
                EnergyShift=atof(str[1].c_str());
            }
            else if( 0 ==strcasecmp(str[2].c_str(), "au") )
            {
                EnergyShift=atof(str[1].c_str());
            }
            else if( 0 ==strcasecmp(str[2].c_str(), "ev") )
            {
                EnergyShift=atof(str[1].c_str())*energy_eV_au;
            }
            else
            {
                printf("Units %s aren't allowed for EnergyShift.\n", str[2].c_str());
                exit(1);
            }
        }
        
        if(!(getline(fp_input, s)))
        {
             printf("Can't find the end of electronicstates section.\n");
             exit(1);
        }
        
    } //end reading electronicstates
}// end of electronicstates

void fmolecule(ifstream& fp_input,string& molecule,string& s,vector<string>& str, size_t pos)
{
    s="";
    while(s==""){getline(fp_input, s);}
    
    while(s.find("}") == -1 || 0 == strcasecmp( s.c_str(), ""))
    {
        str.clear();
        pos =s.length();
        Separate_string(s, str, pos);
        if( 0 ==strcasecmp(str[0].c_str(), "name") )
        {
            molecule=str[1];
        }
        
        if(!(getline(fp_input, s)))
        {
             printf("Can't find the end of molecule section.\n");
             exit(1);
        }
    }// end reading ground state
}// end of ground state block


void fobservables(ifstream& fp_input,vector<bool>& iObservables,string& s,vector<string>& str, size_t pos)
{
    s="";
    while(s==""){getline(fp_input, s);}
    
    while(s.find("}") == -1 || 0 == strcasecmp( s.c_str(), ""))
    {
        str.clear();
        pos =s.length();
        Separate_string(s, str, pos);
        if( 0 ==strcasecmp(str[0].c_str(), "wfd") )
        {
            iObservables[0]=true;
        }
        else if( 0 ==strcasecmp(str[0].c_str(), "amplitudes") )
        {
            iObservables[1]=true;
        }
        else if(0 ==strcasecmp(str[0].c_str(), "XPS") )
        {
            iObservables[2]=true;
        }
        
        if(!(getline(fp_input, s)))
        {
             printf("Can't find the end of observables section.\n");
             exit(1);
        }
    }// end reading observables
}// end of fobservables


void fcontinuum(ifstream& fp_input, vec1s& path,string& s,vector<string>& str, size_t pos,vec1C& contstate,int& ind)
{
    int NPump = 0;
    int NProbe = 0;
    s="";
    while(s==""){getline(fp_input, s);}
    
    while(s.find("}") == -1 || 0 == strcasecmp( s.c_str(), ""))
    {
        str.clear();
        pos =s.length();
        Separate_string(s, str, pos);
        if( 0 ==strcasecmp(str[0].c_str(), "NCoupledPump") )
        {
            NPump=atof(str[1].c_str());
        }
        else if( 0 ==strcasecmp(str[0].c_str(), "WhichStatPump") )
        {
            for(int i=0;i<NPump;i++){
                (contstate[ind].StatCouplPump).push_back(atof(str[1+i].c_str()));
            }
        }
        else if( 0 ==strcasecmp(str[0].c_str(), "NCoupledProbe") )
        {
            NProbe=atof(str[1].c_str());
        }
        else if( 0 ==strcasecmp(str[0].c_str(), "WhichStatProbe") )
        {
            for(int i=0;i<NProbe;i++){
                (contstate[ind].StatCouplProbe).push_back(atof(str[1+i].c_str()));
            }
        }
        else if( 0 ==strcasecmp(str[0].c_str(), "Emax") )
        {
            contstate[ind].Emax = atof(str[1].c_str());

            if(0 ==strcasecmp(str[2].c_str(), "eV"))
            {
                contstate[ind].Emax *= energy_eV_au;
            }
            else if(0 ==strcasecmp(str[2].c_str(), "au") || (0 ==strcasecmp(str[2].c_str(), "atomicunit")))
            {
                continue;
            }
            else
            {
                printf("Units %s aren't allowed for Emax.\n", str[2].c_str());
                exit(1);
            }
        }//end Emax
        else if( 0 ==strcasecmp(str[0].c_str(), "Emin") )
        {
            contstate[ind].Emin = atof(str[1].c_str());

            if(0 ==strcasecmp(str[2].c_str(), "eV"))
            {
                contstate[ind].Emin *= energy_eV_au;
            }
            else if(0 ==strcasecmp(str[2].c_str(), "au") || (0 ==strcasecmp(str[2].c_str(), "atomicunit")))
            {
                continue;
            }
            else
            {
                printf("Units %s aren't allowed for Emin.\n", str[2].c_str());
                exit(1);
            }
        }//end Emin
        else if( 0 ==strcasecmp(str[0].c_str(), "dE") )
        {
            contstate[ind].dE = atof(str[1].c_str());

            if(0 ==strcasecmp(str[2].c_str(), "eV"))
            {
                contstate[ind].dE *= energy_eV_au;
            }
            else if(0 ==strcasecmp(str[2].c_str(), "au") || (0 ==strcasecmp(str[2].c_str(), "atomicunit")))
            {
                continue;
            }
            else
            {
                printf("Units %s aren't allowed for dE.\n", str[2].c_str());
                exit(1);
            }
        }//end dE
       else if( 0 ==strcasecmp(str[0].c_str(), "BE") )
        {   
            int i;
            for(i=0;i<(NPump+NProbe);i++){
                (contstate[ind].BE).push_back(atof(str[i+1].c_str()));
                cout << i+1 << endl;
            }
            if(str.size() < NPump+NProbe)
            {
                continue;
            }
            else if(0 ==strcasecmp(str[i+1].c_str(), "eV"))
            {
                for(int i=0;i<(NPump+NProbe);i++) contstate[ind].BE[i] *= energy_eV_au;
            }
            else if(0 ==strcasecmp(str[i+1].c_str(), "au") || (0 ==strcasecmp(str[i+1].c_str(), "atomicunit")))
            {
                continue;
            }
            else
            {
                printf("Units %s aren't allowed for BE.\n", str[i+1].c_str());
                exit(1);
            }
        }//end E_BE
        else if( 0 ==strcasecmp(str[0].c_str(), "Lmax") )
        {
            contstate[ind].Lmax=atof(str[1].c_str());
        }
        else if( 0 ==strcasecmp(str[0].c_str(), "Mmax") )
        {
            contstate[ind].Mmax=atof(str[1].c_str());
        }
        else if(0 ==strcasecmp(str[0].c_str(), "PESpath")){
            (contstate[ind]).PESpath = str[1];
        }
        else if(0 ==strcasecmp(str[0].c_str(), "DIPpath")){
            for(int i=0;i<(NPump+NProbe);i++){
                (contstate[ind].DIPpath).push_back(str[i+1]);
                std::cout << "states: " << i << endl;
            }
        }
        else if(0 ==strcasecmp(str[0].c_str(), "Gamma")){
            contstate[ind].Gamma = atof(str[1].c_str());

            if(0 ==strcasecmp(str[2].c_str(), "eV"))
            {
                contstate[ind].Gamma *= energy_eV_au;
            }
            else if(0 ==strcasecmp(str[2].c_str(), "au") || (0 ==strcasecmp(str[2].c_str(), "atomicunit")))
            {
                continue;
            }
            else
            {
                printf("Units %s aren't allowed for Gamma.\n", str[2].c_str());
                exit(1);
            }
        }
        if(!(getline(fp_input, s)))
        {
             printf("Can't find the end of continuumstates section.\n");
             exit(1);
        }
    }//end reading tdse
    
} //end of continuum block

void Read_Input(ifstream& fp_input,double& dt, string& iInitialState, int& NR,double& tf, int& NEi,vec1s& path,Laser& pulse1, double& sigma1, bool& gaussian1,Laser& pulse2, double& sigma2, bool& gaussian2,double& DELAY,double& EnergyShift,string& molecule,vector<bool>& iObservables,vec1C& ArrayCont,int& NContStat)
//(ifstream& fp_input,Laser& pulse1, Laser& pulse2, Coord_B& u1, Coord_B& u2, bool& gaussian1, bool& gaussian2, double& sigma1, double& sigma2, double& DELAY,//lasers
//bool& iTightBinding, string& TBtype, double& dt,  //tdse options
//bool& iWFDs, bool& iCurrent, bool& iTAbs, bool& iTAbsK   //observables
{
    NContStat=-1;
    string s;

    printf("Reading input...");
    while(getline(fp_input,s))
    {
        if(s!="")
        {
            vector<string> str;
            size_t pos = s.find("{");

            if (pos == -1)
            {
               printf("Can't find the character { in the string \" %s \". \nPlease, review your input.\n", s.c_str());
               exit(1);
            }
            else  Separate_string(s, str, pos);
            
            if( 0 == strcasecmp(str[0].c_str(), "tdse"))
            {
                ftdse(fp_input,dt,iInitialState,path,NR,tf,s,str,pos);
            }
            else if( 0 == strcasecmp(str[0].c_str(), "electronicstates"))
            {
                felectronicstates(fp_input,NEi,path,EnergyShift,s,str,pos);
            }
            else if( 0 == strcasecmp(str[0].c_str(), "laserpump"))
            {
                flaser(fp_input,pulse1,sigma1,gaussian1,DELAY,s,str,pos);
            }
            else if( 0 == strcasecmp(str[0].c_str(), "laserprobe"))
            {
                flaser(fp_input,pulse2,sigma2,gaussian2,DELAY,s,str,pos);
            }
            else if( 0 == strcasecmp(str[0].c_str(), "molecule"))
            {
                fmolecule(fp_input,molecule,s,str,pos);
            }
            else if( 0 == strcasecmp(str[0].c_str(), "observables"))
            {
                fobservables(fp_input,iObservables,s,str,pos);
            }
            else if(0 == strcasecmp(str[0].c_str(), "continuumstate"))
            {
                NContStat++;
                ArrayCont.push_back(Continuum(NContStat));
                fcontinuum(fp_input,path,s,str,pos,ArrayCont,NContStat);
                cout << "Ncont stats : " << NContStat << endl;
            }
        }
    } //end while
    printf("  done \n");
} //end read_input
#include "simtrack.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//================================  
//
//  Global  constants 
// 
//================================

//  1) coordinate system:
//        right handed, x^ cross y^ =  s^, bending angle >0, bend particle to -x direction.
//  2) canonical variables: 
//        (x, px, y, py, -c*dt, pt); 
//        z= -cdt, time-of-flight times -c;
//        pt_= dE/(P_0 c)= (E- E_0)/ (P_0 c), P_0 and E_0 for reference particle
//  3) if ultra relativistic or H_expand== true:
//        (x, px, y, py, -c*dt, pt) ==> (x, px, y, py, -c*dt, delta)
//        and small angle or para-axial approximation is assumed.

enum ps_index { x_ = 0, px_ = 1, y_ = 2, py_ = 3,  z_ = 4, pt_= 5 }; 

double PI = 3.14159265358979323;
double speed_light = 299792458.;
double Cr= 8.84627632682e-5;     // unit: m*GeV^{-3}, constant for electron radiation
double Cf= (55./24./sqrt(3.) ) * 2.8179403e-15*1.054572e-34/9.1093822e-31;  // constant for electron radiation fluctuation
double tune_low = 0.660, tune_high = 0.70;  //default tune search window

double Fdrift1 = 1.e0 / (2.e0 * (2.e0 - pow(2.e0, 1.0/3.0 )));
double Fdrift2 = 0.5e0 - Fdrift1;
double Fkick1  = 2.e0 * Fdrift1;
double Fkick2  = 1.e0 - 2e0 * Fkick1;

//------for RHIC good enough
//double BLslice = 2.5 ;    // sbend magnets, useful for H_expand and  spin tracking
//double QLslice = 0.25;    // multi magnets:   for RHIC: 0.25m, for AGS: 0.01m,  for eRHIC: 0.01m
//double GLslice = 0.25;    // combined sbends: for RHIC: 0.25m, for AGS: 0.01m,  for eRHIC: 0.01m

//-----for erhic and  RCS:  test
//double BLslice = 1.0;      // sbend magnets, useful for H_expand and  spin tracking
//double QLslice = 0.2 ;     // multi magnets:   for RHIC: 0.25m, for AGS: 0.01m,  for eRHIC: 0.01m
//double GLslice = 0.2 ;     // combined sbends: for RHIC: 0.25m, for AGS: 0.01m,  for eRHIC: 0.01m

//-----for  erhic and ESR
double BLslice = 0.25;      // sbend magnets, useful for H_expand and  spin tracking
double QLslice = 0.05;     // multi magnets:   for RHIC: 0.25m, for AGS: 0.01m,  for eRHIC: 0.01m
double GLslice = 0.05;     // combined sbends: for RHIC: 0.25m, for AGS: 0.01m,  for eRHIC: 0.01m


//-----random number generator seeds
double  seed = 1234.  ;
gsl_rng *gsl_r1, *gsl_r2;

GlobalVariables GP;

void Print_GlobalVariables()
{
  cout<<"-------------------------------"<<endl;
  cout<<"GP.A       :     "<<GP.A<<endl;
  cout<<"GP.Q       :     "<<GP.Q<<endl;
  cout<<"GP.m0      :     "<<GP.m0<<endl;
  cout<<"GP.m       :     "<<GP.m<<endl; 
  cout<<"GP.q       :     "<<GP.q<<endl;  
  cout<<"GP.gamma   :     "<<GP.gamma<<endl;
  cout<<"GP.beta    :     "<<GP.beta<<endl;
  cout<<"GP.energy0 :     "<<GP.energy0<<endl;
  cout<<"GP.energy  :     "<<GP.energy<<endl;
  cout<<"GP.p       :     "<<GP.p<<endl;  
  cout<<"GP.cp      :     "<<GP.cp<<endl;
  cout<<"GP.brho    :     "<<GP.brho<<endl;
  cout<<"GP.G       :     "<<GP.G<<endl;
  cout<<"GP.Gr      :     "<<GP.Gr<<endl;
  cout<<"GP.step_deltap : "<<GP.step_deltap<<endl;
  cout<<"GP.step_deltaz : "<<GP.step_deltaz<<endl;  
  cout<<"GP.harm    :     "<<GP.harm<<endl;
  cout<<"GP.circumference:"<<GP.circumference<<endl;
  cout<<"GP.dgamma  :     "<<GP.dgamma<<endl;
  cout<<"GP.U0rad   :     "<<GP.U0rad<<endl; 
  cout<<"GP.gammat  :     "<<GP.gammat<<endl;
  cout<<"GP.H_expand :    "<<GP.H_expand<<endl;
  cout<<"GP.radiate :     "<<GP.radiate<<endl;
  //cout<<"GP.twiss_6d:     "<<GP.twiss_6d<<endl;
  cout<<"GP.quad_fringe:  "<<GP.quad_fringe<<endl;
  cout<<"-------------------------------"<<endl;
}

void Set_RefPartEnergy(double gamma, double E0, double Q, double A)
//Here: E0 is the TOTAL  REST energy for the reference particle as a whole, in unit of MeV,
//      Q in unit of proton charge,  
//      A is the number of nucleons.
//      for proton:   Q=1, A=1;
//      for electron, Q=-1, A=1
//      for Au79+,    Q=79, A=197
//Following parameters are all for the reference charged particle as a whole:
{
  GP.Q       = Q;
  GP.A       = A;                          // GP.A only used for longitudinal area, parameter etc.
  GP.m0      = E0 * 1.0e6 * 1.6021765e-19 / 299792458. / 299792458. ;  //   in Kg
  GP.m       = gamma * GP.m0;              //   in Kg
  GP.q       = Q * 1.6021765e-19;          //   in C
  GP.gamma   = gamma;
  GP.beta    = sqrt(1.0- 1.0/gamma/gamma);
  GP.energy0 = E0;                         // total rest energy for the whole reference particle, in MeV
  GP.energy  = gamma * GP.energy0;         // total energy for the whole refernce particle, in MeV
  GP.cp      = GP.energy * GP.beta;        // total cp for the whole reference particle  in MeV /c 
  GP.p       = GP.energy * GP.beta * 1.0e6 *  1.6021765e-19 / 299792458. ;   // in SI unit
  GP.brho    = 3.33564095198 * GP.beta * (GP.energy/1000)/ GP.Q ;   //  for the whole reference particle
}

void Set_Spin_G(double g)
//needed before spin tracking, GP.Gr is not used 
{
  GP.G  = g;
  GP.Gr = g * GP.gamma;   
}

void Cal_ParticleDelta(double pt, double &delta1, double &gamma1, double &beta1)
//  from pt to delta1, gamma1, beta1
{
  delta1= sqrt( 1.0 + 2*pt/GP.beta+pt*pt) -1.0;
  gamma1 =GP.gamma+sqrt(GP.gamma*GP.gamma-1)*pt; 
  beta1 = sqrt(1.0-1.0/gamma1/gamma1);
}

double DeltaToPt(double delta)
//  from delta to pt
{
  return (-2.0/GP.beta + sqrt( 4.0/GP.beta/GP.beta - 4.0*(1.0-(1+delta)*(1+delta)) ) )/2.; 
}

double PtToDelta(double pt)
//  from pt to delta
{
return sqrt( 1.0 + 2*pt/GP.beta+pt*pt) -1.0;
}

void  Cal_Momenta_Velocity(double x[6], double p[3], double v[3], double & gamma1, double & delta1, double & beta1)
//----applicable for a straight section, from x[]=(x,px,y,py, -cdt, pt )
{
  int i;
  double c =2.99792458e8;
  double p0, p1;

  delta1= sqrt( 1.0 + 2*x[5]/GP.beta+x[5]*x[5]) -1.0;
  gamma1= GP.gamma+sqrt(GP.gamma*GP.gamma-1)*x[5]; 
  beta1 = sqrt(1.0-1.0/gamma1/gamma1);
    
  p1=(1+ delta1) * GP.p;
  p[0]= x[1] * GP.p;
  p[1]= x[3] * GP.p;
  p[2]= sqrt(p1*p1 - p[0]*p[0] - p[1]*p[1] );
  for(i=0;i<3;i++) v[i]= p[i] / gamma1 / GP.m0;
}

void  Cal_Momenta_Velocity1(double x[6], double p[3], double v[3], double & gamma1, double & delta1, double & beta1)
//----applicable for a straight section, from x[]=(x,x',y,y', -cdt, pt )
{
  int i;
  double q= 1.0* 1.6021773e-19;
  double c =2.99792458e8;
  double ps, p0, p1;

  delta1= sqrt( 1.0 + 2*x[5]/GP.beta+x[5]*x[5]) -1.0;
  gamma1= GP.gamma+sqrt(GP.gamma*GP.gamma-1)*x[5]; 
  beta1 = sqrt(1.0-1.0/gamma1/gamma1);
  
  p1=(1+ delta1) * GP.p;
  ps = sqrt(1+x[1]*x[1]+x[3]*x[3]);
  p[0]= p1 * x[1] / ps;
  p[1]= p1 * x[3] / ps;
  p[2]= p1 * 1.0  / ps;
  for(i=0;i<3;i++) v[i]= p[i] / gamma1 / GP.m0;
}

double  gamma_to_velocity(double gamma)
{
  double beta;
  beta = sqrt( 1.0 - 1.0/gamma/gamma );
  return  2.99792458e8 * beta;
}

double  velocity_to_gamma(double v[3])
{
  double  beta1, v1;

  v1=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  beta1 = v1 / 2.99792458e8 ;
  return 1.0/sqrt(1-beta1*beta1);
}

//===========================================
//
//          supporting functions
//
//===========================================


void SplitString(const std::string& input, int& itemCount, std::vector<std::string>& items) {
    std::istringstream iss(input);
    std::string word;
    items.clear(); // Make sure it's empty before filling

    while (iss >> word) {
        items.push_back(word);
    }

    itemCount = static_cast<int>(items.size());
}

int  StringToInteger(const std::string& input ) {
    std::stringstream ss(input);
    int output;
    ss >> output;

    // Optional: Check if conversion was successful
    if (ss.fail()) {
        std::cerr << "Conversion failed: invalid integer format." << std::endl;
        output = 0; // Or handle error as needed
    }
    return output;
}

double  StringToDouble(const std::string& input ) {
    std::stringstream ss(input);
    double output;
    ss >> output;

    // Optional: Check if conversion was successful
    if (ss.fail()) {
        std::cerr << "Conversion failed: invalid integer format." << std::endl;
        output = 0; // Or handle error as needed
    }
    return output;
}

bool  StringInStringVector(const std::vector<std::string>& vec, const std::string& target) {
    return std::find(vec.begin(), vec.end(), target) != vec.end();
}

int fac( int i ){
  if( i==0 ) {
    return 1;}
  else if ( i==1) {
    return 1;
  }
  else if( i==2 ) {
    return 2;
  }
  else if( i==3 ) {
    return 6;
  }
  else if( i==4 ) {
    return 24;
  }
  else if( i==5 ) {
    return 120;
  }
  else if( i==6 ) {
    return 720;
  }
  else if( i==7 ) {
    return 5040;
  }
  else if( i==8 ) {
    return 40320;
  }
  else if( i==9 ) {
    return 362880;
  }
  else if( i==10) {
    return 3628800;
  }
  else if( i==11) {
    return 39916800;
  }
  else{
    cout<<"factorial input exceed 11. error exit!"<<endl;
    exit(1);
  }
  
}

double rnd(double & r)
{
  int m;
  double s,u,v,p;
  s=65536.0; u=2053.0; v=13849.0;
  m=int(r/s);  r=r-m*s;  r=u*r+v; 
  m=(int)(r/s); r=r-m*s; p=r/s;
  return(p);
}

double gaussian(double u,double g, double & r)
{ int i,m;
 double s,w,v,t;
 s=65536.0; w=2053.0; v=13849.0;
 t=0.0;
 for (i=1; i<=12; i++)
   { r=r*w+v; m=(int)(r/s);
   r=r-m*s; t=t+(r)/s;
   }
 t=u+g*(t-6.0);
 return(t);
}

double guassian_with_cut(double cut )
{
  double  a;
  cut= abs( cut ) ;
  do{
    a = gsl_ran_gaussian(gsl_r1, 1.0);
      } while (a > cut or a < -cut );
  return a;
}

void Init_GSL_Seeds()
{
  const gsl_rng_type * T1;
  gsl_rng_env_setup();
  T1 = gsl_rng_default;
  gsl_r1 = gsl_rng_alloc (T1);
  gsl_rng_set(gsl_r1, 9348.);

  const gsl_rng_type * T2;
  gsl_rng_env_setup();
  T2 = gsl_rng_default;
  gsl_r2 = gsl_rng_alloc (T2);
  gsl_rng_set(gsl_r2, 1324.);
}



double splint(double XA[], double YA[], int N, double X)
{
  int  KLO=0, KHI=N-1, K;
  double A, B, H;

  while( ( KHI-KLO) >  1 ){
    K=(KHI+KLO)/2;
    if(XA[K] > X ){
      KHI=K;
    }
    else{
      KLO=K;
    }
  }
  
  H=XA[KHI]-XA[KLO];
  if(H == 0. ) { cout<<" Bad XA input. XAs must be distinct."<<endl; exit(0); }
  A=(XA[KHI]-X)/H;
  B=(X-XA[KLO])/H;

  //cout<<XA[KLO] <<" "<<X<<" "<<XA[KHI] <<endl;
  //cout<<YA[KLO] <<" "<<A*YA[KLO]+B*YA[KHI]<<" "<<YA[KHI] <<endl;
  
  return A*YA[KLO]+B*YA[KHI];
}

double Interpolate(double x[], double y[], int n, double t)
{
int i,j,k,m;
    double z,s;
    z=0.0;
    if (n<1) return(z);
    if (n==1) { z=y[0]; return(z);}
    if (n==2)
      { z=(y[0]*(t-x[1])-y[1]*(t-x[0]))/(x[0]-x[1]);
        return(z);
      }
    if (t<=x[1]) { k=0; m=2;}
    else if (t>=x[n-2]) { k=n-3; m=n-1;}
    else
      { k=1; m=n;
        while (m-k!=1)
          { i=(k+m)/2;
            if (t<x[i-1]) m=i;
            else k=i;
          }
        k=k-1; m=m-1;
        if (fabs(t-x[k])<fabs(t-x[m])) k=k-1;
        else m=m+1;
      }
    z=0.0;
    for (i=k;i<=m;i++)
      { s=1.0;
        for (j=k;j<=m;j++)
          if (j!=i) s=s*(t-x[j])/(x[i]-x[j]);
        z=z+s*y[i];
      }
    return(z);
}

double Interpolate1(double x[], double y[], int n, double t)
{
  int  i;
  int  flag = 0;
  
  for(i=0; i< n; i++){
    if( ( t - x[i] ) * (t - x[i+1] ) <  0  ){
      flag =1; break;
    }
  }

  if(flag == 1 ) {
    return  y[i]  +  (  y[i+1] - y[i] ) * (t - x[i] )/  ( x[i+1] - x[i]  )  ;
  }
  else{
    cout<<"out of interpolation region. "<<endl;
    exit(0);
  }
  
}

void Generate_SectionMap_4D(double matrix[],   double betax1, double alfax1, double betax2, double alfax2, 
			                       double betay1, double alfay1, double betay2, double alfay2,
                                               double dphix,   double dphiy)
//   dphix, dphiy:  phase advances
{
  int i;
  
  for(i=0;i<16;i++) matrix[i]=0; 
  
  matrix[0*4+0]= sqrt(betax2/betax1) * (  cos(dphix)  + alfax1 * sin(dphix) ) ;  						      
  matrix[0*4+1]= sqrt(betax1*betax2) * sin(dphix);  									      
  matrix[1*4+0]=-( 1+ alfax1 * alfax2 ) * sin(dphix)/sqrt(betax1 *betax2) + (alfax1 - alfax2)*cos(dphix)/sqrt(betax1*betax2);    
  matrix[1*4+1]= sqrt(betax1/betax2) * ( cos(dphix) - alfax2*sin(dphix) );                                                      

  matrix[2*4+2]= sqrt(betay2/betay1) * (  cos(dphiy)  + alfay1 * sin(dphiy) ) ;  						      
  matrix[2*4+3]= sqrt(betay1*betay2) * sin(dphiy);  									      
  matrix[3*4+2]=-( 1+ alfay1 * alfay2 )* sin(dphiy)/sqrt(betay1 *betay2) + (alfay1 - alfay2)*cos(dphiy)/sqrt(betay1*betay2);    
  matrix[3*4+3]= sqrt(betay1/betay2) * ( cos(dphiy) - alfay2*sin(dphiy) );                                                      

}

void Generate_SectionMap_6D(double matrix[],   double betax1, double alfax1, double betax2, double alfax2, 
			                       double betay1, double alfay1, double betay2, double alfay2,
                                               double dphix,   double dphiy)
{
  int i;
  
  for(i=0;i<36;i++) matrix[i]=0; 
  
  matrix[0*6+0]= sqrt(betax2/betax1) * (  cos(dphix)  + alfax1 * sin(dphix) ) ;  						      
  matrix[0*6+1]= sqrt(betax1*betax2) * sin(dphix);  									      
  matrix[1*6+0]=-( 1+ alfax1 * alfax2 ) * sin(dphix)/sqrt(betax1 *betax2) + (alfax1 - alfax2)*cos(dphix)/sqrt(betax1*betax2);    
  matrix[1*6+1]= sqrt(betax1/betax2) * ( cos(dphix) - alfax2*sin(dphix) );                                                      

  matrix[2*6+2]= sqrt(betay2/betay1) * (  cos(dphiy)  + alfay1 * sin(dphiy) ) ;  						      
  matrix[2*6+3]= sqrt(betay1*betay2) * sin(dphiy);  									      
  matrix[3*6+2]=-( 1+ alfay1 * alfay2 )* sin(dphiy)/sqrt(betay1 *betay2) + (alfay1 - alfay2)*cos(dphiy)/sqrt(betay1*betay2);    
  matrix[3*6+3]= sqrt(betay1/betay2) * ( cos(dphiy) - alfay2*sin(dphiy) );                                                      

  matrix[4*6+4]= 1.0; 
  matrix[5*6+5]= 1.0;

}

void Generate_OneTurnMap_4D(double matrix[],  double mux, double betax, double alfax, 
			                      double muy, double betay, double alfay)
// mux, muy:  working point (  mod: 2PI) 
{
  int i;

  for(i=0;i<16;i++) matrix[i]=0;
  
  matrix[0*4+0]= cos(mux*2*PI) + alfax * sin(mux*2*PI);
  matrix[0*4+1]=                 betax * sin(mux*2*PI);
  matrix[1*4+0]= -(1+ alfax*alfax)*sin(mux*2*PI)/betax;
  matrix[1*4+1]= cos(mux*2*PI) - alfax * sin(mux*2*PI);

  matrix[2*4+2]= cos(muy*2*PI) + alfay * sin(muy*2*PI);
  matrix[2*4+3]=                 betay * sin(muy*2*PI);
  matrix[3*4+2]= -(1+ alfay*alfay)*sin(muy*2*PI)/betay;
  matrix[3*4+3]= cos(muy*2*PI) - alfay * sin(muy*2*PI);

}

void Generate_OneTurnMap_6D(double m66[],  double mux, double betax, double alfax, 
                                           double muy, double betay, double alfay, 
                                           double muz, double betaz, double alfaz )
{
  int i;

  for(i=0;i<36;i++) m66[i]=0;
  
  m66[0*6+0]= cos(mux*2*PI) + alfax * sin(mux*2*PI);
  m66[0*6+1]=                 betax * sin(mux*2*PI);
  m66[1*6+0]= -(1+ alfax*alfax)*sin(mux*2*PI)/betax;
  m66[1*6+1]= cos(mux*2*PI) - alfax * sin(mux*2*PI);

  m66[2*6+2]= cos(muy*2*PI) + alfay * sin(muy*2*PI);
  m66[2*6+3]=                 betay * sin(muy*2*PI);
  m66[3*6+2]= -(1+ alfay*alfay)*sin(muy*2*PI)/betay;
  m66[3*6+3]= cos(muy*2*PI) - alfay * sin(muy*2*PI);

  m66[4*6+4]= cos(muz*2*PI) ;
  m66[4*6+5]= betaz * sin(muz*2*PI);
  m66[5*6+4]=-sin(muz*2*PI)/betaz;
  m66[5*6+5]= cos(muz*2*PI) ;

}

void Matrix_Transfer(double matrix[36], double x[6])
{
  int j,k;
  double x0[6];

  for(j=0;j<6;j++) x0[j]=x[j];
  
  for(j=0;j<6;j++) {
    x[j]=0.;
    for(k=0;k<6;k++) x[j] += matrix[j*6+k] * x0[k];
  }
}

void Matrix_Transfer(double matrix[], double x[], int  Np)
{
  int i,j,k;
  double x0[6];

  for(i=0;i<Np;i++) {
    for(j=0;j<6;j++) x0[j]=x[6*i + j];
    for(j=0;j<6;j++) {
      x[6*i+j]=0.;
      for(k=0;k<6;k++) x[6*i+j] += matrix[j*6+k] * x0[k];
    }
  }

}

void Transfer_Ring(int Npart, double x[], double mux, double chromx, double betax, double alfax, double muy, double chromy, double betay, double alfay, double muz, double betaz, double alfaz)
{
  int j;
  int  k;
  double m66[36];
  double xtrack[6];
  double qx, qy;  
  
  for(j=0;j<Npart;j++){
    for(k=0;k<6;k++) xtrack[k]=x[j*6+k];      
    qx=mux  + chromx* xtrack[5];   
    qy=muy  + chromy* xtrack[5];  
    Generate_OneTurnMap_6D(m66, qx, betax, alfax, qy, betay, alfay, muz, betaz, alfaz);
    Matrix_Transfer(m66, xtrack);
    for(k=0;k<6;k++) x[j*6+k]=xtrack[k];      
  } 
}

/*---

void Transfer_Ring2(int Npart, double x[], double mux, double chromx, double betax, double alfax, double muy, double chromy, double betay, double alfay, double muz, double betaz, double alfaz)
{
  int j;
  double sintmp_z=sin(muz*2*PI), costmp_z=cos(muz*2*PI);
  double m11z= costmp_z,  m12z= betaz * sintmp_z, m21z=-sintmp_z/betaz,  m22z= costmp_z;
  
  #pragma omp parallel for
  for(j=0;j<Npart;j++){
    int  k;
    double m11, m12, m21, m22;
    double sintmp, costmp;
    double xtrack[6];
    double qx, qy;
    
    for(k=0;k<6;k++) xtrack[k]=x[j*6+k];      
    qx=mux  + chromx* xtrack[5];   
    qy=muy  + chromy* xtrack[5];

    sintmp= sin(qx*2*PI);
    costmp= cos(qx*2*PI);
    m11= costmp + alfax * sintmp;
    m12= betax * sintmp;
    m21= -(1+ alfax*alfax)*sintmp/betax;
    m22= costmp - alfax * sintmp;
    x[j*6+0] = m11 * xtrack[0] + m12 * xtrack[1]; 
    x[j*6+1] = m21 * xtrack[0] + m22 * xtrack[1];      
    
    sintmp= sin(qy*2*PI);
    costmp= cos(qy*2*PI);   
    m11= costmp + alfay * sintmp;
    m12= betay * sintmp;
    m21= -(1+ alfay*alfay)*sintmp/betay;
    m22= costmp - alfay * sintmp;
    x[j*6+2] = m11 * xtrack[2] + m12 * xtrack[3]; 
    x[j*6+3] = m21 * xtrack[2] + m22 * xtrack[3];

    x[j*6+4] = m11z * xtrack[4] + m12z * xtrack[5]; 
    x[j*6+5] = m21z * xtrack[4] + m22z * xtrack[5];
  }
  
}

void Transfer_Ring2c(int Npart, double x[], double mux, double chromx, double betax, double alfax, double muy, double chromy, double betay, double alfay, double muz, double betaz, double alfaz)
{
  int j;
  double sintmp_z=sin(muz*2*PI), costmp_z=cos(muz*2*PI);
  double m11z= costmp_z,  m12z= betaz * sintmp_z, m21z=-sintmp_z/betaz,  m22z= costmp_z;
  
  #pragma omp parallel for
  for(j=0;j<Npart;j++){
    int  k;
    double m11, m12, m21, m22;
    double sintmp, costmp;
    double tempx0, tempy0;
    double xtrack[6];
    double qx, qy;

    //-----one-turn map  for on-momentum particle
    for(k=0;k<6;k++) xtrack[k]=x[j*6+k];      
    qx=mux;
    qy=muy;

    sintmp= sin(qx*2*PI);
    costmp= cos(qx*2*PI);
    m11= costmp + alfax * sintmp;
    m12= betax * sintmp;
    m21= -(1+ alfax*alfax)*sintmp/betax;
    m22= costmp - alfax * sintmp;
    x[j*6+0] = m11 * xtrack[0] + m12 * xtrack[1]; 
    x[j*6+1] = m21 * xtrack[0] + m22 * xtrack[1];      
    
    sintmp= sin(qy*2*PI);
    costmp= cos(qy*2*PI);   
    m11= costmp + alfay * sintmp;
    m12= betay * sintmp;
    m21= -(1+ alfay*alfay)*sintmp/betay;
    m22= costmp - alfay * sintmp;
    x[j*6+2] = m11 * xtrack[2] + m12 * xtrack[3]; 
    x[j*6+3] = m21 * xtrack[2] + m22 * xtrack[3];

    x[j*6+4] = m11z * xtrack[4] + m12z * xtrack[5]; 
    x[j*6+5] = m21z * xtrack[4] + m22z * xtrack[5];

    //----apply a  thin kick for symplectic chrom implementation

    x[j*6+1]  =   x[j*6+1]  - 4 * PI * chromx * x[j*6+5] * x[j*6+0] / betax ;
    x[j*6+3]  =   x[j*6+3]  - 4 * PI * chromy * x[j*6+5] * x[j*6+2] / betay ;
    x[j*6+4]  =   x[j*6+4]   +  ( 2 * PI * chromx * x[j*6+0] * x[j*6+0] / betax
			        + 2 * PI * chromy * x[j*6+2] * x[j*6+2] / betay ) ;
  }
  
}

---*/

void Cal_Mean(int Npart, double PartDist[], double mean[6] )
{
  int i;
  double sum[6];
  
  for(i=0;i<6;i++) {
    sum[0] = 0;
    sum[1] = 0;
    sum[2] = 0;
    sum[3] = 0; 
    sum[4] = 0;
    sum[5] = 0;     
  }

  for(i=0;i<Npart;i++) {
    sum[0] = sum[0] + PartDist[i*6+0] ; 
    sum[1] = sum[1] + PartDist[i*6+1] ;
    sum[2] = sum[2] + PartDist[i*6+2] ;
    sum[3] = sum[3] + PartDist[i*6+3] ; 
    sum[4] = sum[4] + PartDist[i*6+4] ;
    sum[5] = sum[5] + PartDist[i*6+5] ;     
  }
  
  for(i=0;i<6;i++) {
    sum[0] = sum[0]/Npart;
    sum[1] = sum[1]/Npart;
    sum[2] = sum[2]/Npart;
    sum[3] = sum[3]/Npart; 
    sum[4] = sum[4]/Npart;
    sum[5] = sum[5]/Npart;     
  } 
  
}

void Cal_Max_Min_Mean_RMS( double x[], int n, double & max, double & min,  double & mean,  double & rms)
{
  int i;

  max=x[0];
  for(i=0;i<n;i++){
    if(x[i]>max) max = x[i];
  }

  min=x[0];
  for(i=0;i<n;i++){
    if(x[i]<min) min = x[i];
  }
  
  mean=0;
  for(i=0;i<n;i++){
    mean +=  x[i];
  }
  mean= mean*1.0 / n;
  
  rms=0;
  for(i=0;i<n;i++){
    rms +=  (x[i] - mean ) * (x[i] - mean );
  }
  rms= sqrt(rms*1.0 / (n-1));

}

void Cal_Mean_RMS( double x[], int n, double & mean, double & rms )
{
  int i;  
  mean=0;
  for(i=0;i<n;i++){
    mean +=  x[i];
  }
  mean= mean*1.0 / n;
  
  rms=0;
  for(i=0;i<n;i++){
    rms +=  (x[i] - mean ) * (x[i] - mean );
  }
  rms= sqrt(rms*1.0 / (n-1));
}

void Cal_Mean_RMS(int Npart, double PartDist[], double mean[2], double rms[2])
{
  int i;
  double sumx=0, sumy=0, sumx2=0, sumy2=0;

  for(i=0;i<Npart;i++) {
    int  j= i*6 ;
    double tempx = PartDist[j] ;
    double tempy = PartDist[j+2] ;
    sumx = sumx + tempx;
    sumy = sumy + tempy;
    sumx2 = sumx2 + tempx*tempx;
    sumy2 = sumy2 + tempy*tempy;
  }
  
  mean[0]= sumx / Npart;
  mean[1]= sumy / Npart;
  rms[0]= sqrt( sumx2/Npart-mean[0]*mean[0]);
  rms[1]= sqrt( sumy2/Npart-mean[1]*mean[1]);
}

void Cal_Mean_RMS_6D(int Npart,  double PartDist[],  int id[],  double mean[], double rms[])
{
  int i;
  double sumx2=0, sumy2=0, sumz2=0, sumx=0, sumy=0, sumz=0;
  double tempx, tempy, tempz;

  for(i=0;i<Npart;i++) {
    tempx = PartDist[id[i]*6+0];
    tempy = PartDist[id[i]*6+2];
    tempz = PartDist[id[i]*6+4];
    sumx = sumx + tempx;
    sumy = sumy + tempy;
    sumz = sumz + tempz;
    sumx2 = sumx2 + tempx*tempx;
    sumy2 = sumy2 + tempy*tempy;
    sumz2 = sumz2 + tempz*tempz;
  }
  mean[0]= sumx / Npart;
  mean[1]= sumy / Npart;
  mean[2]= sumz / Npart;
  rms[0]= sqrt( sumx2-sumx*sumx)/ Npart;
  rms[1]= sqrt( sumy2-sumy*sumy)/ Npart;
  rms[2]= sqrt( sumz2-sumz*sumz)/ Npart;
}


void Cal_Twiss_Emit( double x[], double px[], double n, double & beta, double & alpha, double & gamma, double & emit)
{
  int i;
  double sum;
  double xmean, pxmean;
  double sigma_xx,sigma_xpx,sigma_pxpx;

  sum=0.;
  for(i=0;i<n;i++){
    sum=sum+x[i];
  }
  xmean=sum/n;
  sum=0.;
  for(i=0;i<n;i++){
    sum=sum+(x[i]-xmean)*(x[i]-xmean);
  }
  sigma_xx=sum/n;

  sum=0.;
  for(i=0;i<n;i++){
    sum=sum+px[i];
  }
  pxmean=sum/n;
  sum=0.;
  for(i=0;i<n;i++){
    sum=sum+(px[i]-pxmean)*(px[i]-pxmean);
  }
  sigma_pxpx=sum/n;
  
  sum=0.;
  for(i=0;i<n;i++){
    sum=sum+(x[i]-xmean)*(px[i]-pxmean);
  }
  sigma_xpx=sum/n;

  //---here emit is for  a Gaussian distribution
  emit=sqrt(sigma_xx * sigma_pxpx-sigma_xpx * sigma_xpx);
  beta=sigma_xx/emit;
  gamma=sigma_pxpx/emit;
  alpha=-sigma_xpx/emit;

  for(i=0;i<1024;i++){
    x[i] =x[i]-xmean;
    px[i]=px[i]-pxmean;
  }

  //---here emit is for a single particle, emit=2J
  if(false){
    sum=0;
    for(i=0;i<n;i++){
      sum = sum+ gamma * x[i]*x[i] + 2.0* alpha * x[i]* px[i]+ beta * px[i]*px[i];
    }
    emit=sum/n;
  }
} 

void Cal_Emit(int Npart, double PartDist[], double emit[2], double beta[2], double alfa[2], double gama[2])
//  gamma  for energy  and gama for Twiss
{
 
  int i;
  double sum;
  double *z = new double [Npart],  *pz = new double [Npart];
  double zmean, pzmean;
  double sigma_zz,sigma_zpz,sigma_pzpz;

  //----horizontal plane

  for(i=0;i<Npart;i++){
    z[i] = PartDist[i*6+0];
    pz[i]= PartDist[i*6+1];
  }

  sum=0.;
  for(i=0;i<Npart;i++){
    sum=sum+z[i];
  }
  zmean=sum/Npart;
  sum=0.;
  for(i=0;i<Npart;i++){
    sum=sum+(z[i]-zmean)*(z[i]-zmean);
  }
  sigma_zz=sum/Npart;

  sum=0.;
  for(i=0;i<Npart;i++){
    sum=sum+pz[i];
  }
  pzmean=sum/Npart;
  sum=0.;
  for(i=0;i<Npart;i++){
    sum=sum+(pz[i]-pzmean)*(pz[i]-pzmean);
  }
  sigma_pzpz=sum/Npart;
  
  sum=0.;
  for(i=0;i<Npart;i++){
    sum=sum+(z[i]-zmean)*(pz[i]-pzmean);
  }
  sigma_zpz=sum/Npart;

  emit[0]= sqrt(sigma_zz * sigma_pzpz-sigma_zpz * sigma_zpz);
  beta[0]= sigma_zz/emit[0];
  alfa[0]=-sigma_zpz/emit[0];
  gama[0]=sigma_pzpz/emit[0];

 //----vertical plane

  for(i=0;i<Npart;i++){
    z[i] = PartDist[i*6+2];
    pz[i]= PartDist[i*6+3];
  }

  sum=0.;
  for(i=0;i<Npart;i++){
    sum=sum+z[i];
  }
  zmean=sum/Npart;
  sum=0.;
  for(i=0;i<Npart;i++){
    sum=sum+(z[i]-zmean)*(z[i]-zmean);
  }
  sigma_zz=sum/Npart;

  sum=0.;
  for(i=0;i<Npart;i++){
    sum=sum+pz[i];
  }
  pzmean=sum/Npart;
  sum=0.;
  for(i=0;i<Npart;i++){
    sum=sum+(pz[i]-pzmean)*(pz[i]-pzmean);
  }
  sigma_pzpz=sum/Npart;
  
  sum=0.;
  for(i=0;i<Npart;i++){
    sum=sum+(z[i]-zmean)*(pz[i]-pzmean);
  }
  sigma_zpz=sum/Npart;

  emit[1]=sqrt(sigma_zz * sigma_pzpz-sigma_zpz * sigma_zpz);
  beta[1]=sigma_zz/emit[1];
  alfa[1]=-sigma_zpz/emit[1];
  gama[1]=sigma_pzpz/emit[1];
  
}

void Cal_Emit_3D(int Npart, double PartDist[], double emit[], double beta[], double alfa[], double gama[])
//  gamma  for energy  and gama for Twiss
{
 
  int i;
  double sum;
  double *z = new double [Npart],  *pz = new double [Npart];
  double zmean, pzmean;
  double sigma_zz,sigma_zpz,sigma_pzpz;

  //----horizontal plane

  for(i=0;i<Npart;i++){
    z[i] = PartDist[i*6+0];
    pz[i]= PartDist[i*6+1];
  }

  sum=0.;
  for(i=0;i<Npart;i++){
    sum=sum+z[i];
  }
  zmean=sum/Npart;
  sum=0.;
  for(i=0;i<Npart;i++){
    sum=sum+(z[i]-zmean)*(z[i]-zmean);
  }
  sigma_zz=sum/Npart;

  sum=0.;
  for(i=0;i<Npart;i++){
    sum=sum+pz[i];
  }
  pzmean=sum/Npart;
  sum=0.;
  for(i=0;i<Npart;i++){
    sum=sum+(pz[i]-pzmean)*(pz[i]-pzmean);
  }
  sigma_pzpz=sum/Npart;
  
  sum=0.;
  for(i=0;i<Npart;i++){
    sum=sum+(z[i]-zmean)*(pz[i]-pzmean);
  }
  sigma_zpz=sum/Npart;

  emit[0]= sqrt(sigma_zz * sigma_pzpz-sigma_zpz * sigma_zpz);
  beta[0]= sigma_zz/emit[0];
  alfa[0]=-sigma_zpz/emit[0];
  gama[0]=sigma_pzpz/emit[0];

 //----vertical plane

  for(i=0;i<Npart;i++){
    z[i] = PartDist[i*6+2];
    pz[i]= PartDist[i*6+3];
  }

  sum=0.;
  for(i=0;i<Npart;i++){
    sum=sum+z[i];
  }
  zmean=sum/Npart;
  sum=0.;
  for(i=0;i<Npart;i++){
    sum=sum+(z[i]-zmean)*(z[i]-zmean);
  }
  sigma_zz=sum/Npart;

  sum=0.;
  for(i=0;i<Npart;i++){
    sum=sum+pz[i];
  }
  pzmean=sum/Npart;
  sum=0.;
  for(i=0;i<Npart;i++){
    sum=sum+(pz[i]-pzmean)*(pz[i]-pzmean);
  }
  sigma_pzpz=sum/Npart;
  
  sum=0.;
  for(i=0;i<Npart;i++){
    sum=sum+(z[i]-zmean)*(pz[i]-pzmean);
  }
  sigma_zpz=sum/Npart;

  emit[1]=sqrt(sigma_zz * sigma_pzpz-sigma_zpz * sigma_zpz);
  beta[1]=sigma_zz/emit[1];
  alfa[1]=-sigma_zpz/emit[1];
  gama[1]=sigma_pzpz/emit[1];

 //----longitudinal  plane

  for(i=0;i<Npart;i++){
    z[i] = PartDist[i*6+4];
    pz[i]= PartDist[i*6+5];
  }

  sum=0.;
  for(i=0;i<Npart;i++){
    sum=sum+z[i];
  }
  zmean=sum/Npart;
  sum=0.;
  for(i=0;i<Npart;i++){
    sum=sum+(z[i]-zmean)*(z[i]-zmean);
  }
  sigma_zz=sum/Npart;

  sum=0.;
  for(i=0;i<Npart;i++){
    sum=sum+pz[i];
  }
  pzmean=sum/Npart;
  sum=0.;
  for(i=0;i<Npart;i++){
    sum=sum+(pz[i]-pzmean)*(pz[i]-pzmean);
  }
  sigma_pzpz=sum/Npart;
  
  sum=0.;
  for(i=0;i<Npart;i++){
    sum=sum+(z[i]-zmean)*(pz[i]-pzmean);
  }
  sigma_zpz=sum/Npart;

  emit[2]=sqrt(sigma_zz * sigma_pzpz-sigma_zpz * sigma_zpz);
  beta[2]=sigma_zz/emit[2];
  alfa[2]=-sigma_zpz/emit[2];
  gama[2]=sigma_pzpz/emit[2];  
  
}

void Cal_Num_Bin( double x[], int m, double y[], int n, double xmax, double xmin)
{
  int i,j;
  double step = (xmax-xmin)/n;

  for(i=0;i<n;i++) y[i]=0;
 
  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      if( x[i] > xmin + j*step and x[i] < xmin + (j+1)*step ) y[j] ++;  
    }
  }
}

void Cal_Num_Bin_Weight( double x[], int m, double weight[], double y[], int n, double xmax, double xmin)
{
  int i,j;
  double step = (xmax-xmin)/n;

  for(i=0;i<n;i++) y[i]=0;
 
  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      if( x[i] > xmin + j*step and x[i] < xmin + (j+1)*step ) y[j] = y[j] + weight[i];  
    }
  }
}



void  Cal_Max_6D(int Mp2, double x2[],  double maxi[])
{
  int  i;
  double  temp;

  
  temp=0;  
  for(i=0;i<Mp2;i++){
    if(abs( x2[i*6+0] )  > temp ) temp= abs( x2[i*6+0] );
  }
  maxi[0]=temp;

  temp=0;  
  for(i=0;i<Mp2;i++){
    if( abs( x2[i*6+2] )  > temp ) temp= abs( x2[i*6+2] );
  }
  maxi[1]=temp;  

  temp=0;  
  for(i=0;i<Mp2;i++){
    if(abs( x2[i*6+4])  > temp ) temp= abs( x2[i*6+4] );
  }
  maxi[2]=temp;

  //cout<< maxi[0]<<"   "<< maxi[1]<<"  "<< maxi[2]<<endl;

}

void llsq( int n, double x[], double y[], double & a, double  & b )
{
  double bot;
  int i;
  double top;
  double xbar;
  double ybar;
/*
  Special case.
*/
  if ( n == 1 )
  {
    a = 0.0;
    b = y[0];
    return;
  }
/*
  Average X and Y.
*/
  xbar = 0.0;
  ybar = 0.0;
  for ( i = 0; i < n; i++ )
  {
    xbar = xbar + x[i];
    ybar = ybar + y[i];
  }
  xbar = xbar / ( double ) n;
  ybar = ybar / ( double ) n;
/*
  Compute Beta.
*/
  top = 0.0;
  bot = 0.0;
  for ( i = 0; i < n; i++ )
  {
    top = top + ( x[i] - xbar ) * ( y[i] - ybar );
    bot = bot + ( x[i] - xbar ) * ( x[i] - xbar );
  }
  
  a = top / bot;
  b = ybar - a * xbar;

  return;
}

void llsq_err( int n, double x[], double y[], double & a, double  & b ,  double  & ar,  double  & br, double & chi2, double  & r2)
{
  int i;
  double temp, ssrot;
  double S, Sx, Sy, Sxx, Syy, Sxy, Delta;

  S=0; Sx=0; Sy=0; Sxx=0; Syy=0;Sxy=0; 
  for(i=0;i<n; i++){
    S= S+1;
    Sx=Sx + x[i];
    Sy=Sy + y[i];
    Sxx = Sxx  + x[i]*x[i];
    Syy = Syy  + y[i]*y[i];    
    Sxy = Sxy  + x[i]*y[i];
  }

  Delta = S * Sxx - Sx * Sx;
  a = (S * Sxy  - Sx *Sy ) / Delta; 
  b = (Sxx * Sy - Sx *Sxy ) / Delta;

  chi2=0;
  ssrot = 0;
  for(i=0;i<n;i++){
    temp=  a * x[i] + b  - y[i] ;
    chi2 = chi2 + temp*temp;
    temp=  y[i] - Sy/n ;
    ssrot= ssrot + temp*temp;
  }

  ar =  sqrt( S / Delta );
  br =  sqrt( Sxx / Delta );
  
  r2 =  1 - chi2 / ssrot  ;

  temp = sqrt( chi2 / (n-2) );
  ar= ar * temp;
  br= br * temp;
    
  return;
}

void pfit(double x[],double y[], int n, double a[], int m)
{ 
  int m1,i,j,l,ii,k,im,ix[21];
  double h[21],ha,hh,y1,y2,h1,h2,d,hm;
  for (i=0; i<=m; i++) a[i]=0.0;
  if (m>=n) m=n-1;
  if (m>=20) m=19;
  m1=m+1;
  ha=0.0;
  ix[0]=0; ix[m]=n-1;
  l=(n-1)/m; j=l;
  for (i=1; i<=m-1; i++)
    { ix[i]=j; j=j+l;}
  while (1==1)
    { hh=1.0;
    for (i=0; i<=m; i++)
      { a[i]=y[ix[i]]; h[i]=-hh; hh=-hh;}
    for (j=1; j<=m; j++)
      { ii=m1; y2=a[ii-1]; h2=h[ii-1];
      for (i=j; i<=m; i++)
	{ d=x[ix[ii-1]]-x[ix[m1-i-1]];
	y1=a[m-i+j-1];
	h1=h[m-i+j-1];
	a[ii-1]=(y2-y1)/d;
	h[ii-1]=(h2-h1)/d;
	ii=m-i+j; y2=y1; h2=h1;
	}
      }
    hh=-a[m]/h[m];
    for (i=0; i<=m; i++)
      a[i]=a[i]+h[i]*hh;
    for (j=1; j<=m-1; j++)
      { ii=m-j; d=x[ix[ii-1]];
      y2=a[ii-1];
      for (k=m1-j; k<=m; k++)
	{ y1=a[k-1]; a[ii-1]=y2-d*y1;
	y2=y1; ii=k;
	}
      }
    hm=fabs(hh);
    if (hm<=ha) { a[m]=-hm; return;}
    a[m]=hm; ha=hm; im=ix[0]; h1=hh;
    j=0;
    for (i=0; i<=n-1; i++)
      { if (i==ix[j])
	{ if (j<m) j=j+1;}
      else
	{ h2=a[m-1];
	for (k=m-2; k>=0; k--)
	  h2=h2*x[i]+a[k];
	h2=h2-y[i];
	if (fabs(h2)>hm)
	  { hm=fabs(h2); h1=h2; im=i;}
	}
      }
    if (im==ix[0]) return;
    i=0;l=1;
    while (l==1)
      { l=0;
      if (im>=ix[i])
	{ i=i+1;
	if (i<=m) l=1;
	}
      }
    if (i>m) i=m;
    if (i==(i/2)*2) h2=-hh;
    else h2=hh;
    if (h1*h2>=0.0) ix[i]=im;
    else
      { if (im<ix[0])
	{ for (j=m-1; j>=0; j--)
	  ix[j+1]=ix[j];
	ix[0]=im;
	}
      else
	{ if (im>ix[m])
	  { for (j=1; j<=m; j++)
	    ix[j-1]=ix[j];
	  ix[m]=im;
	  }
	else ix[i-1]=im;
	}
      }
    }
}

void  Linear_Fit(double x[],double y[], int n, double & term0, double & term1, double & term2, double & term3)
{
  double coeff[8];
  pfit(x, y,n,coeff,7);
  term0  = coeff[0];
  term1  = coeff[1];
  term2  = coeff[2];
  term3  = coeff[3]; 
}

void Generate_Eye_Matrix( double eye[], int m )
{
  int i;
  for(i=0;i<m*m;i++) eye[i]=0.;
  for(i=0;i<m;i++)   eye[i*m+i]=1.;
}

void Print_Matrix(double  mat[], int m, int n) // m lines  * n columns
{
  int i, j;
  cout<<" ...The matrix is : "<<endl;
  for (i=0;i<m;i++) {
    for(j=0;j<n;j++) cout<<setw(12)<<mat[i*n+j]<<"  ";
    cout<<endl;
  }
  
}

void mat_vect_mult(double mat[], int m, double xin[], double xout[])
// matrix (m,m), xin(m,1), xout(m,1) 
{
  int i,j,k;

  for(i=0;i<m;i++){
    xout[i]=0;
    for(j=0;j<m;j++){
      xout[i]=xout[i]+ mat[i*m+j] * xin[j];
    }
  }
  
}


void mat_vect_mult(double mat[], int m,  int n, double xin[], double xout[])
// matrix (m,n), xin(n,1), xout(m,1) 
{
  int i,j,k;

  for(i=0;i<m;i++){
    xout[i]=0;
    for(j=0;j<n;j++){
      xout[i]=xout[i]+ mat[i*n+j] * xin[j];
    }
  }
  
}

void mat_vect_mult(double m66[],  double xin[6], double xout[6])
{
  int i,j,k;
  
  for(i=0;i<6;i++){
    xout[i]=0;
    for(j=0;j<6;j++){
      xout[i]=xout[i]+ m66[i*6+j] * xin[j];
    }
  }
  
}

void mat_mult(double a[],  double b[], double c[], int m, int n, int k)
// A[m][n] * B[n][k] ->  C[m][k]
{
int i,j,l,u;
    for (i=0; i<=m-1; i++)
    for (j=0; j<=k-1; j++)
      { u=i*k+j; c[u]=0.0;
        for (l=0; l<=n-1; l++)
          c[u]=c[u]+a[i*n+l]*b[l*k+j];
      }
    return;
}

double mat_det(double a[],int n)
{
  int i,j,k,is=0,js=0,l,u,v;
  double f,det,q,d;
  f=1.0; det=1.0;
  for (k=0; k<=n-2; k++)
    { q=0.0;
    for (i=k; i<=n-1; i++)
      for (j=k; j<=n-1; j++)
	{ l=i*n+j; d=fabs(a[l]);
	if (d>q) { q=d; is=i; js=j;}
	}
    if (q+1.0==1.0)
      { det=0.0; return(det);}
    if (is!=k)
      { f=-f;
      for (j=k; j<=n-1; j++)
	{ u=k*n+j; v=is*n+j;
	d=a[u]; a[u]=a[v]; a[v]=d;
	}
      }
    if (js!=k)
      { f=-f;
      for (i=k; i<=n-1; i++)
	{ u=i*n+js; v=i*n+k;
	d=a[u]; a[u]=a[v]; a[v]=d;
	}
      }
    l=k*n+k;
    det=det*a[l];
    for (i=k+1; i<=n-1; i++)
      { d=a[i*n+k]/a[l];
      for (j=k+1; j<=n-1; j++)
	{ u=i*n+j;
	a[u]=a[u]-d*a[k*n+j];
	}
      }
    }
  det=f*det*a[n*n-1];
  return(det);
}

void   mat_transpose(double a[], double at[], int m, int n)
// a[m][n] --> at[n][m]
{
  int i,j;
  
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      at[i*m+j] = a[j*n+i];
    }
  }

}

int mat22_inv(double a[])
{
  double det,temp;

  det = a[0]*a[3] -a[1]*a[2];
  if( det + 1.0  == 1.0 ){
    cout<<"Failed in 2*2 Matrix Inverting. "<<endl;
    cout<<a[0]<<" "<<a[1]<<endl;
    cout<<a[2]<<" "<<a[3]<<endl;
    return 0;
  }
  else{
    temp= a[0];
    a[0]= a[3] /det;
    a[3]= temp/det;
    a[1]= -a[1]/det;
    a[2]= -a[2]/det;
    return 1;
  }
}  

int mat_inv(double a[],int n)
{
  int *is,*js,i,j,k,l,u,v;
  double d,p;
  is=new int[n];
  js=new int[n];
  for (k=0; k<=n-1; k++)
    { d=0.0;
    for (i=k; i<=n-1; i++)
      for (j=k; j<=n-1; j++)
	{ l=i*n+j; p=abs(a[l]);
	if (p>d) { d=p; is[k]=i; js[k]=j;}
	}
    if (d+1.0==1.0)
      { free(is); free(js); printf("Matrix inverting failed. exit. \n");
      return(0);
      }
    if (is[k]!=k)
      for (j=0; j<=n-1; j++)
	{ u=k*n+j; v=is[k]*n+j;
	p=a[u]; a[u]=a[v]; a[v]=p;
	}
    if (js[k]!=k)
      for (i=0; i<=n-1; i++)
	{ u=i*n+k; v=i*n+js[k];
	p=a[u]; a[u]=a[v]; a[v]=p;
	}
    l=k*n+k;
    a[l]=1.0/a[l];
    for (j=0; j<=n-1; j++)
      if (j!=k)
	{ u=k*n+j; a[u]=a[u]*a[l];}
    for (i=0; i<=n-1; i++)
      if (i!=k)
	for (j=0; j<=n-1; j++)
	  if (j!=k)
	    { u=i*n+j;
	    a[u]=a[u]-a[i*n+k]*a[k*n+j];
	    }
    for (i=0; i<=n-1; i++)
      if (i!=k)
	{ u=i*n+k; a[u]=-a[u]*a[l];}
    }
  for (k=n-1; k>=0; k--)
    { if (js[k]!=k)
      for (j=0; j<=n-1; j++)
	{ u=k*n+j; v=js[k]*n+j;
	p=a[u]; a[u]=a[v]; a[v]=p;
	}
    if (is[k]!=k)
      for (i=0; i<=n-1; i++)
	{ u=i*n+k; v=i*n+is[k];
	p=a[u]; a[u]=a[v]; a[v]=p;
	}
    }
  free(is); free(js);
  return(1);
}

void mat_change_hessenberg(double a[], int n)
{ 
  int i,j,k,u,v;
  double d,t;
  for (k=1; k<=n-2; k++)
    { d=0.0;
    for (j=k; j<=n-1; j++)
      { u=j*n+k-1; t=a[u];
      if (fabs(t)>fabs(d))
	{ d=t; i=j;}
      }
    if (fabs(d)+1.0!=1.0)
      { if (i!=k)
	{ for (j=k-1; j<=n-1; j++)
	  { u=i*n+j; v=k*n+j;
	  t=a[u]; a[u]=a[v]; a[v]=t;
	  }
	for (j=0; j<=n-1; j++)
	  { u=j*n+i; v=j*n+k;
	  t=a[u]; a[u]=a[v]; a[v]=t;
	  }
	}
      for (i=k+1; i<=n-1; i++)
	{ u=i*n+k-1; t=a[u]/d; a[u]=0.0;
	for (j=k; j<=n-1; j++)
	  { v=i*n+j;
	  a[v]=a[v]-t*a[k*n+j];
	  }
	for (j=0; j<=n-1; j++)
	  { v=j*n+k;
	  a[v]=a[v]+t*a[j*n+i];
	  }
	}
      }
    }
  return;
}

int mat_root_hessenberg(double a[],int n,double u[],double v[],double eps,int jt)
{ 
  int m,it,i,j,k,l,ii,jj,kk,ll;
  double b,c,w,g,xy,p,q,r,x,s,e,f,z,y;
  it=0; m=n;
  while (m!=0)
    { l=m-1;
    while ((l>0)&&(fabs(a[l*n+l-1])>eps*
		   (fabs(a[(l-1)*n+l-1])+fabs(a[l*n+l])))) l=l-1;
    ii=(m-1)*n+m-1; jj=(m-1)*n+m-2;
    kk=(m-2)*n+m-1; ll=(m-2)*n+m-2;
    if (l==m-1)
      { u[m-1]=a[(m-1)*n+m-1]; v[m-1]=0.0;
      m=m-1; it=0;
      }
    else if (l==m-2)
      { b=-(a[ii]+a[ll]);
      c=a[ii]*a[ll]-a[jj]*a[kk];
      w=b*b-4.0*c;
      y=sqrt(fabs(w));
      if (w>0.0)
	{ xy=1.0;
	if (b<0.0) xy=-1.0;
	u[m-1]=(-b-xy*y)/2.0;
	u[m-2]=c/u[m-1];
	v[m-1]=0.0; v[m-2]=0.0;
	}
      else
	{ u[m-1]=-b/2.0; u[m-2]=u[m-1];
	v[m-1]=y/2.0; v[m-2]=-v[m-1];
	}
      m=m-2; it=0;
      }
    else
      { if (it>=jt)
	{ printf("fail\n");
	return(-1);
	}
      it=it+1;
      for (j=l+2; j<=m-1; j++)
	a[j*n+j-2]=0.0;
      for (j=l+3; j<=m-1; j++)
	a[j*n+j-3]=0.0;
      for (k=l; k<=m-2; k++)
	{ if (k!=l)
	  { p=a[k*n+k-1]; q=a[(k+1)*n+k-1];
	  r=0.0;
	  if (k!=m-2) r=a[(k+2)*n+k-1];
	  }
	else
	  { x=a[ii]+a[ll];
	  y=a[ll]*a[ii]-a[kk]*a[jj];
	  ii=l*n+l; jj=l*n+l+1;
	  kk=(l+1)*n+l; ll=(l+1)*n+l+1;
	  p=a[ii]*(a[ii]-x)+a[jj]*a[kk]+y;
	  q=a[kk]*(a[ii]+a[ll]-x);
	  r=a[kk]*a[(l+2)*n+l+1];
	  }
	if ((fabs(p)+fabs(q)+fabs(r))!=0.0)
	  { xy=1.0;
	  if (p<0.0) xy=-1.0;
	  s=xy*sqrt(p*p+q*q+r*r);
	  if (k!=l) a[k*n+k-1]=-s;
	  e=-q/s; f=-r/s; x=-p/s;
	  y=-x-f*r/(p+s);
	  g=e*r/(p+s);
	  z=-x-e*q/(p+s);
	  for (j=k; j<=m-1; j++)
	    { ii=k*n+j; jj=(k+1)*n+j;
	    p=x*a[ii]+e*a[jj];
	    q=e*a[ii]+y*a[jj];
	    r=f*a[ii]+g*a[jj];
	    if (k!=m-2)
	      { kk=(k+2)*n+j;
	      p=p+f*a[kk];
	      q=q+g*a[kk];
	      r=r+z*a[kk]; a[kk]=r;
	      }
	    a[jj]=q; a[ii]=p;
	    }
	  j=k+3;
	  if (j>=m-1) j=m-1;
	  for (i=l; i<=j; i++)
	    { ii=i*n+k; jj=i*n+k+1;
	    p=x*a[ii]+e*a[jj];
	    q=e*a[ii]+y*a[jj];
	    r=f*a[ii]+g*a[jj];
	    if (k!=m-2)
	      { kk=i*n+k+2;
	      p=p+f*a[kk];
	      q=q+g*a[kk];
	      r=r+z*a[kk]; a[kk]=r;
	      }
	    a[jj]=q; a[ii]=p;
	    }
	  }
	}
      }
    }
  return(1);
}

template<class T> void vector_cross_product(T a[], T b[], T c[])
{
  c[0]=      a[1]*b[2] - a[2]*b[1];
  c[1]=-1.0* a[0]*b[2] + a[2]*b[0];
  c[2]=      a[0]*b[1] - a[1]*b[0];
}

void LinearEquations(double a1, double b1, double c1, double a2, double b2, double c2, double & x, double & y  )
//     a1*x +b1*y= c1
//     a2*x +b2*y= c2
{
  if(a1*b2-a2*b1 != 0. ){
    x=(c1*b2-c2*b1) / ( a1*b2-a2*b1);
    y=(c1*a2-c2*a1) / ( a2*b1-a1*b2);
  }
}

//---counter-clockwise rotation (x0, y0) vetctor to (x1, y1) vetcor
void rotation(double & x1, double & y1, double x0, double y0, double dtheta)
{
  x1= x0*cos(dtheta) - y0*sin(dtheta);
  y1= x0*sin(dtheta) + y0*cos(dtheta);  
}

void  EigenSolver(double Matrix[4][4] , double wr[4], double  wi[4], double vr[4][4], double vi[4][4]) 
{
  int i, j;
  double data[16];

  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      data[i*4+j]= Matrix[i][j];
  
  gsl_matrix_view m  = gsl_matrix_view_array (data, 4, 4);
  gsl_vector_complex *eval = gsl_vector_complex_alloc (4);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc (4, 4);
  gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (4);
  gsl_eigen_nonsymmv (&m.matrix, eval, evec, w);
  gsl_eigen_nonsymmv_free (w);
  gsl_eigen_nonsymmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);
  
  for (i = 0; i < 4; i++)
    {
      gsl_complex eval_i = gsl_vector_complex_get (eval, i);
      gsl_vector_complex_view evec_i = gsl_matrix_complex_column (evec, i);
      wr[i]=  	GSL_REAL(eval_i);
      wi[i]=  	GSL_IMAG(eval_i);
      
      for (j = 0; j < 4; j++)
	{
	  gsl_complex z = gsl_vector_complex_get(&evec_i.vector, j);	
	  vr[i][j]=GSL_REAL(z);
	  vi[i][j]=GSL_IMAG(z);
	}
    }
  
  gsl_vector_complex_free(eval);
  gsl_matrix_complex_free(evec);
}

void  EigenSolver_6D(double Matrix[6][6] , double wr[6], double  wi[6], double vr[6][6], double vi[6][6]) 
{
  int i, j;
  double data[36];

  for(i=0;i<6;i++)
    for(j=0;j<6;j++)
      data[i*6+j]= Matrix[i][j];
  
  gsl_matrix_view m  = gsl_matrix_view_array (data, 6, 6);
  gsl_vector_complex *eval = gsl_vector_complex_alloc (6);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc (6, 6);
  gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (6);
  gsl_eigen_nonsymmv (&m.matrix, eval, evec, w);
  gsl_eigen_nonsymmv_free (w);
  gsl_eigen_nonsymmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);
  
  for (i = 0; i < 6; i++)
    {
      gsl_complex eval_i = gsl_vector_complex_get (eval, i);
      gsl_vector_complex_view evec_i = gsl_matrix_complex_column (evec, i);
      wr[i]=  	GSL_REAL(eval_i);
      wi[i]=  	GSL_IMAG(eval_i);
      
      for (j = 0; j < 6; j++)
	{
	  gsl_complex z = gsl_vector_complex_get(&evec_i.vector, j);	
	  vr[i][j]=GSL_REAL(z);
	  vi[i][j]=GSL_IMAG(z);
	}
    }
  
  gsl_vector_complex_free(eval);
  gsl_matrix_complex_free(evec);
}

void brmul(double a[], double b[] , int m, int n, int k, double c[] )
{ 
  int i,j,l,u;
  for (i=0; i<=m-1; i++)
    for (j=0; j<=k-1; j++)
      { u=i*k+j; c[u]=0.0;
        for (l=0; l<=n-1; l++)
          c[u]=c[u]+a[i*n+l]*b[l*k+j];
      }
  return;
}

static void ppp(double *a,double *e,double *s,double *v,int m,int n)
{ 
  int i,j,p,q;
  double d;
  if (m>=n) i=n;
 else i=m;
 for (j=1; j<=i-1; j++)
   { a[(j-1)*n+j-1]=s[j-1];
   a[(j-1)*n+j]=e[j-1];
      }
    a[(i-1)*n+i-1]=s[i-1];
    if (m<n) a[(i-1)*n+i]=e[i-1];
    for (i=1; i<=n-1; i++)
    for (j=i+1; j<=n; j++)
      { p=(i-1)*n+j-1; q=(j-1)*n+i-1;
        d=v[p]; v[p]=v[q]; v[q]=d;
      }
    return;
  }

static void sss(double fg[2],double cs[2])
{ double r,d;
 if ((fabs(fg[0])+fabs(fg[1]))==0.0)
      { cs[0]=1.0; cs[1]=0.0; d=0.0;}
    else 
      { d=sqrt(fg[0]*fg[0]+fg[1]*fg[1]);
        if (fabs(fg[0])>fabs(fg[1]))
          { d=fabs(d);
            if (fg[0]<0.0) d=-d;
          }
        if (fabs(fg[1])>=fabs(fg[0]))
          { d=fabs(d);
            if (fg[1]<0.0) d=-d;
          }
        cs[0]=fg[0]/d; cs[1]=fg[1]/d;
      }
    r=1.0;
    if (fabs(fg[0])>fabs(fg[1])) r=cs[1];
    else
      if (cs[0]!=0.0) r=1.0/cs[0];
    fg[0]=d; fg[1]=r;
    return;
  }

int bmuav(double a[],int m,int n,double *u,double *v,double eps,int ka)
//  A = U S VT,  ka= max{m,n}+1
//  double *v actually for VT instead of V !
{ 
 int i,j,k,l,it,ll,kk,ix,iy,mm,nn,iz,m1,ks;
 double d,dd,t,sm,sm1,em1,sk,ek,b,c,shh,fg[2],cs[2];
 double *s,*e,*w;
  
 s= new double[ka];
 e= new double[ka];
 w= new double[ka];

 it=60; k=n;
 if (m-1<n) k=m-1;
 l=m;
 if (n-2<m) l=n-2;
 if (l<0) l=0;
 ll=k;
 if (l>k) ll=l;
 if (ll>=1)
   { for (kk=1; kk<=ll; kk++)
     { if (kk<=k)
       { d=0.0;
       for (i=kk; i<=m; i++)
	 { ix=(i-1)*n+kk-1; d=d+a[ix]*a[ix];}
       s[kk-1]=sqrt(d);
       if (s[kk-1]!=0.0)
	 { ix=(kk-1)*n+kk-1;
	 if (a[ix]!=0.0)
	   { s[kk-1]=fabs(s[kk-1]);
	   if (a[ix]<0.0) s[kk-1]=-s[kk-1];
	   }
	 for (i=kk; i<=m; i++)
	   { iy=(i-1)*n+kk-1;
	   a[iy]=a[iy]/s[kk-1];
	   }
	 a[ix]=1.0+a[ix];
	 }
       s[kk-1]=-s[kk-1];
       }
     if (n>=kk+1)
       { for (j=kk+1; j<=n; j++)
	 { if ((kk<=k)&&(s[kk-1]!=0.0))
	   { d=0.0;
	   for (i=kk; i<=m; i++)
	     { ix=(i-1)*n+kk-1;
	     iy=(i-1)*n+j-1;
	     d=d+a[ix]*a[iy];
	     }
	   d=-d/a[(kk-1)*n+kk-1];
	   for (i=kk; i<=m; i++)
	     { ix=(i-1)*n+j-1;
	     iy=(i-1)*n+kk-1;
	     a[ix]=a[ix]+d*a[iy];
	     }
	   }
                    e[j-1]=a[(kk-1)*n+j-1];
	 }
              }
            if (kk<=k)
              { for (i=kk; i<=m; i++)
                  { ix=(i-1)*m+kk-1; iy=(i-1)*n+kk-1;
                    u[ix]=a[iy];
                  }
              }
            if (kk<=l)
              { d=0.0;
                for (i=kk+1; i<=n; i++)
                  d=d+e[i-1]*e[i-1];
                e[kk-1]=sqrt(d);
                if (e[kk-1]!=0.0)
                  { if (e[kk]!=0.0)
                      { e[kk-1]=fabs(e[kk-1]);
                        if (e[kk]<0.0) e[kk-1]=-e[kk-1];
                      }
                    for (i=kk+1; i<=n; i++)
                      e[i-1]=e[i-1]/e[kk-1];
                    e[kk]=1.0+e[kk];
                  }
                e[kk-1]=-e[kk-1];
                if ((kk+1<=m)&&(e[kk-1]!=0.0))
                  { for (i=kk+1; i<=m; i++) w[i-1]=0.0;
                    for (j=kk+1; j<=n; j++)
                      for (i=kk+1; i<=m; i++)
                        w[i-1]=w[i-1]+e[j-1]*a[(i-1)*n+j-1];
                    for (j=kk+1; j<=n; j++)
                      for (i=kk+1; i<=m; i++)
                        { ix=(i-1)*n+j-1;
                          a[ix]=a[ix]-w[i-1]*e[j-1]/e[kk];
                        }
                  }
                for (i=kk+1; i<=n; i++)
                  v[(i-1)*n+kk-1]=e[i-1];
              }
          }
      }
    mm=n;
    if (m+1<n) mm=m+1;
    if (k<n) s[k]=a[k*n+k];
    if (m<mm) s[mm-1]=0.0;
    if (l+1<mm) e[l]=a[l*n+mm-1];
    e[mm-1]=0.0;
    nn=m;
    if (m>n) nn=n;
    if (nn>=k+1)
      { for (j=k+1; j<=nn; j++)
          { for (i=1; i<=m; i++)
              u[(i-1)*m+j-1]=0.0;
            u[(j-1)*m+j-1]=1.0;
          }
      }
    if (k>=1)
      { for (ll=1; ll<=k; ll++)
          { kk=k-ll+1; iz=(kk-1)*m+kk-1;
            if (s[kk-1]!=0.0)
              { if (nn>=kk+1)
                  for (j=kk+1; j<=nn; j++)
                    { d=0.0;
                      for (i=kk; i<=m; i++)
                        { ix=(i-1)*m+kk-1;
                          iy=(i-1)*m+j-1;
                          d=d+u[ix]*u[iy]/u[iz];
                        }
                      d=-d;
                      for (i=kk; i<=m; i++)
                        { ix=(i-1)*m+j-1;
                          iy=(i-1)*m+kk-1;
                          u[ix]=u[ix]+d*u[iy];
                        }
                    }
                  for (i=kk; i<=m; i++)
                    { ix=(i-1)*m+kk-1; u[ix]=-u[ix];}
                  u[iz]=1.0+u[iz];
                  if (kk-1>=1)
                    for (i=1; i<=kk-1; i++)
                      u[(i-1)*m+kk-1]=0.0;
              }
            else
              { for (i=1; i<=m; i++)
                  u[(i-1)*m+kk-1]=0.0;
                u[(kk-1)*m+kk-1]=1.0;
              }
          }
      }
    for (ll=1; ll<=n; ll++)
      { kk=n-ll+1; iz=kk*n+kk-1;
        if ((kk<=l)&&(e[kk-1]!=0.0))
          { for (j=kk+1; j<=n; j++)
              { d=0.0;
                for (i=kk+1; i<=n; i++)
                  { ix=(i-1)*n+kk-1; iy=(i-1)*n+j-1;
                    d=d+v[ix]*v[iy]/v[iz];
                  }
                d=-d;
                for (i=kk+1; i<=n; i++)
                  { ix=(i-1)*n+j-1; iy=(i-1)*n+kk-1;
                    v[ix]=v[ix]+d*v[iy];
                  }
              }
          }
        for (i=1; i<=n; i++)
          v[(i-1)*n+kk-1]=0.0;
        v[iz-n]=1.0;
      }
    for (i=1; i<=m; i++)
    for (j=1; j<=n; j++)
      a[(i-1)*n+j-1]=0.0;
    m1=mm; it=60;
    while (1==1)
      { if (mm==0)
          { ppp(a,e,s,v,m,n);
	    free(s); free(e); free(w); return(1);
          }
        if (it==0)
          { ppp(a,e,s,v,m,n);
            free(s); free(e); free(w); return(-1);
          }
        kk=mm-1;
	while ((kk!=0)&&(fabs(e[kk-1])!=0.0))
          { d=fabs(s[kk-1])+fabs(s[kk]);
            dd=fabs(e[kk-1]);
            if (dd>eps*d) kk=kk-1;
            else e[kk-1]=0.0;
          }
        if (kk==mm-1)
          { kk=kk+1;
            if (s[kk-1]<0.0)
              { s[kk-1]=-s[kk-1];
                for (i=1; i<=n; i++)
                  { ix=(i-1)*n+kk-1; v[ix]=-v[ix];}
              }
            while ((kk!=m1)&&(s[kk-1]<s[kk]))
              { d=s[kk-1]; s[kk-1]=s[kk]; s[kk]=d;
                if (kk<n)
                  for (i=1; i<=n; i++)
                    { ix=(i-1)*n+kk-1; iy=(i-1)*n+kk;
                      d=v[ix]; v[ix]=v[iy]; v[iy]=d;
                    }
                if (kk<m)
                  for (i=1; i<=m; i++)
                    { ix=(i-1)*m+kk-1; iy=(i-1)*m+kk;
                      d=u[ix]; u[ix]=u[iy]; u[iy]=d;
                    }
                kk=kk+1;
              }
            it=60;
            mm=mm-1;
          }
        else
          { ks=mm;
            while ((ks>kk)&&(fabs(s[ks-1])!=0.0))
              { d=0.0;
                if (ks!=mm) d=d+fabs(e[ks-1]);
                if (ks!=kk+1) d=d+fabs(e[ks-2]);
                dd=fabs(s[ks-1]);
                if (dd>eps*d) ks=ks-1;
                else s[ks-1]=0.0;
              }
            if (ks==kk)
              { kk=kk+1;
                d=fabs(s[mm-1]);
                t=fabs(s[mm-2]);
                if (t>d) d=t;
                t=fabs(e[mm-2]);
                if (t>d) d=t;
                t=fabs(s[kk-1]);
                if (t>d) d=t;
                t=fabs(e[kk-1]);
                if (t>d) d=t;
                sm=s[mm-1]/d; sm1=s[mm-2]/d;
                em1=e[mm-2]/d;
                sk=s[kk-1]/d; ek=e[kk-1]/d;
                b=((sm1+sm)*(sm1-sm)+em1*em1)/2.0;
                c=sm*em1; c=c*c; shh=0.0;
                if ((b!=0.0)||(c!=0.0))
                  { shh=sqrt(b*b+c);
                    if (b<0.0) shh=-shh;
                    shh=c/(b+shh);
                  }
                fg[0]=(sk+sm)*(sk-sm)-shh;
                fg[1]=sk*ek;
                for (i=kk; i<=mm-1; i++)
                  { sss(fg,cs);
                    if (i!=kk) e[i-2]=fg[0];
                    fg[0]=cs[0]*s[i-1]+cs[1]*e[i-1];
                    e[i-1]=cs[0]*e[i-1]-cs[1]*s[i-1];
                    fg[1]=cs[1]*s[i];
                    s[i]=cs[0]*s[i];
                    if ((cs[0]!=1.0)||(cs[1]!=0.0))
                      for (j=1; j<=n; j++)
                        { ix=(j-1)*n+i-1;
                          iy=(j-1)*n+i;
                          d=cs[0]*v[ix]+cs[1]*v[iy];
                          v[iy]=-cs[1]*v[ix]+cs[0]*v[iy];
                          v[ix]=d;
                        }
                    sss(fg,cs);
                    s[i-1]=fg[0];
                    fg[0]=cs[0]*e[i-1]+cs[1]*s[i];
                    s[i]=-cs[1]*e[i-1]+cs[0]*s[i];
                    fg[1]=cs[1]*e[i];
                    e[i]=cs[0]*e[i];
                    if (i<m)
                      if ((cs[0]!=1.0)||(cs[1]!=0.0))
                        for (j=1; j<=m; j++)
                          { ix=(j-1)*m+i-1;
                            iy=(j-1)*m+i;
                            d=cs[0]*u[ix]+cs[1]*u[iy];
                            u[iy]=-cs[1]*u[ix]+cs[0]*u[iy];
                            u[ix]=d;
                          }
                  }
                e[mm-2]=fg[0];
                it=it-1;
              }
            else
              { if (ks==mm)
                  { kk=kk+1;
                    fg[1]=e[mm-2]; e[mm-2]=0.0;
                    for (ll=kk; ll<=mm-1; ll++)
                      { i=mm+kk-ll-1;
                        fg[0]=s[i-1];
                        sss(fg,cs);
                        s[i-1]=fg[0];
                        if (i!=kk)
                          { fg[1]=-cs[1]*e[i-2];
                            e[i-2]=cs[0]*e[i-2];
                          }
                        if ((cs[0]!=1.0)||(cs[1]!=0.0))
                          for (j=1; j<=n; j++)
                            { ix=(j-1)*n+i-1;
                              iy=(j-1)*n+mm-1;
                              d=cs[0]*v[ix]+cs[1]*v[iy];
                              v[iy]=-cs[1]*v[ix]+cs[0]*v[iy];
                              v[ix]=d;
                            }
                      }
                  }
                else
                  { kk=ks+1;
                    fg[1]=e[kk-2];
                    e[kk-2]=0.0;
                    for (i=kk; i<=mm; i++)
                      { fg[0]=s[i-1];
                        sss(fg,cs);
                        s[i-1]=fg[0];
                        fg[1]=-cs[1]*e[i-1];
                        e[i-1]=cs[0]*e[i-1];
                        if ((cs[0]!=1.0)||(cs[1]!=0.0))
                          for (j=1; j<=m; j++)
                            { ix=(j-1)*m+i-1;
                              iy=(j-1)*m+kk-2;
                              d=cs[0]*u[ix]+cs[1]*u[iy];
                              u[iy]=-cs[1]*u[ix]+cs[0]*u[iy];
                              u[ix]=d;
                            }
                      }
                  }
              }
          }
      }
}

void mat_inv_svd(int m, int  n, double A[])
{
  int i,j,p, l;
  double  U[m*m],  VT[n*n], Sigma[m*n];
  double  UT[m*m], V[n*n],  SigmaT[n*m], AT[n*m];
  double  eps, cut, cut_scale=1000;

  p=m;
  if(n>m) p=n;
  eps=1.0e-6;  // maximum of m and n for svd input

  l=m;
  if(n<m) l=n; // maximum number of singular values
  
  for(i=0;i<m*m;i++)   U[i]=0.;
  for(i=0;i<n*n;i++)   V[i]=0.;
  for(i=0;i<n*m;i++)   Sigma[i]=0.; 
  
  i=bmuav(A,m,n,U, VT, eps, p+1);
  if( i<= 0) { cout<<" SVD failed"<<endl;  exit(0); }

  for(i=0;i<m*n;i++) Sigma[i]=A[i];
  cut=Sigma[0]/cut_scale;

  for(i=0;i<l;i++) {
    if( Sigma[i*n+i] < cut) {
      Sigma[i*n+i] = 0.;
    }
    else{
      Sigma[i*n+i] = 1.0/ Sigma[i*n+i];
    }
  }
  
  mat_transpose(Sigma, SigmaT, m,n);
  mat_transpose(U, UT, m,m);
  mat_transpose(VT,V, n,n);

  mat_mult(V, SigmaT, AT, n, n, m);
  mat_mult(AT, UT, SigmaT, n, m, m);

  for(i=0;i<m*n;i++) A[i] = SigmaT[i];
  
}

void fft(int m, double*x, double*y)
{
   short int dir=1;
   int n,i,i1,j,k,i2,l,l1,l2;
   double c1,c2,tx,ty,t1,t2,u1,u2,z;

   n = 1;
   for (i=0;i<m;i++) 
      n *= 2;

   i2 = n >> 1;
   j = 0;
   for (i=0;i<n-1;i++) {
      if (i < j) {
         tx = x[i];
         ty = y[i];
         x[i] = x[j];
         y[i] = y[j];
         x[j] = tx;
         y[j] = ty;
      }
      k = i2;
      while (k <= j) {
         j -= k;
         k >>= 1;
      }
      j += k;
   }

   c1 = -1.0; 
   c2 = 0.0;
   l2 = 1;
   for (l=0;l<m;l++) {
      l1 = l2;
      l2 <<= 1;
      u1 = 1.0; 
      u2 = 0.0;
      for (j=0;j<l1;j++) {
         for (i=j;i<n;i+=l2) {
            i1 = i + l1;
            t1 = u1 * x[i1] - u2 * y[i1];
            t2 = u1 * y[i1] + u2 * x[i1];
            x[i1] = x[i] - t1; 
            y[i1] = y[i] - t2;
            x[i] += t1;
            y[i] += t2;
         }
         z =  u1 * c1 - u2 * c2;
         u2 = u1 * c2 + u2 * c1;
         u1 = z;
      }
      c2 = sqrt((1.0 - c1) / 2.0);
      if (dir == 1) 
         c2 = -c2;
      c1 = sqrt((1.0 + c1) / 2.0);
   }

   if (dir == 1) {
      for (i=0;i<n;i++) {
         x[i] /= n;
         y[i] /= n;
      }
   }
}

void ifft(int m, double*x, double*y)
{
   short int dir=-1;
   int n,i,i1,j,k,i2,l,l1,l2;
   double c1,c2,tx,ty,t1,t2,u1,u2,z;

   n = 1;
   for (i=0;i<m;i++) 
      n *= 2;

   i2 = n >> 1;
   j = 0;
   for (i=0;i<n-1;i++) {
      if (i < j) {
         tx = x[i];
         ty = y[i];
         x[i] = x[j];
         y[i] = y[j];
         x[j] = tx;
         y[j] = ty;
      }
      k = i2;
      while (k <= j) {
         j -= k;
         k >>= 1;
      }
      j += k;
   }

   c1 = -1.0; 
   c2 = 0.0;
   l2 = 1;
   for (l=0;l<m;l++) {
      l1 = l2;
      l2 <<= 1;
      u1 = 1.0; 
      u2 = 0.0;
      for (j=0;j<l1;j++) {
         for (i=j;i<n;i+=l2) {
            i1 = i + l1;
            t1 = u1 * x[i1] - u2 * y[i1];
            t2 = u1 * y[i1] + u2 * x[i1];
            x[i1] = x[i] - t1; 
            y[i1] = y[i] - t2;
            x[i] += t1;
            y[i] += t2;
         }
         z =  u1 * c1 - u2 * c2;
         u2 = u1 * c2 + u2 * c1;
         u1 = z;
      }
      c2 = sqrt((1.0 - c1) / 2.0);
      if (dir == 1) 
         c2 = -c2;
      c1 = sqrt((1.0 + c1) / 2.0);
   }

   if (dir == 1) {
      for (i=0;i<n;i++) {
         x[i] /= n;
         y[i] /= n;
      }
   }
}

void FineTuneFinder(int n, double x[], double tune_lower, double tune_upper, double & peaktune)
{
  int i;
  double tune, tune_step = 0.001;
  double sumr, sumi,amp; 
  double peakamp;
  
  do{

    tune=tune_lower;
    peaktune=tune_lower;
    peakamp=0;
    do{
      sumr=0.; sumi=0;
      for(i=0;i<n;i++){
	sumr = sumr + x[i] * cos(-2*PI*i*tune);
	sumi = sumi + x[i] * sin(-2*PI*i*tune) ;
      } 
      amp = sumr*sumr+ sumi*sumi;
      if(amp  > peakamp ) {
	peaktune = tune;
	peakamp =  amp ;
      }
      tune=tune + tune_step;
    } while (tune < tune_upper);

    tune_lower=peaktune - tune_step;
    tune_upper=peaktune + tune_step;
    tune_step =tune_step/10.;
  }while(tune_step > 1.0e-12);
  
}

void FineTuneFinder(int n, double x[], double px[], double tune_lower, double tune_upper, double & peaktune)
{
  int i;
  double tune, tune_step = 0.001;
  double sumr, sumi,  amp; 
  double peakamp;
  
  do{

    tune=tune_lower;
    peaktune=tune_lower;
    peakamp=0;
    do{
      sumr=0; sumi=0.;
      for(i=0;i<n;i++){
	sumr = sumr + x[i] * cos(-2*PI*i*tune) -  px[i] *  sin(-2*PI*i*tune) ;
	sumi = sumi + x[i] * sin(-2*PI*i*tune) +  px[i] *  cos(-2*PI*i*tune) ;
      } 
      amp = sumr*sumr + sumi*sumi;
      if(amp  > peakamp ) {
	peaktune = tune;
	peakamp =  amp ;
      }
      tune=tune + tune_step;
    } while (tune < tune_upper);

    tune_lower=peaktune - tune_step;
    tune_upper=peaktune + tune_step;
    tune_step =tune_step/10.;
  }while(tune_step > 1.0e-12);
  
}


//----4TH ORDER RUNGE-KUTTA METHOD : ONLY WITH B FIELD, in Cartersian coordinate system   r=(x,y,z)

void Track_RungeKutta_BField_1step(double m0, double q, double gamma, double r0[3], double v0[3], void (*GetB)(double r[3], double B[]), double dt)
{
  int i,j;
  double k1[3], l1[3], k2[3], l2[3], k3[3], l3[3], k4[3], l4[3];
  double scale = q / m0 / gamma ;
  double rtemp[3], vtemp[3], Btemp[3], atemp[3];

  GetB(r0,Btemp);
  vector_cross_product(v0,  Btemp,  atemp);
  for(i=0;i<3;i++)  k1[i]= dt * scale * atemp[i];
  for(i=0;i<3;i++)  l1[i]= dt * v0[i];
  
  for(i=0;i<3;i++) vtemp[i] = v0[i] + 0.5 * k1[i];
  for(i=0;i<3;i++) rtemp[i] = r0[i] + 0.5 * l1[i];
  GetB(rtemp,Btemp);
  vector_cross_product(vtemp,  Btemp, atemp);
  for(i=0;i<3;i++)  k2[i]= dt * scale * atemp[i];
  for(i=0;i<3;i++)  l2[i]= dt *( v0[i] + 0.5*k1[i] );

  for(i=0;i<3;i++) vtemp[i] = v0[i] + 0.5 * k2[i];
  for(i=0;i<3;i++) rtemp[i] = r0[i] + 0.5 * l2[i];
  GetB(rtemp,Btemp);
  vector_cross_product(vtemp,  Btemp, atemp);
  for(i=0;i<3;i++)  k3[i]= dt * scale * atemp[i];
  for(i=0;i<3;i++)  l3[i]= dt * ( v0[i] +  0.5*k2[i] );

  for(i=0;i<3;i++) vtemp[i] = v0[i] +  k3[i];
  for(i=0;i<3;i++) rtemp[i] = r0[i] +  l3[i];
  GetB(rtemp,Btemp);
  vector_cross_product(vtemp,  Btemp, atemp);
  for(i=0;i<3;i++)  k4[i]= dt * scale * atemp[i];
  for(i=0;i<3;i++)  l4[i]= dt * (v0[i] +  k3[i]);

  for(i=0;i<3;i++)
    v0[i]= v0[i] + (k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;
  for(i=0;i<3;i++)
    r0[i]= r0[i] + (l1[i]+2*l2[i]+2*l3[i]+l4[i])/6.0;  

}

//----4TH ORDER RUNGE-KUTTA METHOD : WITH BOTH E ADN B FIELDS, in Cartesian coordinate system

void Track_RungeKutta_EMField_1step(double m0, double q, double gamma0, double r0[3], double p0[3], void (*GetB)(double r[3], double B[]), double dt)
{
  int k;
  double v0[3];
  double k1[3], l1[3], k2[3], l2[3], k3[3], l3[3], k4[3], l4[3];
  double scale = q / m0 ;
  double rtemp[3], ptemp[3];
  double gtemp, vtemp[3], Btemp[3], Etemp[3], atemp[3];

  for(k=0;k<3;k++)  v0[k] = p0[k] / m0 / gamma0;
  GetB(r0, Btemp);
  for(k=0;k<3;k++)  Etemp[k] = 0.;
  vector_cross_product(v0,  Btemp,  atemp);
  for(k=0;k<3;k++)  atemp[k] =  (atemp[k]  + Etemp[k] )  * q ;
  for(k=0;k<3;k++)  k1[k]= dt * atemp[k] ;
  for(k=0;k<3;k++)  l1[k]= dt * v0[k];
  
  for(k=0;k<3;k++)  ptemp[k] = p0[k] + 0.5 * k1[k];
  for(k=0;k<3;k++)  rtemp[k] = r0[k] + 0.5 * l1[k];
  GetB(rtemp,Btemp);
  for(k=0;k<3;k++)  Etemp[k] = 0.;
  gtemp= sqrt(  (ptemp[0]*ptemp[0]+ptemp[1]*ptemp[1]+ptemp[2]*ptemp[2]) /  m0 / m0 / 299792458. / 299792458.  + 1.0 ) ;
  for(k=0;k<3;k++)  vtemp[k] = ptemp[k]/m0/gtemp;
  vector_cross_product(vtemp,  Btemp,  atemp);  
  for(k=0;k<3;k++)  atemp[k] =  (atemp[k]  + Etemp[k] )  * q ;
  for(k=0;k<3;k++)  k2[k]= dt * atemp[k];
  for(k=0;k<3;k++)  l2[k]= dt * vtemp[k];

  for(k=0;k<3;k++)  ptemp[k] = p0[k] + 0.5 * k2[k];
  for(k=0;k<3;k++)  rtemp[k] = r0[k] + 0.5 * l2[k];
  GetB(rtemp,Btemp);
  for(k=0;k<3;k++)  Etemp[k] = 0.;
  gtemp= sqrt(  (ptemp[0]*ptemp[0]+ptemp[1]*ptemp[1]+ptemp[2]*ptemp[2]) /  m0 / m0 / 299792458. / 299792458.  + 1.0 ) ;
  for(k=0;k<3;k++)  vtemp[k] = ptemp[k]/m0/gtemp;
  vector_cross_product(vtemp,  Btemp,  atemp);  
  for(k=0;k<3;k++)  atemp[k] =  (atemp[k]  + Etemp[k] )  * q ;
  for(k=0;k<3;k++)  k3[k]= dt * atemp[k];
  for(k=0;k<3;k++)  l3[k]= dt * vtemp[k];

  for(k=0;k<3;k++)  ptemp[k] = p0[k] +  k2[k];
  for(k=0;k<3;k++)  rtemp[k] = r0[k] +  l2[k];
  GetB(rtemp,Btemp);
  for(k=0;k<3;k++)  Etemp[k] = 0.;
  gtemp= sqrt(  (ptemp[0]*ptemp[0]+ptemp[1]*ptemp[1]+ptemp[2]*ptemp[2]) /  m0 / m0 / 299792458. / 299792458.  + 1.0 ) ;
  for(k=0;k<3;k++)  vtemp[k] = ptemp[k]/m0/gtemp;
  vector_cross_product(vtemp,  Btemp,  atemp);  
  for(k=0;k<3;k++)  atemp[k] =  (atemp[k]  + Etemp[k] )  * q ;
  for(k=0;k<3;k++)  k4[k]= dt * atemp[k];
  for(k=0;k<3;k++)  l4[k]= dt * vtemp[k];  
  
  for(k=0;k<3;k++) p0[k]= p0[k] + (k1[k]+2*k2[k]+2*k3[k]+k4[k])/6.0;
  for(k=0;k<3;k++) r0[k]= r0[k] + (l1[k]+2*l2[k]+2*l3[k]+l4[k])/6.0;
  
}

//----4TH ORDER RUNGE-KUTTA METHOD : WITH BOTH E ADN B FIELDS, in Cartesian coordinate system

void Track_RungeKutta_Spin_1step(double m0, double q, double gamma, double r0[3], double p0[3], double s0[3], void (*GetB)(double r[3], double B[]), double dt)
{
  int k;
  double v0[3];
  double k1[3], l1[3], k2[3], l2[3], k3[3], l3[3], k4[3], l4[3], m1[3], m2[3], m3[3], m4[3];
  double scale = q / m0 ;
  double rtemp[3], ptemp[3], stemp[3];
  double gtemp, vtemp[3], Btemp[3], Etemp[3], atemp[3];
  double u[3], vabs, udotB, Bpar[3], htemp[3];

  //----first step
  
  for(k=0;k<3;k++)  v0[k] = p0[k] / gamma /  m0;
  GetB(r0, Btemp);
  for(k=0;k<3;k++)  Etemp[k] = 0.;
  vector_cross_product(v0,  Btemp,  atemp);
  for(k=0;k<3;k++)  atemp[k] =  (atemp[k]  + Etemp[k] )  * q ;
  for(k=0;k<3;k++)  k1[k]= dt * atemp[k] ;
  for(k=0;k<3;k++)  l1[k]= dt * v0[k];

  vabs= sqrt(v0[0]*v0[0] + v0[1]*v0[1] + v0[2]*v0[2] );
  for(k=0;k<3;k++)  u[k]= v0[k]/vabs;
  udotB =  u[0] * Btemp[0] +  u[1] * Btemp[1] +  u[2] * Btemp[2];
  for(k=0;k<3;k++) Bpar[k] =  udotB * u[k];
  for(k=0;k<3;k++) htemp[k] = (1+GP.G * gamma ) * Btemp[k]  - GP.G * (gamma -1.0 ) * Bpar[k];
  for(k=0;k<3;k++) htemp[k] = htemp[k] *q / gamma / m0;
  vector_cross_product(s0,  htemp,  atemp);
  for(k=0;k<3;k++) m1[k]= dt * atemp[k];

  //--second step
					
  for(k=0;k<3;k++)  ptemp[k] = p0[k] + 0.5 * k1[k];
  for(k=0;k<3;k++)  rtemp[k] = r0[k] + 0.5 * l1[k];
  for(k=0;k<3;k++)  stemp[k] = s0[k] + 0.5 * m1[k];  
  GetB(rtemp,Btemp);
  for(k=0;k<3;k++)  Etemp[k] = 0.;
  gtemp= sqrt(  (ptemp[0]*ptemp[0]+ptemp[1]*ptemp[1]+ptemp[2]*ptemp[2]) /  m0 / m0 / 299792458. / 299792458.  + 1.0 ) ;
  for(k=0;k<3;k++)  vtemp[k] = ptemp[k]/m0/gtemp;
  vector_cross_product(vtemp,  Btemp,  atemp);  
  for(k=0;k<3;k++)  atemp[k] =  (atemp[k]  + Etemp[k] )  * q ;
  for(k=0;k<3;k++)  k2[k]= dt * atemp[k];
  for(k=0;k<3;k++)  l2[k]= dt * vtemp[k];

  vabs= sqrt(vtemp[0]*vtemp[0] + vtemp[1]*vtemp[1] + vtemp[2]*vtemp[2] );
  for(k=0;k<3;k++)  u[k]= vtemp[k]/vabs;
  udotB =  u[0] * Btemp[0] +  u[1] * Btemp[1] +  u[2] * Btemp[2];
  for(k=0;k<3;k++) Bpar[k] =  udotB * u[k];
  for(k=0;k<3;k++) htemp[k] = (1+GP.G * gtemp ) * Btemp[k]  - GP.G * (gtemp -1 ) * Bpar[k];
  for(k=0;k<3;k++) htemp[k] = htemp[k] *q / gtemp / m0;
  vector_cross_product(stemp,  htemp,  atemp); 
  for(k=0;k<3;k++) m2[k]= dt * atemp[k];

  //----third step
  
  for(k=0;k<3;k++)  ptemp[k] = p0[k] + 0.5 * k2[k];
  for(k=0;k<3;k++)  rtemp[k] = r0[k] + 0.5 * l2[k];
  for(k=0;k<3;k++)  stemp[k] = s0[k] + 0.5 * m2[k];  
  GetB(rtemp,Btemp);
  for(k=0;k<3;k++)  Etemp[k] = 0.;
  gtemp= sqrt(  (ptemp[0]*ptemp[0]+ptemp[1]*ptemp[1]+ptemp[2]*ptemp[2]) /  m0 / m0 / 299792458. / 299792458.  + 1.0 ) ;
  for(k=0;k<3;k++)  vtemp[k] = ptemp[k]/m0/gtemp;
  vector_cross_product(vtemp,  Btemp,  atemp);  
  for(k=0;k<3;k++)  atemp[k] =  (atemp[k]  + Etemp[k] )  * q ;
  for(k=0;k<3;k++)  k3[k]= dt * atemp[k];
  for(k=0;k<3;k++)  l3[k]= dt * vtemp[k];

  vabs= sqrt(vtemp[0]*vtemp[0] + vtemp[1]*vtemp[1] + vtemp[2]*vtemp[2] );
  for(k=0;k<3;k++)  u[k]= vtemp[k]/vabs;
  udotB =  u[0] * Btemp[0] +  u[1] * Btemp[1] +  u[2] * Btemp[2];
  for(k=0;k<3;k++) Bpar[k] =  udotB * u[k];
  for(k=0;k<3;k++) htemp[k] = (1+GP.G * gtemp ) * Btemp[k]  - GP.G * (gtemp -1 ) *Bpar[k];
  for(k=0;k<3;k++) htemp[k] = htemp[k] *q / gtemp / m0;
  vector_cross_product(stemp,  htemp,  atemp);  
  for(k=0;k<3;k++) m3[k]= dt * atemp[k];

  //--fourth step 

  for(k=0;k<3;k++)  ptemp[k] = p0[k] +  k2[k];
  for(k=0;k<3;k++)  rtemp[k] = r0[k] +  l2[k];
  for(k=0;k<3;k++)  stemp[k] = s0[k] +  m2[k];  
  GetB(rtemp,Btemp);
  for(k=0;k<3;k++)  Etemp[k] = 0.;
  gtemp= sqrt(  (ptemp[0]*ptemp[0]+ptemp[1]*ptemp[1]+ptemp[2]*ptemp[2]) /  m0 / m0 / 299792458. / 299792458.  + 1.0 ) ;
  for(k=0;k<3;k++)  vtemp[k] = ptemp[k]/m0/gtemp;
  vector_cross_product(vtemp,  Btemp,  atemp);  
  for(k=0;k<3;k++)  atemp[k] =  (atemp[k]  + Etemp[k] )  * q ;
  for(k=0;k<3;k++)  k4[k]= dt * atemp[k];
  for(k=0;k<3;k++)  l4[k]= dt * vtemp[k];

  vabs= sqrt(vtemp[0]*vtemp[0] + vtemp[1]*vtemp[1] + vtemp[2]*vtemp[2] );
  for(k=0;k<3;k++)  u[k]= vtemp[k]/vabs;
  udotB =  u[0] * Btemp[0] +  u[1] * Btemp[1] +  u[2] * Btemp[2];
  for(k=0;k<3;k++) Bpar[k] =  udotB * u[k];
  for(k=0;k<3;k++) htemp[k] = (1+GP.G * gtemp ) * Btemp[k]  - GP.G * (gtemp -1 ) *Bpar[k];
  for(k=0;k<3;k++) htemp[k] = htemp[k] *q / gtemp / m0;
  vector_cross_product(stemp,  htemp,  atemp);    
  for(k=0;k<3;k++) m4[k]= dt * atemp[k];

  //-----finish one setp integration

  for(k=0;k<3;k++) p0[k]= p0[k] + (k1[k]+2*k2[k]+2*k3[k]+k4[k])/6.0;
  for(k=0;k<3;k++) r0[k]= r0[k] + (l1[k]+2*l2[k]+2*l3[k]+l4[k])/6.0;
  for(k=0;k<3;k++) s0[k]= s0[k] + (m1[k]+2*m2[k]+2*m3[k]+m4[k])/6.0;

}

int Powerof2(int n,int *m,int *twopm)
{
   if (n <= 1) {
      *m = 0;
      *twopm = 1;
      return(false);
   }

   *m = 1;
   *twopm = 2;
   do {
      (*m)++;
      (*twopm) *= 2;
   } while (2*(*twopm) <= n);

   if (*twopm != n)
      return(false);
   else
      return(true);
}

int FFT1D(double x[],double y[], int m, int dir)
{
   long nn,i,i1,j,k,i2,l,l1,l2;
   double c1,c2,tx,ty,t1,t2,u1,u2,z;

   nn = 1;
   for (i=0;i<m;i++)
      nn *= 2;

   i2 = nn >> 1;
   j = 0;
   for (i=0;i<nn-1;i++) {
      if (i < j) {
         tx = x[i];
         ty = y[i];
         x[i] = x[j];
         y[i] = y[j];
         x[j] = tx;
         y[j] = ty;
      }
      k = i2;
      while (k <= j) {
         j -= k;
         k >>= 1;
      }
      j += k;
   }

   c1 = -1.0;
   c2 = 0.0;
   l2 = 1;
   for (l=0;l<m;l++) {
      l1 = l2;
      l2 <<= 1;
      u1 = 1.0;
      u2 = 0.0;
      for (j=0;j<l1;j++) {
         for (i=j;i<nn;i+=l2) {
            i1 = i + l1;
            t1 = u1 * x[i1] - u2 * y[i1];
            t2 = u1 * y[i1] + u2 * x[i1];
            x[i1] = x[i] - t1;
            y[i1] = y[i] - t2;
            x[i] += t1;
            y[i] += t2;
         }
         z =  u1 * c1 - u2 * c2;
         u2 = u1 * c2 + u2 * c1;
         u1 = z;
      }
      c2 = sqrt((1.0 - c1) / 2.0);
      if (dir == 1)
         c2 = -c2;
      c1 = sqrt((1.0 + c1) / 2.0);
   }

   if (dir == 1) {
      for (i=0;i<nn;i++) {
         x[i] /= (double)nn;
         y[i] /= (double)nn;
      }
   }

   return(true);
}

int FFT2D(double x[], double y[], int nx, int ny,int dir)
{
   int i,j;
   int m, twopm;
   double *real, *imag;

   int* p1 = new int;

   real = new double[nx*ny];
   imag = new double [nx*ny];
   if (real == NULL || imag == NULL)
      return(false);
   if (!Powerof2(nx,&m,&twopm) || twopm != nx)
      return(false);
   for (j=0;j<ny;j++) {
      for (i=0;i<nx;i++) {
         real[i] = x[i*nx+j];
         imag[i] = y[i*nx+j];
      }
      FFT1D(real,imag,m,dir);
      for (i=0;i<nx;i++) {
         x[i*nx+j] = real[i];
         y[i*nx+j] = imag[i];
      }
   }
   delete[] real;
   delete[] imag;

   real = new double[nx*ny];
   imag = new double [nx*ny];
   if (real == NULL || imag == NULL)
      return(false);
   if (!Powerof2(ny,&m,&twopm) || twopm != ny)
      return(false);
   for (i=0;i<nx;i++) {
      for (j=0;j<ny;j++) {
         real[j] = x[i*nx+j];
         imag[j] = y[i*nx+j];
      }
      FFT1D(real,imag,m,dir);
      for (j=0;j<ny;j++) {
         x[i*nx+j] = real[j];
         y[i*nx+j] = imag[j];
      }
   }
   delete[] real;
   delete[] imag;

   return(true);
}

void  Cal_Grid_Weight0(double x1, double y1, double x0, double y0, int Nx, int Ny, double hx,  double hy, int xmesh[], int ymesh[], double xweight[], double yweight[])
{
  int i,j;
  double dx=x1-x0, dy=y1-y0;
  double rx, rx2, ry, ry2;
  double flag1=0, flag2=0;
  
  for(i=1;i<Nx-1;i++){
    if( dx >= i*hx and dx < i*hx + 0.5*hx  ){ 
      xmesh[0]=i; xmesh[1]=i+1; xmesh[2]=i-1; 
      rx = (dx-i*hx)/hx;  rx2=rx*rx ;
      flag1=1;
      break;
    }
    if(dx >= i*hx + 0.5*hx and dx< (i+1)*hx ){ 
      xmesh[0]=i+1; xmesh[1]=i; xmesh[2]=i+2; 
      rx = ((i+1)*hx-dx)/hx;    rx2=rx*rx ;
      flag1=1;
      break;
    }
  }
  xweight[0]= 0.75 - rx2; 
  xweight[1]= 0.5*(0.25+rx+rx2); 
  xweight[2]= 0.5*(0.25-rx+rx2); 
  
  for(i=1;i<Ny-1;i++){
    if( dy >= i*hy and dy < i*hy + 0.5*hy  ){ 
      ymesh[0]=i; ymesh[1]=i+1; ymesh[2]=i-1; 
      ry = (dy-i*hy)/hy;  ry2=ry*ry ;
      flag2=1;
      break;
    }
    if(dy >= i*hy + 0.5*hy and dy< (i+1)*hy ){ 
      ymesh[0]=i+1; ymesh[1]=i; ymesh[2]=i+2; 
      ry = ((i+1)*hy-dy)/hy;    ry2=ry*ry ;
      flag2=1;
      break;
    }
  }
  yweight[0]= 0.75 - ry2; 
  yweight[1]= 0.5*(0.25+ry+ry2); 
  yweight[2]= 0.5*(0.25-ry+ry2);

  if(flag1*flag2==0) {
    cout<<" One partilce out of grids."<<endl;
    for(i=0;i<3;i++) {
      xweight[i]=0;
      yweight[i]=0;
      xmesh[0]=0; xmesh[1]=1; xmesh[2]=2;
      ymesh[0]=0; ymesh[1]=1; ymesh[2]=2;     
    }
  }
  
}

void  Cal_Grid_Weight(double x1, double y1, double x0, double y0, int Nx, int Ny, double hx,  double hy, int xmesh[], int ymesh[], double xweight[], double yweight[])
{
  int i,j;
  double dx=x1-x0, dy=y1-y0;
  double fx, fy;
  int    ifx, ify;
  double rx, rx2, ry, ry2;
  double flag1=0, flag2=0;
  
  fx= dx/hx;   ifx=int(fx);
  fy= dy/hy;   ify=int(fy);
  
  if( (fx-ifx) <= 0.5 ){    
    xmesh[0]=ifx; xmesh[1]=ifx+1; xmesh[2]=ifx-1; 
    rx = fx - ifx;  rx2 = rx * rx;
    }
  else{
    xmesh[0]=ifx+1; xmesh[1]=ifx; xmesh[2]=ifx+2; 
    rx = ifx + 1 - fx;    rx2=rx*rx ;
  }
  xweight[0]= 0.75 - rx2; 
  xweight[1]= 0.5*(0.25+rx+rx2); 
  xweight[2]= 0.5*(0.25-rx+rx2);  

  if( (fy-ify) <= 0.5 ){    
    ymesh[0]=ify; ymesh[1]=ify+1; ymesh[2]=ify-1; 
    ry = fy - ify;  ry2 = ry * ry;
    }
  else{
    ymesh[0]=ify+1; ymesh[1]=ify; ymesh[2]=ify+2; 
    ry = ify + 1 - fy;    ry2=ry*ry ;
  }
  yweight[0]= 0.75 - ry2; 
  yweight[1]= 0.5*(0.25+ry+ry2); 
  yweight[2]= 0.5*(0.25-ry+ry2); 
  
}

void  Charge2Grid(double x, double y, double x0, double y0, int Nx, int Ny, double hx,  double hy, double src[])
{
  int i,j;
  int xmesh[3], ymesh[3];
  double xweight[3], yweight[3];

  Cal_Grid_Weight(x,y,x0,y0,Nx,Ny,hx,hy,xmesh,ymesh,xweight,yweight);
  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      src[ xmesh[i] *Nx + ymesh[j] ]  =  src[ xmesh[i] *Nx + ymesh[j] ] + xweight[i] * yweight[j] ;
    }
  }

}

double int_green(double x, double y)  
{
     double rslt;
     double x2, y2;

     if( 0 == x || 0 == y ){
          rslt = 0;
     }else{
          x2 = x*x;
          y2 = y*y;
          rslt = -3*x*y + x2*atan(y/x) + y2*atan(x/y) + x*y*log(x2+y2);
     }
     return rslt;
}

double calc_egreen(double x1, double x2, double y1, double y2, double hx, double hy)
{
  int i,j;
  double tempx, tempy;

  int   sign_x1y1, sign_x1y2, sign_x2y1, sign_x2y2;
  double val_x1y1, val_x1y2, val_x2y1, val_x2y2;
  double rslt;
  

  if(x1*x2*y1*y2 != 0.) {
    
    sign_x1y1 = 1;
    sign_x1y2 = 1;
    sign_x2y1 = 1;
    sign_x2y2 = 1;
    
    if( x1 > 0 ){
      sign_x1y1 *= -1;
      sign_x1y2 *= -1;
    }
    else if( x1 < 0 && x2 < 0 ){
      sign_x2y1 *= -1;
      sign_x2y2 *= -1;
    }
    
    if( y1 > 0 ){
      sign_x1y1 *= -1;
      sign_x2y1 *= -1;
    }
    else if( y1 < 0 && y2 < 0 ){
      sign_x1y2 *= -1;
      sign_x2y2 *= -1;
    }
    
    x1 = fabs(x1); x2 = fabs(x2);  y1 = fabs(y1);    y2 = fabs(y2);
    val_x1y1 = int_green( x1, y1 );   val_x1y2 = int_green( x1, y2 );
    val_x2y1 = int_green( x2, y1 );   val_x2y2 = int_green( x2, y2 );
    rslt = sign_x1y1*val_x1y1 + sign_x1y2*val_x1y2 + sign_x2y1*val_x2y1 + sign_x2y2*val_x2y2;
    
    return -1.0*rslt/hx/hy/2;
  }
  else{
    rslt=0.;
    for(i=0;i<20;i++){
      for(j=0;j<20;j++){ 
	tempx=x1+ (i-1)*hx/20 + hx/40;
	tempy=y1+ (j-1)*hy/20 + hy/40;
        if(tempx !=0. and tempy !=0. ){
	  rslt=rslt  - 0.5 * log(tempx*tempx + tempy*tempy) ; 
	}
      }
    }
    return rslt/400.;
  }
  
}

void Cal_Effective_Green_Function(double xoffset, double yoffset, int Nx, int Ny, double hx, double hy, double grn_c_r[], double grn_c_i[] )
{
  int i,j;
  double x1, y1,x2,y2;
  
  for(i=0;i<2*Nx*2*Ny;i++) {
    grn_c_r[i]=0.;
    grn_c_i[i]=0.;
  }

  for(i=0;i<2*Nx;i++){
    for(j=0;j<2*Ny;j++){
      
      if( (i>=0 and i<=Nx-1) and   (j>=0 and j<=Ny-1) ){ 
	x1=xoffset + i*hx - 0.5*hx;   x2=x1 + hx ; 
	y1=yoffset + j*hy - 0.5*hy;   y2=y1 + hy;
	grn_c_r[i*2*Nx+j]=calc_egreen(x1, x2,  y1, y2, hx, hy);
      }
      else if((i>=0 and i<=Nx-1) and   (j>Ny-1 and j<=2*Ny-1)  ){ //before :(i>=0 and i<=Nx-1) and   (j>Ny and j<=2*Ny-1) ,changed on June 18, 2024
	x1=xoffset + i*hx - 0.5*hx;        x2=x1 + hx;
	y1=yoffset - (2*Ny-j)*hy-0.5*hy ;  y2=y1+ hy;
	grn_c_r[i*2*Nx+j]=calc_egreen(x1, x2,  y1, y2, hx, hy);
      }
      else if((i>Nx-1 and i<=2*Nx-1) and  (j>=0 and j<=Ny-1) ){ // before: (i>Nx and i<=2*Nx-1) and  (j>=0 and j<=Ny-1),changed on June 18, 2024
       	x1=xoffset - (2*Nx-i)*hx - 0.5*hx;  x2=x1 + hx ; 
	y1=yoffset + j*hy - 0.5*hy;         y2=y1+ hy;
	grn_c_r[i*2*Nx+j]=calc_egreen(x1, x2,  y1, y2, hx, hy);
      }
      else if((i>Nx-1 and i<=2*Nx-1) and  (j>Ny-1 and j<=2*Ny-1) ){ // before: (i>Nx and i<=2*Nx-1) and  (j>Ny and j<=2*Ny-1), changed on June 18, 2024
	x1=xoffset - (2*Nx-i)*hx - 0.5*hx;  x2=x1 + hx;
	y1=yoffset - (2*Ny-j)*hy - 0.5*hy ; y2=y1 + hy;
	grn_c_r[i*2*Nx+j]=calc_egreen(x1, x2,  y1, y2, hx, hy);
      }
      else{
      }
      
    }
  }
  
  FFT2D(grn_c_r, grn_c_i, 2*Nx, 2*Ny,1);
}

void Cal_Potential_Direct(double xoffset, double yoffset, int Nx, int Ny, double hx, double hy, double src[],  double phi[])
{
  int m,n,i,j;
  double x1,y1, sum;
  double eps0=8.854187817e-12;
  
  for(m=0;m<Nx;m++){
    for(n=0;n<Ny;n++){
 
      sum=0;
      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	  if( src[i*Nx+j] != 0. ){
	    x1=xoffset + (m-i)*hx ;
	    y1=yoffset + (n-j)*hy ;
	    if (x1 == 0. and  y1 == 0. ) {
              sum=sum+ calc_egreen(-hx/2, hx/2, -hy/2, hy/2, hx, hy) * src[i*Nx+j];
	    }
	    else{
	      //sum=sum - 0.5  * log( x1*x1  + y1*y1)  * src[i*Nx+j];
	      sum=sum+ calc_egreen(x1-hx/2, x1+hx/2, y1-hy/2, y1+hy/2, hx, hy) * src[i*Nx+j];
	    }
	  }
	}
      }
      phi[m*Nx+n] = sum/2/3.14159265/eps0;

    }
  }
  
}
  
void Cal_Potential_Convolution(int Nx, int Ny, double hx, double hy, double src[], double grn_c_r[], double grn_c_i[], double phi[])
{
  int i,j;
  double src_c_r[2*Nx*2*Ny],  src_c_i[2*Nx*2*Ny], phi_c_r[2*Nx*2*Ny], phi_c_i[2*Nx*2*Ny];
  double eps0=8.854187817e-12;

  for(i=0;i<2*Nx*2*Ny;i++) {
    src_c_r[i]=0.; src_c_i[i]=0; phi_c_r[i]=0.; phi_c_i[i]=0.;
  }
  
  for(i=0;i<2*Nx;i++){
    for(j=0;j<2*Ny;j++){
      if(  (i>=0 and i<=Nx-1)  and  (j>=0 and j<=Ny-1 ) ) { 
	src_c_r[i*2*Nx +j]=src[i*Nx +j];
      }
      else{
	src_c_r[i*2*Nx +j]=0.;
      }
    }
  }
  
  FFT2D(src_c_r, src_c_i, 2*Nx, 2*Ny,1);
  for(i=0;i<2*Nx*2*Ny;i++){
    phi_c_r[i]=src_c_r[i]*grn_c_r[i]-src_c_i[i]*grn_c_i[i];
    phi_c_i[i]=src_c_r[i]*grn_c_i[i]+src_c_i[i]*grn_c_r[i]; 
  }
  FFT2D(phi_c_r,phi_c_i, 2*Nx, 2*Ny,-1);  
  
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      phi[i*Nx+j] = phi_c_r[i*2*Nx+j]*2*Nx*2*Ny/2/3.14159265/eps0;
    }
  }
  
}

void Cal_Potential_Convolution2(int Nx, int Ny, double hx, double hy, double xoffset, double yoffset, double src[], double phi[])
{
  double  grn_c_r[2*Nx*2*Ny], grn_c_i[2*Nx*2*Ny];
  Cal_Effective_Green_Function(xoffset, yoffset, Nx, Ny, hx, hy, grn_c_r, grn_c_i);
  Cal_Potential_Convolution(Nx, Ny, hx, hy, src, grn_c_r, grn_c_i, phi);
}

void Cal_Electric_Field(int Nx, int Ny, double hx, double hy, double phi[], double Ex[], double  Ey[])
{
  int i, j;
  
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      Ex[i*Nx+j]=0.;   Ey[i*Nx+j]=0.;
    }
  }
  
  for(i=1;i<Nx-1;i++){
    for(j=1;j<Ny-1;j++){
      Ex[i*Nx+j] =-(  ( phi[(i+1)*Nx+j+1] - phi[(i-1)*Nx+j+1]) 
                 + 4 *( phi[(i+1)*Nx+j]   - phi[(i-1)*Nx+j]  )  
                    + ( phi[(i+1)*Nx+j-1] - phi[(i-1)*Nx+j-1])  )/12./hx;  
      Ey[i*Nx+j] =-(  ( phi[(i+1)*Nx+j+1] - phi[(i+1)*Nx+j-1]) 
                 + 4 *( phi[i*Nx+j+1]     - phi[i*Nx+j-1]    )  
                    + ( phi[(i-1)*Nx+j+1] - phi[(i-1)*Nx+j-1] ) )/12./hy; 
    }
  }  
    
}

void Interpolate_Electric_Field(int Nx, int Ny, double hx, double hy, double x0, double y0,  double Ex[], double  Ey[], double x, double y, double & ex, double & ey)
{
  int i,j;
  int xmesh[3], ymesh[3];
  double xweight[3], yweight[3];
  
  Cal_Grid_Weight(x,y, x0,y0,Nx,Ny,hx,hy,xmesh,ymesh,xweight,yweight);
  ex=0.;  ey=0.;
  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
     ex= ex+ Ex[ xmesh[i] *Nx + ymesh[j] ] * xweight[i] * yweight[j] ;
     ey= ey+ Ey[ xmesh[i] *Nx + ymesh[j] ] * xweight[i] * yweight[j] ;
    }
  }
   
}

void Cal_Electric_PIC(int Nx, int Ny, double x0, double y0, double xoffset, double yoffset, double hx, double hy, double src[], double Ex[], double  Ey[])
{
  double  grn_c_r[2*Nx*2*Ny], grn_c_i[2*Nx*2*Ny], phi[Nx*Ny];

  Cal_Effective_Green_Function(xoffset, yoffset, Nx, Ny, hx, hy, grn_c_r, grn_c_i);
  Cal_Potential_Convolution(Nx, Ny, hx, hy, src, grn_c_r, grn_c_i, phi);
  Cal_Electric_Field(Nx, Ny, hx, hy, phi, Ex, Ey);
}

void Interpolate_Electric_Field_1Potential(int Nx, int Ny, double hx, double hy, double x0, double y0,   double phi[] , double x, double y, double & ex, double & ey)
{
  int i,j;
  int    xmesh[3], ymesh[3];
  int    index0, indey0;
  double xweight[3], yweight[3];
  double phi25[25], Ex25[25], Ey25[25];
  
  Cal_Grid_Weight(x,y,x0,y0,Nx,Ny,hx,hy,xmesh,ymesh,xweight,yweight);
  
  index0=10000;
  indey0=10000;
  for(i=0;i<3;i++){
    if( xmesh[i] < index0 ) index0 =  xmesh[i];
    if( ymesh[i] < indey0 ) indey0 =  ymesh[i];    
  } 

  for(i=0;i<25;i++) {
    phi25[i]=0; Ex25[i]=0; Ey25[i]=0;
  }
  
  for(i=0;i<5;i++){
    for(j=0;j<5;j++){
      phi25[5*i+j]=phi[Nx*(index0-1+i)+(indey0-1+j)];
    }
  }
  
  for(i=1;i<4;i++){
    for(j=1;j<4;j++){
      Ex25[i*5+j] =-(   ( phi25[(i+1)*5+j+1] - phi25[(i-1)*5+j+1]) 
			+ 4 *( phi25[(i+1)*5+j]   - phi25[(i-1)*5+j]  )  
			+ ( phi25[(i+1)*5+j-1] - phi25[(i-1)*5+j-1])  )/12./hx;  
      
      Ey25[i*5+j] =-(  ( phi25[(i+1)*5+j+1] -  phi25[(i+1)*5+j-1]) 
		       + 4 *( phi25[i*5+j+1]     - phi25[i*5+j-1]    )  
		       + ( phi25[(i-1)*5+j+1] - phi25[(i-1)*5+j-1] ) )/12./hy; 
    }
  }   
  
  ex=0.;  ey=0.;
  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      ex= ex+ Ex25[ (xmesh[i]-index0+1) *5 + (ymesh[j]-indey0+1)  ] * xweight[i] * yweight[j] ;
      ey= ey+ Ey25[ (xmesh[i]-index0+1) *5 + (ymesh[j]-indey0+1)  ] * xweight[i] * yweight[j] ;
    }
  }
   
}

void Interpolate_Electric_Field_2Potential(int Nx, int Ny, double hx, double hy, double x0, double y0, double sfront, double sback,  double phi1[], double  phi2[],  double x, double y, double sp, double & ex, double & ey)
{
  int i,j;
  int xmesh[3], ymesh[3];
  int index0, indey0;
  double xweight[3], yweight[3];
  double phi_front, phi_back;
  double phi25[25], Ex25[25], Ey25[25];
  
  Cal_Grid_Weight(x,y,x0,y0,Nx,Ny,hx,hy,xmesh,ymesh,xweight,yweight);
  
  index0=10000;
  indey0=10000;
  for(i=0;i<3;i++){
    if( xmesh[i] < index0 ) index0 =  xmesh[i];
    if( ymesh[i] < indey0 ) indey0 =  ymesh[i];    
  } 

  for(i=0;i<25;i++) {
    phi25[i]=0; Ex25[i]=0; Ey25[i]=0;
  }
  
  if(sfront == sback ){
    for(i=0;i<5;i++){
      for(j=0;j<5;j++){
	phi25[5*i+j]=phi1[Nx*(index0-1+i)+(indey0-1+j)];
      }
    }
  }
  else{
    for(i=0;i<5;i++){
      for(j=0;j<5;j++){
	phi_front= phi1[Nx*(index0-1+i)+(indey0-1+j)];
	phi_back = phi2[Nx*(index0-1+i)+(indey0-1+j)];
	phi25[5*i+j]=phi_back + (phi_front-phi_back)*(sp-sback)/(sfront-sback) ;
      }
    }
  }
    
  for(i=1;i<4;i++){
    for(j=1;j<4;j++){
      Ex25[i*5+j] =-(   ( phi25[(i+1)*5+j+1] - phi25[(i-1)*5+j+1]) 
                      + 4 *( phi25[(i+1)*5+j]   - phi25[(i-1)*5+j]  )  
                      + ( phi25[(i+1)*5+j-1] - phi25[(i-1)*5+j-1])  )/12./hx;  

      Ey25[i*5+j] =-(  ( phi25[(i+1)*5+j+1] -  phi25[(i+1)*5+j-1]) 
                     + 4 *( phi25[i*5+j+1]     - phi25[i*5+j-1]    )  
                     + ( phi25[(i-1)*5+j+1] - phi25[(i-1)*5+j-1] ) )/12./hy; 
    }
  }   
  
  ex=0.;  ey=0.;
  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      ex= ex+ Ex25[ (xmesh[i]-index0+1) *5 + (ymesh[j]-indey0+1)  ] * xweight[i] * yweight[j] ;
      ey= ey+ Ey25[ (xmesh[i]-index0+1) *5 + (ymesh[j]-indey0+1)  ] * xweight[i] * yweight[j] ;
    }
  }
   
}

void Zboundary_slices_Gaussian(double zb[], int Nslice, double sigmal)
{
  int i;

  zb[0]=  sigmal * 20;
  zb[Nslice]= -sigmal * 20;
  
  for(i=1;i<=Nslice-1;i++){
    zb[i]= -gsl_cdf_ugaussian_Pinv(  i*1.0/Nslice   ) * sigmal;
  }

}

tps::tps() 
{
  int i;
  for (i=0;i<7;i++) sample[i]=0.;  
}

tps::tps(double d )
{
  int i;
  for (i=0;i<7;i++) sample[i]=0;
  sample[0]=d;  
}

tps::tps(double d[7] )
{
  int i;
  for (i=0;i<7;i++) sample[i]=d[i];  
}

double & tps::operator[](int i)
{ 
  return sample[i]; 
}

tps operator-(tps x)
{
  int i;
  tps z;
  
  for (i = 0; i < 7; i++)  z[i] = -1.0 * x[i];
  return z;
} 

ostream &operator<<(ostream &stream, tps x)
{
  int i;
  for(i=0;i<7;i++) stream<<setw(12)<<setprecision(9)<<scientific<<x[i]<<"  ";
  return stream;
}

tps operator+(tps x, tps y)
{
  int i;
  tps z;
  
  for (i = 0; i < 7; i++)  z[i] = x[i] + y[i];
  return z;
}

tps operator-(tps x, tps y)
{
  int i;
  tps z;
  
  for (i = 0; i < 7; i++)  z[i] = x[i] - y[i];
  return z;
}

tps operator*(tps x, tps y)
{
  int i;
  tps z;
  
  z[0] = x[0] * y[0];
  for (i = 1; i < 7; i++)
    z[i] = x[0] * y[i] + x[i] * y[0];
  return z;
}

tps DAinv(tps x)
{
  int i;
  double a, temp;
  tps z;      
  
  z[0] = 1.0 / x[0];
  temp = x[0];
  a = -1.0 / (temp * temp);
  for (i = 1; i < 7; i++) z[i] = a * x[i];
  return z;
}

tps operator/(tps x, tps y)
{
  return x*DAinv(y);
}

tps sqr(tps x)
{
  return x*x;
}

tps sqrt(tps x)
{
  int i;
  double a;
  tps  z;
  
  a = sqrt(x[0]);
  z[0] = a;
  a = 0.5 / a;
  for (i = 1; i < 7; i++) z[i] = a * x[i];
  return z;
}

tps sin(tps x)
{
  int i;
  double a;
  tps z;
  
  z[0] = sin(x[0]); 
  a =    cos(x[0]);
  for (i = 1; i < 7; i++) z[i] = a * x[i];
  return z;
}

tps asin(tps x)
{
  int  i;
  double a;
  tps z;
  
  a = x[0];
  z[0] = asin(a);
  a = 1. / sqrt(1. - a * a);
  for (i = 1; i < 7; i++)  z[i] = a * x[i];
  return z;
}

tps cos(tps x)
{
  int i;
  double a;
  tps z;
  
  z[0] = cos(x[0]); 
  a =   -sin(x[0]);
  for (i = 1; i < 7; i++) z[i] = a * x[i];
  return z;
}

tps tan(tps x)
{
  return sin(x)/cos(x);
}

tps atan(tps x)
{
  int  i;
  double a;
  tps z;
  
  a = x[0];
  z[0] = atan(a);
  a = 1 / (1 + a * a);
  for (i = 1; i < 7; i++)  z[i] = a * x[i];
  return z;
}

tps sinh(tps x)
{
  int i;
  double a;
  tps z;
  
  z[0] = sinh(x[0]); 
  a = cosh(x[0]);
  for (i = 1; i < 7; i++) z[i] = a * x[i];
  return z;
}

tps cosh(tps x)
{
  int i;
  double a;
  tps z;
  
  z[0] = cosh(x[0]); 
  a = sinh(x[0]);
  for (i = 1; i <7; i++) z[i] = a * x[i];
  return z;
}

tps exp(tps x)
{
  int i;
  double a;
  tps z;
  
  a = exp(x[0]);
  z[0] = a;
  for (i = 1; i <7; i++) z[i] = a * x[i];
  return z;
}

tps ln(tps x)
{
  int i;
  tps z;
  
  z[0] = log(x[0]);
  for (i = 1; i < 7; i++) z[i] = x[i] / x[0];
  return z;
} 

linmap:: linmap()
{
  int i,j;
  for(i=0;i<6;i++)
    for(j=0;j<7;j++) map0[i][j]=0.;
}

void linmap::identity()
{
  int i,j;
  for(i=0;i<6;i++)
    for(j=0;j<7;j++) map0[i][j]=0.;
  map0[0][1]=1;
  map0[1][2]=1;
  map0[2][3]=1;
  map0[3][4]=1;
  map0[4][5]=1;
  map0[5][6]=1;
}

linmap::linmap(double x[6])
{
  int i,j;
  for(i=0;i<6;i++)
    for(j=0;j<7;j++) map0[i][j]=0.;
  for(i=0;i<6;i++) map0[i][0]=x[i];
}

tps & linmap::operator[](int i) 
{ 
  return map0[i]; 
}

linmap operator+(linmap x, linmap y)
{
  int i,j;
  linmap z;
  for(i=0;i<6;i++)
    for(j=0;j<7;j++) z[i][j]=x[i][j]+y[i][j];
  return z;
}

void linmap::print()
{
  int i;
  for(i=0;i<6;i++)  cout<<map0[i]<<endl;  
}

void Getmat(linmap map0, double x[36])
{
  int i,j;
  for (i=0;i<6;i++)
    for(j=0;j<6;j++) x[i*6+j]=map0[i][j+1];
}

void Getpos(linmap map0, double x[6])
{
  int i;
  for (i=0;i<6;i++)  x[i]=map0[i][0];
}

//=========================================
//
//         particle transfer functions
//
//=========================================

template<class T>  void  GtoL1(T x[], double DX, double DY, double DT)
// at entrance: only take  DX, DY, DT alignment erros
// DT  is rotation angle with  s as axis
{
  int i;
  T xtemp[6];
  
  if ( abs(DX) +  abs(DY)  > 1.0e-10 ) {
    x[0] = x[0] - DX;   x[2] = x[2] - DY;
  }
  if (  abs(DT)  > 1.0e-10 ) {
    double cosT=cos(DT), sinT=sin(DT);
    for (i=0; i<6;i++) xtemp[i]=x[i];
    x[0] = cosT * xtemp[0] + sinT * xtemp[2];
    x[1] = cosT * xtemp[1] + sinT * xtemp[3];
    x[2] = cosT * xtemp[2] - sinT * xtemp[0];
    x[3] = cosT * xtemp[3] - sinT * xtemp[1];
  }
}

template<class T>  void LtoG1(T x[], double DX, double DY, double DT)
// at exit: only take  DX, DY, DT alignment erros
// DT  is rotation angle with  s as axis
{
  int i;
  T xtemp[6];

  if (  abs(DT)  > 1.0e-10 ) {
    double cosT=cos(DT), sinT=sin(DT);
    for(i=0; i<6;i++) xtemp[i]=x[i];
    x[0] = cosT * xtemp[0] - sinT * xtemp[2];
    x[1] = cosT * xtemp[1] - sinT * xtemp[3];
    x[2] = sinT * xtemp[0] + cosT * xtemp[2];
    x[3] = sinT * xtemp[1] + cosT * xtemp[3];
  }
  if ( abs(DX) +  abs(DY)  > 1.0e-10 ) {
    x[0] = x[0] + DX  ;  x[2] = x[2]+ DY;
  }
}

//-----more general misalignment treatment: just like madx notation

void  Cal_W_Matrix(double theta, double phi, double psi, double W[]) 
// calculate rotation matrix  in 3-D
// equivalent to  subrotine sumtrx() in madx 
{
  double  cosphi, cospsi, costhe, sinphi, sinpsi, sinthe;

  costhe = cos(theta);   cosphi = cos(phi);   cospsi = cos(psi);
  sinthe = sin(theta);   sinphi = sin(phi);   sinpsi = sin(psi);
  
  W[0] = + costhe * cospsi - sinthe * sinphi * sinpsi;
  W[1] = - costhe * sinpsi - sinthe * sinphi * cospsi;
  W[2] =                     sinthe * cosphi;
  W[3] =                              cosphi * sinpsi;
  W[4] =                              cosphi * cospsi;
  W[5] =                              sinphi;
  W[6] = - sinthe * cospsi - costhe * sinphi * sinpsi;
  W[7] = + sinthe * sinpsi - costhe * sinphi * cospsi;  
  W[8] =                     costhe * cosphi;
}

void Cal_VE_WE_Matrx(double length, double angle, double tilt, double VE[], double WE[])
//   applied  to straight and sbends, not to reference coordinate changing elements such as  TRANS, YROTAT, SROTAT.
//   equivalent to subroutine suelem()  in madx, psi->tilt, angle->theta
//   VE[] and WE[] are displacement and rotations of exit w.r.t entrance of an element 
{
  int i;
  double temp;
  double cospsi, costhe, sinpsi, sinthe;
  
  //--- straight element
  if(  abs(angle) < 1.0e-12 ) {
    
    for(i=0;i<3;i++)  VE[i]=0.;
    VE[2]=length;
    
    for(i=0;i<9;i++)  WE[i] =0.;
    for(i=0;i<3;i++)  WE[i*3+i] = 1.0;
    
  }

  //----sbend element : no tilt, only with horizontal bending 'angle'
  if( abs(angle) > 1.0e-12 ) {

    tilt = 0. ;
 
    cospsi = cos(tilt);  sinpsi = sin(tilt);
    costhe = cos(angle); sinthe = sin(angle);

    /*-------
    VE[0] = length * (costhe-1.0)/angle ;
    VE[1] = 0; 
    VE[2] = length * sin(angle)/angle;

    WE[0] =  costhe;
    WE[1] =  0.;
    WE[2] =  sinthe;
    WE[3] =  0.;
    WE[4] =  1.;
    WE[5] =  0.;
    WE[6] =  -sinthe;   
    WE[7] =  0.;
    WE[8] =  costhe;
    -----*/
    
    temp  = length * (costhe-1.0)/angle ;
    VE[0] = temp * cospsi;
    VE[1] = temp * sinpsi;
    VE[2] = length * sin(angle)/angle; 

    WE[0*3+0] =  costhe * cospsi*cospsi + sinpsi*sinpsi;  //  costhe
    WE[1*3+0] =  (costhe - 1.0) * cospsi * sinpsi;        //  0 
    WE[2*3+0] =  sinthe * cospsi;                         //  sinthe
    WE[0*3+1] =  WE[1*3+0];                               //   0
    WE[1*3+1] =  costhe * sinpsi*sinpsi + cospsi*cospsi;  //   1
    WE[2*3+1] =  sinthe * sinpsi;                         //   0
    WE[0*3+2] = -WE[2*3+0];                               //  -sinthe
    WE[1*3+2] = -WE[2*3+1];                               //   0
    WE[2*3+2] =  costhe;                                  //   costhe
  }

}

template<class T>  void  GtoL2(T x[], double  length, double angle, double dx, double dy, double ds, double theta, double phi, double psi)
//  at the entrance of  an element
//  equivalent  to subroutine tmali1() in madx
{
  int    i,j;
  double s2, W[9], RM[36];
  T      xtemp[6];

  Cal_W_Matrix(theta, phi, psi, W);
  
  s2=( W[0*3+2] * dx  + W[1*3+2] * dy  +   W[2*3+2] * ds ) / W[ 2*3+2] ;  

  for(i=0;i<36;i++)  RM[i]=0.;
  for(i=0;i<6; i++)  RM[i*6+i]=1.0;
  
  RM[1*6+1] = W[0*3+0];
  RM[1*6+3] = W[1*3+0];
  RM[1*6+5] = W[2*3+0] / GP.beta;
  
  RM[3*6+1] = W[0*3+1];
  RM[3*6+3] = W[1*3+1];
  RM[3*6+5] = W[2*3+1] / GP.beta;
  
  RM[0*6+0] = W[1*3+1] / W[2*3+2];
  RM[0*6+1] = RM[0*6+0] * s2;
  RM[0*6+2] = -W[0*3+1] / W[2*3+2];
  RM[0*6+3] = RM[0*6+2] * s2;
  
  RM[2*6+0] = -W[1*3+0] / W[2*3+2];
  RM[2*6+1] = RM[2*6+0] * s2;
  RM[2*6+2] =  W[0*3+0] / W[2*3+2];
  RM[2*6+3] = RM[2*6+2] * s2;
  
  RM[4*6+0] = W[0*3+2] / (W[2*3+2] * GP.beta);
  RM[4*6+1] = RM[4*6+0] * s2;
  RM[4*6+2] = W[1*3+2] / (W[2*3+2] * GP.beta);
  RM[4*6+3] = RM[4*6+2] * s2;
  RM[4*6+5] = -s2 / (GP.beta * GP.gamma * GP.beta * GP.gamma); 

  for(i=0;i<6;i++){
    xtemp[i] = 0.;
    for(j=0;j<6;j++){
      xtemp[i] =  xtemp[i] + RM[i*6+j]*x[j];
    }
  }

  x[0]=xtemp[0] - (W[1*3+1] * dx  - W[0*3+1] * dy ) / W[ 2*3+2] ;
  x[1]=xtemp[1] +  W[2*3+0] ;
  x[2]=xtemp[2] - (W[0*3+0] * dy  - W[1*3+0] * dx ) / W[ 2*3+2] ;
  x[3]=xtemp[3] +  W[2*3+1] ;
  x[4]=xtemp[4] - s2/GP.beta;   
  x[5]=xtemp[5];

}

template<class T>  void LtoG2(T x[], double length, double angle, double dx, double dy, double ds, double theta, double phi, double psi)
//  at the  exit of  an element
//  equivalent  to subroutine tmali2() in madx  
{
  int i,j;
  double V[3], W[9];
  double VE[3],WE[9];
  double V1[3], M1[9], M2[9];
  double RM[36];
  double s2;
  double tilt = 0. ;
  T      xtemp[6];
  
  Cal_W_Matrix(theta, phi, psi, W);
  Cal_VE_WE_Matrx(length, angle, tilt, VE, WE);
  
  V[0] = dx + W[0*3+0]*VE[0]+W[0*3+1]*VE[1]+W[0*3+2]*VE[2]-VE[0];
  V[1] = dy + W[1*3+0]*VE[0]+W[1*3+1]*VE[1]+W[1*3+2]*VE[2]-VE[1];  
  V[2] = ds + W[2*3+0]*VE[0]+W[2*3+1]*VE[1]+W[2*3+2]*VE[2]-VE[2];

  //V = matmul(transpose(WE),V);
  mat_transpose(WE, M1, 3, 3);
  for(i=0;i<3;i++){
    V1[i]=0;
    for(j=0;j<3;j++){
      V1[i]=V1[i]+ M1[i*3+j]*V[j];
    }
  }
  for(i=0;i<3;i++) V[i]=V1[i];
  
  //W = matmul(matmul(transpose(WE),W),WE);
  //mat_transpose(WE, M1, 3, 3);
  mat_mult(M1, W, M2, 3, 3, 3);
  mat_mult(M2, WE, W, 3, 3, 3);
  s2 = -(W[0*3+2] * V[0] + W[1*3+2] * V[1] + W[2*3+2] * V[2]) / W[2*3+2]; 

  for(i=0;i<36;i++)  RM[i]=0.;
  for(i=0;i<6; i++)  RM[i*6+i]=1.0;

  RM[0*6+0] = W[0*3+0];
  RM[2*6+0] = W[1*3+0];     
  RM[4*6+0] = W[2*3+0] / GP.beta;
  
  RM[0*6+2] = W[0*3+1];      
  RM[2*6+2] = W[1*3+1];
  RM[4*6+2] = W[2*3+1] / GP.beta;
  
  RM[1*6+1] = W[1*3+1] / W[2*3+2];
  RM[0*6+1] = RM[1*6+1] * s2;
  RM[3*6+1] = -W[0*3+1] / W[2*3+2];
  RM[2*6+1] = RM[3*6+1] * s2;
  
  RM[1*6+3] = -W[1*3+0] / W[2*3+2];
  RM[0*6+3] = RM[1*6+3] * s2;
  RM[3*6+3] =  W[0*3+0] / W[2*3+2];
  RM[2*6+3] = RM[3*6+3] * s2;
  
  RM[1*6+5] =  W[0*3+2] / (W[2*3+2] * GP.beta);
  RM[0*6+5] =  RM[1*6+5] * s2;  
  RM[3*6+5] =  W[1*3+2] / (W[2*3+2] * GP.beta);
  RM[2*6+5] =  RM[3*6+5] * s2;
  RM[4*6+5] = -s2 / (GP.beta * GP.gamma * GP.beta * GP.gamma);

  xtemp[0]=x[0] + (W[1*3+1] * V[0]  - W[0*3+1] * V[1] ) / W[ 2*3+2] ;
  xtemp[1]=x[1] -  W[2*3+0] ;
  xtemp[2]=x[2] + (W[0*3+0] * V[1]  - W[1*3+0] * V[0] ) / W[ 2*3+2] ;  
  xtemp[3]=x[3] -  W[2*3+1] ;
  xtemp[4]=x[4] -  s2/GP.beta;  
  xtemp[5]=x[5];

  for(i=0;i<6;i++){
    x[i] = 0.;
    for(j=0;j<6;j++){
      x[i] =  x[i] + RM[i*6+j]*xtemp[j];
    }
  }

}

template<class T>  void  GtoL(T x[], double length, double angle, double dx, double dy, double ds, double theta, double phi, double psi )
// at entrance
{
  double temp;

  //GtoL2(x, length, angle, dx, dy, ds, theta, phi, psi);
  
  temp = abs(dx) + abs(dy) + abs(ds) + abs(theta) + abs(phi) + abs(psi);
  if( temp  < 1.0e-9 ){                   //   no misalignment  
    return ;
  }
  
  if( abs(theta) + abs(phi) < 1.0e-9) {   //  no X   or  Y rotation
    GtoL1(x, dx, dy, psi);
  }
  else{                                   //  with X   or Y-rotation
    GtoL2(x ,length, angle, dx, dy, ds, theta, phi, psi);
  }

}

template<class T>  void LtoG(T x[], double length, double angle, double dx, double dy, double ds, double theta, double phi, double psi)
//  at exit
{
  double temp;
  
  //LtoG2(x, length, angle, dx, dy, ds, theta, phi, psi);
 
  temp = abs(dx) + abs(dy) + abs(ds) + abs(theta) + abs(phi) + abs(psi);
  if( temp  < 1.0e-9 ){                   //   no misalignment  
    return ;
  }
  
  if(  abs(theta) + abs(phi)  < 1.0e-9) {  // no X  or  Y rotation
    LtoG1(x, dx, dy, psi);
  }
  else{                                   // with X or Y-rotation 
    LtoG2(x ,length, angle, dx, dy, ds, theta, phi, psi);
  }
  
}

template <class T>  void DRIFT_Pass(T x[], double L)
{
  if( GP.H_expand == true){ 
    T u;
    u= L/(1.0+x[5]); 
    x[x_]=x[x_]+x[px_]*u;
    x[y_]=x[y_]+x[py_]*u;
    x[z_]=x[z_]-(x[px_]*x[px_]+x[py_]*x[py_])*u/2.0/(1+x[5]);
  }
  else{
    T pz, betaz;
    T delta1,gamma1,beta1;

    delta1= sqrt(1.0 + 2*x[pt_]/GP.beta+x[pt_]*x[pt_]) -1.0;
    gamma1 =GP.gamma+sqrt(GP.gamma*GP.gamma-1.)*x[pt_];
    beta1 = sqrt(1.0-1.0/gamma1/gamma1);
   
    pz = sqrt( (1.0+delta1)*(1.0+delta1)-x[px_]*x[px_]-x[py_]*x[py_]);
    betaz=pz*beta1/(1.0+delta1);
    x[x_]=x[x_]+x[px_]*L/pz; 
    x[y_]=x[y_]+x[py_]*L/pz;
    x[z_]=x[z_]-(L/betaz-L/GP.beta);
  }
}

template <class T>  void drift_polar_pass(T x[], double L, double href)
{
  int i;
  double angle=L*href, sinangle=sin(angle), cosangle=cos(angle), tanangle=tan(angle);
  T delat, pz, temp, x0[6];
  T delta1,gamma1,beta1;
  
  for(i=0;i<6;i++) x0[i]=x[i];
  
  delta1= sqrt(1.0 + 2*x[pt_]/GP.beta+x[pt_]*x[pt_]) -1.0;
  gamma1 =GP.gamma+sqrt(GP.gamma*GP.gamma-1.)*x[pt_];
  beta1 = sqrt(1.0-1.0/gamma1/gamma1);
  
  pz= sqrt( (1+delta1)*(1+delta1)-x0[px_]*x0[px_]-x0[py_]*x0[py_]);
  temp = 1-x0[px_]*tanangle/pz ;

  x[x_] =x0[x_] /cosangle /temp;
  x[px_]=x0[px_]*cosangle + pz * sinangle;
  x[y_] =x0[y_] + x0[py_]*x0[x_]*tanangle / pz / temp;
  x[z_]=x[z_] - ( (L+(1+delta1)*x0[x_]*tanangle/pz/temp )/beta1 - L/GP.beta);
}

template <class T>  void bend_kick_pass(T x[], double L, double href  )
{
  x[px_]=x[px_]+(href*x[pt_]-href*href*x[x_])*L;
  x[z_]=x[z_]-href*x[x_]*L;
}

template <class T>  void general_bend_kick_pass(T x[], double L, double href, double hreal  )
{
  x[px_]=x[px_] + ( -(hreal-href) + href*x[pt_] - href*hreal*x[x_] ) * L ;
  x[z_]=x[z_]-href*x[x_]*L;
}

template <class T>  void bend_kick_pass_exact(T x[], double L, double href  )  //  not  called  anywhere
{
  int i;
  T x0[9];
  T delta1,gamma1,beta1,pz;
  double temp1, temp2;
  T  temp;

  delta1= sqrt(1.0 + 2*x[pt_]/GP.beta+x[pt_]*x[pt_]) -1.0;
  gamma1 =GP.gamma+sqrt(GP.gamma*GP.gamma-1.)*x[pt_];
  beta1 = sqrt(1.0-1.0/gamma1/gamma1);
  
  temp1=href*L/2;
  x[px_]=x[px_] - href*temp1*x[x_] + temp1*delta1;
  x[z_] =x[z_] -  temp1*x[x_]/beta1; 
  
  for(i=0;i<9;i++) x0[i]=x[i];
  temp1= href*L;  temp2 = temp1*temp1 ; 
  temp= -temp2*x0[px_]*x0[px_] + 2*temp1*(1+delta1)*x0[px_];
  
  pz=sqrt((1+delta1)*(1+delta1)-x0[px_]*x0[px_]-x0[py_]*x0[py_]);
  x[px_]=(x[px_] + temp1 *(1+delta1)*( sqrt(1.0-(x0[px_]*x0[px_]+x0[py_]*x0[py_]-temp)/(1+delta1)/(1+delta1))-1) )/ (1+ temp2);
  x[x_]=x0[x_]+temp1*x0[x_]*x0[px_]/pz;
  x[y_]=x0[y_]+temp1*x0[x_]*x0[py_]/pz;
  x[z_]=x0[z_]-temp1*x0[x_]*((1+delta1)/pz-1)/beta1;
  
  temp1=href*L/2;
  x[px_]=x[px_] - href*temp1*x[x_] + temp1*delta1;
  x[z_] =x[z_] -  temp1*x[x_]/beta1; 
}

template <class T>  void sbend_exact_pass(T x[], double L, double Angle, double cosAngle, double sinAngle)
{
  int i;
  double href=Angle/L;
  T x0[6],pz0,pz1,delta,gamma, beta, derivative;
  T temp, arcsin1, arcsin2;

  for(i=0;i<6;i++) x0[i]=x[i];
  
  delta= sqrt(1.0 + 2*x[pt_]/GP.beta+x[pt_]*x[pt_]) -1.0;
  gamma =GP.gamma+sqrt(GP.gamma*GP.gamma-1)*x[pt_];
  beta = sqrt(1.0-1.0/gamma/gamma);
  
  pz0=sqrt( (1.+delta) * (1.+delta) -x0[1]*x0[1] -x0[3]*x0[3]);
  x[1]=x0[1]*cosAngle + (pz0-1.-href*x0[0])*sinAngle;
  
  pz1=sqrt( (1.+delta) * (1.+delta) -x[1]*x[1] -x[3]*x[3]);
  derivative=-1.0*x0[1]*sinAngle*href + (pz0-1.-x[0]*href)*cosAngle*href ;
  x[0]=pz1/href-derivative/href/href-1.0/href;
  
  temp=sqrt(  (1.+delta) * (1.+delta) -x0[3]*x0[3] );
  arcsin1=asin( x0[1]/temp);  arcsin2=asin( x[1]/temp);
  x[2]=x0[2]+ x0[3]*L  + x0[3]* ( arcsin1 -arcsin2)/href;
 
  temp= (1+delta)*L+(1+delta)*(arcsin1 -arcsin2)/href ;
  x[4]=x[4] - (temp/beta - L/GP.beta); 
}


//template <class T>  void sbend_exact_pass_v0(T x[], double L, double Angle)
void sbend_exact_pass_v0(double x[], double L, double Angle)
//-----sbend solution from geometric analysis by  Yun Luo
{
  int i;
  //T x1[6];
  //T angle1, rho1, rho;
  //T xo1, yo1, xp0, yp0, vx0, vy0, vt, vz0, pz0;
  //T pos_x0, pos_y0, pos_x1, pos_y1, pos_o1x, pos_o1y, theta;
  //T dt1, dt0;
  //T vx1, vz1, xp1,yp1,temp1,temp2;
  
  //T pz, betaz;
  //T delta1,gamma1,beta1;

  double  x1[6];
  double  angle1, rho1, rho;
  double  xo1, yo1, xp0, yp0, vx0, vy0, vt, vz0, pz0;
  double  pos_x0, pos_y0, pos_x1, pos_y1, pos_o1x, pos_o1y, theta;
  double  dt1, dt0;
  double  vx1, vz1, xp1,yp1,temp1,temp2;
  
  double  pz, betaz;
  double  delta1,gamma1,beta1;

  
  delta1= sqrt(1.0 + 2*x[pt_]/GP.beta+x[pt_]*x[pt_]) -1.0;
  gamma1 =GP.gamma+sqrt(GP.gamma*GP.gamma-1)*x[pt_];
  beta1 = sqrt(1.0-1.0/gamma1/gamma1);
  
  if(Angle > 0.) {
    pz0=sqrt( (1+delta1)*(1+delta1)-x[px_]*x[px_]-x[py_]*x[py_]);
    xp0=x[px_]/pz0;  yp0=x[py_]/pz0; 
    temp1=sqrt(xp0*xp0+yp0*yp0+1.0);
    vx0=beta1*xp0/temp1;  vy0=beta1*yp0/temp1; 
    vz0=beta1/temp1;      vt=sqrt(vx0*vx0+vz0*vz0);  
    
    rho=L/Angle;
    rho1=(1.+delta1)*rho*sqrt(xp0*xp0+1)/temp1;
    temp1=sin(Angle);
    temp2=sqrt(1.0- temp1*temp1);
    pos_x0=(rho + x[x_])*temp2;
    pos_y0=(rho + x[x_])*temp1;

    theta=atan(vx0/vz0);
    theta=theta+Angle;
    pos_o1x=pos_x0+rho1*cos(PI+theta);
    pos_o1y=pos_y0+rho1*sin(PI+theta);
    pos_x1=sqrt(rho1*rho1-pos_o1y*pos_o1y)+ pos_o1x;
    pos_y1=0.;
    
    temp1=sqrt((pos_x1-pos_x0)*(pos_x1-pos_x0) +(pos_y1-pos_y0)*(pos_y1-pos_y0));
    angle1=2.*asin(temp1/rho1/2.);
    dt1= angle1*rho1/vt;
    dt0= L/GP.beta;
    
    x[x_]=pos_x1-rho;
    x[y_]=x[y_]+vy0*dt1;
    x[z_]=x[z_]-(dt1-dt0);
    
    theta=asin(pos_o1y/rho1);
    temp1=sin(theta); temp2 =sqrt(1.-temp1*temp1);
    vx1=-vt*temp1;    vz1=vt*temp2; 
    xp1=-temp1/temp2; yp1=vy0/vz1;
    temp1=(1+delta1)*(1+delta1)/(1. + 1./(xp1*xp1+yp1*yp1));
    x[1]=xp1*sqrt((1+delta1)*(1+delta1)-temp1);
    x[3]=yp1*sqrt((1+delta1)*(1+delta1)-temp1); 
  }
  else if(Angle<0.){
    Angle=-Angle;
    pz0=sqrt( (1+delta1)*(1+delta1)-x[px_]*x[px_]-x[py_]*x[py_]);
    xp0=x[px_]/pz0;  yp0=x[py_]/pz0; 
    temp1=sqrt(xp0*xp0+yp0*yp0+1.0);
    vx0=beta1*xp0/temp1;  vy0=beta1*yp0/temp1; 
    vz0=beta1/temp1;      vt=sqrt(vx0*vx0+vz0*vz0);  
    
    rho=L/Angle;
    rho1=(1.+delta1)*rho*sqrt(xp0*xp0+1)/temp1;
    temp1=sin(Angle);
    temp2=sqrt(1.0- temp1*temp1);
    pos_x0=(rho - x[x_])*temp2;
    pos_y0=(rho - x[x_])*temp1;
    
    theta=atan(vx0/vz0);
    theta=Angle-theta;
    pos_o1x=pos_x0+rho1*cos(PI+theta);
    pos_o1y=pos_y0+rho1*sin(PI+theta);
    pos_x1=sqrt(rho1*rho1-pos_o1y*pos_o1y)+ pos_o1x;
    pos_y1=0.;
    
    temp1=sqrt((pos_x1-pos_x0)*(pos_x1-pos_x0) + pos_y0*pos_y0);
    angle1=2.*asin(temp1/rho1/2.);
    dt1= angle1*rho1/vt;
    dt0= L/GP.beta;
    
    x[x_]=-pos_x1+rho;
    x[y_]=x[y_]+vy0*dt1;
    x[z_]=x[z_]-(dt1-dt0);
    
    theta=asin(pos_o1y/rho1);
    temp1=sin(theta); temp2 =sqrt(1.-temp1*temp1);
    vx1=vt*temp1;     vz1=vt*temp2; 
    xp1=temp1/temp2;  yp1=vy0/vz1;
    temp1=(1+delta1)*(1+delta1)/(1. + 1./(xp1*xp1+yp1*yp1));
    x[px_]=xp1*sqrt((1+delta1)*(1+delta1)-temp1);
    x[py_]=yp1*sqrt((1+delta1)*(1+delta1)-temp1); 
  }
  else{
    DRIFT_Pass(x,L);
  }
}

template <class T>  void general_sbend_exact_pass(T x[], double L, double Angle, double K0L, double cosAngle, double sinAngle)
// K0L  = q * By * L / P_0 , Byis for the real dipole field
// L and Angle are for the frame
{
  int i;
  double rhoc=L/Angle;
  double b0  =K0L/L;

  T x0[6],pz0,pz1,delta,gamma, beta, derivative;
  T temp, arcsin1, arcsin2;

  for(i=0;i<6;i++) x0[i]=x[i];
  
  delta= sqrt(1.0 + 2*x[5]/GP.beta+x[5]*x[5]) - 1.0;
  gamma =GP.gamma+sqrt(GP.gamma*GP.gamma-1)*x[5];
  beta = sqrt(1.0-1.0/gamma/gamma);
  
  pz0=sqrt( (1.+delta) * (1.+delta) -x0[1]*x0[1] -x0[3]*x0[3]);
  x[1]=x0[1]*cosAngle + (pz0-b0*rhoc -b0* x0[0])*sinAngle;
  
  pz1=sqrt( (1.+delta) * (1.+delta) -x[1]*x[1] -x[3]*x[3]);
  derivative=-x0[1]*sinAngle/rhoc + (pz0-b0*rhoc -b0* x0[0])*cosAngle/rhoc ;
  x[0]=pz1/b0-derivative*rhoc/b0-rhoc;
  
  temp=sqrt(  (1.+delta) * (1.+delta) -x0[3]*x0[3] );
  arcsin1=asin( x0[1]/temp);  arcsin2=asin( x[1]/temp);
  x[2]=x0[2]+ x0[3]*L/b0/rhoc  + x0[3]* ( arcsin1 -arcsin2)/b0;
 
  temp= (1+delta)*L/b0/rhoc+(1+delta)*(arcsin1 -arcsin2)/b0;
  x[4]=x[4] - (temp/beta - L/GP.beta); 
}

void Cal_gsbend_l_angle(double x[6], double L, double Angle, double K0L, double & actual_l, double & actual_angle)
// return actual l and angle for a GSBEND
// x0[6]  is closed orbit through GSBEND
// K0L  = q * By * L / P_0 , is for the real dipole field, which may not equal to reference frame's .
// L and Angle are for the frame, not about real particle
{
  int i;
  double rhoc=L/Angle;
  double b0  =K0L/L;

  double x0[6],pz0,pz1,delta,gamma, beta, derivative;
  double temp, arcsin1, arcsin2;

  double  cosAngle= cos(Angle), sinAngle = sin(Angle);

  for(i=0;i<6;i++) x0[i]=x[i];
  
  delta= sqrt(1.0 + 2*x[5]/GP.beta+x[5]*x[5]) - 1.0;
  gamma =GP.gamma+sqrt(GP.gamma*GP.gamma-1)*x[5];
  beta = sqrt(1.0-1.0/gamma/gamma);
  
  pz0=sqrt( (1.+delta) * (1.+delta) -x0[1]*x0[1] -x0[3]*x0[3]);
  x[1]=x0[1]*cosAngle + (pz0-b0*rhoc -b0* x0[0])*sinAngle;
  
  pz1=sqrt( (1.+delta) * (1.+delta) -x[1]*x[1] -x[3]*x[3]);
  derivative=-x0[1]*sinAngle/rhoc + (pz0-b0*rhoc -b0* x0[0])*cosAngle/rhoc ;
  x[0]=pz1/b0-derivative*rhoc/b0-rhoc;
  
  temp=sqrt(  (1.+delta) * (1.+delta) -x0[3]*x0[3] );
  arcsin1=asin( x0[1]/temp);  arcsin2=asin( x[1]/temp);
  x[2]=x0[2]+ x0[3]*L/b0/rhoc  + x0[3]* ( arcsin1 -arcsin2)/b0;
 
  temp= (1+delta)*L/b0/rhoc+(1+delta)*(arcsin1 -arcsin2)/b0;
  x[4]=x[4] - (temp/beta - L/GP.beta);
  
  actual_l=temp;
  actual_angle= K0L  * actual_l / L ;
  
}

template <class T>  void quad_kick_pass(T x[], double k1l, double k1sl)
{
  T preal, pimag;
  T  By, Bx;

  preal = x[0];
  pimag = x[2];
  Bx    = k1sl * preal + k1l  * pimag;
  By    = k1l  * preal - k1sl * pimag;
  x[px_]= x[px_]-By;
  x[py_]= x[py_]+Bx;
 
}

template <class T>  void sext_kick_pass(T x[], double k2l, double k2sl)
{
  T preal, pimag;
  T  By, Bx;
  
  preal =  x[0]*x[0]- x[2]*x[2];
  pimag =  2. * x[0] * x[2];
  Bx    = ( k2sl * preal + k2l  * pimag )/2.;
  By    = ( k2l  * preal - k2sl * pimag )/2.;
  x[px_]= x[px_]-By;
  x[py_]= x[py_]+Bx;
}

template <class T> void oct_kick_pass(T x[], double k3l, double k3sl)
{
  T preal, pimag;
  T  By, Bx;
  
  preal = x[0]*x[0]*x[0] - 3 * x[0] * x[2] * x[2];
  pimag = 3*x[0]*x[0]*x[2] - x[2]*x[2]*x[2];
  Bx    = (k3sl * preal + k3l  * pimag)/6.;
  By    = (k3l  * preal - k3sl * pimag)/6.;  
  x[px_]= x[px_]-By;
  x[py_]= x[py_]+Bx;
}

template <class T>  void mult_kick_pass(T x[], int Norder, double KNL[11], double KNSL[11])
{
  int i;
  int fac=1;
  T   Xn, Yn, Xn0, Yn0;
  T   By, Bx;

  By=KNL[0];
  Bx=KNSL[0];
  Xn=1.;
  Yn=0.;

  for(i=1;i<Norder+1;i++){
    Xn0=Xn;
    Yn0=Yn;
    Xn=Xn0*x[x_]-Yn0*x[y_];
    Yn=Xn0*x[y_]+Yn0*x[x_];
    fac=fac*i;
    if ( KNL[i] != 0. || KNSL[i] !=0. ) {
      By=By+(KNL[i]*Xn-KNSL[i]*Yn)/fac;
      Bx=Bx+(KNL[i]*Yn+KNSL[i]*Xn)/fac;
    }
  }
  x[px_]=x[px_]-By;
  x[py_]=x[py_]+Bx;
}

template <class T> void bend_mult_kick_pass(T x[6], double L, double href, int Norder, double KNL[11], double KNSL[11])
{
  int i;
  int fac=1;
  T   Xn, Yn, Xn0, Yn0;
  T   By, Bx;

  By=KNL[0];
  Bx=KNSL[0];
  Xn=1.;
  Yn=0.;

  for(i=1;i<Norder+1;i++){
    Xn0=Xn;
    Yn0=Yn;
    Xn=Xn0*x[x_]-Yn0*x[y_];
    Yn=Xn0*x[y_]+Yn0*x[x_];
    fac=fac*i;
    if ( KNL[i] != 0. || KNSL[i] !=0. ) {
      By=By+(KNL[i]*Xn-KNSL[i]*Yn)/fac;
      Bx=Bx+(KNL[i]*Yn+KNSL[i]*Xn)/fac;
    }
  }
  x[px_]=x[px_]-By  + (href*x[pt_]-href*href*x[x_])*L;;
  x[py_]=x[py_]+Bx;
  x[z_]=x[z_]-href*x[x_]*L;
}

template <class T> void general_bend_mult_kick_pass(T x[6], double L, double href, double hreal, int Norder, double KNL[11], double KNSL[11])
{
  int i;
  int fac=1;
  T   Xn, Yn, Xn0, Yn0;
  T   By, Bx;

  By=KNL[0];
  Bx=KNSL[0];
  Xn=1.;
  Yn=0.;

  for(i=1;i<Norder+1;i++){
    Xn0=Xn;
    Yn0=Yn;
    Xn=Xn0*x[x_]-Yn0*x[y_];
    Yn=Xn0*x[y_]+Yn0*x[x_];
    fac=fac*i;
    if ( KNL[i] != 0. || KNSL[i] !=0. ) {
      By=By+(KNL[i]*Xn-KNSL[i]*Yn)/fac;
      Bx=Bx+(KNL[i]*Yn+KNSL[i]*Xn)/fac;
    }
  }
  x[px_]=x[px_]-By  + ( -(hreal-href) + href*x[pt_] - href*hreal*x[x_] ) * L ; 
  x[py_]=x[py_]+Bx;
  x[z_]=x[z_]-href*x[x_]*L;
}

template <class T> void cal_Bfield(T x[], int Nint, int Norder, double Angle, double KNL[11], double KNSL[11],T BLbrho[])
{
  int i;
  int fac=1;
  T  Xn, Yn, Xn0, Yn0;
  T  By, Bx;

  By=KNL[0];
  Bx=KNSL[0];
  Xn=1.;
  Yn=0.;

  for(i=1;i<Norder+1;i++){
    Xn0=Xn;
    Yn0=Yn;
    Xn=Xn0*x[x_]-Yn0*x[y_];
    Yn=Xn0*x[y_]+Yn0*x[x_];
    fac=fac*i;
    if ( KNL[i] != 0. || KNSL[i] !=0. ) {
      By=By+(KNL[i]*Xn-KNSL[i]*Yn)/fac;
      Bx=Bx+(KNL[i]*Yn+KNSL[i]*Xn)/fac;
    }
  }
  
  BLbrho[0]=Bx/Nint;
  BLbrho[1]=(Angle + By)/Nint;
  BLbrho[2]=0.;
}

template <class T> void radiate(T x[], double L, double href, T BLbrho[])
// calculate radiation effect for a short integration step L, no  quantumn fluctuation
// L:  integration step for reference orbit,
// herf: curvature for the reference orbit
// BL/(Brho)_0, BL for the particle, (Brho)_0=GP.p/Gp.q for reference particle.
{
  int i;
  T x0[6];
  T delta1, gamma1, beta1, delta2;
  T pz, xp, yp, temp1,temp2,ratio;
  T u[3], BLper[3], BLper2;
  
  for (i=0;i<6;i++) x0[i]=x[i];

  delta1= sqrt( 1.0 + 2*x[pt_]/GP.beta+x[pt_]*x[pt_]) -1.0;
  gamma1 =GP.gamma+sqrt(GP.gamma*GP.gamma-1)*x[pt_]; 
  beta1 = sqrt(1.0-1.0/gamma1/gamma1);

  pz=sqrt( (1+delta1)*(1+delta1)-x[px_]*x[px_]-x[py_]*x[py_]);
  xp= (1+ x[x_]*href) * x[px_] / pz ; 
  yp= (1+ x[x_]*href) * x[py_] / pz ; 
  temp1=sqrt( (1+ x[x_]*href) * (1+ x[x_]*href) + xp*xp  + yp*yp );
  u[0]=xp/temp1;   u[1] = yp/temp1;   u[2] = (1+ x[x_]*href)/temp1;
  
  vector_cross_product(BLbrho, u, BLper);
  
  BLper2=BLper[0]*BLper[0] + BLper[1]*BLper[1]  +  BLper[2]*BLper[2];
  //temp2= Cr * (1+delta1) *(1+delta1) * pow( GP.energy/1000.,3) * pow( GP.beta, 3) / beta1 /beta1 / 2 / PI; 
  temp2= Cr * (1+delta1) *(1+delta1) * pow( GP.energy/1000.,3) / 2 / PI;
  //x[pt_]=x[pt_]  - temp2 * BLper2 * temp1 / beta1 / L ;
  x[pt_]=x[pt_]  - temp2 * BLper2 * temp1 / L ;
  
  delta2=sqrt( 1.0 + 2*x[pt_]/GP.beta + x[pt_]*x[pt_]) -1.0;
  ratio = (1+ delta2)/(1+delta1);
  x[px_]=x0[px_] * ratio; 
  x[py_]=x0[py_] * ratio;
}

template <class T> void radiate1(T  x[], double L, double href, T BLbrho[])
// calculate radiation effect for a short integration step L
// L:  integration step for reference orbit,
// herf: curvature for the reference orbit
// BL/(Brho)_0, BL for the particle, (Brho)_0=GP.p/Gp.q for reference particle.
// including statistics treatment of quantumn fluctuation
{
  int i;
  T x0[6];
  T delta1, gamma1, beta1, delta2;
  T pz, xp, yp, temp1,temp2,ratio;
  T u[3], BLper[3], BLper2;
  T g1, dtdl, ptloss, ptfluct;
  double randgaus ;
  
  for (i=0;i<6;i++) x0[i]=x[i];

  delta1= sqrt( 1.0 + 2*x[pt_]/GP.beta+x[pt_]*x[pt_]) -1.0;
  gamma1 =GP.gamma+sqrt(GP.gamma*GP.gamma-1)*x[pt_]; 
  beta1 = sqrt(1.0-1.0/gamma1/gamma1);

  pz=sqrt( (1+delta1)*(1+delta1)-x[px_]*x[px_]-x[py_]*x[py_]);
  xp= (1+ x[x_]*href) * x[px_] / pz ; 
  yp= (1+ x[x_]*href) * x[py_] / pz ; 
  temp1=sqrt( (1+ x[x_]*href) * (1+ x[x_]*href) + xp*xp  + yp*yp );
  u[0]=xp/temp1;   u[1] = yp/temp1;   u[2] = (1+ x[x_]*href)/temp1;
  
  vector_cross_product(BLbrho, u, BLper);
  BLper2=BLper[0]*BLper[0] + BLper[1]*BLper[1]  +  BLper[2]*BLper[2];

  //temp2= Cr * (1+delta1) *(1+delta1) * pow( GP.energy/1000.,3) * pow( GP.beta, 3) / beta1 /beta1 / 2 / PI; 
  temp2= Cr * (1+delta1) *(1+delta1) * pow( GP.energy/1000.,3) / 2 / PI;
  //x[pt_]=x[pt_]  - temp2 * BLper2 * temp1 / beta1 / L ;
  //x[pt_]=x[pt_]  - temp2 * BLper2 * temp1 / L ;
  ptloss =   temp2 * BLper2 * temp1 / L;
  
  dtdl = temp1/beta1/speed_light;
  g1=sqrt( BLper2 / L / L) / (1+delta1);   //---g1 =  1/ rho1, rho1  is real radius, instead of reference frame's rho
  
  randgaus = gsl_ran_gaussian ( gsl_r2, 1.0);  // gaussian(0., 1.0, seed) is a bad generator ;
  ptfluct= sqrt( Cf * gamma1*gamma1*gamma1*gamma1*gamma1 *g1*g1*g1* dtdl* L ) * randgaus; //  gaussian(0., 1.0, seed) ;
  
  x[pt_]=x[pt_]  - ptloss  + ptfluct ;  
  //x[pt_]=x[pt_]  - ptloss;

  //cout<<ptloss <<"  "<<ptfluct<<"   "<<1/g1<<"  "<<randg<<" "<<dtdl<<" "<<L<<endl;
  //exit(0);
  
  delta2=sqrt( 1.0 + 2*x[pt_]/GP.beta + x[pt_]*x[pt_]) -1.0;
  ratio = (1.0 + delta2)/(1+delta1);
  x[px_]=x0[px_] * ratio; 
  x[py_]=x0[py_] * ratio;
}

void radiate2(double x[], double L, double href, double BLbrho[])
// href = 1/ rho_0, rho_0 is for reference orbit, BL/(Brho)_0, BL for the particle, L is for reference orbit.  
// this function to simulate photon emission
// including Monte Carlo  simulation of photon emission
{

  int i;
  double alpha = 1.0/137.0359997;
  
  double x0[6];
  double delta1, gamma1, beta1, delta2;
  double pz, xp, yp, temp1,temp2,ratio;
  double u[3], BLper[3], BLper2;
  double rho1, Nr, Ng, Etot;

  int    nt = 43;
  double xt[]={   .00123e0, .0123e0, .0265e0, .0571e0, .1228e0, .1544e0, 
                  0.194e0, .221e0, .2619e0, .2893e0, .327e0,.36914e0, 
                   .4074e0,.4417e0,.47273e0,.5209e0,.5791e0,.62385e0,.6517e0,
                   .6905e0,  .7165e0,     .7494e0,     .7726e0, 
                   .8104e0,.8401e0, .85954e0, .8835e0,.892e0,.89665e0,.8997e0,
                   .90397e0,.9067e0,.9132e0, .93444e0, .9561e0,  .9661e0,
	  	 .9768e0, .98734e0, .9930e0, .99776e0, .999574e0, .9999156e0, 1.e0};
    
  double  yt[]={ 1.e-9, 1.e-6, 1.e-5, 1.e-4, 1.e-3, 2.e-3, 4.e-3, 6.e-3, 
                   .01e0, .0136e0, 
                   .02e0, .0292e0, .04e0,  .052e0, .065e0, .09e0, .13e0,.17e0,
                   .2e0,.25e0,  .29e0, .35e0, 
                   .4e0, .5e0, .6e0, .68e0, .8e0, .85e0, .88e0, .9e0,.93e0,
                   .95e0, 1.e0, 1.2e0, 1.5e0,1.7e0, 
                   2.e0, 2.5e0, 3.e0, 4.e0,5.5e0,7.e0, 1.e1};
  
  double  facg=12. * sqrt(3.) / ( pow(2., 1.0/3.0) * 5.0 * 2.6789385347077475);
  double  dxspli=xt[nt-1];
  double  rtemp, Etemp;
  double  Ec;
  
  for (i=0;i<6;i++) x0[i]=x[i];

  delta1= sqrt( 1.0 + 2*x[5]/GP.beta+x[5]*x[5]) -1.0;
  gamma1 =GP.gamma+sqrt(GP.gamma*GP.gamma-1)*x[5]; 
  beta1 = sqrt(1.0-1.0/gamma1/gamma1);

  pz=sqrt( (1+delta1)*(1+delta1)-x[1]*x[1]-x[3]*x[3]);
  xp= (1+ x[0]*href) * x[1] / pz ; 
  yp= (1+ x[0]*href) * x[1] / pz ; 
  temp1=sqrt( (1+ x[0]*href) * (1+ x[0]*href) + xp*xp  + yp*yp );
  u[0]=xp/temp1;   u[1] = yp/temp1;   u[2] = (1+ x[0]*href)/temp1;
  
  vector_cross_product(BLbrho, u, BLper);
  BLper2=BLper[0]*BLper[0] + BLper[1]*BLper[1]  +  BLper[2]*BLper[2];

  rho1= (1+x[5]) / sqrt( BLper2 ) *   L ; 
  Ec=1.5 * 6.582119e-16 * 2.99792458e8 * gamma1 * gamma1 * gamma1 / rho1 / 1.0e6;
  Nr = 5. * alpha * gamma1 / 2. / sqrt(3.) / rho1 * L * temp1 /beta1 ;  
  Ng = gsl_ran_poisson(gsl_r1, Nr) ;

  Etot = 0. ;
  for(i=0;i<Ng; i++){
    rtemp = gsl_rng_uniform (gsl_r2);
    if(rtemp < 0.26){
      Etemp = pow(rtemp/facg, 3.0);
    }
    else{
      Etemp = splint(xt, yt, nt, rtemp);
    }
    Etot= Etot + Etemp;
  }
  Etot = Etot * Ec;

  //cout<<" Ec  Nr  Ng  Etot  "<< Ec <<"  "<<Nr<<"  "<<Ng<<"  "<<Etot<<endl;
  
  x[5]=x[5]  - Etot/GP.energy;
  delta2=sqrt( 1.0 + 2*x[5]/GP.beta + x[5]*x[5]) -1.0;
  ratio = (1+ delta2)/(1+delta1);
  x[1]=x0[1] * ratio; 
  x[3]=x0[3] * ratio;
}

//-------first  order approach by E. Forest, K. Ohmi

template<class T> void get_Axy_wiggler(T AxoBrho[], T AyoBrho[], double z, T x[], double b0, double kx, double kz, double phiz0)

{
  int     i;
  double  ky;
  T       cx, cz, sx, sz, chy, shy;
  
  ky=sqrt(kx*kx+kz*kz);
  for (i = 0; i <4; i++) {
    AxoBrho[i] = 0.0; AyoBrho[i] = 0.0;
  }
  cx  = cos(kx*x[0]);     sx  = sin(kx*x[0]);
  chy = cosh(ky*x[2]);    shy = sinh(ky*x[2]);
  cz  = cos(kz*z+phiz0);   sz  = sin(kz*z+phiz0);
  
  AxoBrho[0] = (b0/kz)*cx*chy*sz;
  AyoBrho[0] = (b0*kx/ky/kz)*sx*shy*sz;
  
  AxoBrho[1] = -(b0*kx/kz)*sx*chy*sz;
  AyoBrho[1] =  (b0*kx*kx/ky/kz)*cx*shy*sz;
  
  AxoBrho[2] = (b0*ky/kz)*cx*shy*sz;
  AyoBrho[2] = (b0*kx/kz)*sx*chy*sz;
  
  AxoBrho[3] = b0*cx*chy*cz;
  AyoBrho[3] = (b0*kx/ky)*sx*shy*cz;
  
}

template<class T> void Wiggler_Pass_Forest(T x[],  double l, int nint, double b0, double kx, double kz, double phiz0) 
{
  int i,j;
  double h = l/nint, z;
  T      AxoBrho[4], AyoBrho[4], psi, hodp;
  T      det, t1, t2, x01, x02, y01, y02;
  T      x0[6];
  T      BL[3];
  
  for (i = 0; i < nint; i++) {
    
    z  = h* i ; 
    get_Axy_wiggler(AxoBrho, AyoBrho, z, x, b0, kx, kz, phiz0);

    for(j=0;j<4;j++){
      AxoBrho[j]=AxoBrho[j]/GP.brho;
      AyoBrho[j]=AyoBrho[j]/GP.brho;
    }

    for(j=0;j<6;j++) x0[j] = x[j];
    
    psi = 1.0 + x0[5]; hodp = h/psi;
    det = (  (1- AxoBrho[1]*hodp)*(1- AyoBrho[2]*hodp) - AxoBrho[2] * AyoBrho[1] *  hodp  *  hodp ) ;

    x01= AxoBrho[0] * AxoBrho[1] ; x02= AxoBrho[0] * AxoBrho[2];
    y01= AyoBrho[0] * AyoBrho[1] ; y02= AyoBrho[0] * AyoBrho[2];
    
    
    x[1] = ( x0[1] - ( x01  + y01 )* hodp ) * ( 1- AyoBrho[2] * hodp)
        +  ( x0[3] - ( x02  + y02 )* hodp ) * AyoBrho[1] * hodp ;
    x[1] = x[1]/ det; 

    x[3] = ( x0[3] - (  x02 +  y02 )* hodp ) * ( 1- AxoBrho[1] * hodp)
        +  ( x0[1] - (  x01 +  y01 )* hodp ) * AxoBrho[2] * hodp ;
    x[3] = x[3]/ det;

    x[0]= x0[0] + hodp * ( x[1] - AxoBrho[0] );
    x[2]= x0[2] + hodp * ( x[3] - AyoBrho[0] );

    t1= (x[1]-AxoBrho[0])/psi; t2= (x[3]-AyoBrho[0])/psi;
    x[4] = x[4] - h*( t1*t1+t2*t2)/2.0;

    //cout<<x[0]<<" "<<x[1]<<" "<<x[2]<<" "<<x[3]<<" "<<z+h<<" "<<x[4]<<" "<<x[5]<<endl;
    
    if(GP.radiate  == true  ){
      BL[0] = -AyoBrho[3] * GP.brho * h;
      BL[1] =  AxoBrho[3] * GP.brho * h;
      BL[2] =  (AyoBrho[1] - AxoBrho[2]) * GP.brho * h;
      radiate(x, h, 0, BL);  
    }
    
  }
  
}

//------second order apporach by Y. Wu, E. Forest, D. Robin

template<class T>  void get_Ax_Wu(T & AxoBrho, T & IntpAxpydx, double z, T x[], double b0, double kx, double kz, double phiz0)
{
  double ky;
  T cx, sx,chy, shy, cz, sz;
  
  ky=sqrt(kx*kx+kz*kz);
  AxoBrho = 0.0; IntpAxpydx = 0.;

  cx  = cos(kx*x[0]);     sx  = sin(kx*x[0]);
  chy = cosh(ky*x[2]);    shy = sinh(ky*x[2]);
  cz  = cos(kz*z+phiz0);   sz  = sin(kz*z+phiz0);

  AxoBrho = (b0/kz)*cx*chy*sz;
  if( kx ==0.){
    IntpAxpydx = (b0*ky/kz)*x[0]*shy*sz;
  }
  else{
    IntpAxpydx = (b0*ky/kx/kz)*sx*shy*sz;
  }
  
}

template <class T>  void get_Ay_Wu(T & AyoBrho, T & IntpAypxdy, double z, T x[], double b0, double kx, double kz, double phiz0)
{
  double ky;
  T cx, sx,chy, shy, cz, sz;
  
  ky=sqrt(kx*kx+kz*kz);
  AyoBrho = 0.0; IntpAypxdy = 0.;

  cx  = cos(kx*x[0]);     sx  = sin(kx*x[0]);
  chy = cosh(ky*x[2]);    shy = sinh(ky*x[2]);
  cz  = cos(kz*z+phiz0);   sz  = sin(kz*z+phiz0);

  AyoBrho = (b0*kx/ky/kz)*sx*shy*sz;
  IntpAypxdy = (b0*kx*kx/ky/ky/kz)*cx*chy*sz;

}


template <class T> void get_wiggler_B(T x[], T B[], double b0, double kx, double kz, double phiz0)
{
  int  i;
  double ky;
  T cx, sx,chy, shy, cz, sz;
  
  ky=sqrt(kx*kx+kz*kz);
  for(i=0;i<3;i++) B[i]=0.;

  cx  = cos(kx*x[0]);     sx  = sin(kx*x[0]);
  chy = cosh(ky*x[2]);    shy = sinh(ky*x[2]);
  cz  = cos(kz*x[4]+phiz0);   sz  = sin(kz*x[4]+phiz0);

  B[0] =  (kx/ky)*b0*sx*shy*cz;
  B[1] =  b0*cx*chy*cz;
  B[2] = -(kz*b0/ky)*cx*shy*sz;

}

template <class T> void Wiggler_Pass_Wu(T x[],  double l, int nint, double b0, double kx, double kz, double phiz0) 
{
  int     i,j;
  double  h, z;
  T       hd, AxoBrho, AyoBrho, IntpAxpydx,  IntpAypxdy, B[3];
  
  h = l / nint; z = 0.0;

  for (i = 0; i < nint; i++) {
    hd = h/(1.0+x[5]);

    // 1: K1 * h /2 
    z += 0.5*h;

    // 2: K2 * h /2
    get_Ay_Wu(AyoBrho, IntpAypxdy, z, x, b0, kx, kz, phiz0);
    x[1] = x[1] -  IntpAypxdy/GP.brho; x[3] = x[3] - AyoBrho/GP.brho;

    x[2] = x[2] + 0.5*hd*x[3];
    x[4] = x[4] - hd*x[3]*x[3]/(1.0+x[5])/(1.0+x[5])/4;

    get_Ay_Wu(AyoBrho, IntpAypxdy, z, x, b0, kx, kz, phiz0);
    x[1] = x[1] + IntpAypxdy/GP.brho; x[3] = x[3] + AyoBrho/GP.brho;
    
    // 3: k3 * h
    get_Ax_Wu(AxoBrho, IntpAxpydx, z, x, b0, kx, kz, phiz0);
    x[1] =x[1] - AxoBrho/GP.brho; x[3] = x[3] - IntpAxpydx/GP.brho;

    x[0] = x[0] +  hd*x[1];
    x[4] = x[4] - 0.5*hd*x[1]*x[1]/(1.0+x[5]);

    get_Ax_Wu(AxoBrho, IntpAxpydx, z, x, b0, kx, kz, phiz0);
    x[1] =x[1] + AxoBrho/GP.brho; x[3] = x[3] + IntpAxpydx/GP.brho;

    // 4: K2 * h /2
    get_Ay_Wu(AyoBrho, IntpAypxdy, z, x, b0, kx, kz, phiz0);
    x[1] = x[1] -  IntpAypxdy/GP.brho; x[3] = x[3] - AyoBrho/GP.brho;

    x[2] = x[2] + 0.5*hd*x[3];
    x[4] = x[4] - hd*x[3]*x[3]/(1.0+x[5])/(1.0+x[5])/4;

    get_Ay_Wu(AyoBrho, IntpAypxdy, z, x, b0, kx, kz, phiz0);
    x[1] = x[1] + IntpAypxdy/GP.brho; x[3] = x[3] + AyoBrho/GP.brho;
        
    // 5: K1 * h /2 
    z += 0.5*h;

    if (GP.radiate == true) {
      get_wiggler_B(x, B, b0, kx, kz, phiz0);
      for(j=0;j<3;j++) B[j] =B[j]*h;
      radiate(x, h, 0.0, B);
    }

    //cout<<x[0]<<" "<<x[1]<<" "<<x[2]<<" "<<x[3]<<" "<<z<<" "<<x[4]<<" "<<x[5]<<endl;
    
  }
}

void spin(double x[], double angle, double href, double BLbrho[])
// spin tracking for a short tracking step, only spin or x[6-8] will be updated
// angle, href=1/rho: for reference orbit; href used for calculating v direction only 
// BLbrho:  (BL)/(Brho)_0, or (Brho)_0= GP.p/GP.q for reference particle
// (BL), B-->test particle, L-->reference orbit path length
// Gr :  for the test particle, not for reference particle
{
  int i,j;
  double u[3];
  double Bpar[3], udotB;
  double temp;
  double delta1=sqrt( 1.0 + 2.*x[pt_]/GP.beta+x[pt_]*x[pt_])-1.0;  
  double gamma1=GP.gamma+sqrt(GP.gamma*GP.gamma-1.)*x[pt_]; 
  double Gr   = GP.G * gamma1;  
  double F[3], Fabs, uF[3];
  double theta, costheta, sintheta;
  double rot[9];
  double spin0[3], spin[3];
  
  if(  BLbrho[0]* BLbrho[0] +  BLbrho[1]* BLbrho[1] +  BLbrho[2]* BLbrho[2] !=0. ){

    for(i=0;i<3;i++) spin0[i]=x[6+i];

    temp=sqrt( (1+delta1) * (1+delta1)-x[px_]*x[px_]-x[py_]*x[py_] );
    u[0] = x[px_] / temp;
    u[1] = x[py_] / temp;
    u[2] = 1.0 ;
    temp=sqrt( u[0]*u[0] +  u[1]*u[1] +  u[2]*u[2] );
    for(i=0;i<3;i++) u[i]=u[i]/temp;
    temp=(1 + href * x[0])* temp;  
 
    udotB =  u[0] * BLbrho[0] +  u[1] * BLbrho[1] +  u[2] * BLbrho[2];
    for(i=0;i<3;i++) Bpar[i] =  udotB * u[i];
    for(i=0;i<3;i++) F[i] = ( 1+ Gr ) * BLbrho[i]  + GP.G*(1-gamma1)*Bpar[i] ;
    for(i=0;i<3;i++) F[i] = F[i] * temp / ( 1+ delta1) ; 
    F[1]=F[1] - angle;
    Fabs=sqrt(F[0]* F[0] + F[1]*F[1]  +  F[2]* F[2] );
    for(i=0;i<3;i++) uF[i] = F[i]/Fabs; 
    theta=-Fabs; costheta = cos(theta); sintheta= sin(theta);
 
    rot[0] = uF[0]*uF[0]*(1-costheta) + costheta ;  
    rot[1] = uF[0]*uF[1]*(1-costheta) - uF[2]*sintheta;  
    rot[2] = uF[0]*uF[2]*(1-costheta) + uF[1]*sintheta;  
    rot[3] = uF[0]*uF[1]*(1-costheta) + uF[2]*sintheta;  
    rot[4] = uF[1]*uF[1]*(1-costheta) + costheta ;   
    rot[5] = uF[1]*uF[2]*(1-costheta) - uF[0]*sintheta;  
    rot[6] = uF[0]*uF[2]*(1-costheta) - uF[1]*sintheta; 
    rot[7] = uF[1]*uF[2]*(1-costheta) + uF[0]*sintheta;  
    rot[8] = uF[2]*uF[2]*(1-costheta) + costheta;
    for(i=0;i<3;i++){
      spin[i]=0.;
      for(j=0;j<3;j++) spin[i] = spin[i] + rot[i*3 + j] * spin0[j];
    }

    for(i=0;i<3;i++) x[6+i]= spin[i];
  }
}

template <class T>  void SBEND_DAPass(T x[], double L, int Nint, double Angle, double E1, double E2) 
{
  int i,j;
  double href=Angle/L;
  double Lint=L/Nint;
  double edgefocus;
  double cosAngle, sinAngle;
  T x1[6], BLbrho[3];
  
  if(Angle==0.){
    DRIFT_Pass(x,L); return ;
  }
  
  edgefocus = tan(E1)*href;
  x[1] = x[1]+ edgefocus*x[0];
  x[3] = x[3]- edgefocus*x[2];

  if( GP.H_expand == true  and GP.radiate == false){
    for(i=0;i<Nint;i++){
      DRIFT_Pass(x,Fdrift1*Lint);
      bend_kick_pass(x, Fkick1*Lint, href);
      DRIFT_Pass(x,Fdrift2*Lint);
      bend_kick_pass(x, Fkick2*Lint, href);
      DRIFT_Pass(x,Fdrift2*Lint);
      bend_kick_pass(x, Fkick1*Lint, href);
      DRIFT_Pass(x,Fdrift1*Lint); 
    }
  }
  
  if( GP.H_expand == true  and GP.radiate == true){
    BLbrho[0] = 0.;  BLbrho[1] = Angle/Nint;  BLbrho[2] = 0.;  
    for(i=0;i<Nint;i++){
      DRIFT_Pass(x,Fdrift1*Lint);
      bend_kick_pass(x, Fkick1*Lint, href);
      DRIFT_Pass(x,Fdrift2*Lint);
      bend_kick_pass(x, Fkick2*Lint, href);
      DRIFT_Pass(x,Fdrift2*Lint);
      bend_kick_pass(x, Fkick1*Lint, href);
      DRIFT_Pass(x,Fdrift1*Lint); 
      radiate(x, Lint, href, BLbrho);
    }
  }
  
  if(GP.H_expand == false  and GP.radiate == false ){
    sbend_exact_pass(x,L,Angle,cos(Angle),sin(Angle));
  }

  if(GP.H_expand == false  and GP.radiate == true ){
    BLbrho[0] = 0.;   BLbrho[1] = Angle/Nint;  BLbrho[2] = 0.; 
    cosAngle=cos(Angle/Nint);   sinAngle=sin(Angle/Nint);
    for(i=0;i<Nint;i++){
      sbend_exact_pass(x,Lint,Angle/Nint,cosAngle,sinAngle);
      radiate(x, Lint, href, BLbrho);
    }
  }
  
  edgefocus = tan(E2)*href;
  x[1] = x[1]+edgefocus*x[0];
  x[3] = x[3]-edgefocus*x[2];
}

template <class T>  void SBEND_Pass(T x[], double L, int Nint, double Angle, double E1, double E2) 
{
  int i,j;
  double href=Angle/L;
  double Lint=L/Nint;
  double edgefocus;
  double cosAngle, sinAngle;
  T x1[6], BLbrho[3];
  
  if(Angle==0.){
    DRIFT_Pass(x,L); return ;
  }
  
  edgefocus = tan(E1)*href;
  x[1] = x[1] + edgefocus*x[0];
  x[3] = x[3] - edgefocus*x[2];

  if( GP.H_expand == true  and GP.radiate == false){
    for(i=0;i<Nint;i++){
      DRIFT_Pass(x,Fdrift1*Lint);
      bend_kick_pass(x, Fkick1*Lint, href);
      DRIFT_Pass(x,Fdrift2*Lint);
      bend_kick_pass(x, Fkick2*Lint, href);
      DRIFT_Pass(x,Fdrift2*Lint);
      bend_kick_pass(x, Fkick1*Lint, href);
      DRIFT_Pass(x,Fdrift1*Lint); 
    }
  }
  
  if( GP.H_expand == true  and GP.radiate == true){
    BLbrho[0] = 0.;  BLbrho[1] = Angle/Nint;  BLbrho[2] = 0.;  
    for(i=0;i<Nint;i++){
      DRIFT_Pass(x,Fdrift1*Lint);
      bend_kick_pass(x, Fkick1*Lint, href);
      DRIFT_Pass(x,Fdrift2*Lint);
      bend_kick_pass(x, Fkick2*Lint, href);
      DRIFT_Pass(x,Fdrift2*Lint);
      bend_kick_pass(x, Fkick1*Lint, href);
      DRIFT_Pass(x,Fdrift1*Lint); 
      radiate1(x, Lint, href, BLbrho);
    }
  }
  
  if(GP.H_expand == false  and GP.radiate == false ){
    sbend_exact_pass(x,L,Angle,cos(Angle),sin(Angle));
  }

  if(GP.H_expand == false  and GP.radiate == true ){
    BLbrho[0] = 0.;   BLbrho[1] = Angle/Nint;  BLbrho[2] = 0.; 
    cosAngle=cos(Angle/Nint);   sinAngle=sin(Angle/Nint);
    for(i=0;i<Nint;i++){
      sbend_exact_pass(x,Lint,Angle/Nint,cosAngle,sinAngle);
      radiate1(x, Lint, href, BLbrho);
    }
  }
  
  edgefocus = tan(E2)*href;
  x[1] = x[1]+edgefocus*x[0];
  x[3] = x[3]-edgefocus*x[2];
}

template <class T>  void GSBEND_DAPass(T x[], double L, int Nint, double Angle, double K0L, double E1, double E2) 
{
  int i,j;
  double href=Angle/L;
  double hreal=K0L/L;
  
  double Lint=L/Nint;
  double edgefocus;
  double cosAngle, sinAngle;
  T x1[6], BLbrho[3];
  
  if(Angle==0.){
    DRIFT_Pass(x,L); return ;
  }
  
  edgefocus = tan(E1)*hreal;
  x[1] = x[1]+ edgefocus*x[0];
  x[3] = x[3]- edgefocus*x[2];

  if( GP.H_expand == true  and GP.radiate == false){
    for(i=0;i<Nint;i++){
      DRIFT_Pass(x,Fdrift1*Lint);
      general_bend_kick_pass(x, Fkick1*Lint, href, hreal);
      DRIFT_Pass(x,Fdrift2*Lint);
      general_bend_kick_pass(x, Fkick2*Lint, href, hreal);
      DRIFT_Pass(x,Fdrift2*Lint);
      general_bend_kick_pass(x, Fkick1*Lint, href, hreal);
      DRIFT_Pass(x,Fdrift1*Lint); 
    }
  }
  
  if( GP.H_expand == true  and GP.radiate == true){
    BLbrho[0] = 0.;  BLbrho[1] = K0L/Nint;  BLbrho[2] = 0.;
    for(i=0;i<Nint;i++){
      DRIFT_Pass(x,Fdrift1*Lint);
      general_bend_kick_pass(x, Fkick1*Lint, href, hreal);
      DRIFT_Pass(x,Fdrift2*Lint);
      general_bend_kick_pass(x, Fkick2*Lint, href, hreal);
      DRIFT_Pass(x,Fdrift2*Lint);
      general_bend_kick_pass(x, Fkick1*Lint, href, hreal);
      DRIFT_Pass(x,Fdrift1*Lint); 
      radiate(x, Lint, href, BLbrho);
    }
  }
  
  if(GP.H_expand == false  and GP.radiate == false ){
    general_sbend_exact_pass(x,L,Angle,K0L,cos(Angle),sin(Angle));
  }

  if(GP.H_expand == false  and GP.radiate == true ){
    BLbrho[0] = 0.;   BLbrho[1] = K0L/Nint;  BLbrho[2] = 0.; 
    cosAngle=cos(Angle/Nint);   sinAngle=sin(Angle/Nint);
    for(i=0;i<Nint;i++){
      general_sbend_exact_pass(x, Lint, Angle/Nint, K0L/Nint, cosAngle, sinAngle);
      radiate(x, Lint, href, BLbrho);
    }
  }
  
  edgefocus = tan(E2)*hreal;
  x[1] = x[1]+edgefocus*x[0];
  x[3] = x[3]-edgefocus*x[2];

}

template <class T>  void GSBEND_Pass(T x[], double L, int Nint, double Angle, double K0L, double E1, double E2) 
{
  int i,j;
  double href=Angle/L;
  double hreal=K0L/L;
  
  double Lint=L/Nint;
  double edgefocus;
  double cosAngle, sinAngle;
  T x1[6], BLbrho[3];
  
  if(Angle==0.){
    DRIFT_Pass(x,L); return ;
  }
  
  edgefocus = tan(E1)*hreal;
  x[1] = x[1]+ edgefocus*x[0];
  x[3] = x[3]- edgefocus*x[2];

  if( GP.H_expand == true  and GP.radiate == false){
    for(i=0;i<Nint;i++){
      DRIFT_Pass(x,Fdrift1*Lint);
      general_bend_kick_pass(x, Fkick1*Lint, href, hreal);
      DRIFT_Pass(x,Fdrift2*Lint);
      general_bend_kick_pass(x, Fkick2*Lint, href, hreal);
      DRIFT_Pass(x,Fdrift2*Lint);
      general_bend_kick_pass(x, Fkick1*Lint, href, hreal);
      DRIFT_Pass(x,Fdrift1*Lint); 
    }
  }
  
  if( GP.H_expand == true  and GP.radiate == true){
    BLbrho[0] = 0.;  BLbrho[1] = K0L/Nint;  BLbrho[2] = 0.;
    for(i=0;i<Nint;i++){
      DRIFT_Pass(x,Fdrift1*Lint);
      general_bend_kick_pass(x, Fkick1*Lint, href, hreal);
      DRIFT_Pass(x,Fdrift2*Lint);
      general_bend_kick_pass(x, Fkick2*Lint, href, hreal);
      DRIFT_Pass(x,Fdrift2*Lint);
      general_bend_kick_pass(x, Fkick1*Lint, href, hreal);
      DRIFT_Pass(x,Fdrift1*Lint); 
      radiate2(x, Lint, href, BLbrho);
    }
  }
  
  if(GP.H_expand == false  and GP.radiate == false ){
    general_sbend_exact_pass(x,L,Angle,K0L,cos(Angle),sin(Angle));
  }

  if(GP.H_expand == false  and GP.radiate == true ){
    BLbrho[0] = 0.;   BLbrho[1] = K0L/Nint;  BLbrho[2] = 0.; 
    cosAngle=cos(Angle/Nint);   sinAngle=sin(Angle/Nint);
    for(i=0;i<Nint;i++){
      general_sbend_exact_pass(x, Lint, Angle/Nint, K0L/Nint, cosAngle, sinAngle);
      radiate2(x, Lint, href, BLbrho);
    }
  }
  
  edgefocus = tan(E2)*hreal;
  x[1] = x[1]+edgefocus*x[0];
  x[3] = x[3]-edgefocus*x[2];

}

template <class T>  void QUAD_Pass(T x[], double L, int Nint, double k1l, double k1sl)
{
  int  i, j;
  
  if(k1l==0. and  k1sl ==0.){
    DRIFT_Pass(x,L);  return;
  }

  if(L==0.) 
    {
      quad_kick_pass(x, k1l, k1sl);  return;
    }

  if( L != 0.  and  GP.quad_fringe == true ){
    T x0[6];
    double k= k1l / L;
    
    for(i=0;i<6;i++) x0[i]=x[i];

    //T temp = k * ( 3 *x0[0]*x0[0] - 3 * x0[2] *x0[2] ) / 12;
    //T D    = 1- temp*temp;
    //x[0] = x0[0] + k * ( x0[0]*x0[0]*x0[0] + 3*x0[0]*x0[2]*x0[2])/12 ;
    //x[1] = x0[1] - k * ( ( x0[0]*x0[0] + x0[2]*x0[2] ) * x0[1] - 2*x0[0]*x0[2]*x0[3])/4;
    //x[1] = x[1] / D;
    //x[2] = x0[2] - k * ( 3 *x0[0]*x0[0]*x0[2] +  x0[2]*x0[2]*x0[2] ) / 12;
    //x[3] = x0[3] + k * ( ( x0[0]*x0[0] + x0[2]*x0[2] ) * x0[3] - 2*x0[0]*x0[2]*x0[1])/4;
    //x[3] = x[3] / D;

    x[0] = x0[0] + k * ( x0[0]*x0[0]*x0[0] + 3*x0[0]*x0[2]*x0[2])/12/(1+x0[5]) ;
    x[1] = x0[1] - k * ( ( x0[0]*x0[0] + x0[2]*x0[2] ) * x0[1] - 2*x0[0]*x0[2]*x0[3])/4/(1+x0[5]) ;
    x[2] = x0[2] - k * ( 3 *x0[0]*x0[0]*x0[2] +  x0[2]*x0[2]*x0[2] ) / 12 /(1+x0[5])  ;
    x[3] = x0[3] + k * ( ( x0[0]*x0[0] + x0[2]*x0[2] ) * x0[3] - 2*x0[0]*x0[2]*x0[1])/4/(1+x0[5]);
    x[4] = x0[4] + k * (  x0[2]*x0[2]*x0[2]*x0[3] - x0[0]*x0[0]*x0[0]*x0[1]
			  + 3*x0[0]*x0[0]*x0[2]*x0[3] - 3*x0[2]*x0[2]*x0[0]*x0[1] )/12/(1 + x0[5])/(1 + x0[5]);
  }

  double Lint=L/Nint;
  double k1l_kick1,k1sl_kick1;
  double k1l_kick2,k1sl_kick2;
  
  k1l_kick1 =Fkick1*k1l/Nint;
  k1sl_kick1=Fkick1*k1sl/Nint;
  k1l_kick2 =Fkick2*k1l/Nint;
  k1sl_kick2=Fkick2*k1sl/Nint;
  
  for(i=0;i<Nint;i++){
    DRIFT_Pass(x,Fdrift1*Lint);
    quad_kick_pass(x, k1l_kick1, k1sl_kick1);
    DRIFT_Pass(x,Fdrift2*Lint);
    quad_kick_pass(x, k1l_kick2, k1sl_kick2);
    DRIFT_Pass(x,Fdrift2*Lint);
    quad_kick_pass(x, k1l_kick1, k1sl_kick1);
    DRIFT_Pass(x,Fdrift1*Lint);
  }

  if( L != 0.  and  GP.quad_fringe == true ){
    T x0[6];
    double  k= k1l / L;
    for(i=0;i<6;i++) x0[i]=x[i];
    
    //T temp = k * ( 3 *x0[0]*x0[0] - 3 * x0[2] *x0[2] ) / 12;
    //T D    = 1- temp*temp;
    //x[0] = x0[0] - k * ( x0[0]*x0[0]*x0[0] + 3*x0[0]*x0[2]*x0[2])/12 ;
    //x[1] = x0[1] + k * ( ( x0[0]*x0[0] + x0[2]*x0[2] ) * x0[1] - 2*x0[0]*x0[2]*x0[3])/4;
    //x[1] = x[1] / D;
    //x[2] = x0[2] + k * ( 3 *x0[0]*x0[0]*x0[2] +  x0[2]*x0[2]*x0[2] ) / 12;
    //x[3] = x0[3] - k * ( ( x0[0]*x0[0] + x0[2]*x0[2] ) * x0[3] - 2*x0[0]*x0[2]*x0[1])/4;
    //x[3] = x[3] / D;

    x[0] = x0[0] - k * ( x0[0]*x0[0]*x0[0] + 3*x0[0]*x0[2]*x0[2])/12/(1+x0[5]) ;
    x[1] = x0[1] + k * ( ( x0[0]*x0[0] + x0[2]*x0[2] ) * x0[1] - 2*x0[0]*x0[2]*x0[3])/4/(1+x0[5]) ;
    x[2] = x0[2] + k * ( 3 *x0[0]*x0[0]*x0[2] +  x0[2]*x0[2]*x0[2] ) / 12 /(1+x0[5])  ;
    x[3] = x0[3] - k * ( ( x0[0]*x0[0] + x0[2]*x0[2] ) * x0[3] - 2*x0[0]*x0[2]*x0[1])/4/(1+x0[5]);
    x[4] = x0[4] - k * (  x0[2]*x0[2]*x0[2]*x0[3] - x0[0]*x0[0]*x0[0]*x0[1]
			  + 3*x0[0]*x0[0]*x0[2]*x0[3] - 3*x0[2]*x0[2]*x0[0]*x0[1] )/12/(1 + x0[5])/(1 + x0[5]);   
  }
  
}

template <class T>  void QUAD_Pass_Radiate(T x[], double L, int Nint, double k1l, double k1sl)
{
  if(k1l==0. and  k1sl ==0.){
    DRIFT_Pass(x,L);  return;
  }

  if(L==0.) 
    {
      quad_kick_pass(x, k1l, k1sl);  return;
    }
  
  int i,j;
  double Lint=L/Nint;
  double k1l_kick1,k1sl_kick1;
  double k1l_kick2,k1sl_kick2;
  T      x1[6], preal, pimag, BLbrho[3];

  k1l_kick1 =Fkick1*k1l/Nint;
  k1sl_kick1=Fkick1*k1sl/Nint;
  k1l_kick2 =Fkick2*k1l/Nint;
  k1sl_kick2=Fkick2*k1sl/Nint;
  
  for(i=0;i<Nint;i++){
    for(j=0;j<6;j++) x1[j]=x[j];
    
    DRIFT_Pass(x,Fdrift1*Lint);
    quad_kick_pass(x, k1l_kick1, k1sl_kick1);
    DRIFT_Pass(x,Fdrift2*Lint);
    quad_kick_pass(x, k1l_kick2, k1sl_kick2);
    DRIFT_Pass(x,Fdrift2*Lint);
    quad_kick_pass(x, k1l_kick1, k1sl_kick1);
    DRIFT_Pass(x,Fdrift1*Lint);
    
    for(j=0;j<6;j++) x1[j]=(x[j]+x1[j])/2.;
    preal = x1[0];   pimag = x1[2];
    BLbrho[0]=(k1sl * preal + k1l  * pimag )/Nint;
    BLbrho[1]=(k1l  * preal - k1sl * pimag )/Nint;
    BLbrho[2]= 0.  ;
    radiate(x,Lint,0, BLbrho);
  } 
  
}

template <class T>  void SEXT_Pass(T x[], double L, int Nint, double k2l, double k2sl)
{
  if(k2l==0. and k2sl==0.){
    DRIFT_Pass(x,L);  return;
  }
  
  if(L==0.) {
    sext_kick_pass(x, k2l, k2sl); return;
  }
  
  int i,j;
  double Lint=L/Nint;
  double k2l_kick1,k2sl_kick1;
  double k2l_kick2,k2sl_kick2;
  
  k2l_kick1 =Fkick1*k2l/Nint;
  k2sl_kick1=Fkick1*k2sl/Nint;
  k2l_kick2 =Fkick2*k2l/Nint;
  k2sl_kick2=Fkick2*k2sl/Nint;
  
  for(i=0;i<Nint;i++){
    DRIFT_Pass(x,Fdrift1*Lint);
    sext_kick_pass(x, k2l_kick1, k2sl_kick1);
    DRIFT_Pass(x,Fdrift2*Lint);
    sext_kick_pass(x, k2l_kick2, k2sl_kick2);
    DRIFT_Pass(x,Fdrift2*Lint);
    sext_kick_pass(x, k2l_kick1, k2sl_kick1);
    DRIFT_Pass(x,Fdrift1*Lint);
  }
  
}

template <class T>  void SEXT_Pass_Radiate(T x[], double L, int Nint, double k2l, double k2sl)
{
  if(k2l==0. and k2sl==0.){
    DRIFT_Pass(x,L);  return;
  }
  
  if(L==0.) {
    sext_kick_pass(x, k2l, k2sl); return;
  }
  
  int i,j;
  double Lint=L/Nint;
  double k2l_kick1,k2sl_kick1;
  double k2l_kick2,k2sl_kick2;
  T      x1[6], preal, pimag, BLbrho[3];
  
  k2l_kick1 =Fkick1*k2l/Nint;
  k2sl_kick1=Fkick1*k2sl/Nint;
  k2l_kick2 =Fkick2*k2l/Nint;
  k2sl_kick2=Fkick2*k2sl/Nint;
  
  for(i=0;i<Nint;i++){
    for(j=0;j<6;j++) x1[j]=x[j];
    
    DRIFT_Pass(x,Fdrift1*Lint);
    sext_kick_pass(x, k2l_kick1, k2sl_kick1);
    DRIFT_Pass(x,Fdrift2*Lint);
    sext_kick_pass(x, k2l_kick2, k2sl_kick2);
    DRIFT_Pass(x,Fdrift2*Lint);
    sext_kick_pass(x, k2l_kick1, k2sl_kick1);
    DRIFT_Pass(x,Fdrift1*Lint);
    
    for(j=0;j<6;j++) x1[j]=(x[j]+x1[j])/2.;
    preal =  x1[0]*x1[0]- x1[2]*x1[2];  pimag =  2. * x1[0] * x1[2];
    BLbrho[0]= ( k2sl * preal + k2l  * pimag )/2./Nint;
    BLbrho[1]= ( k2l  * preal - k2sl * pimag )/2./Nint;
    BLbrho[2]= 0.  ; 
    radiate(x,Lint,0, BLbrho);
  }
}

template <class T> void OCT_Pass(T x[], double L, int Nint, double k3l, double k3sl)
{
  if(k3l == 0. and k3sl ==0.){
    DRIFT_Pass(x,L);  return;
  }
  
  if(L==0.) {
    oct_kick_pass(x, k3l, k3sl); return ;
  }

  int i,j;
  double Lint=L/Nint;
  double k3l_kick1,k3sl_kick1;
  double k3l_kick2,k3sl_kick2;

  k3l_kick1 =Fkick1*k3l/Nint;
  k3sl_kick1=Fkick1*k3sl/Nint;
  k3l_kick2 =Fkick2*k3l/Nint;
  k3sl_kick2=Fkick2*k3sl/Nint;

  for(i=0;i<Nint;i++){
    DRIFT_Pass(x,Fdrift1*Lint);
    oct_kick_pass(x, k3l_kick1, k3sl_kick1);
    DRIFT_Pass(x,Fdrift2*Lint);
    oct_kick_pass(x, k3l_kick2, k3sl_kick2);
    DRIFT_Pass(x,Fdrift2*Lint);
    oct_kick_pass(x, k3l_kick1, k3sl_kick1);
    DRIFT_Pass(x,Fdrift1*Lint);
  }

}

template <class T> void OCT_Pass_Radiate(T x[], double L, int Nint, double k3l, double k3sl)
{
  if(k3l == 0. and k3sl ==0.){
    DRIFT_Pass(x,L);  return;
  }
  
  if(L==0.) {
    oct_kick_pass(x, k3l, k3sl); return ;
  }

  int i,j;
  double Lint=L/Nint;
  double k3l_kick1,k3sl_kick1;
  double k3l_kick2,k3sl_kick2;
  T      x1[6], preal, pimag, BLbrho[3];

  k3l_kick1 =Fkick1*k3l/Nint;
  k3sl_kick1=Fkick1*k3sl/Nint;
  k3l_kick2 =Fkick2*k3l/Nint;
  k3sl_kick2=Fkick2*k3sl/Nint;

  for(i=0;i<Nint;i++){
    for(j=0;j<6;j++) x1[j]=x[j];
    
    DRIFT_Pass(x,Fdrift1*Lint);
    oct_kick_pass(x, k3l_kick1, k3sl_kick1);
    DRIFT_Pass(x,Fdrift2*Lint);
    oct_kick_pass(x, k3l_kick2, k3sl_kick2);
    DRIFT_Pass(x,Fdrift2*Lint);
    oct_kick_pass(x, k3l_kick1, k3sl_kick1);
    DRIFT_Pass(x,Fdrift1*Lint);
    
    for(j=0;j<6;j++) x1[j]=(x[j]+x1[j])/2.;
    preal = x1[0]*x1[0]*x1[0] - 3 * x1[0] * x1[2]* x1[2];
    pimag = 3*x1[0]*x1[0]*x1[2] - x1[2]*x1[2]*x1[2];
    BLbrho[0]= (k3sl * preal + k3l  * pimag)/6/Nint ;
    BLbrho[1]= (k3l  * preal - k3sl * pimag)/6/Nint;
    BLbrho[2]= 0.  ; 
    radiate(x,Lint,0, BLbrho); 
  }

}

template <class T> void MULT_Pass(T x[], double L, int Nint, int Norder, double KNL[11], double KNSL[11]) 
{
  if(L==0.) {
    mult_kick_pass(x, Norder, KNL, KNSL); return ;
  }
  
  int i,j;
  double Lint=L/Nint;
  double knl_kick1[11],knsl_kick1[11];
  double knl_kick2[11],knsl_kick2[11];

  for(i=0;i<11;i++) {
    knl_kick1[i] =Fkick1*KNL[i]/Nint;
    knsl_kick1[i]=Fkick1*KNSL[i]/Nint;
  }
  for(i=0;i<11;i++) {
    knl_kick2[i] =Fkick2*KNL[i]/Nint;
    knsl_kick2[i]=Fkick2*KNSL[i]/Nint;
  }
  
  for(i=0;i<Nint;i++){
    DRIFT_Pass(x,Fdrift1*Lint);
    mult_kick_pass(x, Norder,knl_kick1, knsl_kick1);
    DRIFT_Pass(x,Fdrift2*Lint);
    mult_kick_pass(x, Norder,knl_kick2, knsl_kick2);
    DRIFT_Pass(x,Fdrift2*Lint);
    mult_kick_pass(x, Norder,knl_kick1, knsl_kick1);
    DRIFT_Pass(x,Fdrift1*Lint);
  }

} 

template <class T> void MULT_Pass_Radiate(T x[], double L, int Nint, int Norder, double KNL[11], double KNSL[11]) 
{
  if(L==0.) {
    mult_kick_pass(x, Norder, KNL, KNSL); return ;
  }
  
  int i,j;
  double Lint=L/Nint;
  double knl_kick1[11],knsl_kick1[11];
  double knl_kick2[11],knsl_kick2[11];
  T      x1[6], BLbrho[3];

  for(i=0;i<11;i++) {
    knl_kick1[i] =Fkick1*KNL[i]/Nint;
    knsl_kick1[i]=Fkick1*KNSL[i]/Nint;
  }
  for(i=0;i<11;i++) {
    knl_kick2[i] =Fkick2*KNL[i]/Nint;
    knsl_kick2[i]=Fkick2*KNSL[i]/Nint;
  }
  
  for(i=0;i<Nint;i++){
    for(j=0;j<6;j++) x1[j]=x[j];
    
    DRIFT_Pass(x,Fdrift1*Lint);
    mult_kick_pass(x, Norder,knl_kick1, knsl_kick1);
    DRIFT_Pass(x,Fdrift2*Lint);
    mult_kick_pass(x, Norder,knl_kick2, knsl_kick2);
    DRIFT_Pass(x,Fdrift2*Lint);
    mult_kick_pass(x, Norder,knl_kick1, knsl_kick1);
    DRIFT_Pass(x,Fdrift1*Lint);
    
    for(j=0;j<6;j++) x1[j]=(x[j]+x1[j])/2.;
    cal_Bfield(x1, Nint, Norder, 0, KNL, KNSL, BLbrho);
    radiate(x,Lint,0, BLbrho); 
  }

}

template <class T> void GMULT_Pass(T x[], double L, int Nint,  int Norder, double Angle, double E1, double E2, double KNL[11], double KNSL[11])
{
  if(L==0.) {
    mult_kick_pass(x, Norder, KNL, KNSL); return;
  }
  
  int i,j;
  double href=Angle/L;
  double Lint=L/Nint;
  double edgefocus;
  double knl_kick1[11],knsl_kick1[11];
  double knl_kick2[11],knsl_kick2[11];
  double cosAngle, sinAngle, cosAngle1, sinAngle1;
  T      x1[6], BLbrho[3];

  edgefocus = tan(E1)*href;
  x[1] = x[1]+edgefocus*x[0];
  x[3] = x[3]-edgefocus*x[2];
  
  if(GP.H_expand == true and GP.radiate == false){
    for(j=0;j<11;j++) {
      knl_kick1[j] =Fkick1*KNL[j]/Nint;
      knsl_kick1[j]=Fkick1*KNSL[j]/Nint;
    }
    for(j=0;j<11;j++) {
      knl_kick2[j] =Fkick2*KNL[j]/Nint;
      knsl_kick2[j]=Fkick2*KNSL[j]/Nint;
    }
    for(i=0;i<Nint;i++){
      DRIFT_Pass(x,Fdrift1*Lint);
      bend_mult_kick_pass(x, Fkick1*Lint, href, Norder,knl_kick1, knsl_kick1);
      DRIFT_Pass(x,Fdrift2*Lint);
      bend_mult_kick_pass(x, Fkick2*Lint, href, Norder,knl_kick2, knsl_kick2);
      DRIFT_Pass(x,Fdrift2*Lint);
      bend_mult_kick_pass(x, Fkick1*Lint, href, Norder,knl_kick1, knsl_kick1);
      DRIFT_Pass(x,Fdrift1*Lint);
    }
  }

  if(GP.H_expand == true  and GP.radiate == true){
    for(j=0;j<11;j++) {
      knl_kick1[j] =Fkick1*KNL[j]/Nint;
      knsl_kick1[j]=Fkick1*KNSL[j]/Nint;
    }
    for(j=0;j<11;j++) {
      knl_kick2[j] =Fkick2*KNL[j]/Nint;
      knsl_kick2[j]=Fkick2*KNSL[j]/Nint;
    }
    for(i=0;i<Nint;i++){
      for(j=0;j<6;j++) x1[j]=x[j];

      DRIFT_Pass(x,Fdrift1*Lint);
      bend_mult_kick_pass(x, Fkick1*Lint, href, Norder,knl_kick1, knsl_kick1);
      DRIFT_Pass(x,Fdrift2*Lint);
      bend_mult_kick_pass(x, Fkick2*Lint, href, Norder,knl_kick2, knsl_kick2);
      DRIFT_Pass(x,Fdrift2*Lint);
      bend_mult_kick_pass(x, Fkick1*Lint, href, Norder,knl_kick1, knsl_kick1);
      DRIFT_Pass(x,Fdrift1*Lint);
      
      for(j=0;j<6;j++) x1[j]=(x[j]+x1[j])/2.;
      cal_Bfield(x1, Nint, Norder, Angle, KNL, KNSL, BLbrho);
      radiate(x,Lint,href,BLbrho);
    }
  }

  if(GP.H_expand == false and GP.radiate == false){
    for(j=0;j<11;j++) {
      knl_kick1[j] =Fkick1*KNL[j]/Nint;
      knsl_kick1[j]=Fkick1*KNSL[j]/Nint;
    }
    for(j=0;j<11;j++) {
      knl_kick2[j] =Fkick2*KNL[j]/Nint;
      knsl_kick2[j]=Fkick2*KNSL[j]/Nint;
    }
    cosAngle1=cos(Fdrift1*Angle/Nint);
    sinAngle1=sin(Fdrift1*Angle/Nint);
    cosAngle=cos(Fdrift2*Angle/Nint);
    sinAngle=sin(Fdrift2*Angle/Nint);
    for(i=0;i<Nint;i++){
      sbend_exact_pass(x,Fdrift1*Lint,Fdrift1*Angle/Nint,cosAngle1,sinAngle1);
      mult_kick_pass(x, Norder,knl_kick1, knsl_kick1);
      sbend_exact_pass(x,Fdrift2*Lint,Fdrift2*Angle/Nint,cosAngle,sinAngle);
      mult_kick_pass(x, Norder,knl_kick2, knsl_kick2);
      sbend_exact_pass(x,Fdrift2*Lint,Fdrift2*Angle/Nint,cosAngle,sinAngle);
      mult_kick_pass(x, Norder,knl_kick1, knsl_kick1);
      sbend_exact_pass(x,Fdrift1*Lint,Fdrift1*Angle/Nint,cosAngle1,sinAngle1);
    }
  }

  if(GP.H_expand == false and GP.radiate == true){
    for(j=0;j<11;j++) {
      knl_kick1[j] =Fkick1*KNL[j]/Nint;
      knsl_kick1[j]=Fkick1*KNSL[j]/Nint;
    }
    for(j=0;j<11;j++) {
      knl_kick2[j] =Fkick2*KNL[j]/Nint;
      knsl_kick2[j]=Fkick2*KNSL[j]/Nint;
    }
    cosAngle1=cos(Fdrift1*Angle/Nint);
    sinAngle1=sin(Fdrift1*Angle/Nint);
    cosAngle=cos(Fdrift2*Angle/Nint);
    sinAngle=sin(Fdrift2*Angle/Nint);
    for(i=0;i<Nint;i++){
      for(j=0;j<6;j++) x1[j]=x[j];

      sbend_exact_pass(x,Fdrift1*Lint,Fdrift1*Angle/Nint,cosAngle1,sinAngle1);
      mult_kick_pass(x, Norder,knl_kick1, knsl_kick1);
      sbend_exact_pass(x,Fdrift2*Lint,Fdrift2*Angle/Nint,cosAngle,sinAngle);
      mult_kick_pass(x, Norder,knl_kick2, knsl_kick2);
      sbend_exact_pass(x,Fdrift2*Lint,Fdrift2*Angle/Nint,cosAngle,sinAngle);
      mult_kick_pass(x, Norder,knl_kick1, knsl_kick1);
      sbend_exact_pass(x,Fdrift1*Lint,Fdrift1*Angle/Nint,cosAngle1,sinAngle1);
      
      for(j=0;j<6;j++) x1[j]=(x[j]+x1[j])/2.;
      cal_Bfield(x1, Nint, Norder, Angle, KNL, KNSL, BLbrho);
      radiate(x,Lint,href,BLbrho);
    }
  }
  
  edgefocus = tan(E2)*href;
  x[1] = x[1]+edgefocus*x[0];
  x[3] = x[3]-edgefocus*x[2];
}


template <class T> void SBENDMULT_Pass(T x[], double L, int Nint,  int Norder, double Angle, double E1, double E2, double KNL[11], double KNSL[11])
{
  if(L==0.) {
    mult_kick_pass(x, Norder, KNL, KNSL); return;
  }
  
  int i,j;
  double href=Angle/L;
  double Lint=L/Nint;
  double edgefocus;
  double knl_kick1[11],knsl_kick1[11];
  double knl_kick2[11],knsl_kick2[11];
  double cosAngle, sinAngle, cosAngle1, sinAngle1;
  T      x1[6], BLbrho[3];

  edgefocus = tan(E1)*href;
  x[1] = x[1]+edgefocus*x[0];
  x[3] = x[3]-edgefocus*x[2];
  
  if(GP.H_expand == true and GP.radiate == false){
    for(j=0;j<11;j++) {
      knl_kick1[j] =Fkick1*KNL[j]/Nint;
      knsl_kick1[j]=Fkick1*KNSL[j]/Nint;
    }
    for(j=0;j<11;j++) {
      knl_kick2[j] =Fkick2*KNL[j]/Nint;
      knsl_kick2[j]=Fkick2*KNSL[j]/Nint;
    }
    for(i=0;i<Nint;i++){
      DRIFT_Pass(x,Fdrift1*Lint);
      bend_mult_kick_pass(x, Fkick1*Lint, href, Norder,knl_kick1, knsl_kick1);
      DRIFT_Pass(x,Fdrift2*Lint);
      bend_mult_kick_pass(x, Fkick2*Lint, href, Norder,knl_kick2, knsl_kick2);
      DRIFT_Pass(x,Fdrift2*Lint);
      bend_mult_kick_pass(x, Fkick1*Lint, href, Norder,knl_kick1, knsl_kick1);
      DRIFT_Pass(x,Fdrift1*Lint);
    }
  }

  if(GP.H_expand == true  and GP.radiate == true){
    for(j=0;j<11;j++) {
      knl_kick1[j] =Fkick1*KNL[j]/Nint;
      knsl_kick1[j]=Fkick1*KNSL[j]/Nint;
    }
    for(j=0;j<11;j++) {
      knl_kick2[j] =Fkick2*KNL[j]/Nint;
      knsl_kick2[j]=Fkick2*KNSL[j]/Nint;
    }
    for(i=0;i<Nint;i++){
      for(j=0;j<6;j++) x1[j]=x[j];

      DRIFT_Pass(x,Fdrift1*Lint);
      bend_mult_kick_pass(x, Fkick1*Lint, href, Norder,knl_kick1, knsl_kick1);
      DRIFT_Pass(x,Fdrift2*Lint);
      bend_mult_kick_pass(x, Fkick2*Lint, href, Norder,knl_kick2, knsl_kick2);
      DRIFT_Pass(x,Fdrift2*Lint);
      bend_mult_kick_pass(x, Fkick1*Lint, href, Norder,knl_kick1, knsl_kick1);
      DRIFT_Pass(x,Fdrift1*Lint);
      
      for(j=0;j<6;j++) x1[j]=(x[j]+x1[j])/2.;
      cal_Bfield(x1, Nint, Norder, Angle, KNL, KNSL, BLbrho);
      radiate(x,Lint,href,BLbrho);
    }
  }

  if(GP.H_expand == false and GP.radiate == false){
    for(j=0;j<11;j++) {
      knl_kick1[j] =Fkick1*KNL[j]/Nint;
      knsl_kick1[j]=Fkick1*KNSL[j]/Nint;
    }
    for(j=0;j<11;j++) {
      knl_kick2[j] =Fkick2*KNL[j]/Nint;
      knsl_kick2[j]=Fkick2*KNSL[j]/Nint;
    }
    cosAngle1=cos(Fdrift1*Angle/Nint);
    sinAngle1=sin(Fdrift1*Angle/Nint);
    cosAngle=cos(Fdrift2*Angle/Nint);
    sinAngle=sin(Fdrift2*Angle/Nint);
    for(i=0;i<Nint;i++){
      sbend_exact_pass(x,Fdrift1*Lint,Fdrift1*Angle/Nint,cosAngle1,sinAngle1);
      mult_kick_pass(x, Norder,knl_kick1, knsl_kick1);
      sbend_exact_pass(x,Fdrift2*Lint,Fdrift2*Angle/Nint,cosAngle,sinAngle);
      mult_kick_pass(x, Norder,knl_kick2, knsl_kick2);
      sbend_exact_pass(x,Fdrift2*Lint,Fdrift2*Angle/Nint,cosAngle,sinAngle);
      mult_kick_pass(x, Norder,knl_kick1, knsl_kick1);
      sbend_exact_pass(x,Fdrift1*Lint,Fdrift1*Angle/Nint,cosAngle1,sinAngle1);
    }
  }

  if(GP.H_expand == false and GP.radiate == true){
    for(j=0;j<11;j++) {
      knl_kick1[j] =Fkick1*KNL[j]/Nint;
      knsl_kick1[j]=Fkick1*KNSL[j]/Nint;
    }
    for(j=0;j<11;j++) {
      knl_kick2[j] =Fkick2*KNL[j]/Nint;
      knsl_kick2[j]=Fkick2*KNSL[j]/Nint;
    }
    cosAngle1=cos(Fdrift1*Angle/Nint);
    sinAngle1=sin(Fdrift1*Angle/Nint);
    cosAngle=cos(Fdrift2*Angle/Nint);
    sinAngle=sin(Fdrift2*Angle/Nint);
    for(i=0;i<Nint;i++){
      for(j=0;j<6;j++) x1[j]=x[j];

      sbend_exact_pass(x,Fdrift1*Lint,Fdrift1*Angle/Nint,cosAngle1,sinAngle1);
      mult_kick_pass(x, Norder,knl_kick1, knsl_kick1);
      sbend_exact_pass(x,Fdrift2*Lint,Fdrift2*Angle/Nint,cosAngle,sinAngle);
      mult_kick_pass(x, Norder,knl_kick2, knsl_kick2);
      sbend_exact_pass(x,Fdrift2*Lint,Fdrift2*Angle/Nint,cosAngle,sinAngle);
      mult_kick_pass(x, Norder,knl_kick1, knsl_kick1);
      sbend_exact_pass(x,Fdrift1*Lint,Fdrift1*Angle/Nint,cosAngle1,sinAngle1);
      
      for(j=0;j<6;j++) x1[j]=(x[j]+x1[j])/2.;
      cal_Bfield(x1, Nint, Norder, Angle, KNL, KNSL, BLbrho);
      radiate(x,Lint,href,BLbrho);
    }
  }
  
  edgefocus = tan(E2)*href;
  x[1] = x[1]+edgefocus*x[0];
  x[3] = x[3]-edgefocus*x[2];
}

template <class T> void SMULT_Pass(T x[], double L, int Nint,  int Norder, double Angle, double E1, double E2, double KNL[11], double KNSL[11])
{
  if(L==0.) 
    {
      mult_kick_pass(x, Norder, KNL, KNSL); return;
    }

  int i,j;
  double href=Angle/L;
  double Lint=L/Nint;
  double edgefocus;
  double knl_kick1[11],knsl_kick1[11];
  double knl_kick2[11],knsl_kick2[11];
  double cosAngle, sinAngle, cosAngle1, sinAngle1;
  T      x1[6], BLbrho[3];

  edgefocus = tan(E1)*href;
  x[1] = x[1]+edgefocus*x[0];
  x[3] = x[3]-edgefocus*x[2];
  
  if(GP.H_expand == true and GP.radiate == false){
    for(j=0;j<11;j++) {
      knl_kick1[j] =Fkick1*KNL[j]/Nint;
      knsl_kick1[j]=Fkick1*KNSL[j]/Nint;
    }
    for(j=0;j<11;j++) {
      knl_kick2[j] =Fkick2*KNL[j]/Nint;
      knsl_kick2[j]=Fkick2*KNSL[j]/Nint;
    }
    for(i=0;i<Nint;i++){
      DRIFT_Pass(x,Fdrift1*Lint);
      bend_mult_kick_pass(x, Fkick1*Lint, href, Norder, knl_kick1, knsl_kick1);
      DRIFT_Pass(x,Fdrift2*Lint);
      bend_mult_kick_pass(x, Fkick2*Lint, href, Norder, knl_kick2, knsl_kick2);
      DRIFT_Pass(x,Fdrift2*Lint);
      bend_mult_kick_pass(x, Fkick1*Lint, href, Norder, knl_kick1, knsl_kick1);
      DRIFT_Pass(x,Fdrift1*Lint);
    }
  }

  if(GP.H_expand == true  and GP.radiate == true){
    for(j=0;j<11;j++) {
      knl_kick1[j] =Fkick1*KNL[j]/Nint;
      knsl_kick1[j]=Fkick1*KNSL[j]/Nint;
    }
    for(j=0;j<11;j++) {
      knl_kick2[j] =Fkick2*KNL[j]/Nint;
      knsl_kick2[j]=Fkick2*KNSL[j]/Nint;
    }
    for(i=0;i<Nint;i++){
      for(j=0;j<6;j++) x1[j]=x[j];

      DRIFT_Pass(x,Fdrift1*Lint);
      bend_mult_kick_pass(x, Fkick1*Lint, href, Norder,knl_kick1, knsl_kick1);
      DRIFT_Pass(x,Fdrift2*Lint);
      bend_mult_kick_pass(x, Fkick2*Lint, href, Norder,knl_kick2, knsl_kick2);
      DRIFT_Pass(x,Fdrift2*Lint);
      bend_mult_kick_pass(x, Fkick1*Lint, href, Norder,knl_kick1, knsl_kick1);
      DRIFT_Pass(x,Fdrift1*Lint);
      
      for(j=0;j<6;j++) x1[j]=(x[j]+x1[j])/2.;
      cal_Bfield(x1, Nint, Norder, Angle, KNL, KNSL, BLbrho);
      radiate(x,Lint,href,BLbrho);
    }
  }

  if(GP.H_expand == false and GP.radiate == false){
    for(j=0;j<11;j++) {
      knl_kick1[j] =Fkick1*KNL[j]/Nint;
      knsl_kick1[j]=Fkick1*KNSL[j]/Nint;
    }
    for(j=0;j<11;j++) {
      knl_kick2[j] =Fkick2*KNL[j]/Nint;
      knsl_kick2[j]=Fkick2*KNSL[j]/Nint;
    }
    cosAngle1=cos(Fdrift1*Angle/Nint);
    sinAngle1=sin(Fdrift1*Angle/Nint);
    cosAngle=cos(Fdrift2*Angle/Nint);
    sinAngle=sin(Fdrift2*Angle/Nint);
    for(i=0;i<Nint;i++){
      sbend_exact_pass(x,Fdrift1*Lint,Fdrift1*Angle/Nint,cosAngle1,sinAngle1);
      mult_kick_pass(x, Norder,knl_kick1, knsl_kick1);
      sbend_exact_pass(x,Fdrift2*Lint,Fdrift2*Angle/Nint,cosAngle,sinAngle);
      mult_kick_pass(x, Norder,knl_kick2, knsl_kick2);
      sbend_exact_pass(x,Fdrift2*Lint,Fdrift2*Angle/Nint,cosAngle,sinAngle);
      mult_kick_pass(x, Norder,knl_kick1, knsl_kick1);
      sbend_exact_pass(x,Fdrift1*Lint,Fdrift1*Angle/Nint,cosAngle1,sinAngle1);
    }
  }

  if(GP.H_expand == false and GP.radiate == true){
    for(j=0;j<11;j++) {
      knl_kick1[j] =Fkick1*KNL[j]/Nint;
      knsl_kick1[j]=Fkick1*KNSL[j]/Nint;
    }
    for(j=0;j<11;j++) {
      knl_kick2[j] =Fkick2*KNL[j]/Nint;
      knsl_kick2[j]=Fkick2*KNSL[j]/Nint;
    }
    cosAngle1=cos(Fdrift1*Angle/Nint);
    sinAngle1=sin(Fdrift1*Angle/Nint);
    cosAngle=cos(Fdrift2*Angle/Nint);
    sinAngle=sin(Fdrift2*Angle/Nint);
    for(i=0;i<Nint;i++){
      for(j=0;j<6;j++) x1[j]=x[j];

      sbend_exact_pass(x,Fdrift1*Lint,Fdrift1*Angle/Nint,cosAngle1,sinAngle1);
      mult_kick_pass(x, Norder,knl_kick1, knsl_kick1);
      sbend_exact_pass(x,Fdrift2*Lint,Fdrift2*Angle/Nint,cosAngle,sinAngle);
      mult_kick_pass(x, Norder,knl_kick2, knsl_kick2);
      sbend_exact_pass(x,Fdrift2*Lint,Fdrift2*Angle/Nint,cosAngle,sinAngle);
      mult_kick_pass(x, Norder,knl_kick1, knsl_kick1);
      sbend_exact_pass(x,Fdrift1*Lint,Fdrift1*Angle/Nint,cosAngle1,sinAngle1);
      
      for(j=0;j<6;j++) x1[j]=(x[j]+x1[j])/2.;
      cal_Bfield(x1, Nint, Norder, Angle, KNL, KNSL, BLbrho);
      radiate(x,Lint,href,BLbrho);
    }
  }
  
  edgefocus = tan(E2)*href;
  x[1] = x[1]+edgefocus*x[0];
  x[3] = x[3]-edgefocus*x[2];
}

template <class T> void GSBENDMULT_Pass(T x[], double L, int Nint,  int Norder, double Angle, double K0L, double E1, double E2, double KNL[11], double KNSL[11])
{
  if(L==0.) {
    mult_kick_pass(x, Norder, KNL, KNSL); return;
  }
  
  int i,j;
  double href=Angle/L;
  double hreal=K0L/L;
  
  double Lint=L/Nint;
  double edgefocus;
  double knl_kick1[11],knsl_kick1[11];
  double knl_kick2[11],knsl_kick2[11];
  double cosAngle, sinAngle, cosAngle1, sinAngle1;
  T      x1[6], BLbrho[3];

  edgefocus = tan(E1)*hreal;
  x[1] = x[1]+edgefocus*x[0];
  x[3] = x[3]-edgefocus*x[2];
  
  if(GP.H_expand == true and GP.radiate == false){
    for(j=0;j<11;j++) {
      knl_kick1[j] =Fkick1*KNL[j]/Nint;
      knsl_kick1[j]=Fkick1*KNSL[j]/Nint;
    }
    for(j=0;j<11;j++) {
      knl_kick2[j] =Fkick2*KNL[j]/Nint;
      knsl_kick2[j]=Fkick2*KNSL[j]/Nint;
    }
    for(i=0;i<Nint;i++){
      DRIFT_Pass(x,Fdrift1*Lint);
      general_bend_mult_kick_pass(x, Fkick1*Lint, href, hreal, Norder, knl_kick1, knsl_kick1);
      DRIFT_Pass(x,Fdrift2*Lint);
      general_bend_mult_kick_pass(x, Fkick2*Lint, href, hreal, Norder, knl_kick2, knsl_kick2);
      DRIFT_Pass(x,Fdrift2*Lint);
      general_bend_mult_kick_pass(x, Fkick1*Lint, href, hreal, Norder, knl_kick1, knsl_kick1);
      DRIFT_Pass(x,Fdrift1*Lint);
    }
  }

  if(GP.H_expand == true  and GP.radiate == true){
    for(j=0;j<11;j++) {
      knl_kick1[j] =Fkick1*KNL[j]/Nint;
      knsl_kick1[j]=Fkick1*KNSL[j]/Nint;
    }
    for(j=0;j<11;j++) {
      knl_kick2[j] =Fkick2*KNL[j]/Nint;
      knsl_kick2[j]=Fkick2*KNSL[j]/Nint;
    }
    for(i=0;i<Nint;i++){
      for(j=0;j<6;j++) x1[j]=x[j];

      DRIFT_Pass(x,Fdrift1*Lint);
      general_bend_mult_kick_pass(x, Fkick1*Lint, href, hreal, Norder,knl_kick1, knsl_kick1); 
      DRIFT_Pass(x,Fdrift2*Lint);
      general_bend_mult_kick_pass(x, Fkick2*Lint, href, hreal, Norder,knl_kick2, knsl_kick2);
      DRIFT_Pass(x,Fdrift2*Lint);
      general_bend_mult_kick_pass(x, Fkick1*Lint, href, hreal, Norder,knl_kick1, knsl_kick1);
      DRIFT_Pass(x,Fdrift1*Lint);
      
      for(j=0;j<6;j++) x1[j]=(x[j]+x1[j])/2.;
      //cal_Bfield(x1, Nint, Norder, Angle, KNL, KNSL, BLbrho); // only dipole component radiates, need to check K0L effect
      BLbrho[0] = 0.;   BLbrho[1] = K0L/Nint;  BLbrho[2] = 0.; 
      radiate(x,Lint,href,BLbrho);
    }
  }

  if(GP.H_expand == false and GP.radiate == false){
    for(j=0;j<11;j++) {
      knl_kick1[j] =Fkick1*KNL[j]/Nint;
      knsl_kick1[j]=Fkick1*KNSL[j]/Nint;
    }
    for(j=0;j<11;j++) {
      knl_kick2[j] =Fkick2*KNL[j]/Nint;
      knsl_kick2[j]=Fkick2*KNSL[j]/Nint;
    }
    cosAngle1=cos(Fdrift1*Angle/Nint);
    sinAngle1=sin(Fdrift1*Angle/Nint);
    cosAngle=cos(Fdrift2*Angle/Nint);
    sinAngle=sin(Fdrift2*Angle/Nint);
    for(i=0;i<Nint;i++){
      general_sbend_exact_pass(x,Fdrift1*Lint,Fdrift1*Angle/Nint,Fdrift1*K0L/Nint, cosAngle1, sinAngle1);
      mult_kick_pass(x, Norder,knl_kick1, knsl_kick1);
      general_sbend_exact_pass(x,Fdrift2*Lint,Fdrift2*Angle/Nint,Fdrift2*K0L/Nint, cosAngle, sinAngle);
      mult_kick_pass(x, Norder,knl_kick2, knsl_kick2);
      general_sbend_exact_pass(x,Fdrift2*Lint,Fdrift2*Angle/Nint,Fdrift2*K0L/Nint, cosAngle, sinAngle);
      mult_kick_pass(x, Norder,knl_kick1, knsl_kick1);
      general_sbend_exact_pass(x,Fdrift1*Lint,Fdrift1*Angle/Nint,Fdrift1*K0L/Nint, cosAngle1, sinAngle1);
    }
  }

  if(GP.H_expand == false and GP.radiate == true){
    for(j=0;j<11;j++) {
      knl_kick1[j] =Fkick1*KNL[j]/Nint;
      knsl_kick1[j]=Fkick1*KNSL[j]/Nint;
    }
    for(j=0;j<11;j++) {
      knl_kick2[j] =Fkick2*KNL[j]/Nint;
      knsl_kick2[j]=Fkick2*KNSL[j]/Nint;
    }
    cosAngle1=cos(Fdrift1*Angle/Nint);
    sinAngle1=sin(Fdrift1*Angle/Nint);
    cosAngle=cos(Fdrift2*Angle/Nint);
    sinAngle=sin(Fdrift2*Angle/Nint);
    for(i=0;i<Nint;i++){
      for(j=0;j<6;j++) x1[j]=x[j];

      general_sbend_exact_pass(x,Fdrift1*Lint,Fdrift1*Angle/Nint,Fdrift1*K0L/Nint, cosAngle1, sinAngle1);
      mult_kick_pass(x, Norder,knl_kick1, knsl_kick1);
      general_sbend_exact_pass(x,Fdrift2*Lint,Fdrift2*Angle/Nint,Fdrift2*K0L/Nint, cosAngle, sinAngle);
      mult_kick_pass(x, Norder,knl_kick2, knsl_kick2);
      general_sbend_exact_pass(x,Fdrift2*Lint,Fdrift2*Angle/Nint,Fdrift2*K0L/Nint, cosAngle, sinAngle);
      mult_kick_pass(x, Norder,knl_kick1, knsl_kick1);
      general_sbend_exact_pass(x,Fdrift1*Lint,Fdrift1*Angle/Nint,Fdrift1*K0L/Nint, cosAngle1, sinAngle1);      
      
      for(j=0;j<6;j++) x1[j]=(x[j]+x1[j])/2.;
      //cal_Bfield(x1, Nint, Norder, Angle, KNL, KNSL, BLbrho);  needd  to  considerd K0L
      BLbrho[0] = 0.;  BLbrho[1] = K0L/Nint;  BLbrho[2] = 0.; 
      radiate(x,Lint,href,BLbrho);
    }
  }
  
  edgefocus = tan(E2)*hreal;
  x[1] = x[1]+edgefocus*x[0];
  x[3] = x[3]-edgefocus*x[2];
}

template <class T> void SOLEN_Pass_Old(T x[], double L, int Nint, double KS)
{
  int i, j, k;	  
  double g=KS/2, theta=KS*L/2, costheta= cos( theta ), sintheta=  sin( theta ) ;
  double cctheta = costheta * costheta, sstheta=sintheta * sintheta, cstheta=sintheta * costheta;
  double M[36];
  T      xtemp[6];
  
  if( KS == 0. ) {
    DRIFT_Pass(x, L);}
  else{
    //----calculate the matrix
    for(i=0;i<36;i++) M[i]= 0.0;

    
    M[0*6+0] = costheta * costheta;
    M[0*6+1] = sintheta * costheta / g; 
    M[0*6+2] = sintheta * costheta;
    M[0*6+3] = sintheta *  sintheta / g;
    
    M[1*6+0] = -g * sintheta * costheta ;
    M[1*6+1] =      costheta * costheta ;
    M[1*6+2] = -g*  sintheta * sintheta ;
    M[1*6+3] =      sintheta *  costheta;
    
    M[2*6+0] = -sintheta * costheta     ;
    M[2*6+1] = -sintheta * sintheta / g ;
    M[2*6+2] =  costheta * costheta     ;
    M[2*6+3] =  sintheta * costheta / g ; 
    
    M[3*6+0] = g *  sintheta * sintheta ;
    M[3*6+1] =     -sintheta * costheta ;
    M[3*6+2] = -g * sintheta * costheta ;
    M[3*6+3] =      costheta * costheta ;
    
    M[4*6+4] = 1.0;
    M[5*6+5] = 1.0;
    
    /*---
    M[0*6+0] = cctheta;
    M[0*6+1] = cstheta / g; 
    M[0*6+2] = cstheta;
    M[0*6+3] = sstheta / g;
    
    M[1*6+0] = -g * cstheta ;
    M[1*6+1] =      cctheta ;
    M[1*6+2] = -g * sstheta ;
    M[1*6+3] =      cstheta;
    
    M[2*6+0] = -cstheta     ;
    M[2*6+1] = -sstheta / g ;
    M[2*6+2] =  cctheta     ;
    M[2*6+3] =  cstheta / g ; 
    
    M[3*6+0] =  g * sstheta ;
    M[3*6+1] =     -cstheta ;
    M[3*6+2] = -g * cstheta ;
    M[3*6+3] =      cctheta ;
    
    M[4*6+4] = 1.0;
    M[5*6+5] = 1.0;
    ---*/

    //---transfer
    x[1]=x[1]/(1+x[5]);  // x'=px / (1+delta), is it needed ?
    x[3]=x[3]/(1+x[5]);
    for(i=0;i<6;i++) xtemp[i]=x[i];
    
    for(k=0;k<6;k++) {
      x[k]=0.;
      for(j=0;j<6;j++) x[k]= x[k]+ M[k*6+j]*xtemp[j];
    }
    x[1]=x[1]*(1+x[5]); // px=x' * (1+delta), is it needed ?
    x[3]=x[3]*(1+x[5]);
    
  }
  
}

template <class T> void SOLEN_Pass(T x[], double L, int Nint, double KS)
{
  int i, j, k;
  T   xtemp[6], ps;
  T   gamma1, delta1, beta1, vz;
  T   Bz= KS * GP.brho, K, Lint=L/Nint, theta, sintheta, costheta, m44[16];

  if(KS == 0 ) {
    DRIFT_Pass(x, L);
    return;
  }
  
  //----switch to (x, x', y, y')
  delta1 = sqrt(1.0 + 2*x[5]/GP.beta+x[5]*x[5]) -1.0;
  gamma1 = GP.gamma+sqrt(GP.gamma*GP.gamma-1.)*x[5];
  beta1  = sqrt(1.0-1.0/gamma1/gamma1);
  
  ps=sqrt( ( 1+delta1 ) * (1+delta1 ) - x[1]*x[1] - x[3]*x[3] ); 
  x[1]=x[1]/ps; 
  x[3]=x[3]/ps;

  //----entrance
  vz =  2.99792458e8 * beta1 / sqrt( 1.0 + x[1]*x[1] +x[3]*x[3] );  
  K  =  GP.q * Bz / 2.0 / gamma1 / GP.m0 / vz;
  x[1] = x[1] + x[2] * K ;
  x[3] = x[3] - x[0] * K ;

  //---body
  theta =  -2.0 * K * L / Nint;

  sintheta = sin(theta); costheta = cos(theta); 
  m44[0*4+0] = 1.0 ;   m44[0*4+1] = Lint*sintheta/theta;      m44[0*4+2] = 0. ;   m44[0*4+3] = -(1-costheta)*Lint/theta;
  m44[1*4+0] = 0. ;    m44[1*4+1] = costheta;                 m44[1*4+2] = 0. ;   m44[1*4+3] = -sintheta;
  m44[2*4+0] = 0. ;    m44[2*4+1] = (1-costheta)*Lint/theta;  m44[2*4+2] = 1. ;   m44[2*4+3] = Lint*sintheta/theta;
  m44[3*4+0] = 0. ;    m44[3*4+1] = sintheta;                 m44[3*4+2] = 0. ;   m44[3*4+3] = costheta;

  for(i=0;i<Nint; i++){
    for(k=0;k<4;k++) xtemp[k]=x[k];
    for(k=0;k<4;k++){
      x[k]=0.;
      for(j=0;j<4;j++) x[k]= x[k]+ m44[k*4+j]*xtemp[j];
    }
  }

  //----exit
  x[1] = x[1] - x[2] * K ;
  x[3] = x[3] + x[0] * K ;
      
  //---switch to (x, px, y, py)
  x[1]=gamma1* GP.m0 * x[1] * vz / GP.p; 
  x[3]=gamma1* GP.m0 * x[3] * vz / GP.p; 
  x[4]=x[4] - 2.99792458e8 * ( L/vz - L/GP.beta/ 2.99792458e8);
 
}

template <class T> void WIGGLER_Pass(T x[], double L, int Nint, double B0, double KX, double KZ, double PHIZ0)
{
  int i, j, k;
  T xtemp[6];
  
  if(B0 == 0 ) {
    DRIFT_Pass(x, L);}
  else{
    wiggler_pass(x, Nint, L, B0, KX, KZ, PHIZ0);
  }
}

template<class T> void MATRIX_Pass( T x[], double L, double XCO_IN[6], double XCO_OUT[6], double M66[36])
{
  int i,j;
  T x1[6];

  for(i=0;i<6;i++) x[i]=x[i]-XCO_IN[i];

  for(i=0;i<6;i++) {
    x1[i]=0.0;
    for(j=0;j<6;j++)  x1[i] = x1[i] + M66[i*6+j] * x[j];
  }   
  for(i=0;i<6;i++)    x[i]=x1[i];
  
  for(i=0;i<6;i++) x[i]=x[i]+XCO_OUT[i];
}

template <class T> void KICK_Pass(T x[], double L, double HKICK, double VKICK) 
{
  DRIFT_Pass(x,L/2.0); 
  x[1]=x[1]+HKICK;     
  x[3]=x[3]+VKICK;
  DRIFT_Pass(x,L/2.0);
}

template <class T>  void   ACMULT_Pass(T  x[], double L, int Norder, double KL, double KSL, int TTURNS, double PHI0)
{
  int i;
  int fac=1;
  T Xn, Yn, Xn0, Yn0;
  double  KNL, KNSL;
  T By, Bx;

  DRIFT_Pass(x, L/2.0);
  
  if(Norder==0) {
    By=KL;
    Bx=KSL;
  }
  else{
    KNL=KL;
    KNSL=KSL;
    By=0.;  Bx=0.;  Xn=1.;  Yn=0.;
    for(i=1;i<Norder+1;i++){
      Xn0=Xn;
      Yn0=Yn;
      Xn=Xn0*x[x_]-Yn0*x[y_];
      Yn=Xn0*x[y_]+Yn0*x[x_];
      fac=fac*i;
    }
    By=By+(KNL*Xn-KNSL*Yn)/fac;
    Bx=Bx+(KNL*Yn+KNSL*Xn)/fac;
  }
  
  x[px_]=x[px_]- By * sin( 2.0 * PI * GP.turn / TTURNS + PHI0);
  x[py_]=x[py_]+ Bx * sin( 2.0 * PI * GP.turn / TTURNS + PHI0);
  
  DRIFT_Pass(x, L/2.0);
}

template <class T> void ACDIP_Pass(T x[], double L, double HKICKMAX, double VKICKMAX, double NUD, double TURNS, double TURNE, double PHID) 
{
  if( GP.turn < TURNS ) {
    DRIFT_Pass(x, L);
  }
  else if ( GP.turn >= TURNS  and GP.turn <= TURNE  ){
    DRIFT_Pass(x, L/2.);
    x[1] += ((GP.turn-TURNS)*1.0 * HKICKMAX /(TURNE-TURNS )) * sin( 2.0* 3.14159265*NUD*(GP.turn-TURNS)+PHID);
    x[3] += ((GP.turn-TURNS)*1.0 * VKICKMAX /(TURNE-TURNS )) * sin( 2.0* 3.14159265*NUD*(GP.turn-TURNS)+PHID);
    DRIFT_Pass(x, L/2.);
  }
  else {
    DRIFT_Pass(x, L/2.);
    x[1] +=  HKICKMAX* sin( 2.0* 3.14159265 * NUD *(GP.turn-TURNS) +PHID );
    x[3] +=  VKICKMAX* sin( 2.0* 3.14159265 * NUD *(GP.turn-TURNS) +PHID );
    DRIFT_Pass(x, L/2.);
  }
}

void  RFCAV_Pass_old(double x[], double L, double VRF, double FRF, double PHASE0)
{
  if(VRF ==0.){
    DRIFT_Pass(x, L); return ;
  }
  
  DRIFT_Pass(x, L/2.);

  if(GP.dgamma == 0. and GP.radiate == false  ){  //   situation: at store for hadrons, PHASE0 = 0 for sure
    if(GP.gamma > GP.gammat ){
      x[pt_] =x[pt_] + abs(VRF*GP.Q)*sin(2.0*PI*FRF*x[z_]/2.99792458e8)/GP.beta/GP.energy;
    }
    else{
      x[pt_] =x[pt_] + abs(VRF*GP.Q)*sin(-2.0*PI*FRF*x[z_]/2.99792458e8  + 0 )/GP.beta/GP.energy;
    }
  }
  else{                            //  situation: on acceleration or with SR radiation with "ONE" RF cavity
    double dE0, dE1;               //  energy increase for reference dE0 and test particle dE1
    double gamma0, gamma1, beta1;  //  all for reference particle: this turn and 1 turn after
    double ratio;                  //  for test particle: to scale px, py
    
    dE0 = GP.dgamma * GP.energy0  + GP.U0rad;  // energy compensation for reference particle in one turn:  acceleration + SR energy loss 
    if(  dE0  >  abs( VRF*GP.Q) ){
      cout<<"RF can't accelerate so fast. exit. "<<endl;
      exit(0);
    }
    else{
      PHASE0 = asin( dE0 / ( VRF*GP.Q))  ;    // PHASE0 always between 0 and PI/2
    }

    gamma0 =GP.gamma;                   //  gamma for reference particle before entering RF
    gamma1 =GP.gamma + GP.dgamma ;      //  gamma for reference particle at exit of RF
    beta1  =sqrt(1.-1.0/gamma1/gamma1); //  beta for reference particle at exit of RF 
    FRF    =GP.harm * (2.99792458e8 * beta1 / GP.circumference);   //  NOTE: we keep harmonic number constant, RF freq. decided on fly
    
    if(gamma1 > GP.gammat ){            // this paragraph  to calulate test particle's energy gain from RF 
      dE1 =  abs(VRF*GP.Q)*sin(-2.0*PI*FRF*x[z_]/2.99792458e8  +  PI-PHASE0);
    }
    else{
      dE1 =  abs(VRF*GP.Q)*sin(-2.0*PI*FRF*x[z_]/2.99792458e8  +  PHASE0);
    }

    x[pt_]= (x[pt_]*GP.energy0*sqrt(gamma0*gamma0-1.0) + GP.U0rad +  dE1 - dE0 )/(GP.energy0*sqrt(gamma1*gamma1-1.0));
    ratio=sqrt(gamma0 * gamma0 - 1.)/sqrt(gamma1 * gamma1 - 1.); 
    x[px_]=x[px_] * ratio;
    x[py_]=x[py_] * ratio;

    Set_RefPartEnergy(gamma1, GP.energy0, GP.Q, GP.A);
  }

  DRIFT_Pass(x, L/2.);
}

void  RFCAV_Pass(double x[], double L, double VRF, double FRF, double PHASE0)
{
  if(VRF ==0.){
    DRIFT_Pass(x, L); return ;
  }
  
  DRIFT_Pass(x, L/2.);
  x[pt_] =x[pt_] + abs(VRF*GP.Q)*sin(-2.0*PI*FRF*x[z_]/2.99792458e8  +  PI-PHASE0 )/GP.beta/GP.energy;
  //x[pt_] =x[pt_] + abs(VRF*GP.Q)*sin(2.0*PI*FRF*x[z_]/2.99792458e8)/GP.beta/GP.energy;
  DRIFT_Pass(x, L/2.);
}

void  CRABRF_Pass(double x[], double L, double VRF, double FRF, double PHASE0)
{
  if(VRF ==0.){
    DRIFT_Pass(x, L); 
  }
  else{
    DRIFT_Pass(x, L/2.);
    x[px_]=x[px_]- (VRF*GP.Q)*sin(2.0*PI*FRF*x[z_]/2.99792458e8 + PHASE0)/GP.beta/GP.energy;
    x[pt_]=x[pt_]- (VRF*GP.Q)*cos(2.0*PI*FRF*x[z_]/2.99792458e8 + PHASE0)*x[0]*(2.0*PI*FRF/2.99792458e8)/GP.beta/GP.energy;
    DRIFT_Pass(x, L/2.);
  }
}

void  CCMULT_Pass(double x[], double L, double VRF, double FRF, double PHASE0, double B2, double A2, double B3, double A3, double B4, double A4, double B5, double A5 )
// B2, A2 as quadrupole components
{
  if(VRF ==0.){
    DRIFT_Pass(x, L); 
  }
  else{
    DRIFT_Pass(x, L/2.);

    //--crab multipole kicks
    double c =2.99792458e8;
    double scale= 2*PI*FRF/c/GP.brho ;
    
    double sinangle  = sin(2.0*PI*FRF*x[z_]/2.99792458e8 + PHASE0);
    double cosangle  = cos(2.0*PI*FRF*x[z_]/2.99792458e8 + PHASE0);
    
    double x2=x[x_]*x[x_], x3=x2*x[x_],  x4=x3*x[x_], x5=x4*x[x_];
    double y2=x[y_]*x[y_], y3=y2*x[y_],  y4=y3*x[y_], y5=y4*x[y_];
    double xy=x[x_]*x[y_], x2y=x2*x[y_], xy2=x[x_]*y2, x2y2=x[x_]*xy2;
    double x3y=x[y_]*x3,  xy3=x[x_]*y3;
    double x2y3=x2y2*x[y_],  x3y2=x2y2*x[x_], x4y=x4*x[y_], xy4=y4*x[x_];
    
    x[px_] = x[px_] + (B2*x[x_] + A2*x[y_])*cosangle / GP.brho;
    x[py_] = x[py_] + (A2*x[x_] - B2*x[y_])*cosangle / GP.brho;
    x[pt_] = x[pt_] + (B2*(x2-y2)*sinangle + A2*xy*sinangle)*scale;

    x[px_] = x[px_] + (B3*(x2-y2) + 2*A3*xy)*cosangle/GP.brho;
    x[py_] = x[py_] + (2*B3*xy + A3*(x2-y2) )*cosangle/GP.brho;  
    x[pt_] = x[pt_] + (B3*(x3-3*xy2)*sinangle + A3*(3*x2y-y3)*sinangle)*scale/3.;

    x[px_] = x[px_] + ( B4*(x3-3*xy2) - A4*(y3-3*x2y) )*cosangle/GP.brho;
    x[py_] = x[py_] + ( B4*(y3-3*x2y) + A4*(x3-3*xy2) )*cosangle/GP.brho;
    x[pt_] = x[pt_] + ( B4* (x4-6*x2*y2 +y4 )*sinangle + A4 *(4*x3y - 4* xy3)*sinangle )*scale/4.;

    x[px_] = x[px_] + ( B5*(x4-6*x2y2+y4) + A5*(4*x3y-4*xy3) )*cosangle/GP.brho;
    x[py_] = x[py_] + ( B5*(4*xy3-4*x3y) + A5*(x4-6*x2y2+y4) )*cosangle/GP.brho;
    x[pt_] = x[pt_] + ( B5* (x5-10*x3y2 +5*xy4 )*sinangle + A5 *(5*x4y - 10*x2y3+y5)*sinangle )*scale/5.;

    DRIFT_Pass(x, L/2.);
  }
}

void  LBT_Pass(double  x[], double theta)
{
  double sintheta=sin(theta), costheta=cos(theta), tantheta=sintheta/costheta;  
  double x0, px0, y0, py0, z0, pz0, ps0, h0;
  double x_star, px_star, y_star, py_star, z_star, pz_star, ps_star, h_star;

  x0= x[0];   px0=x[1];   y0= x[2];
  py0=x[3];   z0= x[4];   pz0=x[5];
  ps0= sqrt( (1+pz0 )*(1+pz0)- px0*px0 - py0*py0 );  
  h0 = (1+pz0)- ps0;

  px_star = ( px0 - tantheta * h0 )  / cos(theta );
  py_star = py0 / costheta;
  pz_star = pz0 - tantheta * px0 + tantheta * tantheta * h0;

  ps_star = sqrt( (1+pz_star )*(1+pz_star)- px_star * px_star - py_star * py_star ); 
  h_star  = (1+pz_star) - ps_star;
  x_star  = tantheta*z0 + ( 1 + (px_star / ps_star) * sintheta ) * x0 ;
  y_star  = y0 + sintheta * (py_star / ps_star) * x0;
  z_star  = z0/costheta + (-h_star  / ps_star) * sintheta * x0;
  
  x[0] = x_star;  x[1] = px_star;  
  x[2] = y_star;  x[3] = py_star;
  x[4] = z_star;  x[5] = pz_star;
}

void ILBT_Pass(double x[], double theta)
{
  double sintheta=sin(theta), costheta=cos(theta), tantheta=sintheta/costheta;  
  double x0, px0, y0, py0, z0, pz0, ps0, h0;
  double x_star, px_star, y_star, py_star, z_star, pz_star;
  double hx_star, hy_star, hz_star, ps_star, h_star;

  x_star= x[0];   px_star=x[1];   
  y_star= x[2];   py_star=x[3];
  z_star= x[4];   pz_star=x[5];

  ps_star = sqrt( (1+pz_star )*(1+pz_star)- px_star * px_star - py_star * py_star );
  h_star  = 1+ pz_star -  ps_star;
  hx_star = px_star / ps_star ;
  hy_star = py_star / ps_star ;
  hz_star = -h_star /  ps_star ;
  
  px0 = (px_star + h_star * sintheta ) * costheta  ;  
  py0 = py_star * costheta;
  pz0 = pz_star  + px_star * sintheta;

  x0  = (x_star - z_star * sintheta ) / ( 1 +  hx_star*sintheta - hz_star *sintheta*sintheta );
  z0  = (z_star - hz_star * sintheta * x0  ) * costheta;
  y0  = y_star - sintheta * hy_star * x0;

  x[0]=x0 ;     x[1]=px0;
  x[2]=y0 ;     x[3]=py0;
  x[4]=z0 ;     x[5]=pz0;
}

template <class T> void DIFFUSE_Pass(T x[], double DIFF_X, double  DIFF_Y, double DIFF_DELTA)
{
  x[1]= x[1]+ DIFF_X * rnd(seed);
  x[3]= x[3]+ DIFF_Y * rnd(seed);
  x[5]= x[5]+ DIFF_DELTA *rnd(seed);
}

template <class T> void COOLING_Pass(T x[], double ALPHA)
{
  int i;
  for(i=0;i<4;i++) x[i]= (1.- ALPHA )* x[i];
}

//----coordinate system  change elements

template<class T> void TRANS_Pass( T x[], double DX, double DY, double DS)
//   translation in field-free region
{
  //if( abs(DS) > 1.0e-12) {
  //  DRIFT_Pass(x, DS);
  //}

  x[0]=x[0]-DX;
  x[2]=x[2]-DY;
  x[4]=x[4]-DS;
  
}

template<class T> void SROTAT_Pass( T x[], double PSI)
// frame s  rotation, folowing madx treatment, same as E. Forest's book
{
  if( PSI == 0. ){
    return;
  }
  else{
    int i;
    T xtemp[6];
    double cosT=cos(PSI), sinT=sin(PSI);
    
    for (i=0; i<6;i++) xtemp[i]=x[i];
    x[0] = cosT * xtemp[0] + sinT * xtemp[2];
    x[1] = cosT * xtemp[1] + sinT * xtemp[3];
    x[2] = cosT * xtemp[2] - sinT * xtemp[0];
    x[3] = cosT * xtemp[3] - sinT * xtemp[1];
  }
}


template<class T> void YROTAT_MADX_Pass( T x[], double THETA)
//  coordinate frame y rotation, y, py will not  change, following  madx treatment
{
  if( abs(THETA) < 1.0e-9 )  {
    return;
  }
  else{
    int i;
    T xtemp[6];
    T delta1= sqrt(1.0 + 2*x[pt_]/GP.beta+x[pt_]*x[pt_]) -1.0;  
    double cosT=cos(THETA), sinT=sin(THETA), tanT=tan(THETA);
    
    for (i=0; i<6;i++) xtemp[i]=x[i];
    x[0] = (1.0/cosT) * xtemp[0];
    x[1] = sinT  + cosT*xtemp[1] + (sinT/GP.beta) * xtemp[5];
    x[4] = x[4] + (-tanT/GP.beta) * xtemp[0];
  }
}

template<class T> void YROTAT_Pass( T x[], double THETA)
//  coordinate frame y rotation,  following  E. Forest's  book
{
  if(  abs(THETA) < 1.0e-9 )  {
    return;
  }
  else{
    int  i;
    T xtemp[6], temp, pz;
    T delta1= sqrt(1.0 + 2*x[pt_]/GP.beta+x[pt_]*x[pt_]) -1.0;
    T gamma1 =GP.gamma+sqrt(GP.gamma*GP.gamma-1.)*x[pt_];
    T beta1 = sqrt(1.0-1.0/gamma1/gamma1);
    double cosT=cos(THETA), sinT=sin(THETA), tanT=tan(THETA);
    
    for (i=0; i<6;i++) xtemp[i]=x[i];
    pz= sqrt( (1+delta1)*(1+delta1) - xtemp[1]*xtemp[1] -xtemp[3]*xtemp[3] );
    temp =  1- xtemp[1]*tanT/pz;
    
    x[0] =  xtemp[0]/cosT/temp;
    x[1] =  xtemp[1]*cosT + pz * sinT;
    x[2] =  xtemp[2] + xtemp[3]*xtemp[0]*tanT/pz/temp;
    x[3] =  xtemp[3];
    x[4] =  xtemp[4] - (1+delta1)*xtemp[0]*tanT/pz/temp/beta1;
    x[5] =  xtemp[5];
  }
  
}

template<class T> void XROTAT_Pass( T x[], double PHI  )
//  coordinate frame x rotation, following  E. Forest's  book
{
  if( abs(PHI) < 1.0e-9 )  {
    return;
  }
  else{
    int  i;
    T xtemp[6], temp, pz;
    T delta1= sqrt(1.0 + 2*x[pt_]/GP.beta+x[pt_]*x[pt_]) -1.0;
    T gamma1 =GP.gamma+sqrt(GP.gamma*GP.gamma-1.)*x[pt_];
    T beta1 = sqrt(1.0-1.0/gamma1/gamma1);
    double cosT=cos(PHI), sinT=sin(PHI), tanT=tan(PHI);
    
    for (i=0; i<6;i++) xtemp[i]=x[i];
    pz= sqrt( (1+delta1)*(1+delta1) - xtemp[1]*xtemp[1] -xtemp[3]*xtemp[3] );
    temp =  1- xtemp[3]*tanT/pz;
    
    x[2] =  xtemp[2]/cosT/temp;
    x[3] =  xtemp[3]*cosT + pz * sinT;
    x[0] =  xtemp[0] + xtemp[1]*xtemp[2]*tanT/pz/temp;
    x[1] =  xtemp[1];
    x[4] =  xtemp[4] - (1+delta1)*xtemp[2]*tanT/pz/temp/beta1;
    x[5] =  xtemp[5];
  }
  
}

template<class T> void PATCH_Pass( T x[], double DX, double DY, double DS, double THETA)
{
  
  TRANS_Pass(x, DX, DY, DS);
  YROTAT_Pass(x, THETA);
  //TRANS_Pass(x, DX, DY, DS);
  // trans first, yrot second, cuase very small x diference between Bmad and SimTRACK, but with larger z difference
  // yrot first, trans second, cause large x difference etween Bmad and SimTRACK, but with very small z difference
}

void cerrf( double xx, double yy, double & wx, double & wy )
{
  int n,nc,nu ;
  double x,y,q,h,xl,xh,yh,tx,ty,tn,sx,sy,saux;
  double rx[33],ry[33];
  double cc=1.12837916709551, zero=0, one=1, two=2, half=.5, 
         xlim=5.33, ylim=4.29, fac1=3.2, fac2=23, fac3=21;

  x = abs(xx);
  y = abs(yy);
  if (y < ylim  and  x  < xlim) 
    {
      q  = (one - y / ylim) * sqrt(one - (x/xlim)*(x/xlim) );
      h  = one / (fac1 * q);
      nc = 7 + int(fac2*q);
      xl = pow(h,(1 - nc));
      xh = y + half/h;
      yh = x;
      nu = 10 + int(fac3*q);
      rx[nu+1] = zero;
      ry[nu+1] = zero;
      for ( n=nu; n>=1; n--) {
	tx = xh + n * rx[n+1];
	ty = yh - n * ry[n+1];
	tn = tx*tx + ty*ty;
	rx[n] = half * tx / tn;
	ry[n] = half * ty / tn;
      }
      sx = zero;
      sy = zero;
      for( n=nc;n>=1; n--) {
	saux = sx + xl;
	sx = rx[n] * saux - ry[n] * sy;
	sy = rx[n] * sy + ry[n] * saux;
	xl = h * xl;
      }
      wx = cc * sx;
      wy = cc * sy;
    }
  else 
    {
      xh = y;
      yh = x;
      rx[1] = zero;
      ry[1] = zero;
      for ( n = 9; n>= 1; n--) {
	tx = xh + n * rx[1];
	ty = yh - n * ry[1];
	tn = tx*tx + ty*ty;
	rx[1] = half * tx / tn;
	ry[1] = half * ty / tn;
      }
      wx = cc * rx[1];
      wy = cc * ry[1];
    }
  if(yy < zero) 
    {
      wx =   two * exp(y*y-x*x) * cos(two*x*y) - wx;
      wy = - two * exp(y*y-x*x) * sin(two*x*y) - wy;
      if(xx >  zero) 	wy = -wy;
    }
  else
    {
      if(xx <  zero) wy = -wy;
    }	
}

void BB4D(double x, double y, double gamma, double N, double sigmax, double sigmay, double & Dpx, double & Dpy)
//apply to p-p collision
{
  double rp=1.53469826769e-18;
  double x1, y1, x2, y2, signx, signy;

  signx = ((x>=0.0)?(1):(-1));
  signy = ((y>=0.0)?(1):(-1));

  if (x==0. && y==0. )
    {
      Dpx=0.;
      Dpy=0.;
    }
  else
    {
      if ( abs(sigmax -sigmay)/sigmax < 1.0e-6 )
	{  
	  double r2= x*x+ y*y;
	  double temp1= 2*N*rp/gamma/r2;
	  double temp2= 1-exp(-r2/2/sigmax/sigmax);
	  Dpx=temp1 * x * temp2;
	  Dpy=temp1 * y * temp2;
	}
      else if ( sigmax > sigmay) 
	{
	  double temp1=(2.*N*rp/gamma)*sqrt(PI/2./(sigmax*sigmax-sigmay*sigmay));
	  double temp2=exp(-x*x/2./sigmax/sigmax - y*y/2./sigmay/sigmay);
	  double temp3=sqrt(2.*(sigmax*sigmax-sigmay*sigmay));
	  double w1_real, w1_imag;
	  double w2_real, w2_imag;
	  x1 = fabs(x)/temp3;
	  y1 = fabs(y)/temp3;
	  x2 = fabs(x)*sigmay/sigmax/temp3;
	  y2 = fabs(y)*sigmax/sigmay/temp3;
	  cerrf( x1, y1, w1_real, w1_imag);
	  cerrf( x2, y2, w2_real, w2_imag);
	  Dpy=temp1 *( w1_real-temp2* w2_real)*signy;
	  Dpx=temp1 *( w1_imag-temp2* w2_imag)*signx;  
	}
      else if ( sigmax < sigmay) 
	{
	  double temp1=(2*N*rp/gamma)*sqrt(PI/2/(sigmay*sigmay-sigmax*sigmax));      
	  double temp2=exp(-x*x/2/sigmax/sigmax - y*y/2/sigmay/sigmay);
	  double temp3=sqrt(2*(sigmay*sigmay-sigmax*sigmax));
	  double w1_real, w1_imag;
	  double w2_real, w2_imag;
	  x1 = fabs(x)/temp3;
	  y1 = fabs(y)/temp3;
	  x2 = fabs(x)*sigmay/sigmax/temp3;
	  y2 = fabs(y)*sigmax/sigmay/temp3;
	  cerrf( y1, x1, w1_real, w1_imag);
	  cerrf( y2, x2, w2_real, w2_imag);
	  Dpy=temp1 *( w1_imag-temp2* w2_imag)*signy; 
	  Dpx=temp1 *( w1_real-temp2* w2_real)*signx;
	}
    }
}

void BB4D(double xi[6], double gamma, double N, double bbscale, double sigmax, double sigmay)
//----apply to p-p collision or for othe collision with bbscale correction
{
  double rp=1.53469826769e-18;
  double x=xi[0], y=xi[2];
  double x1, y1, x2, y2, signx, signy;
  double Dpx, Dpy;

  signx = ((x>=0.0)?(1):(-1));
  signy = ((y>=0.0)?(1):(-1));

  if (x==0. && y==0. )
    {
      Dpx=0.;
      Dpy=0.;
    }
  else
    {
      if ( abs(sigmax -sigmay)/sigmax < 1.0e-6 )
	{  
	  double r2= x*x+ y*y;
	  double temp1= 2*N*rp/gamma/r2;
	  double temp2= 1-exp(-r2/2/sigmax/sigmax);
	  Dpx=temp1 * x * temp2;
	  Dpy=temp1 * y * temp2;
	}
      else if ( sigmax > sigmay) 
	{
	  double temp1=(2.*N*rp/gamma)*sqrt(PI/2./(sigmax*sigmax-sigmay*sigmay));
	  double temp2=exp(-x*x/2./sigmax/sigmax - y*y/2./sigmay/sigmay);
	  double temp3=sqrt(2.*(sigmax*sigmax-sigmay*sigmay));
	  double w1_real, w1_imag;
	  double w2_real, w2_imag;
	  x1 = fabs(x)/temp3;
	  y1 = fabs(y)/temp3;
	  x2 = fabs(x)*sigmay/sigmax/temp3;
	  y2 = fabs(y)*sigmax/sigmay/temp3;
	  cerrf( x1, y1, w1_real, w1_imag);
	  cerrf( x2, y2, w2_real, w2_imag);
	  Dpy=temp1 *( w1_real-temp2* w2_real)*signy;
	  Dpx=temp1 *( w1_imag-temp2* w2_imag)*signx;  
	}
      else if ( sigmax < sigmay) 
	{
	  double temp1=(2*N*rp/gamma)*sqrt(PI/2/(sigmay*sigmay-sigmax*sigmax));      
	  double temp2=exp(-x*x/2/sigmax/sigmax - y*y/2/sigmay/sigmay);
	  double temp3=sqrt(2*(sigmay*sigmay-sigmax*sigmax));
	  double w1_real, w1_imag;
	  double w2_real, w2_imag;
	  x1 = fabs(x)/temp3;
	  y1 = fabs(y)/temp3;
	  x2 = fabs(x)*sigmay/sigmax/temp3;
	  y2 = fabs(y)*sigmax/sigmay/temp3;
	  cerrf( y1, x1, w1_real, w1_imag);
	  cerrf( y2, x2, w2_real, w2_imag);
	  Dpy=temp1 *( w1_imag-temp2* w2_imag)*signy; 
	  Dpx=temp1 *( w1_real-temp2* w2_real)*signx;
	}
    }

  xi[1]=xi[1] + Dpx * bbscale;
  xi[3]=xi[3] + Dpy * bbscale;
   
}

void BB4D(double x, double y, double Q1, double E1, double N2, double Q2, double sigmax, double sigmay, double & Dpx, double & Dpy)
//   general BB kick calculation from  a bi-Gaussian charge distribution, return Dpx, Dpy
//   Q1, Q2:  unit in proton charge, with signs.
//   E1:  total energy of test particle  is E1, in unit of MeV
//   N2, Q2,  sigmax2, sigmay2: RIGID BUNHC'S PARAMETERS
{
  double e0=1.6021765e-19;
  double epsilon0 = 8.854187817e-12;
  double x1, y1, x2, y2, signx, signy;

  signx = ((x>=0.0)?(1):(-1));
  signy = ((y>=0.0)?(1):(-1));

  if (x==0. && y==0. )
    {
      Dpx=0.;
      Dpy=0.;
    }
  else
    {
      if ( abs(sigmax -sigmay)/sigmax < 1.0e-6 )
	{  
	  double r2= x*x+ y*y;
	  double temp1= N2*Q2*Q1*e0/(E1*1e6)/(2.0*PI*epsilon0)/r2;
	  double temp2= 1-exp(-r2/2/sigmax/sigmax);
	  Dpx=temp1 * x * temp2;
	  Dpy=temp1 * y * temp2;
	}
      else if ( sigmax > sigmay) 
	{
	  double temp1=N2*Q2*Q1*e0/(E1*1e6)/(2.0*PI*epsilon0)*sqrt(PI/2./(sigmax*sigmax-sigmay*sigmay));
	  double temp2=exp(-x*x/2./sigmax/sigmax - y*y/2./sigmay/sigmay);
	  double temp3=sqrt(2.*(sigmax*sigmax-sigmay*sigmay));
	  double w1_real, w1_imag;
	  double w2_real, w2_imag;
	  x1 = fabs(x)/temp3;
	  y1 = fabs(y)/temp3;
	  x2 = fabs(x)*sigmay/sigmax/temp3;
	  y2 = fabs(y)*sigmax/sigmay/temp3;
	  cerrf( x1, y1, w1_real, w1_imag);
	  cerrf( x2, y2, w2_real, w2_imag);
	  Dpy=temp1 *( w1_real-temp2* w2_real)*signy;
	  Dpx=temp1 *( w1_imag-temp2* w2_imag)*signx;  
	}
      else if ( sigmax < sigmay) 
	{
	  double temp1=N2*Q2*Q1*e0/(E1*1e6)/(2.0*PI*epsilon0)*sqrt(PI/2/(sigmay*sigmay-sigmax*sigmax));      
	  double temp2=exp(-x*x/2/sigmax/sigmax - y*y/2/sigmay/sigmay);
	  double temp3=sqrt(2*(sigmay*sigmay-sigmax*sigmax));
	  double w1_real, w1_imag;
	  double w2_real, w2_imag;
	  x1 = fabs(x)/temp3;
	  y1 = fabs(y)/temp3;
	  x2 = fabs(x)*sigmay/sigmax/temp3;
	  y2 = fabs(y)*sigmax/sigmay/temp3;
	  cerrf( y1, x1, w1_real, w1_imag);
	  cerrf( y2, x2, w2_real, w2_imag);
	  Dpy=temp1 *( w1_imag-temp2* w2_imag)*signy; 
	  Dpx=temp1 *( w1_real-temp2* w2_real)*signx;
	}
    }
}

void BB4D(double xi[6], double Q1, double E1, double N2, double Q2, double sigmax, double sigmay)
//   general BB kick calculation from a bi-Gaussian charge distribution  
//   Q1, Q2: unit is proton charge, with signs.
//   E1:  total energy of test particle  is A1 * E1, in unit of MeV
//   N2, Q2,  sigmax2, sigmay2:  another bunch paramter
{
  double e0=1.6021765e-19;
  double epsilon0 = 8.854187817e-12;
  double x=xi[0], y=xi[2];
  double x1, y1, x2, y2, signx, signy;
  double Dpx, Dpy;

  signx = ((x>=0.0)?(1):(-1));
  signy = ((y>=0.0)?(1):(-1));

  if (x==0. && y==0. )
    {
      Dpx=0.;
      Dpy=0.;
    }
  else
    {
      if ( abs(sigmax -sigmay)/sigmax < 1.0e-6 )
	{  
	  double r2= x*x+ y*y;
	  double temp1= N2*Q2*Q1*e0/(E1*1e6)/(2.0*PI*epsilon0)/r2;
	  double temp2= 1-exp(-r2/2/sigmax/sigmax);
	  Dpx=temp1 * x * temp2;
	  Dpy=temp1 * y * temp2;
	}
      else if ( sigmax > sigmay) 
	{
	  double temp1=N2*Q2*Q1*e0/(E1*1e6)/(2.0*PI*epsilon0)*sqrt(PI/2./(sigmax*sigmax-sigmay*sigmay));
	  double temp2=exp(-x*x/2./sigmax/sigmax - y*y/2./sigmay/sigmay);
	  double temp3=sqrt(2.*(sigmax*sigmax-sigmay*sigmay));
	  double w1_real, w1_imag;
	  double w2_real, w2_imag;
	  x1 = fabs(x)/temp3;
	  y1 = fabs(y)/temp3;
	  x2 = fabs(x)*sigmay/sigmax/temp3;
	  y2 = fabs(y)*sigmax/sigmay/temp3;
	  cerrf( x1, y1, w1_real, w1_imag);
	  cerrf( x2, y2, w2_real, w2_imag);
	  Dpy=temp1 *( w1_real-temp2* w2_real)*signy;
	  Dpx=temp1 *( w1_imag-temp2* w2_imag)*signx;  
	}
      else if ( sigmax < sigmay) 
	{
	  double temp1=N2*Q2*Q1*e0/(E1*1e6)/(2.0*PI*epsilon0)*sqrt(PI/2/(sigmay*sigmay-sigmax*sigmax));      
	  double temp2=exp(-x*x/2/sigmax/sigmax - y*y/2/sigmay/sigmay);
	  double temp3=sqrt(2*(sigmay*sigmay-sigmax*sigmax));
	  double w1_real, w1_imag;
	  double w2_real, w2_imag;
	  x1 = fabs(x)/temp3;
	  y1 = fabs(y)/temp3;
	  x2 = fabs(x)*sigmay/sigmax/temp3;
	  y2 = fabs(y)*sigmax/sigmay/temp3;
	  cerrf( y1, x1, w1_real, w1_imag);
	  cerrf( y2, x2, w2_real, w2_imag);
	  Dpy=temp1 *( w1_imag-temp2* w2_imag)*signy; 
	  Dpx=temp1 *( w1_real-temp2* w2_real)*signx;
	}
    }
  xi[1]=xi[1] + Dpx;
  xi[3]=xi[3] + Dpy;
}

void BB6D(double x[], double gamma, double Np, double bbscale, double sigma_l, int N_slice, 
                       double emitx_rms,  double betax_star, double alfx_star,
	               double emity_rms,  double betay_star, double alfy_star)
// apply to p-p collision, other type of collision can be adjusted with bbscale.
// emitx_rms, emity_rms:  un-normalized rms emittance,  sigma=SQRT[ emitx_rms * betax ]  
{
  int i,j; 
  double rp= 1.534698e-18;
  double x0[6];
  double Np_slice[N_slice], z_star[N_slice];
  double S, gx_star, gy_star, betax, alfx, betay, alfy;
  double sigmax, sigmay, dsigmax2ds, dsigmay2ds, dUdsigmax2, dUdsigmay2; 
  double X, Y,  Dpx, Dpy, temp1, temp2, temp3, temp4;
  
  //----center of each slice of strong bunch, each slice has not same particle population.

  for(i=0;i<N_slice;i++)  Np_slice[i] =  1.0* Np /  N_slice;

  if(N_slice == 11 ) {
    z_star[0]=  1.7996765362153655;
    z_star[1]=  1.1049618149811977;
    z_star[2]=  0.75072324235797372;
    z_star[3]=  0.47407703558819131;
    z_star[4]=  0.23041113693387344;
    z_star[5]=  0.;
    z_star[6]= -0.23041113693387344;
    z_star[7]= -0.47407703558819131; 
    z_star[8]= -0.75072324235797372;
    z_star[9]= -1.1049618149811977;
    z_star[10]=-1.7996765362153655; 
    for(i=0;i<11;i++) z_star[i]=z_star[i] * sigma_l;
  }
  else{
    double y[N_slice-1];
    for(i=0;i<N_slice-1;i++){
      y[i]= gsl_cdf_ugaussian_Pinv(  (i+1)*1.0/N_slice   );
    }
    z_star[0]= N_slice *( exp(-y[0]*y[0]/2 ) ) / sqrt(2.0*3.14159265);
    z_star[N_slice-1]= N_slice *( 0-exp(-y[0]*y[0]/2 ) )/ sqrt(2.0*3.14159265);
    for(i=1;i<N_slice-1;i++)  z_star[i]= N_slice *( exp(-y[i]*y[i]/2) - exp(-y[i-1]*y[i-1]/2 )  ) / sqrt(2.0*3.14159265);
    for(i=0;i<N_slice;i++) z_star[i]=  z_star[i] * sigma_l;
  }

  //-----calculate the changes in x[6]
  for(i=0; i<N_slice;i++) {

    for (j=0;j<6;j++) x0[j]=x[j];
    S=(x0[4]-z_star[i])/2.0;

    gx_star= (1.0 + alfx_star * alfx_star )/ betax_star;
    gy_star= (1.0 + alfy_star * alfy_star )/ betay_star;
    betax= betax_star + 2* S * alfx_star + S * S *  gx_star;  //  strong slide drift to  -S
    alfx =(alfx_star  + gx_star * S );
    betay= betay_star + 2* S * alfy_star + S * S *  gy_star; 
    alfy =(alfy_star  + gy_star * S );

    sigmax  =  sqrt( emitx_rms * betax ) ;
    sigmay  =  sqrt( emity_rms * betay ) ;
    dsigmax2ds = 2.0* emitx_rms * (alfx_star +  S *  gx_star) ;  //---use S  as  variable 
    dsigmay2ds = 2.0* emity_rms * (alfy_star +  S *  gy_star);   // same as  test particle's direction

    //----calculate the kicks for each slice
    X=x0[0] + x0[1]*S;
    Y=x0[2] + x0[3]*S;
    BB4D(X, Y, gamma, Np_slice[i], sigmax, sigmay, Dpx, Dpy);
    temp1 = 1.0/2./(sigmax*sigmax-sigmay*sigmay) ;
    temp2 = X *  Dpx + Y *  Dpy;
    temp3 = 2.* Np_slice[i] *rp /gamma ;
    temp4 = exp ( - X * X / 2. / sigmax /sigmax -  Y * Y / 2. / sigmay /sigmay ) ;
    
    x[0]=x0[0] - S * Dpx * bbscale ;
    x[1]=x0[1] + Dpx  * bbscale;   
    x[2]=x0[2] - S * Dpy * bbscale ;
    x[3]=x0[3] + Dpy * bbscale;   
    x[4]=x0[4];
    
    if ( abs(sigmax- sigmay)/sigmax < 1.0e-6 ) 
      {
      	x[5]=x0[5]+ 0.5 * Dpx * bbscale * ( x0[1] + 0.5* Dpx * bbscale ) 
                  + 0.5*  Dpy * bbscale * ( x0[3] + 0.5* Dpy * bbscale ) 
                  + (1.0/sigmax/sigmax ) * dsigmax2ds  * ( temp3 * bbscale/ 4. ) * temp4;
      }
    else if (sigmax > sigmay ) 
      {
	dUdsigmax2 =  temp1 * bbscale * ( temp2  + temp3 *( sigmay/ sigmax * temp4 -1. ) );
	dUdsigmay2 = -temp1 * bbscale * ( temp2  + temp3 *( sigmax/ sigmay * temp4 -1. ) );
	x[5]= x0[5] + 0.5 * Dpx * bbscale * ( x0[1] + 0.5* Dpx * bbscale) 
                    + 0.5 * Dpy * bbscale * ( x0[3] + 0.5* Dpy * bbscale) 
                    - ( 0.5 * dsigmax2ds * dUdsigmax2 + 0.5 * dsigmay2ds *  dUdsigmay2 ) ;
      }
    else
      {
	dUdsigmay2 = -temp1 *  bbscale * ( temp2  + temp3 *( sigmax/ sigmay * temp4 -1. ) );
	dUdsigmax2 =  temp1 *  bbscale * ( temp2  + temp3 *( sigmay/ sigmax * temp4 -1. ) );
	x[5]= x0[5] + 0.5 * Dpx * bbscale *  ( x0[1] + 0.5* Dpx* bbscale) 
                    + 0.5 * Dpy * bbscale *  ( x0[3] + 0.5* Dpy *bbscale) 
                    - ( 0.5 * dsigmax2ds * dUdsigmax2 + 0.5 * dsigmay2ds *  dUdsigmay2 ) ;
      }
  }
}

void BB6D_Lumi(double x[], double gamma, double Np, double bbscale, double sigma_l, int N_slice, 
               double emitx_rms,  double betax_star, double alfx_star,
	       double emity_rms,  double betay_star, double alfy_star,
               double & lumi_part)
 // apply to p-p collision, other type of collision can be adjusted by bbscale.
 // emitx_rms, emity_rms:  un-normalized rms emittance,  sigma=SQRT[ emitx_rms * betax ]  
{
  int i,j; 
  double rp= 1.534698e-18;
  double x0[6];
  double Np_slice[N_slice], z_star[N_slice];
  double S, gx_star, gy_star, betax, alfx, betay, alfy;
  double sigmax, sigmay, dsigmax2ds, dsigmay2ds, dUdsigmax2, dUdsigmay2; 
  double X, Y,  Dpx, Dpy, temp1, temp2, temp3, temp4;
  
  //----center of each slice of strong bunch, each slice has not same particle population.

  for(i=0;i<N_slice;i++)  Np_slice[i] =  1.0* Np /  N_slice;

  if(N_slice == 11 ) {
    z_star[0]=  1.7996765362153655;
    z_star[1]=  1.1049618149811977;
    z_star[2]=  0.75072324235797372;
    z_star[3]=  0.47407703558819131;
    z_star[4]=  0.23041113693387344;
    z_star[5]=  0.;
    z_star[6]= -0.23041113693387344;
    z_star[7]= -0.47407703558819131; 
    z_star[8]= -0.75072324235797372;
    z_star[9]= -1.1049618149811977;
    z_star[10]=-1.7996765362153655; 
    for(i=0;i<11;i++) z_star[i]=z_star[i] * sigma_l;
  }
  else{
    double y[N_slice-1];
    for(i=0;i<N_slice-1;i++){
      y[i]= gsl_cdf_ugaussian_Pinv(  (i+1)*1.0/N_slice   );
    }
    z_star[0]= N_slice *( exp(-y[0]*y[0]/2 ) ) / sqrt(2.0*3.14159265);
    z_star[N_slice-1]= N_slice *( 0-exp(-y[0]*y[0]/2 ) )/ sqrt(2.0*3.14159265);
    for(i=1;i<N_slice-1;i++)  z_star[i]= N_slice *( exp(-y[i]*y[i]/2) - exp(-y[i-1]*y[i-1]/2 )  ) / sqrt(2.0*3.14159265);
    for(i=0;i<N_slice;i++) z_star[i]=  z_star[i] * sigma_l;
    //for(i=0;i<N_slice;i++) cout<<  z_star[i] /  sigma_l<<endl; 
  }

  //-----calculate the changes in x[6]
  lumi_part  = 0 ;

  for(i=0; i<N_slice;i++) {

    for (j=0;j<6;j++) x0[j]=x[j];
    S=(x0[4]-z_star[i])/2.0;

    gx_star= (1.0 + alfx_star * alfx_star )/ betax_star;
    gy_star= (1.0 + alfy_star * alfy_star )/ betay_star;
    betax= betax_star - 2* S * alfx_star + S * S *  gx_star; 
    alfx =(alfx_star -gx_star * S );
    betay=  betay_star - 2* S * alfy_star + S * S *  gy_star; 
    alfy =(alfy_star -gy_star * S );

    sigmax  =  sqrt( emitx_rms * betax ) ;
    sigmay  =  sqrt( emity_rms * betay ) ;
    dsigmax2ds =-2.0*alfx* ( emitx_rms );
    dsigmay2ds =-2.0*alfy* ( emity_rms );

    //----calculate the kicks for each slice
    X=x0[0] + x0[1]*S;
    Y=x0[2] + x0[3]*S;
    BB4D(X, Y, gamma, Np_slice[i], sigmax, sigmay, Dpx, Dpy);
    temp1 = 1.0/2./(sigmax*sigmax-sigmay*sigmay) ;
    temp2 = X *  Dpx + Y *  Dpy;
    temp3 = 2.* Np_slice[i] *rp /gamma ;
    temp4 = exp ( - X * X / 2. / sigmax /sigmax -  Y * Y / 2. / sigmay /sigmay ) ;
    
    x[0]=x0[0] - S * Dpx * bbscale ;
    x[1]=x0[1] + Dpx  * bbscale;   
    x[2]=x0[2] - S * Dpy * bbscale ;
    x[3]=x0[3] + Dpy * bbscale;   
    x[4]=x0[4];
    
    if ( abs(sigmax- sigmay)/sigmax < 1.0e-6 ) 
      {
      	x[5]=x0[5]+ 0.5 * Dpx * bbscale * ( x0[1] + 0.5* Dpx * bbscale ) 
                  + 0.5*  Dpy * bbscale * ( x0[3] + 0.5* Dpy * bbscale ) 
                  + (1.0/sigmax/sigmax ) * dsigmax2ds  * ( temp3 * bbscale/ 4. ) * temp4;
      }
    else if (sigmax > sigmay ) 
      {
	dUdsigmax2 =  temp1 * bbscale * ( temp2  + temp3 *( sigmay/ sigmax * temp4 -1. ) );
	dUdsigmay2 = -temp1 * bbscale * ( temp2  + temp3 *( sigmax/ sigmay * temp4 -1. ) );
	x[5]= x0[5] + 0.5 * Dpx * bbscale * ( x0[1] + 0.5* Dpx * bbscale) 
                    + 0.5 * Dpy * bbscale * ( x0[3] + 0.5* Dpy * bbscale) 
                    - ( 0.5 * dsigmax2ds * dUdsigmax2 + 0.5 * dsigmay2ds *  dUdsigmay2 ) ;
      }
    else
      {
	dUdsigmay2 = -temp1 *  bbscale * ( temp2  + temp3 *( sigmax/ sigmay * temp4 -1. ) );
	dUdsigmax2 =  temp1 *  bbscale * ( temp2  + temp3 *( sigmay/ sigmax * temp4 -1. ) );
	x[5]= x0[5] + 0.5 * Dpx * bbscale *  ( x0[1] + 0.5* Dpx* bbscale) 
                    + 0.5 * Dpy * bbscale *  ( x0[3] + 0.5* Dpy *bbscale) 
                    - ( 0.5 * dsigmax2ds * dUdsigmax2 + 0.5 * dsigmay2ds *  dUdsigmay2 ) ;
      }

    //-----luminosity contribution  
    lumi_part  =  lumi_part + 1.0 * Np_slice[i] * exp( - X * X /2 /sigmax /sigmax -  Y * Y /2 /sigmay /sigmay  ) / 2 / PI /sigmax / sigmay ; 

  }
}


void  Lorentz_Transfer(double x[6], double theta)
{
  double sintheta=sin(theta), costheta=cos(theta), tantheta=sintheta/costheta;  
  double x0, px0, y0, py0, z0, pz0, ps0, h0;
  double x_star, px_star, y_star, py_star, z_star, pz_star, ps_star, h_star;

  x0= x[0];   px0=x[1];   y0= x[2];
  py0=x[3];   z0= x[4];   pz0=x[5];
  ps0= sqrt( (1+pz0 )*(1+pz0)- px0*px0 - py0*py0 );  
  h0 = (1+pz0)- ps0;

  px_star = ( px0 - tantheta * h0 )  / cos(theta );
  py_star = py0 / costheta;
  pz_star = pz0 - tantheta * px0 + tantheta * tantheta * h0;

  ps_star = sqrt( (1+pz_star )*(1+pz_star)- px_star * px_star - py_star * py_star ); 
  h_star  = (1+pz_star) - ps_star;
  x_star  = tantheta*z0 + ( 1 + (px_star / ps_star) * sintheta ) * x0 ;
  y_star  = y0 + sintheta * (py_star / ps_star) * x0;
  z_star  = z0/costheta + (-h_star  / ps_star) * sintheta * x0;
  
  x[0] = x_star;  x[1] = px_star;  
  x[2] = y_star;  x[3] = py_star;
  x[4] = z_star;  x[5] = pz_star;
}


void  Lorentz_Transfer_Inverse(double x[6], double theta)
{
  double sintheta=sin(theta), costheta=cos(theta), tantheta=sintheta/costheta;  
  double x0, px0, y0, py0, z0, pz0, ps0, h0;
  double x_star, px_star, y_star, py_star, z_star, pz_star;
  double hx_star, hy_star, hz_star, ps_star, h_star;

  x_star= x[0];   px_star=x[1];   
  y_star= x[2];   py_star=x[3];
  z_star= x[4];   pz_star=x[5];

  ps_star = sqrt( (1+pz_star )*(1+pz_star)- px_star * px_star - py_star * py_star );
  h_star  = 1+ pz_star -  ps_star;
  hx_star = px_star / ps_star ;
  hy_star = py_star / ps_star ;
  hz_star = -h_star /  ps_star ;
  
  px0 = (px_star + h_star * sintheta ) * costheta  ;  
  py0 = py_star * costheta;
  pz0 = pz_star  + px_star * sintheta;

  x0  = (x_star - z_star * sintheta ) / ( 1 +  hx_star*sintheta - hz_star *sintheta*sintheta );
  z0  = (z_star - hz_star * sintheta * x0  ) * costheta;
  y0  = y_star - sintheta * hy_star * x0;

  x[0]=x0 ;     x[1]=px0;
  x[2]=y0 ;     x[3]=py0;
  x[4]=z0 ;     x[5]=pz0;
}

void  Lorentz_Transfer(int  Np, double x[], double theta)
{
  int i;
  
  double sintheta=sin(theta), costheta=cos(theta), tantheta=sintheta/costheta;  
  double x0, px0, y0, py0, z0, pz0, ps0, h0;
  double x_star, px_star, y_star, py_star, z_star, pz_star, ps_star, h_star;

  for(i=0;i<Np;i++){
    
    x0= x[i*6+0];   px0=x[i*6+1];   y0= x[i*6+2];
    py0=x[i*6+3];   z0= x[i*6+4];   pz0=x[i*6+5];
    ps0= sqrt( (1+pz0 )*(1+pz0)- px0*px0 - py0*py0 );  
    h0 = (1+pz0)- ps0;
    
    px_star = ( px0 - tantheta * h0 )  / cos(theta );
    py_star = py0 / costheta;
    pz_star = pz0 - tantheta * px0 + tantheta * tantheta * h0;
    
    ps_star = sqrt( (1+pz_star )*(1+pz_star)- px_star * px_star - py_star * py_star ); 
    h_star  = (1+pz_star) - ps_star;
    x_star  = tantheta*z0 + ( 1 + (px_star / ps_star) * sintheta ) * x0 ;
    y_star  = y0 + sintheta * (py_star / ps_star) * x0;
    z_star  = z0/costheta + (-h_star  / ps_star) * sintheta * x0;
    
    x[i*6+0] = x_star;  x[i*6+1] = px_star;  
    x[i*6+2] = y_star;  x[i*6+3] = py_star;
    x[i*6+4] = z_star;  x[i*6+5] = pz_star;
    
  }
  
}

void  Lorentz_Transfer_Inverse(int Np, double x[6], double theta)
{
  int i;
  
  double sintheta=sin(theta), costheta=cos(theta), tantheta=sintheta/costheta;  
  double x0, px0, y0, py0, z0, pz0, ps0, h0;
  double x_star, px_star, y_star, py_star, z_star, pz_star;
  double hx_star, hy_star, hz_star, ps_star, h_star;

  for(i=0;i<Np;i++){
    
    x_star= x[i*6+0];   px_star=x[i*6+1];   
    y_star= x[i*6+2];   py_star=x[i*6+3];
    z_star= x[i*6+4];   pz_star=x[i*6+5];
    
    ps_star = sqrt( (1+pz_star )*(1+pz_star)- px_star * px_star - py_star * py_star );
    h_star  = 1+ pz_star -  ps_star;
    hx_star = px_star / ps_star ;
    hy_star = py_star / ps_star ;
    hz_star = -h_star /  ps_star ;
    
    px0 = (px_star + h_star * sintheta ) * costheta  ;  
    py0 = py_star * costheta;
    pz0 = pz_star  + px_star * sintheta;
    
    x0  = (x_star - z_star * sintheta ) / ( 1 +  hx_star*sintheta - hz_star *sintheta*sintheta );
    z0  = (z_star - hz_star * sintheta * x0  ) * costheta;
    y0  = y_star - sintheta * hy_star * x0;
    
    x[i*6+0]=x0 ;     x[i*6+1]=px0;
    x[i*6+2]=y0 ;     x[i*6+3]=py0;
    x[i*6+4]=z0 ;     x[i*6+5]=pz0;
    
  }
  
}

void BB6D_Angle(double x[6], double gamma, double Np, double bbscale, double sigma_l, double theta, int N_slice, 
                       double emitx_rms,  double betax_star, double alfx_star,
	               double emity_rms,  double betay_star, double alfy_star)
 // assumed for p-p collision, other type of collision can be adjusted by bbscale.
 // emitx_rms, emity_rms:  un-normalized rms emittance,  sigma=SQRT[ emitx_rms * betax ]  
{
  int i; 
  double rp= 1.534698e-18;
  double sintheta=sin(theta), costheta=cos(theta), tantheta=sintheta/costheta;
  double Np_slice[N_slice], Z0[N_slice], Z0_star[N_slice];
  double S, gx_star, gy_star, betax, alfx, betay, alfy;
  double sigmax, sigmay, dsigmax2ds, dsigmay2ds, dUdsigmax2, dUdsigmay2; 
  double X, Y,  Dpx, Dpy, temp1, temp2, temp3, temp4;
  
  double  x0, y0, z0, px0, py0, pz0, ps0, h0;

  double  x_star, y_star, z_star, px_star, py_star, pz_star;
  double  hx_star, hy_star, hz_star, ps_star, h_star;

  double  x_star_0, y_star_0, z_star_0, px_star_0, py_star_0, pz_star_0;

  //----center of each slice of strong bunch, here each slice has not same particle population.
  for(i=0;i<N_slice;i++)  Np_slice[i] =  1.0* Np /  N_slice;

  if(N_slice == 11 ) {
    Z0[0]=  1.7996765362153655;
    Z0[1]=  1.1049618149811977;
    Z0[2]=  0.75072324235797372;
    Z0[3]=  0.47407703558819131;
    Z0[4]=  0.23041113693387344;
    Z0[5]=  0.;
    Z0[6]= -0.23041113693387344;
    Z0[7]= -0.47407703558819131; 
    Z0[8]= -0.75072324235797372;
    Z0[9]= -1.1049618149811977;
    Z0[10]=-1.7996765362153655; 
    for(i=0;i<11;i++) Z0[i]=Z0[i] * sigma_l;
  }
  else{
    double y[N_slice-1];
    for(i=0;i<N_slice-1;i++){
      y[i]= gsl_cdf_ugaussian_Pinv(  (i+1)*1.0/N_slice   );
    }
    Z0[0]= N_slice *( exp(-y[0]*y[0]/2 ) ) / sqrt(2.0*3.14159265);
    Z0[N_slice-1]= N_slice *( 0-exp(-y[0]*y[0]/2 ) )/ sqrt(2.0*3.14159265);
    for(i=1;i<N_slice-1;i++)  Z0[i]= N_slice *( exp(-y[i]*y[i]/2) - exp(-y[i-1]*y[i-1]/2 )  ) / sqrt(2.0*3.14159265);
    for(i=0;i<N_slice;i++) Z0[i]=  Z0[i] * sigma_l;
  }

  //----Perform Lorentz booster to the test particle's (x,px,y,py,z,delta)

  x0= x[0];   px0=x[1];   y0= x[2];
  py0=x[3];   z0= x[4];   pz0=x[5];
  ps0= sqrt( (1+pz0 )*(1+pz0)- px0*px0 - py0*py0 );  
  h0 = (1+pz0)- ps0;

  px_star = ( px0 - tantheta * h0 )  / cos(theta );
  py_star = py0 / costheta;
  pz_star = pz0 - tantheta * px0 + tantheta * tantheta * h0;

  ps_star = sqrt( (1+pz_star )*(1+pz_star)- px_star * px_star - py_star * py_star ); 
  h_star  = (1+pz_star) - ps_star;
  x_star  = tantheta*z0 + ( 1 + (px_star / ps_star) * sintheta ) * x0 ;
  y_star  = y0 + sintheta * (py_star / ps_star) * x0;
  z_star  = z0/costheta + (-h_star  / ps_star) * sintheta * x0;
  
  //----Perform Lorentz booster to the center of  sliced strong beam
  for(i=0;i<N_slice;i++)
    Z0_star[i]=Z0[i]/costheta; 

  //----Twiss parmeters and emittances in the new frame
  
  gx_star= (1.0 + alfx_star * alfx_star )/ betax_star;
  gy_star= (1.0 + alfy_star * alfy_star )/ betay_star;

  gx_star=gx_star/costheta;
  gy_star=gy_star/costheta;
  betax_star= betax_star *costheta;
  betay_star= betay_star *costheta;
  alfx_star= alfx_star;
  alfy_star= alfy_star;
  emitx_rms = emitx_rms / costheta;
  emity_rms = emity_rms / costheta;

  //-----head-on collision in head-on frame
  for(i=0; i<N_slice;i++) {

    x_star_0 =x_star; 
    px_star_0=px_star;
    y_star_0 =y_star; 
    py_star_0=py_star;
    z_star_0 =z_star; 
    pz_star_0=pz_star;

    S=(z_star_0 - Z0_star[i])/2.0;

    betax= betax_star - 2* S * alfx_star + S * S *  gx_star; 
    alfx =(alfx_star -gx_star * S );
    betay= betay_star - 2* S * alfy_star + S * S *  gy_star; 
    alfy =(alfy_star -gy_star * S );

    sigmax  =  sqrt( emitx_rms * betax ) ;
    sigmay  =  sqrt( emity_rms * betay ) ;
    dsigmax2ds =-2.0*alfx* ( emitx_rms );
    dsigmay2ds =-2.0*alfy* ( emity_rms );
    
    X= x_star_0 + px_star_0 * S - Z0_star[i] * sintheta;
    Y= y_star_0 + py_star_0 * S - 0 ;
    BB4D(X, Y, gamma, Np_slice[i], sigmax, sigmay, Dpx, Dpy);
    temp1 = 1.0/2./(sigmax*sigmax-sigmay*sigmay) ;
    temp2 = X *  Dpx + Y *  Dpy;
    temp3 = 2.* Np_slice[i] *rp /gamma ;
    temp4 = exp ( - X * X / 2. / sigmax /sigmax -  Y * Y / 2. / sigmay /sigmay ) ;
    
    x_star =x_star_0  -  S * Dpx * bbscale ;
    px_star=px_star_0 +  Dpx * bbscale;   
    y_star =y_star_0  -  S * Dpy * bbscale ;
    py_star=py_star_0 +  Dpy * bbscale;   
    z_star =z_star_0;

    if ( abs(sigmax- sigmay)/sigmax < 1.0e-6 ) 
      {
      	pz_star=pz_star_0 + 0.5 * Dpx * bbscale * ( px_star_0 + 0.5* Dpx * bbscale ) 
                          + 0.5*  Dpy * bbscale * ( py_star_0 + 0.5* Dpy * bbscale ) 
                          + (1.0/sigmax/sigmax ) * dsigmax2ds  * ( temp3 * bbscale/ 4. ) * temp4;
      }
    else if (sigmax > sigmay ) 
      {
	dUdsigmax2 =  temp1 * bbscale * ( temp2  + temp3 *( sigmay/ sigmax * temp4 -1. ) );
	dUdsigmay2 = -temp1 * bbscale * ( temp2  + temp3 *( sigmax/ sigmay * temp4 -1. ) );
	pz_star= pz_star_0 + 0.5 * Dpx * bbscale * ( px_star_0 + 0.5* Dpx *  bbscale) 
                           + 0.5 * Dpy * bbscale * ( py_star_0 + 0.5* Dpy *  bbscale) 
                           - ( 0.5 * dsigmax2ds * dUdsigmax2 + 0.5 * dsigmay2ds *  dUdsigmay2 ) ;
      }
    else
      {
	dUdsigmay2 = -temp1 *  bbscale * ( temp2  + temp3 *( sigmax/ sigmay * temp4 -1. ) );
	dUdsigmax2 =  temp1 *  bbscale * ( temp2  + temp3 *( sigmay/ sigmax * temp4 -1. ) );
	pz_star    =  pz_star_0 + 0.5 * Dpx *bbscale * ( px_star_0 + 0.5* Dpx * bbscale) 
                                + 0.5 * Dpy *bbscale * ( py_star_0 + 0.5* Dpy * bbscale) 
                                - ( 0.5 * dsigmax2ds * dUdsigmax2 + 0.5 * dsigmay2ds *  dUdsigmay2 ) ;
      }
   }

  //----perform inverse Lorentz booster to the test particle's coordinates
  ps_star = sqrt( (1+pz_star )*(1+pz_star)- px_star * px_star - py_star * py_star );
  h_star  = 1+ pz_star -  ps_star;
  hx_star = px_star / ps_star ;
  hy_star = py_star / ps_star ;
  hz_star = -h_star /  ps_star ;
  
  px0 = (px_star + h_star * sintheta ) * costheta  ;  
  py0 = py_star * costheta;
  pz0 = pz_star  + px_star * sintheta;

  x0  = (x_star - z_star * sintheta ) / ( 1 +  hx_star*sintheta - hz_star *sintheta*sintheta );
  z0  = (z_star - hz_star * sintheta * x0  ) * costheta;
  y0  = y_star - sintheta * hy_star * x0;

  x[0]=x0 ;  x[1]=px0;
  x[2]=y0 ;  x[3]=py0;
  x[4]=z0 ;  x[5]=pz0;
}

void BB6D_Angle_Lumi(double x[6], double gamma, double Np, double bbscale, double sigma_l, double theta, int N_slice, 
                     double emitx_rms,  double betax_star, double alfx_star,
		     double emity_rms,  double betay_star, double alfy_star,
                     double & lumi_part)
 // assumed for p-p collision, other type of collision can be adjusted by bbscale.  
 // emitx_rms, emity_rms:  un-normalized rms emittance,  sigma=SQRT[ emitx_rms * betax ]  
{
  int i; 
  double rp= 1.534698e-18;
  double sintheta=sin(theta), costheta=cos(theta), tantheta=sintheta/costheta;
  double Np_slice[N_slice], Z0[N_slice], Z0_star[N_slice];
  double S, gx_star, gy_star, betax, alfx, betay, alfy;
  double sigmax, sigmay, dsigmax2ds, dsigmay2ds, dUdsigmax2, dUdsigmay2; 
  double X, Y,  Dpx, Dpy, temp1, temp2, temp3, temp4;
  
  double  x0, y0, z0, px0, py0, pz0, ps0, h0;

  double  x_star, y_star, z_star, px_star, py_star, pz_star;
  double  hx_star, hy_star, hz_star, ps_star, h_star;

  double  x_star_0, y_star_0, z_star_0, px_star_0, py_star_0, pz_star_0;

  //----center of each slice of strong bunch, here each slice has not same particle population.
  for(i=0;i<N_slice;i++)  Np_slice[i] =  1.0* Np /  N_slice;

  if(N_slice == 11 ) {
    Z0[0]=  1.7996765362153655;
    Z0[1]=  1.1049618149811977;
    Z0[2]=  0.75072324235797372;
    Z0[3]=  0.47407703558819131;
    Z0[4]=  0.23041113693387344;
    Z0[5]=  0.;
    Z0[6]= -0.23041113693387344;
    Z0[7]= -0.47407703558819131; 
    Z0[8]= -0.75072324235797372;
    Z0[9]= -1.1049618149811977;
    Z0[10]=-1.7996765362153655; 
    for(i=0;i<11;i++) Z0[i]=Z0[i] * sigma_l;
  }
  else{
    double y[N_slice-1];
    for(i=0;i<N_slice-1;i++){
      y[i]= gsl_cdf_ugaussian_Pinv(  (i+1)*1.0/N_slice   );
    }
    Z0[0]= N_slice *( exp(-y[0]*y[0]/2 ) ) / sqrt(2.0*3.14159265);
    Z0[N_slice-1]= N_slice *( 0-exp(-y[0]*y[0]/2 ) )/ sqrt(2.0*3.14159265);
    for(i=1;i<N_slice-1;i++)  Z0[i]= N_slice *( exp(-y[i]*y[i]/2) - exp(-y[i-1]*y[i-1]/2 )  ) / sqrt(2.0*3.14159265);
    for(i=0;i<N_slice;i++) Z0[i]=  Z0[i] * sigma_l;
  }

  //----Perform Lorentz booster to the test particle's (x,px,y,py,z,delta)
  
  x0= x[0];   px0=x[1];   y0= x[2];
  py0=x[3];   z0= x[4];   pz0=x[5];
  ps0= sqrt( (1+pz0 )*(1+pz0)- px0*px0 - py0*py0 );  
  h0 = (1+pz0)- ps0;

  px_star = ( px0 - tantheta * h0 )  / cos(theta );
  py_star = py0 / costheta;
  pz_star = pz0 - tantheta * px0 + tantheta * tantheta * h0;

  ps_star = sqrt( (1+pz_star )*(1+pz_star)- px_star * px_star - py_star * py_star ); 
  h_star  = (1+pz_star) - ps_star;
  x_star  = tantheta*z0 + ( 1 + (px_star / ps_star) * sintheta ) * x0 ;
  y_star  = y0 + sintheta * (py_star / ps_star) * x0;
  z_star  = z0/costheta + (-h_star  / ps_star) * sintheta * x0;

  //----Perform Lorentz booster to the center of  sliced strong beam
  for(i=0;i<N_slice;i++)
    Z0_star[i]=Z0[i]/costheta; 

  //----Twiss parmeters and emittances in the new frame
  
  gx_star= (1.0 + alfx_star * alfx_star )/ betax_star;
  gy_star= (1.0 + alfy_star * alfy_star )/ betay_star;

  gx_star=gx_star/costheta;
  gy_star=gy_star/costheta;
  betax_star= betax_star *costheta;
  betay_star= betay_star *costheta;
  alfx_star= alfx_star;
  alfy_star= alfy_star;
  emitx_rms = emitx_rms / costheta;
  emity_rms = emity_rms / costheta;

  //-----head-on collision in head-on frame
  lumi_part  = 0 ;

  for(i=0; i<N_slice;i++) {

    x_star_0 =x_star; 
    px_star_0=px_star;
    y_star_0 =y_star; 
    py_star_0=py_star;
    z_star_0 =z_star; 
    pz_star_0=pz_star;

    S=(z_star_0 - Z0_star[i])/2.0;

    betax= betax_star - 2* S * alfx_star + S * S *  gx_star; 
    alfx =(alfx_star -gx_star * S );
    betay= betay_star - 2* S * alfy_star + S * S *  gy_star; 
    alfy =(alfy_star -gy_star * S );

    sigmax  =  sqrt( emitx_rms * betax ) ;
    sigmay  =  sqrt( emity_rms * betay ) ;
    dsigmax2ds =-2.0*alfx* ( emitx_rms );
    dsigmay2ds =-2.0*alfy* ( emity_rms );
    
    X= x_star_0 + px_star_0 * S - Z0_star[i] * sintheta;
    Y= y_star_0 + py_star_0 * S - 0 ;
    BB4D(X, Y, gamma, Np_slice[i], sigmax, sigmay, Dpx, Dpy);
    temp1 = 1.0/2./(sigmax*sigmax-sigmay*sigmay) ;
    temp2 = X *  Dpx + Y *  Dpy;
    temp3 = 2.* Np_slice[i] *rp /gamma ;
    temp4 = exp ( - X * X / 2. / sigmax /sigmax -  Y * Y / 2. / sigmay /sigmay ) ;
    
    x_star =x_star_0  -  S * Dpx * bbscale ;
    px_star=px_star_0 +  Dpx * bbscale;   
    y_star =y_star_0  -  S * Dpy * bbscale ;
    py_star=py_star_0 +  Dpy * bbscale;   
    z_star =z_star_0;

    if ( abs(sigmax- sigmay)/sigmax < 1.0e-6 ) 
      {
      	pz_star=pz_star_0 + 0.5 * Dpx * bbscale * ( px_star_0 + 0.5* Dpx * bbscale ) 
                          + 0.5*  Dpy * bbscale * ( py_star_0 + 0.5* Dpy * bbscale ) 
                          + (1.0/sigmax/sigmax ) * dsigmax2ds  * ( temp3 * bbscale/ 4. ) * temp4;
      }
    else if (sigmax > sigmay ) 
      {
	dUdsigmax2 =  temp1 * bbscale * ( temp2  + temp3 *( sigmay/ sigmax * temp4 -1. ) );
	dUdsigmay2 = -temp1 * bbscale * ( temp2  + temp3 *( sigmax/ sigmay * temp4 -1. ) );
	pz_star= pz_star_0 + 0.5 * Dpx * bbscale * ( px_star_0 + 0.5* Dpx *  bbscale) 
                           + 0.5 * Dpy * bbscale * ( py_star_0 + 0.5* Dpy *  bbscale) 
                           - ( 0.5 * dsigmax2ds * dUdsigmax2 + 0.5 * dsigmay2ds *  dUdsigmay2 ) ;
      }
    else
      {
	dUdsigmay2 = -temp1 *  bbscale * ( temp2  + temp3 *( sigmax/ sigmay * temp4 -1. ) );
	dUdsigmax2 =  temp1 *  bbscale * ( temp2  + temp3 *( sigmay/ sigmax * temp4 -1. ) );
	pz_star    =  pz_star_0 + 0.5 * Dpx *bbscale * ( px_star_0 + 0.5* Dpx * bbscale) 
                                + 0.5 * Dpy *bbscale * ( py_star_0 + 0.5* Dpy * bbscale) 
                                - ( 0.5 * dsigmax2ds * dUdsigmax2 + 0.5 * dsigmay2ds *  dUdsigmay2 ) ;
      }
    
    //-----luminosity contribution  
    lumi_part  =  lumi_part + 1.0 * Np_slice[i] * exp( - X * X /2 /sigmax /sigmax -  Y * Y /2 /sigmay /sigmay  ) / 2 / PI /sigmax / sigmay ; 

   }

  //----perform inverse Lorentz booster to the test particle's coordinates
  ps_star = sqrt( (1+pz_star )*(1+pz_star)- px_star * px_star - py_star * py_star );
  h_star  = 1+ pz_star -  ps_star;
  hx_star = px_star / ps_star ;
  hy_star = py_star / ps_star ;
  hz_star = -h_star / ps_star ;
  
  px0 = (px_star + h_star * sintheta ) * costheta  ;  
  py0 = py_star * costheta;
  pz0 = pz_star  + px_star * sintheta;

  x0  = (x_star - z_star * sintheta ) / ( 1 +  hx_star*sintheta - hz_star *sintheta*sintheta );
  z0  = (z_star - hz_star * sintheta * x0  ) * costheta;
  y0  = y_star - sintheta * hy_star * x0;

  x[0]=x0 ;  x[1]=px0;
  x[2]=y0 ;  x[3]=py0;
  x[4]=z0 ;  x[5]=pz0;
}

void elens_pass_round_Gaussian_topoff(double x[], double gamma, 
	   double Ne, double Le, double beta_e, int N_slice, 
	   double sigmax, double sigmay) 
// I assume sigmax=sigmay, top off [-sigmax*a, sigmax*a ], Ne from perfect Gaussian
{
  int i;
  double Ne_slice= Ne*1.0/ N_slice; 
  double Le_slice= Le*1.0/ N_slice; 
  double a=0.4;                
  double scale;
  double r, rp=1.534698e-18;

  for(i=0;i<N_slice;i++){
    DRIFT_Pass(x, Le_slice/2.);
    r=  sqrt(x[0]*x[0]+x[2]*x[2]);
    if(  r <=  a * sigmax  ) 
      {
	scale = -2 * Ne_slice *(  exp(-a*a/2) * (r/sigmax) * (r/sigmax) /2  ) * rp  / gamma ;
	x[1]=x[1] + scale * x[0]/ r/ r ;
	x[3]=x[3] + scale * x[2]/ r/ r;
      }
    else
      {
	scale = -2 * Ne_slice * ( 1 -exp(-r*r/2/sigmax/sigmax) -  (1 -exp(-a*a/2) ) + exp(-a*a/2) * a * a /2  ) *  rp  / gamma ;
	x[1]=x[1] +  scale * x[0]/ r/r ;
	x[3]=x[3] +  scale * x[2]/ r/r ;
      }
    DRIFT_Pass(x,Le_slice/2.);
  }
}

void elens_pass_round_Gaussian_truncated(double x[], double gamma, 
	   double Ne, double Le, double beta_e, int N_slice, 
	   double sigmax, double sigmay) 
// I assume sigmax=sigmay, Gaussian tail cut off from Nc*sigmax
{
  int i;
  double Le_slice=Le*1.0/ N_slice; 
  double Ne_slice=Ne*1.0/N_slice;
  double Dpx, Dpy;
  double Nc = 100, Nsigma;   
  
  for(i=0;i<N_slice;i++){
    DRIFT_Pass(x, Le_slice/2.);
    Nsigma= sqrt(x[0]*x[0]+ x[2]*x[2]) /  sigmax ;
    if(  Nsigma <=  Nc  ) 
      {
	BB4D(x[0], x[2], gamma, Ne_slice, sigmax, sigmay, Dpx, Dpy);
      }
    else
      {
	BB4D(x[0]* Nc/Nsigma, x[2]* Nc/Nsigma, gamma, Ne_slice, sigmax, sigmay, Dpx, Dpy);      
	Dpx =  Dpx * Nc / Nsigma;
	Dpy =  Dpy * Nc / Nsigma;
      }
    x[1]=x[1]- Dpx * ( 1.+ beta_e);
    x[3]=x[3]- Dpy * ( 1.+ beta_e);
    DRIFT_Pass(x,Le_slice/2.);
  }
}

void elens_pass_round_uniform(double x[], double gamma, 
	   double Ne, double Le, double beta_e, int N_slice, 
	   double sigmax, double sigmay) 
// I assume round uniform distribution in r<=a
{
  int i;
  double Le_slice=Le*1.0/ N_slice; 
  double Ne_slice=Ne*1.0/N_slice;
  double a = 0.31e-3 * 5 ;
  double scale, r2;
  double rp=1.534698e-18;
  
  for(i=0;i<N_slice;i++){
    DRIFT_Pass(x, Le_slice/2.);
    scale=-2*Ne_slice * rp * ( 1 + beta_e ) / GP.gamma;
    r2=(x[0]*x[0]+x[2]*x[2]);
    if(  sqrt(r2) < a  ) 
      {
	x[1]=x[1] + scale * x[0]/ a/a ;
	x[3]=x[3] + scale * x[2]/ a/a ;
      }
    else
      {
	x[1]=x[1] + scale * x[0]/ r2 ;
	x[3]=x[3] + scale * x[2]/ r2 ;
      }
    DRIFT_Pass(x,Le_slice/2.);
  }
}

void BEAMBEAM_Pass(double x[], int TREATMENT, double NP, double BBSCALE, double SIGMAL, int NSLICE, double EMITX,  double BETAX, double ALFAX, double EMITY, double BETAY, double ALFAY)
{
  if ( NP != 0. ) {
      if(int(TREATMENT) == 6){
	BB6D(x, GP.gamma, NP, BBSCALE, SIGMAL, NSLICE, EMITX,  BETAX, ALFAX, EMITY, BETAY, ALFAY); 
      }
      else{
	double Dpx, Dpy;
	BB4D(x[0], x[2], GP.gamma, NP, sqrt( EMITX * BETAX), sqrt( EMITY * BETAY ), Dpx,  Dpy);
	x[1]=x[1]+Dpx * BBSCALE ;
	x[3]=x[3]+Dpy * BBSCALE ;
      }
  } 
}

void LRBB_Pass(double x[], double NP, double BBSCALE, double SEPX, double SEPY, double SIGMAX, double SIGMAY)
{
  x[0]=x[0]+SEPX;  x[2]=x[2]+SEPY;
  BB4D(x, GP.gamma, NP, BBSCALE, SIGMAX, SIGMAY);
  x[0]=x[0]-SEPX;  x[2]=x[2]-SEPY;
}

void CRBB_Pass(double x[], double NP, double THETA, double BBSCALE, double SIGMAL, int NSLICE, double EMITX,  double BETAX, double ALFAX, double EMITY, double BETAY, double ALFAY)
{
  BB6D_Angle(x, GP.gamma, NP, BBSCALE, SIGMAL, THETA, NSLICE, EMITX, BETAX, ALFAX, EMITY, BETAY, ALFAY);
}

void ELENS_Pass(double x[], double Le, double Ne, double bbscale, double beta_e, int N_slice, 
		double sigmax, double sigmay)
{ 
  int i;
  double Le_slice=Le*1.0/ N_slice; 
  double Ne_slice=Ne*1.0/N_slice;
  double Dpx, Dpy;
  
  if(Ne != 0.){
    for(i=0;i<N_slice;i++){
      DRIFT_Pass(x, Le_slice/2.);
      BB4D(x[0], x[2], GP.gamma, Ne_slice, sigmax, sigmay, Dpx, Dpy);
      x[1]=x[1]- Dpx * ( 1.+ beta_e) * bbscale;
      x[3]=x[3]- Dpy * ( 1.+ beta_e) * bbscale;
      DRIFT_Pass(x,Le_slice/2.);
    }
  }
  else{
    DRIFT_Pass(x,Le);
  }
}

void HELENS_Pass(double x[],  double Le, double Ne,  double bbscale, double beta_e, int N_slice,
		 double rout, double rin)
//---written by Xiaofeng Gu
{
  int i;
  double Ne_slice= Ne*1.0/ N_slice;
  double Le_slice= Le*1.0/ N_slice;
  double RINNER = rin;
  double ROUTER = rout;
  double scale;
  double I_e;
  double Theta_max;
  double r, rp=1.534698e-18;

  I_e = 1.6021765e-19*Ne*beta_e*speed_light/Le;
  Theta_max = -2e-7*Le_slice*I_e*(1+beta_e * GP.beta)/beta_e/GP.beta/ROUTER/GP.brho;

  for(i=0;i<N_slice;i++){
    DRIFT_Pass(x, Le_slice/2.);
    
    r=  sqrt(x[0]*x[0]+x[2]*x[2]);
    if(  r <=  RINNER   )
      {
	scale = 0 ;
      }
    else if (  r >  RINNER   && r <=  ROUTER )
      {
	scale = (r*r-RINNER*RINNER)*ROUTER/r/(ROUTER*ROUTER-RINNER*RINNER);
      }
    else if (  r >  ROUTER   )
      {
	scale = ROUTER/r ;
      }
    x[1]=x[1] +  scale * Theta_max * x[0]/r;
    x[3]=x[3] +  scale * Theta_max * x[2]/r;
    
    DRIFT_Pass(x,Le_slice/2.);
  }
  
}

void ERHICBB_Pass(double  x[], double gamma, double Ne, double bbscale)
{
  int i, Nslice=200;
  double z_col[Nslice], x_col[Nslice], y_col[Nslice], sigmax_col[Nslice], sigmay_col[Nslice];
  double x0[6];
  double S, X, Y, xmean, ymean, sigmax, sigmay, Dpx, Dpy;
  double theta=5e-3;  
  FILE *f1;
  double temp1, temp2;

  //----Lorentz transfer
  Lorentz_Transfer(x, theta);

  //----the electron beam sizes at the colliding location
  f1=fopen("./ebeam_info.dat" ,"r");
  for(i=0;i<Nslice;i++){
    fscanf(f1,"%lf  %lf  %lf  %lf  %lf  %lf  %lf", &z_col[i], &x_col[i], &y_col[i], &sigmax_col[i], &sigmay_col[i], &temp1, &temp2);
  }
  fclose(f1); 
  
  //---the colliding location and electron slice sigmas
  for (i=0;i<6;i++) x0[i]=x[i];
  S = x0[4] / 2;
  for(i=0;i<Nslice;i++){
    if( S > z_col[i+1]  and  S < z_col[i] ){
      break;
    }
  }
  if(i==Nslice-1 ) {
    cout<<"ERHICBB: colliding out of the known collision range. exit . "<<endl;
    exit(0);
  }
  xmean=  x_col[i+1] + (x_col[i]-x_col[i+1] ) * ( S - z_col[i+1] ) / ( z_col[i] -  z_col[i+1] );  
  ymean=  y_col[i+1] + (y_col[i]-y_col[i+1] ) * ( S - z_col[i+1] ) / ( z_col[i] -  z_col[i+1] ); 
  sigmax= sigmax_col[i+1] + (  sigmax_col[i]-sigmax_col[i+1] ) * ( S - z_col[i+1] ) / ( z_col[i] -  z_col[i+1] );
  sigmay= sigmay_col[i+1] + (  sigmay_col[i]-sigmay_col[i+1] ) * ( S - z_col[i+1] ) / ( z_col[i] -  z_col[i+1] );
  
  //----calculate the beam-beam given by the electron slice
  X=x0[0] + x0[1]*S;
  Y=x0[2] + x0[3]*S;
  BB4D(X-xmean, Y-ymean, gamma, Ne, sigmax, sigmay, Dpx, Dpy);
  
  x[0]=x0[0] - S * Dpx * bbscale ;
  x[1]=x0[1] + Dpx  * bbscale;   
  x[2]=x0[2] - S * Dpy * bbscale ;
  x[3]=x0[3] + Dpy * bbscale;   
  x[4]=x0[4];
  x[5]=x0[5]; 

  //----Inverse Lorentz Transfer
  Lorentz_Transfer_Inverse(x, theta);
}

void SBEND_sPass(double x[], double L, int Nint, double Angle, double E1, double E2) 
{
  int i,j;
  double href=Angle/L;
  double Lint=L/Nint;
  double edgefocus;
  double BLbrho[3];
  double x0[9],x1[9];
  
  if(Angle==0.){
    DRIFT_Pass(x,L); return;
  }

  edgefocus = tan(E1)*href; 
  x[1] = x[1]+edgefocus*x[0];
  x[3] = x[3]-edgefocus*x[2];

  BLbrho[0] = 0.;  
  BLbrho[1] = Angle/Nint;  
  BLbrho[2] = 0.; 
  if(GP.H_expand == true){
    for(i=0;i<Nint;i++){
      for(j=0;j<9;j++) x0[j]= x[j];
      DRIFT_Pass(x,Fdrift1*Lint);
      bend_kick_pass(x, Fkick1*Lint, href);
      DRIFT_Pass(x,Fdrift2*Lint);
      bend_kick_pass(x, Fkick2*Lint, href);
      DRIFT_Pass(x,Fdrift2*Lint);
      bend_kick_pass(x, Fkick1*Lint, href);
      DRIFT_Pass(x,Fdrift1*Lint); 
      //for(j=0;j<9;j++) x1[j]=(x[j]+x0[j])/2.;
      spin(x0, Angle/Nint, href, BLbrho);
      x[6]= x0[6];  x[7]= x0[7];    x[8]= x0[8];
      if(GP.radiate == true) radiate(x,Lint,href,BLbrho);
     }
  }
  else{
    for(i=0;i<Nint;i++){
      for(j=0;j<9;j++) x0[j]= x[j];
      sbend_exact_pass(x,Lint,Angle/Nint,cos(Angle/Nint),sin(Angle/Nint));
      //for(j=0;j<9;j++) x1[j]=(x[j]+x0[j])/2.;
      spin(x0, Angle/Nint, href, BLbrho);
      x[6]= x0[6];  x[7]= x0[7];    x[8]= x0[8];
      if(GP.radiate == true) radiate(x,Lint,href, BLbrho);
    }
  }
  
  edgefocus = tan(E2)*href; 
  x[1] = x[1]+edgefocus*x[0];
  x[3] = x[3]-edgefocus*x[2];
}

void QUAD_sPass(double x[], double L, int Nint, double k1l, double k1sl)
{
  int i,j;
  double Lint=L/Nint;
  double preal, pimag;
  double BLbrho[3];
  double x0[9],x1[9];
  
  if(k1l==0. and  k1sl ==0.){
    DRIFT_Pass(x,L); return;
  }

  if(L==0.) 
    {
      for(j=0;j<9;j++) x0[j]=x[j];
      quad_kick_pass(x, k1l, k1sl);
      for(j=0;j<9;j++) x1[j]=(x[j]+x0[j])/2.;
      preal = x1[0];
      pimag = x1[2];
      BLbrho[0]= k1sl * preal + k1l  * pimag;
      BLbrho[1]= k1l  * preal - k1sl * pimag;
      BLbrho[2]= 0. ; 
      spin(x0, 0, 0, BLbrho);
      x[6]= x0[6];  x[7]= x0[7];   x[8]= x0[8];
    }
  else 
    {
      double k1l_kick1,k1sl_kick1;
      double k1l_kick2,k1sl_kick2;
      k1l_kick1 =Fkick1*k1l/Nint;
      k1sl_kick1=Fkick1*k1sl/Nint;
      k1l_kick2 =Fkick2*k1l/Nint;
      k1sl_kick2=Fkick2*k1sl/Nint;
      for(i=0;i<Nint;i++){
        for(j=0;j<9;j++) x0[j]=x[j];
	DRIFT_Pass(x,Fdrift1*Lint);
	quad_kick_pass(x, k1l_kick1, k1sl_kick1);
	DRIFT_Pass(x,Fdrift2*Lint);
	quad_kick_pass(x, k1l_kick2, k1sl_kick2);
	DRIFT_Pass(x,Fdrift2*Lint);
	quad_kick_pass(x, k1l_kick1, k1sl_kick1);
	DRIFT_Pass(x,Fdrift1*Lint);
	for(j=0;j<9;j++) x1[j]=(x[j]+x0[j])/2.;
	preal = x1[0];
	pimag = x1[2];
	BLbrho[0]=(k1sl * preal + k1l  * pimag )/Nint;
	BLbrho[1]=(k1l  * preal - k1sl * pimag )/Nint;
	BLbrho[2]= 0.  ;
 	spin(x0, 0, 0, BLbrho);
	x[6]= x0[6];  x[7]= x0[7];  x[8]= x0[8];
      }
    }
} 

void SEXT_sPass(double x[], double L, int Nint, double k2l, double k2sl)
{
  int i,j;
  double Lint=L/Nint;
  double preal, pimag; 
  double BLbrho[3];
  double x0[9], x1[9];

  if(k2l==0. and k2sl==0.){
    DRIFT_Pass(x,L);  return;
  }

  if(L==0.) 
    {
      for(j=0;j<9;j++) x0[j]=x[j];
      sext_kick_pass(x, k2l, k2sl);
      for(j=0;j<9;j++) x1[j]=(x[j]+x0[j])/2.;
      preal =  x1[0]*x1[0]- x1[2]*x1[2];
      pimag =  2. * x1[0] * x1[2];
      BLbrho[0]= ( k2sl * preal + k2l  * pimag )/2.;
      BLbrho[1]= ( k2l  * preal - k2sl * pimag )/2.;
      BLbrho[2]= 0.  ; 
      spin(x0,0,0,BLbrho);
      x[6]= x0[6];  x[7]= x0[7];   x[8]= x0[8];
    }
  else 
    {
      double k2l_kick1,k2sl_kick1;
      double k2l_kick2,k2sl_kick2;
      k2l_kick1 =Fkick1*k2l/Nint;
      k2sl_kick1=Fkick1*k2sl/Nint;
      k2l_kick2 =Fkick2*k2l/Nint;
      k2sl_kick2=Fkick2*k2sl/Nint;
      for(i=0;i<Nint;i++){
	for(j=0;j<9;j++) x0[j]=x[j];
	DRIFT_Pass(x,Fdrift1*Lint);
	sext_kick_pass(x, k2l_kick1, k2sl_kick1);
	DRIFT_Pass(x,Fdrift2*Lint);
	sext_kick_pass(x, k2l_kick2, k2sl_kick2);
	DRIFT_Pass(x,Fdrift2*Lint);
	sext_kick_pass(x, k2l_kick1, k2sl_kick1);
	DRIFT_Pass(x,Fdrift1*Lint);
	for(j=0;j<9;j++) x1[j]=(x[j]+x0[j])/2.;
	preal =  x1[0]*x1[0]- x1[2]*x1[2];
	pimag =  2. * x1[0] * x1[2];
	BLbrho[0]= ( k2sl * preal + k2l  * pimag )/2./Nint;
	BLbrho[1]= ( k2l  * preal - k2sl * pimag )/2./Nint;
	BLbrho[2]= 0.  ; 
	spin(x0,0,0,BLbrho);
	x[6]= x0[6];  x[7]= x0[7];  x[8]= x0[8];
      }
    }
} 

void OCT_sPass(double x[], double L, int Nint, double k3l, double k3sl)
{
  int i,j;
  double Lint=L/Nint;
  double preal, pimag;
  double BLbrho[3];
  double x0[9], x1[9];
  
  if(k3l == 0. and k3sl ==0.){
    DRIFT_Pass(x,L); return;
  }

  if(L==0.) 
    {
      for(j=0;j<9;j++) x0[j]=x[j];
      oct_kick_pass(x, k3l, k3sl);
      for(j=0;j<9;j++) x1[j]=(x[j]+x0[j])/2.;
      preal = x1[0]*x1[0]*x1[0] - 3 * x1[0] * x1[2]* x1[2];
      pimag = 3*x1[0]*x1[0]*x1[2] - x1[2]*x1[2]*x1[2];
      BLbrho[0]= (k3sl * preal + k3l  * pimag)/6/Nint ;
      BLbrho[1]= (k3l  * preal - k3sl * pimag)/6/Nint;
      BLbrho[2]= 0.  ; 
      spin(x0,0,0,BLbrho);
      x[6]= x0[6];  x[7]= x0[7];  x[8]= x0[8];    
    }
  else 
    {
      double k3l_kick1,k3sl_kick1;
      double k3l_kick2,k3sl_kick2;
      k3l_kick1 =Fkick1*k3l/Nint;
      k3sl_kick1=Fkick1*k3sl/Nint;
      k3l_kick2 =Fkick2*k3l/Nint;
      k3sl_kick2=Fkick2*k3sl/Nint;
      for(i=0;i<Nint;i++){
	for(j=0;j<9;j++) x0[j]=x[j];
	DRIFT_Pass(x,Fdrift1*Lint);
	oct_kick_pass(x, k3l_kick1, k3sl_kick1);
	DRIFT_Pass(x,Fdrift2*Lint);
	oct_kick_pass(x, k3l_kick2, k3sl_kick2);
	DRIFT_Pass(x,Fdrift2*Lint);
	oct_kick_pass(x, k3l_kick1, k3sl_kick1);
	DRIFT_Pass(x,Fdrift1*Lint);
	for(j=0;j<9;j++) x1[j]=(x[j]+x0[j])/2.;
	preal = x1[0]*x1[0]*x1[0] - 3 * x1[0] * x1[2]* x1[2];
	pimag = 3*x1[0]*x1[0]*x1[2] - x1[2]*x1[2]*x1[2];
	BLbrho[0]= (k3sl * preal + k3l  * pimag)/6/Nint ;
	BLbrho[1]= (k3l  * preal - k3sl * pimag)/6/Nint;
	BLbrho[2]= 0.  ; 
	spin(x0,0,0,BLbrho);
        x[6]= x0[6];  x[7]= x0[7];  x[8]= x0[8];
      }
    }
}

void MULT_sPass(double x[], double L, int Nint, int Norder, double KNL[11], double KNSL[11]) 
{
  int i,j;
  double Lint=L/Nint;
  double BLbrho[3];
  double x0[9], x1[9];
  
  if(L==0.) 
    {
      cal_Bfield(x, 1, Norder, 0, KNL, KNSL, BLbrho);
      spin(x,0,0,BLbrho);
      mult_kick_pass(x, Norder, KNL, KNSL);
    }
  else 
    {
      double knl_kick1[11],knsl_kick1[11];
      double knl_kick2[11],knsl_kick2[11];
      for(i=0;i<11;i++) {
	knl_kick1[i] =Fkick1*KNL[i]/Nint;
	knsl_kick1[i]=Fkick1*KNSL[i]/Nint;
      }
      for(i=0;i<11;i++) {
	knl_kick2[i] =Fkick2*KNL[i]/Nint;
	knsl_kick2[i]=Fkick2*KNSL[i]/Nint;
      }
      for(i=0;i<Nint;i++){
        for(j=0;j<9;j++) x0[j]=x[j];
	DRIFT_Pass(x,Fdrift1*Lint);
	mult_kick_pass(x, Norder,knl_kick1, knsl_kick1);
	DRIFT_Pass(x,Fdrift2*Lint);
	mult_kick_pass(x, Norder,knl_kick2, knsl_kick2);
	DRIFT_Pass(x,Fdrift2*Lint);
	mult_kick_pass(x, Norder,knl_kick1, knsl_kick1);
	DRIFT_Pass(x,Fdrift1*Lint);
	for(j=0;j<9;j++) x1[j]=(x[j]+x0[j])/2.;
	cal_Bfield(x1, Nint, Norder, 0, KNL, KNSL, BLbrho);
	spin(x0,0,0,BLbrho);
        x[6]= x0[6];  x[7]= x0[7];  x[8]= x0[8];
      }
    }
} 

void GMULT_sPass(double x[], double L, int Nint,  int Norder, double Angle,  double E1, double E2, double KNL[11], double KNSL[11])
{
  int i,j;
  double href=Angle/L;
  double Lint=L/Nint;
  double edgefocus;
  double knl_kick1[11],knsl_kick1[11];
  double knl_kick2[11],knsl_kick2[11];
  double cosAngle, sinAngle, cosAngle1, sinAngle1;
  double BLbrho[3];
  double x0[9],x1[9];
  
  if(L==0.) 
    {
      for(j=0;j<9;j++) x0[j]=x[j];
      mult_kick_pass(x, Norder, KNL, KNSL);
      for(j=0;j<9;j++) x1[j]=(x[j]+x0[j])/2.;
      cal_Bfield(x1, 1, Norder, Angle, KNL, KNSL, BLbrho);
      spin(x0,Angle,0,BLbrho);
      x[6]= x0[6];  x[7]= x0[7];   x[8]= x0[8]; 
      return ;
    }
  else{  

    edgefocus = tan(E1)*href; 
    x[1] = x[1]+edgefocus*x[0];
    x[3] = x[3]-edgefocus*x[2];
    
    if(GP.H_expand == true){  
      
      for(j=0;j<11;j++) {
	knl_kick1[j] =Fkick1*KNL[j]/Nint;
	knsl_kick1[j]=Fkick1*KNSL[j]/Nint;
      }
      for(j=0;j<11;j++) {
	knl_kick2[j] =Fkick2*KNL[j]/Nint;
	knsl_kick2[j]=Fkick2*KNSL[j]/Nint;
      }
      for(i=0;i<Nint;i++){
	for(j=0;j<9;j++) x0[j]=x[j];
	DRIFT_Pass(x,Fdrift1*Lint);
	bend_mult_kick_pass(x, Fkick1*Lint, href, Norder,knl_kick1, knsl_kick1);
	DRIFT_Pass(x,Fdrift2*Lint);
	bend_mult_kick_pass(x, Fkick2*Lint, href, Norder,knl_kick2, knsl_kick2);
	DRIFT_Pass(x,Fdrift2*Lint);
	bend_mult_kick_pass(x, Fkick1*Lint, href, Norder,knl_kick1, knsl_kick1);
	DRIFT_Pass(x,Fdrift1*Lint);
	for(j=0;j<9;j++) x1[j]=(x[j]+x0[j])/2.;
	cal_Bfield(x1, Nint, Norder, Angle, KNL, KNSL, BLbrho);
	spin(x0,Angle/Nint,href,BLbrho);
	x[6]= x0[6];  x[7]= x0[7];    x[8]= x0[8];
	if(GP.radiate == true) radiate(x,Lint,href,BLbrho);
      }

    }
    else{
      
      for(j=0;j<11;j++) {
	knl_kick1[j] =Fkick1*KNL[j]/Nint;
	knsl_kick1[j]=Fkick1*KNSL[j]/Nint;
      }
      for(j=0;j<11;j++) {
	knl_kick2[j] =Fkick2*KNL[j]/Nint;
	knsl_kick2[j]=Fkick2*KNSL[j]/Nint;
      }
      cosAngle1=cos(Fdrift1*Angle/Nint);
      sinAngle1=sin(Fdrift1*Angle/Nint);
      cosAngle=cos(Fdrift2*Angle/Nint);
      sinAngle=sin(Fdrift2*Angle/Nint);
      for(i=0;i<Nint;i++){
	for(j=0;j<9;j++) x0[j]=x[j];
	sbend_exact_pass(x,Fdrift1*Lint,Fdrift1*Angle/Nint,cosAngle1,sinAngle1);
	mult_kick_pass(x, Norder,knl_kick1, knsl_kick1);
	sbend_exact_pass(x,Fdrift2*Lint,Fdrift2*Angle/Nint,cosAngle,sinAngle);
	mult_kick_pass(x, Norder,knl_kick2, knsl_kick2);
	sbend_exact_pass(x,Fdrift2*Lint,Fdrift2*Angle/Nint,cosAngle,sinAngle);
	mult_kick_pass(x, Norder,knl_kick1, knsl_kick1);
	sbend_exact_pass(x,Fdrift1*Lint,Fdrift1*Angle/Nint,cosAngle1,sinAngle1);
	for(j=0;j<9;j++) x1[j]=(x[j]+x0[j])/2.;
	cal_Bfield(x1, Nint, Norder, Angle, KNL, KNSL, BLbrho);
	spin(x0,Angle/Nint,href, BLbrho);
	x[6]= x0[6];  x[7]= x0[7];    x[8]= x0[8];
	if(GP.radiate == true) radiate(x,Lint,href,BLbrho);
      }

    }
    
    edgefocus = tan(E2)*href; 
    x[1] = x[1]+edgefocus*x[0];
    x[3] = x[3]-edgefocus*x[2];

  } 
}

void SBENDMULT_sPass(double x[], double L, int Nint,  int Norder, double Angle,  double E1, double E2, double KNL[11], double KNSL[11])
{
  int i,j;
  double href=Angle/L;
  double Lint=L/Nint;
  double edgefocus;
  double knl_kick1[11],knsl_kick1[11];
  double knl_kick2[11],knsl_kick2[11];
  double cosAngle, sinAngle, cosAngle1, sinAngle1;
  double BLbrho[3];
  double x0[9],x1[9];
  
  if(L==0.) 
    {
      for(j=0;j<9;j++) x0[j]=x[j];
      mult_kick_pass(x, Norder, KNL, KNSL);
      for(j=0;j<9;j++) x1[j]=(x[j]+x0[j])/2.;
      cal_Bfield(x1, 1, Norder, Angle, KNL, KNSL, BLbrho);
      spin(x0,Angle,0,BLbrho);
      x[6]= x0[6];  x[7]= x0[7];   x[8]= x0[8]; 
      return ;
    }
  else{  

    edgefocus = tan(E1)*href; 
    x[1] = x[1]+edgefocus*x[0];
    x[3] = x[3]-edgefocus*x[2];
    
    if(GP.H_expand == true){  
      
      for(j=0;j<11;j++) {
	knl_kick1[j] =Fkick1*KNL[j]/Nint;
	knsl_kick1[j]=Fkick1*KNSL[j]/Nint;
      }
      for(j=0;j<11;j++) {
	knl_kick2[j] =Fkick2*KNL[j]/Nint;
	knsl_kick2[j]=Fkick2*KNSL[j]/Nint;
      }
      for(i=0;i<Nint;i++){
	for(j=0;j<9;j++) x0[j]=x[j];
	DRIFT_Pass(x,Fdrift1*Lint);
	bend_mult_kick_pass(x, Fkick1*Lint, href, Norder,knl_kick1, knsl_kick1);
	DRIFT_Pass(x,Fdrift2*Lint);
	bend_mult_kick_pass(x, Fkick2*Lint, href, Norder,knl_kick2, knsl_kick2);
	DRIFT_Pass(x,Fdrift2*Lint);
	bend_mult_kick_pass(x, Fkick1*Lint, href, Norder,knl_kick1, knsl_kick1);
	DRIFT_Pass(x,Fdrift1*Lint);
	for(j=0;j<9;j++) x1[j]=(x[j]+x0[j])/2.;
	cal_Bfield(x1, Nint, Norder, Angle, KNL, KNSL, BLbrho);
	spin(x0,Angle/Nint,href,BLbrho);
	x[6]= x0[6];  x[7]= x0[7];    x[8]= x0[8];
	if(GP.radiate == true) radiate(x,Lint,href,BLbrho);
      }

    }
    else{
      
      for(j=0;j<11;j++) {
	knl_kick1[j] =Fkick1*KNL[j]/Nint;
	knsl_kick1[j]=Fkick1*KNSL[j]/Nint;
      }
      for(j=0;j<11;j++) {
	knl_kick2[j] =Fkick2*KNL[j]/Nint;
	knsl_kick2[j]=Fkick2*KNSL[j]/Nint;
      }
      cosAngle1=cos(Fdrift1*Angle/Nint);
      sinAngle1=sin(Fdrift1*Angle/Nint);
      cosAngle=cos(Fdrift2*Angle/Nint);
      sinAngle=sin(Fdrift2*Angle/Nint);
      for(i=0;i<Nint;i++){
	for(j=0;j<9;j++) x0[j]=x[j];
	sbend_exact_pass(x,Fdrift1*Lint,Fdrift1*Angle/Nint,cosAngle1,sinAngle1);
	mult_kick_pass(x, Norder,knl_kick1, knsl_kick1);
	sbend_exact_pass(x,Fdrift2*Lint,Fdrift2*Angle/Nint,cosAngle,sinAngle);
	mult_kick_pass(x, Norder,knl_kick2, knsl_kick2);
	sbend_exact_pass(x,Fdrift2*Lint,Fdrift2*Angle/Nint,cosAngle,sinAngle);
	mult_kick_pass(x, Norder,knl_kick1, knsl_kick1);
	sbend_exact_pass(x,Fdrift1*Lint,Fdrift1*Angle/Nint,cosAngle1,sinAngle1);
	for(j=0;j<9;j++) x1[j]=(x[j]+x0[j])/2.;
	cal_Bfield(x1, Nint, Norder, Angle, KNL, KNSL, BLbrho);
	spin(x0,Angle/Nint,href, BLbrho);
	x[6]= x0[6];  x[7]= x0[7];    x[8]= x0[8];
	if(GP.radiate == true) radiate(x,Lint,href,BLbrho);
      }

    }
    
    edgefocus = tan(E2)*href; 
    x[1] = x[1]+edgefocus*x[0];
    x[3] = x[3]-edgefocus*x[2];

  } 
}

void SMULT_sPass(double x[], double L, int Nint,  int Norder, double Angle,  double E1, double E2, double KNL[11], double KNSL[11])
{
  int i,j;
  double href=Angle/L;
  double Lint=L/Nint;
  double edgefocus;
  double knl_kick1[11],knsl_kick1[11];
  double knl_kick2[11],knsl_kick2[11];
  double cosAngle, sinAngle, cosAngle1, sinAngle1;
  double BLbrho[3];
  double x0[9],x1[9];
  
  if(L==0.) 
    {
      for(j=0;j<9;j++) x0[j]=x[j];
      mult_kick_pass(x, Norder, KNL, KNSL);
      for(j=0;j<9;j++) x1[j]=(x[j]+x0[j])/2.;
      cal_Bfield(x1, 1, Norder, Angle, KNL, KNSL, BLbrho);
      spin(x0,Angle,0,BLbrho);
      x[6]= x0[6];  x[7]= x0[7];   x[8]= x0[8]; 
      return ;
    }
  else{  

    edgefocus = tan(E1)*href; 
    x[1] = x[1]+edgefocus*x[0];
    x[3] = x[3]-edgefocus*x[2];
    
    if(GP.H_expand == true){  
      
      for(j=0;j<11;j++) {
	knl_kick1[j] =Fkick1*KNL[j]/Nint;
	knsl_kick1[j]=Fkick1*KNSL[j]/Nint;
      }
      for(j=0;j<11;j++) {
	knl_kick2[j] =Fkick2*KNL[j]/Nint;
	knsl_kick2[j]=Fkick2*KNSL[j]/Nint;
      }
      for(i=0;i<Nint;i++){
	for(j=0;j<9;j++) x0[j]=x[j];
	DRIFT_Pass(x,Fdrift1*Lint);
	bend_mult_kick_pass(x, Fkick1*Lint, href, Norder,knl_kick1, knsl_kick1);
	DRIFT_Pass(x,Fdrift2*Lint);
	bend_mult_kick_pass(x, Fkick2*Lint, href, Norder,knl_kick2, knsl_kick2);
	DRIFT_Pass(x,Fdrift2*Lint);
	bend_mult_kick_pass(x, Fkick1*Lint, href, Norder,knl_kick1, knsl_kick1);
	DRIFT_Pass(x,Fdrift1*Lint);
	for(j=0;j<9;j++) x1[j]=(x[j]+x0[j])/2.;
	cal_Bfield(x1, Nint, Norder, Angle, KNL, KNSL, BLbrho);
	spin(x0,Angle/Nint,href,BLbrho);
	x[6]= x0[6];  x[7]= x0[7];    x[8]= x0[8];
	if(GP.radiate == true) radiate(x,Lint,href,BLbrho);
      }

    }
    else{
      
      for(j=0;j<11;j++) {
	knl_kick1[j] =Fkick1*KNL[j]/Nint;
	knsl_kick1[j]=Fkick1*KNSL[j]/Nint;
      }
      for(j=0;j<11;j++) {
	knl_kick2[j] =Fkick2*KNL[j]/Nint;
	knsl_kick2[j]=Fkick2*KNSL[j]/Nint;
      }
      cosAngle1=cos(Fdrift1*Angle/Nint);
      sinAngle1=sin(Fdrift1*Angle/Nint);
      cosAngle=cos(Fdrift2*Angle/Nint);
      sinAngle=sin(Fdrift2*Angle/Nint);
      for(i=0;i<Nint;i++){
	for(j=0;j<9;j++) x0[j]=x[j];
	sbend_exact_pass(x,Fdrift1*Lint,Fdrift1*Angle/Nint,cosAngle1,sinAngle1);
	mult_kick_pass(x, Norder,knl_kick1, knsl_kick1);
	sbend_exact_pass(x,Fdrift2*Lint,Fdrift2*Angle/Nint,cosAngle,sinAngle);
	mult_kick_pass(x, Norder,knl_kick2, knsl_kick2);
	sbend_exact_pass(x,Fdrift2*Lint,Fdrift2*Angle/Nint,cosAngle,sinAngle);
	mult_kick_pass(x, Norder,knl_kick1, knsl_kick1);
	sbend_exact_pass(x,Fdrift1*Lint,Fdrift1*Angle/Nint,cosAngle1,sinAngle1);
	for(j=0;j<9;j++) x1[j]=(x[j]+x0[j])/2.;
	cal_Bfield(x1, Nint, Norder, Angle, KNL, KNSL, BLbrho);
	spin(x0,Angle/Nint,href, BLbrho);
	x[6]= x0[6];  x[7]= x0[7];    x[8]= x0[8];
	if(GP.radiate == true) radiate(x,Lint,href,BLbrho);
      }

    }
    
    edgefocus = tan(E2)*href; 
    x[1] = x[1]+edgefocus*x[0];
    x[3] = x[3]-edgefocus*x[2];

  } 
}

void GSBENDMULT_sPass(double x[], double L, int Nint,  int Norder, double Angle, double K0L, double E1, double E2, double KNL[11], double KNSL[11])
{
  int i,j;
  double href=Angle/L;
  double Lint=L/Nint;
  double edgefocus;
  double knl_kick1[11],knsl_kick1[11];
  double knl_kick2[11],knsl_kick2[11];
  double cosAngle, sinAngle, cosAngle1, sinAngle1;
  double BLbrho[3];
  double x0[9],x1[9];
  
  if(L==0.) 
    {
      for(j=0;j<9;j++) x0[j]=x[j];
      mult_kick_pass(x, Norder, KNL, KNSL);
      for(j=0;j<9;j++) x1[j]=(x[j]+x0[j])/2.;
      cal_Bfield(x1, 1, Norder, Angle, KNL, KNSL, BLbrho);
      spin(x0,Angle,0,BLbrho);
      x[6]= x0[6];  x[7]= x0[7];   x[8]= x0[8]; 
      return ;
    }
  else{  

    edgefocus = tan(E1)*href; 
    x[1] = x[1]+edgefocus*x[0];
    x[3] = x[3]-edgefocus*x[2];
    
    if(GP.H_expand == true){  
      
      for(j=0;j<11;j++) {
	knl_kick1[j] =Fkick1*KNL[j]/Nint;
	knsl_kick1[j]=Fkick1*KNSL[j]/Nint;
      }
      for(j=0;j<11;j++) {
	knl_kick2[j] =Fkick2*KNL[j]/Nint;
	knsl_kick2[j]=Fkick2*KNSL[j]/Nint;
      }
      for(i=0;i<Nint;i++){
	for(j=0;j<9;j++) x0[j]=x[j];
	DRIFT_Pass(x,Fdrift1*Lint);
	bend_mult_kick_pass(x, Fkick1*Lint, href, Norder,knl_kick1, knsl_kick1);
	DRIFT_Pass(x,Fdrift2*Lint);
	bend_mult_kick_pass(x, Fkick2*Lint, href, Norder,knl_kick2, knsl_kick2);
	DRIFT_Pass(x,Fdrift2*Lint);
	bend_mult_kick_pass(x, Fkick1*Lint, href, Norder,knl_kick1, knsl_kick1);
	DRIFT_Pass(x,Fdrift1*Lint);
	for(j=0;j<9;j++) x1[j]=(x[j]+x0[j])/2.;
	cal_Bfield(x1, Nint, Norder, Angle, KNL, KNSL, BLbrho);
	spin(x0,Angle/Nint,href,BLbrho);
	x[6]= x0[6];  x[7]= x0[7];    x[8]= x0[8];
	if(GP.radiate == true) radiate(x,Lint,href,BLbrho);
      }

    }
    else{
      
      for(j=0;j<11;j++) {
	knl_kick1[j] =Fkick1*KNL[j]/Nint;
	knsl_kick1[j]=Fkick1*KNSL[j]/Nint;
      }
      for(j=0;j<11;j++) {
	knl_kick2[j] =Fkick2*KNL[j]/Nint;
	knsl_kick2[j]=Fkick2*KNSL[j]/Nint;
      }
      cosAngle1=cos(Fdrift1*Angle/Nint);
      sinAngle1=sin(Fdrift1*Angle/Nint);
      cosAngle=cos(Fdrift2*Angle/Nint);
      sinAngle=sin(Fdrift2*Angle/Nint);
      for(i=0;i<Nint;i++){
	for(j=0;j<9;j++) x0[j]=x[j];
	sbend_exact_pass(x,Fdrift1*Lint,Fdrift1*Angle/Nint,cosAngle1,sinAngle1);
	mult_kick_pass(x, Norder,knl_kick1, knsl_kick1);
	sbend_exact_pass(x,Fdrift2*Lint,Fdrift2*Angle/Nint,cosAngle,sinAngle);
	mult_kick_pass(x, Norder,knl_kick2, knsl_kick2);
	sbend_exact_pass(x,Fdrift2*Lint,Fdrift2*Angle/Nint,cosAngle,sinAngle);
	mult_kick_pass(x, Norder,knl_kick1, knsl_kick1);
	sbend_exact_pass(x,Fdrift1*Lint,Fdrift1*Angle/Nint,cosAngle1,sinAngle1);
	for(j=0;j<9;j++) x1[j]=(x[j]+x0[j])/2.;
	cal_Bfield(x1, Nint, Norder, Angle, KNL, KNSL, BLbrho);
	spin(x0,Angle/Nint,href, BLbrho);
	x[6]= x0[6];  x[7]= x0[7];    x[8]= x0[8];
	if(GP.radiate == true) radiate(x,Lint,href,BLbrho);
      }

    }
    
    edgefocus = tan(E2)*href; 
    x[1] = x[1]+edgefocus*x[0];
    x[3] = x[3]-edgefocus*x[2];

  } 
}

void SOLEN_sPass(double x[], double L, int Nint, double KS)
{
  int i, j, k;
  double   ps, gamma1, delta1, beta1, vz;
  double   Bz= KS * GP.brho / L, K, Lint=L/Nint, theta, sintheta, costheta,m44[16];
  double   BLbrho[3], xtemp[9];
  
  if(KS == 0 ) {
    DRIFT_Pass(x, L);
    return;
  }
  
  //----switch from (px, py) to (x', y')
  delta1= sqrt(1.0 + 2*x[5]/GP.beta+x[5]*x[5]) -1.0;
  gamma1 =GP.gamma+sqrt(GP.gamma*GP.gamma-1.)*x[5];
  beta1 = sqrt(1.0-1.0/gamma1/gamma1);
  
  ps=sqrt( ( 1+delta1 ) * (1+delta1 ) - x[1]*x[1] - x[3]*x[3] ); 
  x[1]=x[1]/ps; 
  x[3]=x[3]/ps;

  //----entrance
  vz =  2.99792458e8 * beta1 / sqrt( 1.0 + x[1]*x[1] +x[3]*x[3] );  
  K  =  GP.q * Bz / 2.0 / gamma1 / GP.m0 / vz;
  x[1] = x[1] + x[2] * K ;
  x[3] = x[3] - x[0] * K ;

  //---body
  theta =  -2.0 * K * L / Nint;
  sintheta=sin(theta); costheta=cos(theta);

  sintheta = sin(theta); costheta = cos(theta); 
  m44[0*4+0] = 1.0 ;   m44[0*4+1] = Lint*sintheta/theta;      m44[0*4+2] = 0. ;   m44[0*4+3] = -(1-costheta)*Lint/theta;
  m44[1*4+0] = 0. ;    m44[1*4+1] = costheta;                 m44[1*4+2] = 0. ;   m44[1*4+3] = -sintheta;
  m44[2*4+0] = 0. ;    m44[2*4+1] = (1-costheta)*Lint/theta;  m44[2*4+2] = 1. ;   m44[2*4+3] = Lint*sintheta/theta;
  m44[3*4+0] = 0. ;    m44[3*4+1] = sintheta;                 m44[3*4+2] = 0. ;   m44[3*4+3] = costheta;

  for(i=0;i<Nint; i++){
    for(k=0;k<9;k++) xtemp[k]=x[k];
    for(k=0;k<4;k++){
      x[k]=0.;
      for(j=0;j<4;j++) x[k]= x[k]+ m44[k*4+j]*xtemp[j];
    }
    
    xtemp[1]=gamma1* GP.m0 * xtemp[1] * vz / GP.p; 
    xtemp[3]=gamma1* GP.m0 * xtemp[3] * vz / GP.p; 
    BLbrho[0]=0.;    BLbrho[1]=0.;    BLbrho[2]=Bz*Lint/GP.brho;
    spin(xtemp, 0, 0, BLbrho);
    x[6]= xtemp[6];  x[7]= xtemp[7];    x[8]= xtemp[8];
  }

  //----exit
  x[1] = x[1] - x[2] * K ;
  x[3] = x[3] + x[0] * K ;
      
  //---switch to (x, px, y, py)
  x[1]=gamma1* GP.m0 * x[1] * vz / GP.p; 
  x[3]=gamma1* GP.m0 * x[3] * vz / GP.p; 
  x[4]=x[4] - 2.99792458e8 * ( L/vz - L/GP.beta/ 2.99792458e8);
 
}

void WIGGLER_sPass(double x[], double L, int Nint, double B0, double KX, double KZ, double PHIZ0)
{
  int i, j, k;
  double xtemp[9];
  
  if(B0 == 0 ) {
    DRIFT_Pass(x, L);}
  else{
    DRIFT_Pass(x, L);
  }
  
}

void KICK_sPass(double x[], double L, double HKICK, double VKICK) 
{
  double BLbrho[3];
  
  DRIFT_Pass(x,L/2.0); 
  
  BLbrho[0] =  VKICK;
  BLbrho[1] = -HKICK;
  BLbrho[2] =  0.;
  spin(x,0,0,BLbrho);
  
  x[1]=x[1]+HKICK;     
  x[3]=x[3]+VKICK;
  
  DRIFT_Pass(x,L/2.0);
}

void   ACMULT_sPass(double x[], double L, int Norder, double KL, double KSL, int TTURNS, double PHI0)
{
  int i;
  int fac=1;
  double Xn, Yn, Xn0, Yn0;
  double KNL, KNSL;
  double By, Bx;
  double BLbrho[3];

  DRIFT_Pass(x, L/2.0);
  
  if(Norder==0) {
    By=KL;
    Bx=KSL;
  }
  else{
    KNL=KL;
    KNSL=KSL;
    By=0.;  Bx=0.;  Xn=1.;  Yn=0.;
    for(i=1;i<Norder+1;i++){
      Xn0=Xn;
      Yn0=Yn;
      Xn=Xn0*x[x_]-Yn0*x[y_];
      Yn=Xn0*x[y_]+Yn0*x[x_];
      fac=fac*i;
    }
    By=By+(KNL*Xn-KNSL*Yn)/fac;
    Bx=Bx+(KNL*Yn+KNSL*Xn)/fac;
  }
  
  BLbrho[0] = Bx* sin( 2.0 * PI * GP.turn / TTURNS + PHI0);
  BLbrho[1] = By* sin( 2.0 * PI * GP.turn / TTURNS + PHI0);
  BLbrho[2] = 0.;
  spin(x,0,0,BLbrho);

  x[px_]=x[px_]- By * sin( 2.0 * PI * GP.turn / TTURNS + PHI0);
  x[py_]=x[py_]+ Bx * sin( 2.0 * PI * GP.turn / TTURNS + PHI0);
  
  DRIFT_Pass(x, L/2.0);
}

void ACDIP_sPass(double x[], double L, double HKICKMAX, double VKICKMAX, double NUD, double TURNS, double TURNE, double PHID) 
{
  double BLbrho[3];
  
  if( GP.turn < TURNS ) {
    DRIFT_Pass(x, L);
  }
  else if ( GP.turn >= TURNS  and GP.turn <= TURNE  ){
    DRIFT_Pass(x, L/2.);
    
    BLbrho[0] = ((GP.turn-TURNS)*1.0 * VKICKMAX /(TURNE-TURNS )) * sin( 2.0* 3.14159265*NUD*(GP.turn-TURNS)+PHID);
    BLbrho[1] =-((GP.turn-TURNS)*1.0 * HKICKMAX /(TURNE-TURNS )) * sin( 2.0* 3.14159265*NUD*(GP.turn-TURNS)+PHID); 
    BLbrho[2] = 0.;
    spin(x,0,0,BLbrho);
    
    x[1] += ((GP.turn-TURNS)*1.0 * HKICKMAX /(TURNE-TURNS )) * sin( 2.0* 3.14159265*NUD*(GP.turn-TURNS)+PHID);
    x[3] += ((GP.turn-TURNS)*1.0 * VKICKMAX /(TURNE-TURNS )) * sin( 2.0* 3.14159265*NUD*(GP.turn-TURNS)+PHID);

    DRIFT_Pass(x, L/2.);
  }
  else {
    DRIFT_Pass(x, L/2.);
    
    BLbrho[0] =  VKICKMAX* sin( 2.0* 3.14159265 * NUD *(GP.turn-TURNS) +PHID );
    BLbrho[1] = -HKICKMAX* sin( 2.0* 3.14159265 * NUD *(GP.turn-TURNS) +PHID );
    BLbrho[2] =  0;
    spin(x,0,0,BLbrho);
    
    x[1] +=  HKICKMAX* sin( 2.0* 3.14159265 * NUD *(GP.turn-TURNS) +PHID );
    x[3] +=  VKICKMAX* sin( 2.0* 3.14159265 * NUD *(GP.turn-TURNS) +PHID );
    
    DRIFT_Pass(x, L/2.);
  }
}

void ROTAT_sPass(double x[9], double L, double n[3], double angle)
{
  int i,j;
  double spin0[3], spin[3];  
  double F, uF[3], costheta=cos(-angle), sintheta=sin(-angle), rot[9]; 
 
  DRIFT_Pass(x, L/2.);

  for(i=0;i<3;i++) spin0[i]=x[6+i];
  
  F=sqrt( n[0]*n[0] +  n[1]*n[1] +  n[2]*n[2] );
  for(i=0;i<3;i++) uF[i]= n[i]/F;
  
  rot[0] = uF[0]*uF[0]*(1-costheta) + costheta ;  
  rot[1] = uF[0]*uF[1]*(1-costheta) - uF[2]*sintheta;  
  rot[2] = uF[0]*uF[2]*(1-costheta) + uF[1]*sintheta;  
  rot[3] = uF[0]*uF[1]*(1-costheta) + uF[2]*sintheta;  
  rot[4] = uF[1]*uF[1]*(1-costheta) + costheta ;   
  rot[5] = uF[1]*uF[2]*(1-costheta) - uF[0]*sintheta;  
  rot[6] = uF[0]*uF[2]*(1-costheta) - uF[1]*sintheta; 
  rot[7] = uF[1]*uF[2]*(1-costheta) + uF[0]*sintheta;  
  rot[8] = uF[2]*uF[2]*(1-costheta) + costheta;
  
  for(i=0;i<3;i++){
    spin[i]=0.;
    for(j=0;j<3;j++) spin[i] = spin[i] + rot[i*3 + j] * spin0[j];
  }
  for(i=0;i<3;i++) x[6+i]= spin[i];
  
  DRIFT_Pass(x,L/2.);
}

void SNAKE_sPass(double x[9], double L, double n[3], double angle)
{
  int i,j;
  double spin0[3], spin[3];  
  double F, uF[3], costheta=cos(-angle), sintheta=sin(-angle), rot[9]; 
 
  DRIFT_Pass(x, L/2.);

  for(i=0;i<3;i++) spin0[i]=x[6+i];
  
  F=sqrt( n[0]*n[0] +  n[1]*n[1] +  n[2]*n[2] );
  for(i=0;i<3;i++) uF[i]= n[i]/F;
  
  rot[0] = uF[0]*uF[0]*(1-costheta) + costheta ;  
  rot[1] = uF[0]*uF[1]*(1-costheta) - uF[2]*sintheta;  
  rot[2] = uF[0]*uF[2]*(1-costheta) + uF[1]*sintheta;  
  rot[3] = uF[0]*uF[1]*(1-costheta) + uF[2]*sintheta;  
  rot[4] = uF[1]*uF[1]*(1-costheta) + costheta ;   
  rot[5] = uF[1]*uF[2]*(1-costheta) - uF[0]*sintheta;  
  rot[6] = uF[0]*uF[2]*(1-costheta) - uF[1]*sintheta; 
  rot[7] = uF[1]*uF[2]*(1-costheta) + uF[0]*sintheta;  
  rot[8] = uF[2]*uF[2]*(1-costheta) + costheta;
  
  for(i=0;i<3;i++){
    spin[i]=0.;
    for(j=0;j<3;j++) spin[i] = spin[i] + rot[i*3 + j] * spin0[j];
  }
  for(i=0;i<3;i++) x[6+i]= spin[i];
  
  DRIFT_Pass(x,L/2.);
}


//===========================================
//
//       Element definition
//
//============================================

Element:: Element(string name)
{
  int i;
  NAME=name;
  DX=0.;      // offset of magnet center w.r.t local x-y frame, DX>0 means magnet center moved to positive x
  DY=0.;      // offset of magnet center w.r.t local x-y frame, DY>0 means magnet center moved to positive y
  DS=0.;
  //DT=0.;       replaced by DPSI, tilt angle of magnet w.r.t x-y frame, roll from e_x-->e_y, DT > 0 ;  Normal Q rolled -PI/4 to get a skewQ, k1s=k1
  DTHETA = 0;
  DPHI = 0;
  DPSI = 0;
  for(i=0;i<6;i++)   X[i]=0.;
  for(i=0;i<36;i++)  T[i]=0.;
  for(i=0;i<36;i++)  M[i]=0.;
  for(i=0;i<36;i++)  A[i]=0.;
  Beta1=0.;
  Beta2=0.;
  Beta3=0.;
  Alfa1=0.;
  Alfa2=0.;
  Alfa3=0.;
  r=0;
  c11=0.;
  c12=0.;
  c21=0.;
  c22=0.;
  Etax=0.;
  Etay=0.;
  Etaxp=0.;
  Etayp=0.;
  Ksix=0.;
  Ksiy=0.;
  Ksixp=0.;
  Ksiyp=0.;
  Mu1=0.;
  Mu2=0.;
  Mu3=0.;
  APx=1.;
  APy=1.;
}

//---------------DRIFT-----------------------------------
DRIFT::DRIFT(string name, double l): Element(name)
{ 
  if (l >=0 )
    {
      TYPE=string("DRIFT");
      GROUP=string("");
      L=l;
    }
  else
    {
      cout<<"Error: check signs of L and S of DRFIT "<<name<<endl;
      exit(1); 
    }
}

void DRIFT::SetP(const char *name, double value) 
{
  cout<<"No parameter to be set for DRIFT."<<endl;
  exit(0);
}

double DRIFT::GetP(const char *name) 
{
  cout<<"No parameter to be returned for DRIFT."<<endl;
  exit(0);
}

void DRIFT::Pass(double x[6]){
  DRIFT_Pass(x, L);
}

void DRIFT::DAPass(tps x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  DRIFT_Pass(x, L);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void DRIFT::sPass(double x[9]){
  DRIFT_Pass(x, L);
}

//--------------------------SBEND------------------------------------
SBEND::SBEND(string name, double l, double angle, double e1, double e2): Element(name)
{
  if (l >0 )
    { 
      TYPE=string("SBEND");
      GROUP=string("");
      L=l;
      ANGLE=angle;
      E1=e1;
      E2=e2;
      Nint=int(L/BLslice)+1;
    }
  else
    {
      cout<<"Error: SBEND length should be positive. "<<endl;
      exit(1); 
    }
}

void SBEND::SetP(const char *name, double value)
{
  if (strcmp( name, "ANGLE" ) == 0) 
    {
      ANGLE=value;
    }
  else if (strcmp( name, "E1" ) == 0) 
    {
      E1=value;
    }
  else if (strcmp( name, "E2" ) == 0) 
    {
      E2=value;
    }
  else if (strcmp( name, "Nint" ) == 0) 
    {
      Nint=int(value);
    }
  else 
    {
      cout<<"SBEND does not have  a parameter of "<<name<<endl; 
      exit(0); 
    } 
} 

double SBEND::GetP(const char *name)
{
  if (strcmp( name, "ANGLE" ) == 0) 
    {
      return ANGLE;
    }
  else if (strcmp( name, "E1" ) == 0) 
    {
      return E1;
    }
  else if (strcmp( name, "E2" ) == 0) 
    {
      return E2;
    }
  else if (strcmp( name, "Nint" ) == 0) 
    {
      return Nint;
    }
  else 
    {
      cout<<"SBEND does not have a parameter of "<<name<<endl; 
      exit(0); 
    } 
}

void SBEND::Pass(double x[6]){
  GtoL(x,L,ANGLE,DX,DY,DS,DTHETA,DPHI,DPSI);
  SBEND_Pass(x,L,Nint,ANGLE,E1,E2);
  LtoG(x,L,ANGLE,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void SBEND::DAPass(tps x[6]){
  GtoL(x,L,ANGLE,DX,DY,DS,DTHETA,DPHI,DPSI);
  SBEND_DAPass(x,L,Nint,ANGLE,E1,E2);
  LtoG(x,L,ANGLE,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void SBEND::sPass(double x[9]){
  GtoL(x,L,ANGLE,DX,DY,DS,DTHETA,DPHI,DPSI);
  SBEND_sPass(x,L,Nint,ANGLE,E1,E2);
  LtoG(x,L,ANGLE,DX,DY,DS,DTHETA,DPHI,DPSI);
}

//--------------------------GSBEND------------------------------------
GSBEND::GSBEND(string name, double l, double angle, double k0l, double e1, double e2): Element(name)
{
  if (l >0 )
    { 
      TYPE=string("GSBEND");
      GROUP=string("");
      L=l;
      ANGLE=angle;
      K0L=k0l;
      E1=e1;
      E2=e2;
      Nint=int(L/BLslice)+1;
    }
  else
    {
      cout<<"Error: GSBEND length should be positive. "<<endl;
      exit(1); 
    }
}

void GSBEND::SetP(const char *name, double value)
{
  if (strcmp( name, "ANGLE" ) == 0) 
    {
      ANGLE=value;
    }
  else if (strcmp( name, "K0L" ) == 0) 
    {
      K0L=value;
    }  
  else if (strcmp( name, "E1" ) == 0) 
    {
      E1=value;
    }
  else if (strcmp( name, "E2" ) == 0) 
    {
      E2=value;
    }
  else if (strcmp( name, "Nint" ) == 0) 
    {
      Nint=int(value);
    }
  else 
    {
      cout<<"GSBEND does not have  a parameter of "<<name<<endl; 
      exit(0); 
    } 
} 

double GSBEND::GetP(const char *name)
{
  if (strcmp( name, "ANGLE" ) == 0) 
    {
      return ANGLE;
    }
  else if (strcmp( name, "K0L" ) == 0) 
    {
      return K0L;
    }  
  else if (strcmp( name, "E1" ) == 0) 
    {
      return E1;
    }
  else if (strcmp( name, "E2" ) == 0) 
    {
      return E2;
    }
  else if (strcmp( name, "Nint" ) == 0) 
    {
      return Nint;
    }
  else 
    {
      cout<<"GSBEND does not have a parameter of "<<name<<endl; 
      exit(0); 
    } 
}

void GSBEND::Pass(double x[6]){
  GtoL(x,L,ANGLE,DX,DY,DS,DTHETA,DPHI,DPSI);
  GSBEND_Pass(x,L,Nint,ANGLE,K0L,E1,E2);
  LtoG(x,L,ANGLE,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void GSBEND::DAPass(tps x[6]){
  GtoL(x,L,ANGLE,DX,DY,DS,DTHETA,DPHI,DPSI);
  GSBEND_DAPass(x,L,Nint,ANGLE,K0L,E1,E2);
  LtoG(x,L,ANGLE,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void GSBEND::sPass(double x[9]){
  GtoL(x,L,ANGLE,DX,DY,DS,DTHETA,DPHI,DPSI);
  SBEND_sPass(x,L,Nint,ANGLE,E1,E2);         // 
  LtoG(x,L,ANGLE,DX,DY,DS,DTHETA,DPHI,DPSI); 
}

//---------------------QUAD-----------------------------------------
QUAD::QUAD(string name, double l, double k1l, double k1sl ): Element(name)
{ 
  if (l >=0 )
    {
      TYPE=string("QUAD");
      GROUP=string("");
      L=l;
      K1L=k1l;
      K1SL=k1sl;
      Nint=int(L/QLslice)+1;
      Norder=1;
    }
  else
    {
      cout<<"Error: QUAD length should be non-negative."<<endl;
      exit(1); 
    }
}

void QUAD::SetP(const char *name, double value)
{ 
  if (strcmp( name, "K1L" ) == 0) 
    {
      K1L=value;
    }
  else if (strcmp( name, "K1SL" ) == 0) 
    {
      K1SL=value;
    }
  else if (strcmp( name, "Nint" ) == 0) 
    {
      Nint=int(value);
    }
  else 
    {
      cout<<"QUAD does not have a parameter of "<<name<<endl; 
      exit(0); 
    } 
} 

double QUAD::GetP(const char *name)
{
  if (strcmp( name, "K1L" ) == 0) 
    {
      return K1L;
    }
  else if (strcmp( name, "K1SL" ) == 0) 
    {
      return K1SL;
    }
  else if (strcmp( name, "Nint" ) == 0) 
    {
      return Nint;
    }
  else 
    {
      cout<<"QUAD does not have a parameter of "<<name<<endl; 
      exit(0); 
    } 
}  

void QUAD::Pass(double x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  QUAD_Pass(x,L,Nint,K1L,K1SL);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void QUAD::DAPass(tps x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  QUAD_Pass(x,L,Nint,K1L,K1SL);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void QUAD::sPass(double x[9]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  QUAD_sPass(x,L,Nint,K1L,K1SL);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

//--------------------------------SKEWQ------------------------------
SKEWQ::SKEWQ(string name, double l, double k1sl): Element(name)
{ 
  if (l >=0 ) { 
    TYPE=string("SKEWQ");
    GROUP=string("");
    L=l; 
    K1SL= k1sl; 
    Nint=int(L/QLslice)+1;
    Norder=1;}
  else { 
    cout<<"Error: SKEWQ length should be non-negative."<<endl; 
    exit(1); 
  }
}

void SKEWQ::SetP(const char *name, double value)     
{
  if (strcmp( name, "K1SL" ) == 0) 	
    {  
      K1SL=value;	
    }
  else if (strcmp( name, "Nint" ) == 0) 
    {
      Nint=int(value);
    }
  else  { 
    cout<<"SKEWQ does not have a parameter of "<<name<<endl; 
    exit(0);
  } 
}

double SKEWQ::GetP(const char *name)
{
  if (strcmp( name, "K1SL" ) == 0) 	
    {
      return K1SL;	
    }
  else if (strcmp( name, "Nint" ) == 0) 
    {
      return Nint;
    }
  else  
    {
      cout<<"SKEWQ does not have a parameter of "<<name<<endl; 
      exit(0); 
    }
}  

void SKEWQ::Pass(double x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  QUAD_Pass(x,L,Nint, 0., K1SL);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void SKEWQ::DAPass(tps x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  QUAD_Pass(x,L,Nint, 0., K1SL);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void SKEWQ::sPass(double x[9]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  QUAD_sPass(x,L,Nint, 0., K1SL);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

//----------------------------SEXT-------------------------------------------
SEXT::SEXT(string name, double l, double k2l, double k2sl): Element(name)
{ 
  if (l >=0)
    {
      TYPE=string("SEXT");
      GROUP=string("");
      L=l;
      K2L=k2l;
      K2SL=k2sl;
      Nint=int(L/QLslice)+1;
      Norder=2;
    }
  else
    {
      cout<<"Error: SEXT length should be non-negative."<<endl;
      exit(1); 
    }
}

void SEXT::SetP(const char *name, double value)
{ 
  if (strcmp( name, "K2L" ) == 0) 
    {
      K2L=value;
    }
  else if (strcmp( name, "K2SL" ) == 0) 
    {
      K2SL=value;
    }
  else if (strcmp( name, "Nint" ) == 0) 
    {
      Nint=int(value);
    }
  else 
    {
      cout<<"SEXT does not have a parameter of "<<name<<endl; 
      exit(0); 
    } 
} 

double SEXT::GetP(const char *name)
{
  if (strcmp( name, "K2L" ) == 0) 
    {
      return K2L;
    }
  else if (strcmp( name, "K2SL" ) == 0) 
    {
      return K2SL;
    }
  else if (strcmp( name, "Nint" ) == 0) 
    {
      return Nint;
    }
  else 
    {
      cout<<"SEXT does not have a parameter of "<<name<<endl; 
      exit(0); 
    } 
}  

void SEXT::Pass(double x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  SEXT_Pass(x,L,Nint,K2L,K2SL);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void SEXT::DAPass(tps x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  SEXT_Pass(x,L,Nint,K2L,K2SL);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void SEXT::sPass(double x[9]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  SEXT_sPass(x,L,Nint,K2L,K2SL);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

//-------------------------OCT---------------------------------------------
OCT::OCT(string name, double l, double k3l, double k3sl): Element(name)
{ 
  if (l >=0  )
    {
      TYPE=string("OCT");
      GROUP=string("");
      L=l;
      K3L=k3l;
      K3SL=k3sl;
      Nint=int(L/QLslice)+1;
      Norder=3;
    }
  else
    {
      cout<<"Error: OCT length should be non-negative."<<endl;
      exit(1); 
    }
}

void OCT::SetP(const char *name, double value)
{ 
  if (strcmp( name, "K3L" ) == 0) 
    {
      K3L=value;
    }
  else if (strcmp( name, "K3SL" ) == 0) 
    {
      K3SL=value;
    }
  else if (strcmp( name, "Nint" ) == 0) 
    {
      Nint=int(value);
    }
  else 
    {
      cout<<"OCT does not have a parameter of "<<name<<endl; 
      exit(0); 
    } 
} 

double OCT::GetP(const char *name)
{
  if (strcmp( name, "K3L" ) == 0) 
    {
      return K3L;
    }
  else if (strcmp( name, "K3SL" ) == 0) 
    {
      return K3SL;
    }
  else if (strcmp( name, "Nint" ) == 0) 
    {
      return Nint;
    }
  else 
    {
      cout<<"OCT does not have a parameter of "<<name<<endl; 
      exit(0); 
    } 
}  

void OCT::Pass(double x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  OCT_Pass(x,L,Nint,K3L,K3SL);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void OCT::DAPass(tps x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  OCT_Pass(x,L,Nint,K3L,K3SL);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void OCT::sPass(double x[9]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  OCT_sPass(x,L,Nint,K3L,K3SL);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

//-------------------------------MULT--------------------------------------
MULT::MULT(string name, double l, double knl[11], double knsl[11]): Element(name)
{ 
  int i;
  if (l >=0  )
    {
      TYPE=string("MULT");
      GROUP=string("");
      L=l;
      for(i=0;i<11;i++) {
	KNL[i]=knl[i];  KNSL[i]=knsl[i]; }
      Nint=int(L/QLslice)+1;
      Norder=1;
      for(i=0;i<10;i++) {
	if( KNL[10-i] != 0. ||  KNSL[10-i] != 0. ){
	  Norder=10-i;
	  break;
	}
      }
    }
  else
    {
      cout<<"Error: MULT length should be non-negative."<<endl;
      exit(1); 
    }
}

void MULT::SetP(const char *name, double value)
{ 
  if (     strcmp( name, "K0L" ) == 0) 	{  KNL[0]=value; }
  else if (strcmp( name, "Norder" ) == 0)  {  Norder= int(value); }   
  else if (strcmp( name, "K1L" ) == 0) 	{  KNL[1]=value; }      
  else if (strcmp( name, "K2L" ) == 0) 	{  KNL[2]=value; }
  else if (strcmp( name, "K3L" ) == 0) 	{  KNL[3]=value; }
  else if (strcmp( name, "K4L" ) == 0) 	{  KNL[4]=value; }      
  else if (strcmp( name, "K5L" ) == 0) 	{  KNL[5]=value; }
  else if (strcmp( name, "K6L" ) == 0) 	{  KNL[6]=value; }
  else if (strcmp( name, "K7L" ) == 0) 	{  KNL[7]=value; }      
  else if (strcmp( name, "K8L" ) == 0) 	{  KNL[8]=value; }
  else if (strcmp( name, "K9L" ) == 0) 	{  KNL[9]=value; }      
  else if (strcmp( name, "K10L" ) == 0)    {  KNL[10]=value; }
  else if (strcmp( name, "K0SL" ) == 0)    {  KNSL[0]=value; }
  else if (strcmp( name, "K1SL" ) == 0)    {  KNSL[1]=value; }      
  else if (strcmp( name, "K2SL" ) == 0)    {  KNSL[2]=value; }
  else if (strcmp( name, "K3SL" ) == 0)    {  KNSL[3]=value; }
  else if (strcmp( name, "K4SL" ) == 0)    {  KNSL[4]=value; }      
  else if (strcmp( name, "K5SL" ) == 0)    {  KNSL[5]=value; }
  else if (strcmp( name, "K6SL" ) == 0)    {  KNSL[6]=value; }
  else if (strcmp( name, "K7SL" ) == 0)    {  KNSL[7]=value; }      
  else if (strcmp( name, "K8SL" ) == 0)    {  KNSL[8]=value; }
  else if (strcmp( name, "K9SL" ) == 0)    {  KNSL[9]=value; }      
  else if (strcmp( name, "K10SL" ) == 0)   {  KNSL[10]=value;}
  else if (strcmp( name, "Nint" ) == 0)    {  Nint=int(value);}
  else if (strcmp( name, "Norder" ) == 0)  {  Norder=int(value);}
  else 
    {
      cout<<"MULT does not have a parameter of  "<<name<<endl; 
      exit(0); 
    } 
}

double MULT::GetP(const char *name)
{
  if (     strcmp( name, "K0L" ) == 0)  { return   KNL[0]; }
  else if (strcmp( name, "Norder" ) == 0)  { return   Norder; }    
  else if (strcmp( name, "K1L" ) == 0)  { return   KNL[1]; }      
  else if (strcmp( name, "K2L" ) == 0)  { return   KNL[2]; }
  else if (strcmp( name, "K3L" ) == 0)  { return   KNL[3]; }
  else if (strcmp( name, "K4L" ) == 0)  { return   KNL[4]; }      
  else if (strcmp( name, "K5L" ) == 0)  { return   KNL[5]; }
  else if (strcmp( name, "K6L" ) == 0)  { return   KNL[6]; }
  else if (strcmp( name, "K7L" ) == 0)  { return   KNL[7]; }      
  else if (strcmp( name, "K8L" ) == 0)  { return   KNL[8]; }
  else if (strcmp( name, "K9L" ) == 0)  { return   KNL[9]; }      
  else if (strcmp( name, "K10L" ) == 0) { return   KNL[10]; }
  else if (strcmp( name, "K0SL" ) == 0) { return   KNSL[0]; }
  else if (strcmp( name, "K1SL" ) == 0) { return   KNSL[1]; }      
  else if (strcmp( name, "K2SL" ) == 0) { return   KNSL[2]; }
  else if (strcmp( name, "K3SL" ) == 0) { return   KNSL[3]; }
  else if (strcmp( name, "K4SL" ) == 0) { return   KNSL[4]; }      
  else if (strcmp( name, "K5SL" ) == 0) { return   KNSL[5]; }
  else if (strcmp( name, "K6SL" ) == 0) { return   KNSL[6]; }
  else if (strcmp( name, "K7SL" ) == 0) { return   KNSL[7]; }      
  else if (strcmp( name, "K8SL" ) == 0) { return   KNSL[8]; }
  else if (strcmp( name, "K9SL" ) == 0) { return   KNSL[9]; }      
  else if (strcmp( name, "K10SL" ) == 0){ return   KNSL[10];}
  else if (strcmp( name, "Norder" ) == 0){ return  Norder;}
  else if (strcmp( name, "Nint" ) == 0)  { return  Nint;}
  else 
    {
      cout<<"MULT does not have a parameter of  "<<name<<endl; 
      exit(0); 
    }   
}

void MULT::Pass(double x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  MULT_Pass(x,L,Nint, Norder,KNL,KNSL);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void MULT::DAPass(tps x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  MULT_Pass(x,L,Nint, Norder,KNL,KNSL);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void MULT::sPass(double x[9]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  MULT_sPass(x,L,Nint, Norder,KNL,KNSL);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

//-------------------------------GMULT--------------------------------------
GMULT::GMULT(string name, double l, double angle, double e1, double e2, double knl[11], double knsl[11]): Element(name)
{ 
  int i;
  if (l >=0  )
    {
      TYPE=string("GMULT");
      GROUP=string("");
      L=l;
      ANGLE=angle;
      E1=e1;
      E2=e2;
      Nint=int(L/BLslice)+1;
      for(i=0;i<11;i++) {
	KNL[i]=knl[i];  KNSL[i]=knsl[i]; }
      Norder=1;
      for(i=0;i<10;i++) {
	if( KNL[10-i] != 0. ||  KNSL[10-i] != 0. ){
	  Norder=10-i;
	  break;
	}
      }
    }
  else
    {
      cout<<"Error: GMULT length should be non-negative."<<endl;
      exit(1); 
    }
}

void GMULT::SetP(const char *name, double value)
{ 
  if (strcmp( name, "ANGLE" ) == 0) 
    {
      ANGLE=value;
    }
  else if (strcmp( name, "E1" ) == 0) 
    {
      E1=value;
    }
  else if (strcmp( name, "E2" ) == 0) 
    {
      E2=value;
    }
  else if (strcmp( name, "Nint" ) == 0) 
	 {
	   Nint=int(value);
	 }
  else if (strcmp( name, "K0L" ) == 0) 	{  KNL[0]=value; }
  else if (strcmp( name, "K1L" ) == 0) 	{  KNL[1]=value; }      
  else if (strcmp( name, "K2L" ) == 0) 	{  KNL[2]=value; }
  else if (strcmp( name, "K3L" ) == 0) 	{  KNL[3]=value; }
  else if (strcmp( name, "K4L" ) == 0) 	{  KNL[4]=value; }      
  else if (strcmp( name, "K5L" ) == 0) 	{  KNL[5]=value; }
  else if (strcmp( name, "K6L" ) == 0) 	{  KNL[6]=value; }
  else if (strcmp( name, "K7L" ) == 0) 	{  KNL[7]=value; }      
  else if (strcmp( name, "K8L" ) == 0) 	{  KNL[8]=value; }
  else if (strcmp( name, "K9L" ) == 0) 	{  KNL[9]=value; }      
  else if (strcmp( name, "K10L" ) == 0)    {  KNL[10]=value; }
  else if (strcmp( name, "K0SL" ) == 0)    {  KNSL[0]=value; }
  else if (strcmp( name, "K1SL" ) == 0)    {  KNSL[1]=value; }      
  else if (strcmp( name, "K2SL" ) == 0)    {  KNSL[2]=value; }
  else if (strcmp( name, "K3SL" ) == 0)    {  KNSL[3]=value; }
  else if (strcmp( name, "K4SL" ) == 0)    {  KNSL[4]=value; }      
  else if (strcmp( name, "K5SL" ) == 0)    {  KNSL[5]=value; }
  else if (strcmp( name, "K6SL" ) == 0)    {  KNSL[6]=value; }
  else if (strcmp( name, "K7SL" ) == 0)    {  KNSL[7]=value; }      
  else if (strcmp( name, "K8SL" ) == 0)    {  KNSL[8]=value; }
  else if (strcmp( name, "K9SL" ) == 0)    {  KNSL[9]=value; }      
  else if (strcmp( name, "K10SL" ) == 0)   {  KNSL[10]=value;}
  else if (strcmp( name, "Nint" ) == 0)    {  Nint=int(value);}
  else 
    {
      cout<<"GMULT does not have a parameter of  "<<name<<endl; 
      exit(0); 
    } 
}

double GMULT::GetP(const char *name)
{
  if (strcmp( name, "ANGLE" ) == 0) 
    {
      return ANGLE;
    }
  else if (strcmp( name, "E1" ) == 0) 
    {
      return E1;
    }
  else if (strcmp( name, "E2" ) == 0) 
    {
      return E2;
    }
  else if (strcmp( name, "Nint" ) == 0) 
    {
      return Nint;
    }
  else if (strcmp( name, "K0L" ) == 0)  { return   KNL[0]; }
  else if (strcmp( name, "K1L" ) == 0)  { return   KNL[1]; }      
  else if (strcmp( name, "K2L" ) == 0)  { return   KNL[2]; }
  else if (strcmp( name, "K3L" ) == 0)  { return   KNL[3]; }
  else if (strcmp( name, "K4L" ) == 0)  { return   KNL[4]; }      
  else if (strcmp( name, "K5L" ) == 0)  { return   KNL[5]; }
  else if (strcmp( name, "K6L" ) == 0)  { return   KNL[6]; }
  else if (strcmp( name, "K7L" ) == 0)  { return   KNL[7]; }      
  else if (strcmp( name, "K8L" ) == 0)  { return   KNL[8]; }
  else if (strcmp( name, "K9L" ) == 0)  { return   KNL[9]; }      
  else if (strcmp( name, "K10L" ) == 0) { return   KNL[10]; }
  else if (strcmp( name, "K0SL" ) == 0) { return   KNSL[0]; }
  else if (strcmp( name, "K1SL" ) == 0) { return   KNSL[1]; }      
  else if (strcmp( name, "K2SL" ) == 0) { return   KNSL[2]; }
  else if (strcmp( name, "K3SL" ) == 0) { return   KNSL[3]; }
  else if (strcmp( name, "K4SL" ) == 0) { return   KNSL[4]; }      
  else if (strcmp( name, "K5SL" ) == 0) { return   KNSL[5]; }
  else if (strcmp( name, "K6SL" ) == 0) { return   KNSL[6]; }
  else if (strcmp( name, "K7SL" ) == 0) { return   KNSL[7]; }      
  else if (strcmp( name, "K8SL" ) == 0) { return   KNSL[8]; }
  else if (strcmp( name, "K9SL" ) == 0) { return   KNSL[9]; }      
  else if (strcmp( name, "K10SL" ) == 0){ return   KNSL[10];}
  else if (strcmp( name, "Norder" ) == 0){ return  Norder;}
  else if (strcmp( name, "Nint" ) == 0)  { return  Nint;}
  else 
    {
      cout<<"GMULT does not have a parameter of  "<<name<<endl; 
      exit(0); 
    }   
}

void GMULT::Pass(double x[6]){
  GtoL(x,L,ANGLE,DX,DY,DS,DTHETA,DPHI,DPSI);
  GMULT_Pass(x,L,Nint,Norder,ANGLE,E1,E2,KNL,KNSL);
  LtoG(x,L,ANGLE,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void GMULT::DAPass(tps x[6]){
   GtoL(x,L,ANGLE,DX,DY,DS,DTHETA,DPHI,DPSI);
   GMULT_Pass(x,L,Nint,Norder,ANGLE,E1,E2,KNL,KNSL);
   LtoG(x,L,ANGLE,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void GMULT::sPass(double x[9]){
  GtoL(x,L,ANGLE,DX,DY,DS,DTHETA,DPHI,DPSI);
  GMULT_sPass(x,L,Nint,Norder,ANGLE,E1,E2,KNL,KNSL);
  LtoG(x,L,ANGLE,DX,DY,DS,DTHETA,DPHI,DPSI);
}


//-------------------------------SBENDMULT: same as GMULT--------------------------------------
SBENDMULT::SBENDMULT(string name, double l, double angle, double e1, double e2, double knl[11], double knsl[11]): Element(name)
{ 
  int i;
  if (l >=0  )
    {
      TYPE=string("SBENDMULT");
      GROUP=string("");
      L=l;
      ANGLE=angle;
      E1=e1;
      E2=e2;
      Nint=int(L/BLslice)+1;
      for(i=0;i<11;i++) {
	KNL[i]=knl[i];  KNSL[i]=knsl[i]; }
      Norder=1;
      for(i=0;i<10;i++) {
	if( KNL[10-i] != 0. ||  KNSL[10-i] != 0. ){
	  Norder=10-i;
	  break;
	}
      }
    }
  else
    {
      cout<<"Error: SBENDMULT length should be non-negative."<<endl;
      exit(1); 
    }
}

void SBENDMULT::SetP(const char *name, double value)
{ 
  if (strcmp( name, "ANGLE" ) == 0) 
    {
      ANGLE=value;
    }
  else if (strcmp( name, "E1" ) == 0) 
    {
      E1=value;
    }
  else if (strcmp( name, "E2" ) == 0) 
    {
      E2=value;
    }
  else if (strcmp( name, "K0L" ) == 0) 	{  KNL[0]=value; }
  else if (strcmp( name, "K1L" ) == 0) 	{  KNL[1]=value; }      
  else if (strcmp( name, "K2L" ) == 0) 	{  KNL[2]=value; }
  else if (strcmp( name, "K3L" ) == 0) 	{  KNL[3]=value; }
  else if (strcmp( name, "K4L" ) == 0) 	{  KNL[4]=value; }      
  else if (strcmp( name, "K5L" ) == 0) 	{  KNL[5]=value; }
  else if (strcmp( name, "K6L" ) == 0) 	{  KNL[6]=value; }
  else if (strcmp( name, "K7L" ) == 0) 	{  KNL[7]=value; }      
  else if (strcmp( name, "K8L" ) == 0) 	{  KNL[8]=value; }
  else if (strcmp( name, "K9L" ) == 0) 	{  KNL[9]=value; }      
  else if (strcmp( name, "K10L" ) == 0)    {  KNL[10]=value; }
  else if (strcmp( name, "K0SL" ) == 0)    {  KNSL[0]=value; }
  else if (strcmp( name, "K1SL" ) == 0)    {  KNSL[1]=value; }      
  else if (strcmp( name, "K2SL" ) == 0)    {  KNSL[2]=value; }
  else if (strcmp( name, "K3SL" ) == 0)    {  KNSL[3]=value; }
  else if (strcmp( name, "K4SL" ) == 0)    {  KNSL[4]=value; }      
  else if (strcmp( name, "K5SL" ) == 0)    {  KNSL[5]=value; }
  else if (strcmp( name, "K6SL" ) == 0)    {  KNSL[6]=value; }
  else if (strcmp( name, "K7SL" ) == 0)    {  KNSL[7]=value; }      
  else if (strcmp( name, "K8SL" ) == 0)    {  KNSL[8]=value; }
  else if (strcmp( name, "K9SL" ) == 0)    {  KNSL[9]=value; }      
  else if (strcmp( name, "K10SL" ) == 0)   {  KNSL[10]=value;}
  else if (strcmp( name, "Norder" ) == 0)  {  Norder=int(value);}
  else if (strcmp( name, "Nint" ) == 0)    {  Nint=int(value);}
  else 
    {
      cout<<"SBENDMULT does not have a parameter of  "<<name<<endl; 
      exit(0); 
    } 
}

double SBENDMULT::GetP(const char *name)
{
  if (strcmp( name, "ANGLE" ) == 0) 
    {
      return ANGLE;
    }
  else if (strcmp( name, "E1" ) == 0) 
    {
      return E1;
    }
  else if (strcmp( name, "E2" ) == 0) 
    {
      return E2;
    }
  else if (strcmp( name, "K0L" ) == 0)  { return   KNL[0]; }
  else if (strcmp( name, "K1L" ) == 0)  { return   KNL[1]; }      
  else if (strcmp( name, "K2L" ) == 0)  { return   KNL[2]; }
  else if (strcmp( name, "K3L" ) == 0)  { return   KNL[3]; }
  else if (strcmp( name, "K4L" ) == 0)  { return   KNL[4]; }      
  else if (strcmp( name, "K5L" ) == 0)  { return   KNL[5]; }
  else if (strcmp( name, "K6L" ) == 0)  { return   KNL[6]; }
  else if (strcmp( name, "K7L" ) == 0)  { return   KNL[7]; }      
  else if (strcmp( name, "K8L" ) == 0)  { return   KNL[8]; }
  else if (strcmp( name, "K9L" ) == 0)  { return   KNL[9]; }      
  else if (strcmp( name, "K10L" ) == 0) { return   KNL[10]; }
  else if (strcmp( name, "K0SL" ) == 0) { return   KNSL[0]; }
  else if (strcmp( name, "K1SL" ) == 0) { return   KNSL[1]; }      
  else if (strcmp( name, "K2SL" ) == 0) { return   KNSL[2]; }
  else if (strcmp( name, "K3SL" ) == 0) { return   KNSL[3]; }
  else if (strcmp( name, "K4SL" ) == 0) { return   KNSL[4]; }      
  else if (strcmp( name, "K5SL" ) == 0) { return   KNSL[5]; }
  else if (strcmp( name, "K6SL" ) == 0) { return   KNSL[6]; }
  else if (strcmp( name, "K7SL" ) == 0) { return   KNSL[7]; }      
  else if (strcmp( name, "K8SL" ) == 0) { return   KNSL[8]; }
  else if (strcmp( name, "K9SL" ) == 0) { return   KNSL[9]; }      
  else if (strcmp( name, "K10SL" ) == 0){ return   KNSL[10];}
  else if (strcmp( name, "Norder" ) == 0){ return  Norder;}
  else if (strcmp( name, "Nint" ) == 0)  { return  Nint;}
  else 
    {
      cout<<"SBENDMULT does not have a parameter of  "<<name<<endl; 
      exit(0); 
    }   
}

void SBENDMULT::Pass(double x[6]){
  GtoL(x,L,ANGLE,DX,DY,DS,DTHETA,DPHI,DPSI);
  SBENDMULT_Pass(x,L,Nint,Norder,ANGLE,E1,E2,KNL,KNSL);
  LtoG(x,L,ANGLE,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void SBENDMULT::DAPass(tps x[6]){
   GtoL(x,L,ANGLE,DX,DY,DS,DTHETA,DPHI,DPSI);
   SBENDMULT_Pass(x,L,Nint,Norder,ANGLE,E1,E2,KNL,KNSL);
   LtoG(x,L,ANGLE,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void SBENDMULT::sPass(double x[9]){
  GtoL(x,L,ANGLE,DX,DY,DS,DTHETA,DPHI,DPSI);
  SBENDMULT_sPass(x,L,Nint,Norder,ANGLE,E1,E2,KNL,KNSL);
  LtoG(x,L,ANGLE,DX,DY,DS,DTHETA,DPHI,DPSI);
}

//-------------------------------SMULT--------------------------------------
SMULT::SMULT(string name, double l, double angle, double e1, double e2, double knl[11], double knsl[11]): Element(name)
{ 
  int i;
  if (l >=0  )
    {
      TYPE=string("SMULT");
      GROUP=string("");
      L=l;
      ANGLE=angle;
      E1=e1;
      E2=e2;
      Nint=int(L/GLslice)+1;
      for(i=0;i<11;i++) {
	KNL[i]=knl[i];  KNSL[i]=knsl[i]; }
      Norder=1;
      for(i=0;i<10;i++) {
	if( KNL[10-i] != 0. ||  KNSL[10-i] != 0. ){
	  Norder=10-i;
	  break;
	}
      }
    }
  else
    {
      cout<<"Error: SMULT length should be non-negative."<<endl;
      exit(1); 
    }
}

void SMULT::SetP(const char *name, double value)
{ 
  if (strcmp( name, "ANGLE" ) == 0) 
    {
      ANGLE=value;
    }
  else if (strcmp( name, "E1" ) == 0) 
    {
      E1=value;
    }
  else if (strcmp( name, "E2" ) == 0) 
    {
      E2=value;
    }
  else if (strcmp( name, "Nint" ) == 0) 
    {
      Nint=int(value);
    }
   else if (strcmp( name, "Noder" ) == 0) 
    {
      Nint=int(Norder);
    }  
  else if (strcmp( name, "K0L" ) == 0) 	{  KNL[0]=value; }
  else if (strcmp( name, "K1L" ) == 0) 	{  KNL[1]=value; }      
  else if (strcmp( name, "K2L" ) == 0) 	{  KNL[2]=value; }
  else if (strcmp( name, "K3L" ) == 0) 	{  KNL[3]=value; }
  else if (strcmp( name, "K4L" ) == 0) 	{  KNL[4]=value; }      
  else if (strcmp( name, "K5L" ) == 0) 	{  KNL[5]=value; }
  else if (strcmp( name, "K6L" ) == 0) 	{  KNL[6]=value; }
  else if (strcmp( name, "K7L" ) == 0) 	{  KNL[7]=value; }      
  else if (strcmp( name, "K8L" ) == 0) 	{  KNL[8]=value; }
  else if (strcmp( name, "K9L" ) == 0) 	{  KNL[9]=value; }      
  else if (strcmp( name, "K10L" ) == 0)    {  KNL[10]=value; }
  else if (strcmp( name, "K0SL" ) == 0)    {  KNSL[0]=value; }
  else if (strcmp( name, "K1SL" ) == 0)    {  KNSL[1]=value; }      
  else if (strcmp( name, "K2SL" ) == 0)    {  KNSL[2]=value; }
  else if (strcmp( name, "K3SL" ) == 0)    {  KNSL[3]=value; }
  else if (strcmp( name, "K4SL" ) == 0)    {  KNSL[4]=value; }      
  else if (strcmp( name, "K5SL" ) == 0)    {  KNSL[5]=value; }
  else if (strcmp( name, "K6SL" ) == 0)    {  KNSL[6]=value; }
  else if (strcmp( name, "K7SL" ) == 0)    {  KNSL[7]=value; }      
  else if (strcmp( name, "K8SL" ) == 0)    {  KNSL[8]=value; }
  else if (strcmp( name, "K9SL" ) == 0)    {  KNSL[9]=value; }      
  else if (strcmp( name, "K10SL" ) == 0)   {  KNSL[10]=value;}
  else if (strcmp( name, "Nint" ) == 0)    {  Nint=int(value);}
  else 
    {
      cout<<"SMULT does not have a parameter of  "<<name<<endl; 
      exit(0); 
    } 
}

double SMULT::GetP(const char *name)
{
  if (strcmp( name, "ANGLE" ) == 0) 
    {
      return ANGLE;
    }
  else if (strcmp( name, "E1" ) == 0) 
    {
      return E1;
    }
  else if (strcmp( name, "E2" ) == 0) 
    {
      return E2;
    }
  else if (strcmp( name, "Nint" ) == 0) 
    {
      return Nint;
    }
  else if (strcmp( name, "Norder" ) == 0) 
    {
      return Norder;
    } 
  else if (strcmp( name, "K0L" ) == 0)  { return   KNL[0]; }
  else if (strcmp( name, "K1L" ) == 0)  { return   KNL[1]; }      
  else if (strcmp( name, "K2L" ) == 0)  { return   KNL[2]; }
  else if (strcmp( name, "K3L" ) == 0)  { return   KNL[3]; }
  else if (strcmp( name, "K4L" ) == 0)  { return   KNL[4]; }      
  else if (strcmp( name, "K5L" ) == 0)  { return   KNL[5]; }
  else if (strcmp( name, "K6L" ) == 0)  { return   KNL[6]; }
  else if (strcmp( name, "K7L" ) == 0)  { return   KNL[7]; }      
  else if (strcmp( name, "K8L" ) == 0)  { return   KNL[8]; }
  else if (strcmp( name, "K9L" ) == 0)  { return   KNL[9]; }      
  else if (strcmp( name, "K10L" ) == 0) { return   KNL[10]; }
  else if (strcmp( name, "K0SL" ) == 0) { return   KNSL[0]; }
  else if (strcmp( name, "K1SL" ) == 0) { return   KNSL[1]; }      
  else if (strcmp( name, "K2SL" ) == 0) { return   KNSL[2]; }
  else if (strcmp( name, "K3SL" ) == 0) { return   KNSL[3]; }
  else if (strcmp( name, "K4SL" ) == 0) { return   KNSL[4]; }      
  else if (strcmp( name, "K5SL" ) == 0) { return   KNSL[5]; }
  else if (strcmp( name, "K6SL" ) == 0) { return   KNSL[6]; }
  else if (strcmp( name, "K7SL" ) == 0) { return   KNSL[7]; }      
  else if (strcmp( name, "K8SL" ) == 0) { return   KNSL[8]; }
  else if (strcmp( name, "K9SL" ) == 0) { return   KNSL[9]; }      
  else if (strcmp( name, "K10SL" ) == 0){ return   KNSL[10];}
  else if (strcmp( name, "Norder" ) == 0){ return  Norder;}
  else if (strcmp( name, "Nint" ) == 0)  { return  Nint;}
  else 
    {
      cout<<"SMULT does not have a parameter of  "<<name<<endl; 
      exit(0); 
    }   
}

void SMULT::Pass(double x[6]){
  GtoL(x,L,ANGLE,DX,DY,DS,DTHETA,DPHI,DPSI);
  SMULT_Pass(x,L,Nint,Norder,ANGLE,E1,E2,KNL,KNSL);
  LtoG(x,L,ANGLE,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void SMULT::DAPass(tps x[6]){
   GtoL(x,L,ANGLE,DX,DY,DS,DTHETA,DPHI,DPSI);
   SMULT_Pass(x,L,Nint,Norder,ANGLE,E1,E2,KNL,KNSL);
   LtoG(x,L,ANGLE,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void SMULT::sPass(double x[9]){
  GtoL(x,L,ANGLE,DX,DY,DS,DTHETA,DPHI,DPSI);
  SMULT_sPass(x,L,Nint,Norder,ANGLE,E1,E2,KNL,KNSL);
  LtoG(x,L,ANGLE,DX,DY,DS,DTHETA,DPHI,DPSI);
}

//-------------------------------GSBENDMULT--------------------------------------
GSBENDMULT::GSBENDMULT(string name, double l, double angle, double k0l, double e1, double e2, double knl[11], double knsl[11]): Element(name)
{ 
  int i;
  if (l >=0  )
    {
      TYPE=string("GSBENDMULT");
      GROUP=string("");
      L=l;
      ANGLE=angle;
      K0L=k0l;
      E1=e1;
      E2=e2;
      Nint=int(L/BLslice)+1;
      for(i=0;i<11;i++) {
	KNL[i]=knl[i];  KNSL[i]=knsl[i]; }
      Norder=1;
      for(i=0;i<10;i++) {
	if( KNL[10-i] != 0. ||  KNSL[10-i] != 0. ){
	  Norder=10-i;
	  break;
	}
      }
    }
  else
    {
      cout<<"Error: GSBENDMULT length should be non-negative."<<endl;
      exit(1); 
    }
}

void GSBENDMULT::SetP(const char *name, double value)
{ 
  if (strcmp( name, "ANGLE" ) == 0) 
    {
      ANGLE=value;
    }
  else if (strcmp( name, "K0L" ) == 0) 
    {
      K0L=value;
    }   
  else if (strcmp( name, "E1" ) == 0) 
    {
      E1=value;
    }
  else if (strcmp( name, "E2" ) == 0) 
    {
      E2=value;
    }
  else if (strcmp( name, "Nint" ) == 0) 
    {
      Nint=int(value);
    }
  else if (strcmp( name, "Noder" ) == 0) 
    {
      Nint=int(Norder);
    } 
  else if (strcmp( name, "K0L" ) == 0) 	{  KNL[0]=value; }
  else if (strcmp( name, "K1L" ) == 0) 	{  KNL[1]=value; }      
  else if (strcmp( name, "K2L" ) == 0) 	{  KNL[2]=value; }
  else if (strcmp( name, "K3L" ) == 0) 	{  KNL[3]=value; }
  else if (strcmp( name, "K4L" ) == 0) 	{  KNL[4]=value; }      
  else if (strcmp( name, "K5L" ) == 0) 	{  KNL[5]=value; }
  else if (strcmp( name, "K6L" ) == 0) 	{  KNL[6]=value; }
  else if (strcmp( name, "K7L" ) == 0) 	{  KNL[7]=value; }      
  else if (strcmp( name, "K8L" ) == 0) 	{  KNL[8]=value; }
  else if (strcmp( name, "K9L" ) == 0) 	{  KNL[9]=value; }      
  else if (strcmp( name, "K10L" ) == 0)    {  KNL[10]=value; }
  else if (strcmp( name, "K0SL" ) == 0)    {  KNSL[0]=value; }
  else if (strcmp( name, "K1SL" ) == 0)    {  KNSL[1]=value; }      
  else if (strcmp( name, "K2SL" ) == 0)    {  KNSL[2]=value; }
  else if (strcmp( name, "K3SL" ) == 0)    {  KNSL[3]=value; }
  else if (strcmp( name, "K4SL" ) == 0)    {  KNSL[4]=value; }      
  else if (strcmp( name, "K5SL" ) == 0)    {  KNSL[5]=value; }
  else if (strcmp( name, "K6SL" ) == 0)    {  KNSL[6]=value; }
  else if (strcmp( name, "K7SL" ) == 0)    {  KNSL[7]=value; }      
  else if (strcmp( name, "K8SL" ) == 0)    {  KNSL[8]=value; }
  else if (strcmp( name, "K9SL" ) == 0)    {  KNSL[9]=value; }      
  else if (strcmp( name, "K10SL" ) == 0)   {  KNSL[10]=value;}
  else if (strcmp( name, "Nint" ) == 0)    {  Nint=int(value);}
  else 
    {
      cout<<"GSBENDMULT does not have a parameter of  "<<name<<endl; 
      exit(0); 
    } 
}

double GSBENDMULT::GetP(const char *name)
{
  if (strcmp( name, "ANGLE" ) == 0) 
    {
      return ANGLE;
    }
  else if (strcmp( name, "K0L" ) == 0) 
    {
      return K0L;
    }  
  else if (strcmp( name, "E1" ) == 0) 
    {
      return E1;
    }
  else if (strcmp( name, "E2" ) == 0) 
    {
      return E2;
    }
  else if (strcmp( name, "Nint" ) == 0) 
    {
      return Nint;
    }
  else if (strcmp( name, "Norder" ) == 0) 
    {
      return Norder;
    } 
  else if (strcmp( name, "K0L" ) == 0)  { return   KNL[0]; }
  else if (strcmp( name, "K1L" ) == 0)  { return   KNL[1]; }      
  else if (strcmp( name, "K2L" ) == 0)  { return   KNL[2]; }
  else if (strcmp( name, "K3L" ) == 0)  { return   KNL[3]; }
  else if (strcmp( name, "K4L" ) == 0)  { return   KNL[4]; }      
  else if (strcmp( name, "K5L" ) == 0)  { return   KNL[5]; }
  else if (strcmp( name, "K6L" ) == 0)  { return   KNL[6]; }
  else if (strcmp( name, "K7L" ) == 0)  { return   KNL[7]; }      
  else if (strcmp( name, "K8L" ) == 0)  { return   KNL[8]; }
  else if (strcmp( name, "K9L" ) == 0)  { return   KNL[9]; }      
  else if (strcmp( name, "K10L" ) == 0) { return   KNL[10]; }
  else if (strcmp( name, "K0SL" ) == 0) { return   KNSL[0]; }
  else if (strcmp( name, "K1SL" ) == 0) { return   KNSL[1]; }      
  else if (strcmp( name, "K2SL" ) == 0) { return   KNSL[2]; }
  else if (strcmp( name, "K3SL" ) == 0) { return   KNSL[3]; }
  else if (strcmp( name, "K4SL" ) == 0) { return   KNSL[4]; }      
  else if (strcmp( name, "K5SL" ) == 0) { return   KNSL[5]; }
  else if (strcmp( name, "K6SL" ) == 0) { return   KNSL[6]; }
  else if (strcmp( name, "K7SL" ) == 0) { return   KNSL[7]; }      
  else if (strcmp( name, "K8SL" ) == 0) { return   KNSL[8]; }
  else if (strcmp( name, "K9SL" ) == 0) { return   KNSL[9]; }      
  else if (strcmp( name, "K10SL" ) == 0){ return   KNSL[10];}
  else if (strcmp( name, "Norder" ) == 0){ return  Norder;}
  else if (strcmp( name, "Nint" ) == 0)  { return  Nint;}
  else 
    {
      cout<<"GSBENDMULT does not have a parameter of  "<<name<<endl; 
      exit(0); 
    }   
}

void GSBENDMULT::Pass(double x[6]){
  GtoL(x,L,ANGLE,DX,DY,DS,DTHETA,DPHI,DPSI);
  GSBENDMULT_Pass(x,L,Nint,Norder,ANGLE,K0L,E1,E2,KNL,KNSL);
  LtoG(x,L,ANGLE,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void GSBENDMULT::DAPass(tps x[6]){
   GtoL(x,L,ANGLE,DX,DY,DS,DTHETA,DPHI,DPSI);
   GSBENDMULT_Pass(x,L,Nint,Norder,ANGLE,K0L,E1,E2,KNL,KNSL);
   LtoG(x,L,ANGLE,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void GSBENDMULT::sPass(double x[9]){
  GtoL(x,L,ANGLE,DX,DY,DS,DTHETA,DPHI,DPSI);
  GSBENDMULT_sPass(x,L,Nint,Norder,ANGLE,K0L,E1,E2,KNL,KNSL);
  LtoG(x,L,ANGLE,DX,DY,DS,DTHETA,DPHI,DPSI);
}

//----------------------------SOLEN-----------------------------------
SOLEN::SOLEN(string name, double l, double ks): Element(name)
    { 
      if ( l >0 )  // if ( l >0 && ks !=0.)
	{
	  TYPE=string("SOLEN");
	  GROUP=string("");
	  L=l;
          KS=ks;       //  to be consistent with madx, KS=Bz/GP.brho
          Nint=1;      //  default, can be changed
	}
      else
	{
	  cout<<"Solenoid: L must be positive. "<<endl;
          exit(1); 
	}
    }

void SOLEN::SetP(const char *name, double value) 
    {
      if (strcmp( name, "KS" ) == 0) 
	{
	  KS=value;
	}
      else 
	{
	  cout<<"SOLEN does not have a parameter of "<<name<<endl; 
	  exit(0); 
       } 
    }

double SOLEN::GetP(const char *name) 
    {
     if (strcmp( name, "KS" ) == 0) 
	{
	  return KS;
	}
      else 
	{
	  cout<<"SOLEN does not have a parameter of "<<name<<endl; 
          exit(0); 
	} 
    }

void SOLEN::Pass(double x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  SOLEN_Pass(x, L, Nint, KS);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void SOLEN::DAPass(tps x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  SOLEN_Pass(x, L, Nint, KS);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void SOLEN::sPass(double x[9]){ 
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  SOLEN_sPass(x, L, Nint, KS);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

//----------------------------WIGGLER-----------------------------------
WIGGLER::WIGGLER(string name, double l, int nint, double b0, double kx, double kz, double phiz0) : Element(name)
    { 
      if ( l >0 )
	{
	  TYPE=string("WIGGLER");
	  GROUP=string("");
	  L=l;
	  Nint=nint;
	  B0=b0;
          KX=kx;
	  KZ=kz;
	  PHIZ0=phiz0;
	}
      else
	{
	  cout<<"Wiggler: L must be non-zero."<<endl;
          exit(1); 
	}
    }

void WIGGLER::SetP(const char *name, double value) 
    {
      if (strcmp( name, "B0" ) == 0) {    B0=value;  }
      else if (strcmp( name, "Nint" ) == 0) 	{  Nint=int(value);}      
      else if (strcmp( name, "KX" ) == 0) 	{  KX=value;}
      else if (strcmp( name, "KZ" ) == 0) 	{  KZ=value;}
      else if (strcmp( name, "PHIZ0" ) == 0) 	{  PHIZ0=value; }
      else 
	{
	  cout<<"WIGGLER does not have a parameter of "<<name<<endl; 
	  exit(0); 
	} 
    }

double  WIGGLER:: GetP(const char *name) 
    {
      if (     strcmp( name, "B0" ) == 0) 	{  return B0; }
      else if (strcmp( name, "Nint" ) == 0) 	{  return Nint;}  
      else if (strcmp( name, "KX" ) == 0) 	{  return KX; }
      else if (strcmp( name, "KZ" ) == 0) 	{  return KZ; }
      else if (strcmp( name, "PHIZ0" ) == 0) 	{  return PHIZ0; }     
      else 
	{
	  cout<<"WIGGLER does not have a parameter of  "<<name<<endl; 
	  exit(0); 
	} 
    }

void WIGGLER::Pass(double x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  Wiggler_Pass_Forest(x,L,Nint,B0,KX,KZ,PHIZ0);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void WIGGLER ::DAPass(tps x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  Wiggler_Pass_Wu(x,L,Nint,B0,KX,KZ,PHIZ0);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void WIGGLER::sPass(double x[9]){    //  need  updating spin tracking
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  Wiggler_Pass_Wu(x,L,Nint,B0,KX,KZ,PHIZ0);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

//-----------------------MATRIX-----------------------------------
MATRIX::MATRIX(string name, double l, double xco_in[6], double xco_out[6], double m66[36]): Element(name)
    {
      int i; 
      TYPE=string("MATRIX");
      GROUP=string("");
      L=l;
      for(i=0;i<36;i++) M66[i]=m66[i];
      for(i=0;i<6;i++)  XCO_IN[i]=xco_in[i];
      for(i=0;i<6;i++)  XCO_OUT[i]=xco_out[i];
    }

void MATRIX ::SetP(const char *name, double value) 
    {
      if (     strcmp( name, "M11" ) == 0) 	{  M66[0]=value; }
      else if (strcmp( name, "M12" ) == 0) 	{  M66[1]=value; }      
      else if (strcmp( name, "M13" ) == 0) 	{  M66[2]=value; }
      else if (strcmp( name, "M14" ) == 0) 	{  M66[3]=value; }
      else if (strcmp( name, "M15" ) == 0) 	{  M66[4]=value; }      
      else if (strcmp( name, "M16" ) == 0) 	{  M66[5]=value; }
      else if (strcmp( name, "M21" ) == 0) 	{  M66[6]=value; }   
      else if (strcmp( name, "M22" ) == 0) 	{  M66[7]=value; }      
      else if (strcmp( name, "M23" ) == 0) 	{  M66[8]=value; }
      else if (strcmp( name, "M24" ) == 0) 	{  M66[9]=value; }
      else if (strcmp( name, "M25" ) == 0) 	{  M66[10]=value; }      
      else if (strcmp( name, "M26" ) == 0) 	{  M66[11]=value; }
      else if (strcmp( name, "M31" ) == 0) 	{  M66[12]=value; }   
      else if (strcmp( name, "M32" ) == 0) 	{  M66[13]=value; }      
      else if (strcmp( name, "M33" ) == 0) 	{  M66[14]=value; }
      else if (strcmp( name, "M34" ) == 0) 	{  M66[15]=value; }
      else if (strcmp( name, "M35" ) == 0) 	{  M66[16]=value; }      
      else if (strcmp( name, "M36" ) == 0) 	{  M66[17]=value; }
      else if (strcmp( name, "M41" ) == 0) 	{  M66[18]=value; }   
      else if (strcmp( name, "M42" ) == 0) 	{  M66[19]=value; }      
      else if (strcmp( name, "M43" ) == 0) 	{  M66[20]=value; }
      else if (strcmp( name, "M44" ) == 0) 	{  M66[21]=value; }
      else if (strcmp( name, "M45" ) == 0) 	{  M66[22]=value; }      
      else if (strcmp( name, "M46" ) == 0) 	{  M66[23]=value; }
      else if (strcmp( name, "M51" ) == 0) 	{  M66[24]=value; }   
      else if (strcmp( name, "M52" ) == 0) 	{  M66[25]=value; }      
      else if (strcmp( name, "M53" ) == 0) 	{  M66[26]=value; }
      else if (strcmp( name, "M54" ) == 0) 	{  M66[27]=value; }
      else if (strcmp( name, "M55" ) == 0) 	{  M66[28]=value; }      
      else if (strcmp( name, "M56" ) == 0) 	{  M66[29]=value; }
      else if (strcmp( name, "M61" ) == 0) 	{  M66[30]=value; }   
      else if (strcmp( name, "M62" ) == 0) 	{  M66[31]=value; }      
      else if (strcmp( name, "M63" ) == 0) 	{  M66[32]=value; }
      else if (strcmp( name, "M64" ) == 0) 	{  M66[33]=value; }
      else if (strcmp( name, "M65" ) == 0) 	{  M66[34]=value; }      
      else if (strcmp( name, "M66" ) == 0) 	{  M66[35]=value; }
      else if (strcmp( name, "XCO_IN_X" ) == 0)     {  XCO_IN[0]=value; }
      else if (strcmp( name, "XCO_IN_PX" ) == 0)    {  XCO_IN[1]=value; }
      else if (strcmp( name, "XCO_IN_Y" ) == 0)     {  XCO_IN[2]=value; }
      else if (strcmp( name, "XCO_IN_PY" ) == 0)    {  XCO_IN[3]=value; }
      else if (strcmp( name, "XCO_IN_Z" ) == 0)     {  XCO_IN[4]=value; }
      else if (strcmp( name, "XCO_IN_DELTA" ) == 0) {  XCO_IN[5]=value; }
      else if (strcmp( name, "XCO_OUT_X" ) == 0)    {  XCO_OUT[0]=value; }
      else if (strcmp( name, "XCO_OUT_PX" ) == 0)   {  XCO_OUT[1]=value; }
      else if (strcmp( name, "XCO_OUT_Y" ) == 0)    {  XCO_OUT[2]=value; }
      else if (strcmp( name, "XCO_OUT_PY" ) == 0)   {  XCO_OUT[3]=value; }
      else if (strcmp( name, "XCO_OUT_Z" ) == 0)    {  XCO_OUT[4]=value; }
      else if (strcmp( name, "XCO_OUT_DELTA" ) == 0){  XCO_OUT[5]=value; }
      else 
	{
	  cout<<"Matrix does not have parameter of  "<<name<<endl; 
	  exit(0); 
	} 
    }

double  MATRIX:: GetP(const char *name) 
    {
      if (     strcmp( name, "M11" ) == 0) 	{  return M66[0]; }
      else if (strcmp( name, "M12" ) == 0) 	{  return M66[1]; }      
      else if (strcmp( name, "M13" ) == 0) 	{  return M66[2]; }
      else if (strcmp( name, "M14" ) == 0) 	{  return M66[3]; }
      else if (strcmp( name, "M15" ) == 0) 	{  return M66[4]; }      
      else if (strcmp( name, "M16" ) == 0) 	{  return M66[5]; }
      else if (strcmp( name, "M21" ) == 0) 	{  return M66[6]; }   
      else if (strcmp( name, "M22" ) == 0) 	{  return M66[7]; }      
      else if (strcmp( name, "M23" ) == 0) 	{  return M66[8]; }
      else if (strcmp( name, "M24" ) == 0) 	{  return M66[9]; }
      else if (strcmp( name, "M25" ) == 0) 	{  return M66[10]; }      
      else if (strcmp( name, "M26" ) == 0) 	{  return M66[11]; }
      else if (strcmp( name, "M31" ) == 0) 	{  return M66[12]; }   
      else if (strcmp( name, "M32" ) == 0) 	{  return M66[13]; }      
      else if (strcmp( name, "M33" ) == 0) 	{  return M66[14]; }
      else if (strcmp( name, "M34" ) == 0) 	{  return M66[15]; }
      else if (strcmp( name, "M35" ) == 0) 	{  return M66[16]; }      
      else if (strcmp( name, "M36" ) == 0) 	{  return M66[17]; }
      else if (strcmp( name, "M41" ) == 0) 	{  return M66[18]; }   
      else if (strcmp( name, "M42" ) == 0) 	{  return M66[19]; }      
      else if (strcmp( name, "M43" ) == 0) 	{  return M66[20]; }
      else if (strcmp( name, "M44" ) == 0) 	{  return M66[21]; }
      else if (strcmp( name, "M45" ) == 0) 	{  return M66[22]; }      
      else if (strcmp( name, "M46" ) == 0) 	{  return M66[23]; }
      else if (strcmp( name, "M51" ) == 0) 	{  return M66[24]; }   
      else if (strcmp( name, "M52" ) == 0) 	{  return M66[25]; }      
      else if (strcmp( name, "M53" ) == 0) 	{  return M66[26]; }
      else if (strcmp( name, "M54" ) == 0) 	{  return M66[27]; }
      else if (strcmp( name, "M55" ) == 0) 	{  return M66[28]; }      
      else if (strcmp( name, "M56" ) == 0) 	{  return M66[29]; }
      else if (strcmp( name, "M61" ) == 0) 	{  return M66[30]; }   
      else if (strcmp( name, "M62" ) == 0) 	{  return M66[31]; }      
      else if (strcmp( name, "M63" ) == 0) 	{  return M66[32]; }
      else if (strcmp( name, "M64" ) == 0) 	{  return M66[33]; }
      else if (strcmp( name, "M65" ) == 0) 	{  return M66[34]; }      
      else if (strcmp( name, "M66" ) == 0) 	{  return M66[35]; }
      else if (strcmp( name, "XCO_IN_X" ) == 0)     {  return XCO_IN[0]; }
      else if (strcmp( name, "XCO_IN_PX" ) == 0)    {  return XCO_IN[1]; }
      else if (strcmp( name, "XCO_IN_Y" ) == 0)     {  return XCO_IN[2]; }
      else if (strcmp( name, "XCO_IN_PY" ) == 0)    {  return XCO_IN[3]; }
      else if (strcmp( name, "XCO_IN_Z" ) == 0)     {  return XCO_IN[4]; }
      else if (strcmp( name, "XCO_IN_DELTA" ) == 0) {  return XCO_IN[5]; }
      else if (strcmp( name, "XCO_OUT_X" ) == 0)    {  return XCO_OUT[0]; }
      else if (strcmp( name, "XCO_OUT_PX" ) == 0)   {  return XCO_OUT[1]; }
      else if (strcmp( name, "XCO_OUT_Y" ) == 0)    {  return XCO_OUT[2]; }
      else if (strcmp( name, "XCO_OUT_PY" ) == 0)   {  return XCO_OUT[3]; }
      else if (strcmp( name, "XCO_OUT_Z" ) == 0)    {  return XCO_OUT[4]; }
      else if (strcmp( name, "XCO_OUT_DELTA" ) == 0){  return XCO_OUT[5]; }
      else 
	{
	  cout<<"Matrix does not have a parameter of  "<<name<<endl; 
	  exit(0); 
	} 
    }

void MATRIX::Pass(double x[6]){
  MATRIX_Pass(x,L,XCO_IN, XCO_OUT, M66);
}

void MATRIX ::DAPass(tps x[6]){
  MATRIX_Pass(x,L,XCO_IN, XCO_OUT, M66); 
}

void MATRIX::sPass(double x[9]){
  MATRIX_Pass(x,L,XCO_IN, XCO_OUT, M66);  
}

//---------------------------KICKER-----------------------------------------
KICKER::KICKER(string name, double l, double hkick, double vkick): Element(name)
    { 
      if (l >=0 )
	{
	  TYPE=string("KICKER");
	  GROUP=string("");
	  L=l;
          HKICK=hkick;
          VKICK=vkick;
	}
      else
	{
	  cout<<"Error: KICKER length should be non-negative."<<endl;
          exit(1); 
	}
    }

void KICKER::SetP(const char *name, double value)
    { 
      if (strcmp( name, "HKICK" ) == 0) 
	{
	  HKICK=value;
	}
      else if (strcmp( name, "VKICK" ) == 0) 
	{
	  VKICK=value;
	}
      else 
	{
	  cout<<"KICKRE does not have a parameter of "<<name<<endl; 
          exit(0); 
	} 
    } 

double KICKER::GetP(const char *name)
    {
      if (strcmp( name, "HKICK" ) == 0) 
	{
	  return HKICK;
	}
      else if (strcmp( name, "VKICK" ) == 0) 
	{
	  return VKICK;
	}
      else 
	{
	  cout<<"KICKRE does not have a parameter of "<<name<<endl; 
          exit(0); 
	} 
    }  

void KICKER::Pass(double x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  KICK_Pass(x,L,HKICK, VKICK);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void KICKER::DAPass(tps x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  KICK_Pass(x,L,HKICK, VKICK);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void KICKER::sPass(double x[9]){       
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  KICK_sPass(x,L,HKICK, VKICK);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

//----------------------------HKICKER--------------------------------------
HKICKER::HKICKER(string name, double l, double hkick): Element(name)
    { 
      if (l >=0  )
	{
	  TYPE=string("HKICKER");
	  GROUP=string("");
	  L=l;
          HKICK=hkick;
	}
      else
	{
	  cout<<"Error: HKICKER length should be non-negative."<<endl;
          exit(1); 
	}
    }

void HKICKER::SetP(const char *name, double value)
    { 
      if (strcmp( name, "HKICK" ) == 0) 
	{
	  HKICK=value;
	}
      else 
	{
	  cout<<"HKICKER does not have a parameter of "<<name<<endl; 
          exit(0); 
	} 
    } 

double HKICKER::GetP(const char *name)
    {
      if (strcmp( name, "HKICK" ) == 0) 
	{
	  return HKICK;
	}
      else 
	{
	  cout<<"HKICKER does not have a parameter of "<<name<<endl; 
          exit(0); 
	} 
    }  

void HKICKER::Pass(double x[6]){ 
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  KICK_Pass(x,L, HKICK, 0.);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void HKICKER::DAPass(tps x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  KICK_Pass(x,L, HKICK, 0.);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void HKICKER::sPass(double x[9]){   
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  KICK_sPass(x,L, HKICK, 0.);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

//------------------------------VKICKER------------------------------------------
VKICKER::VKICKER(string name, double l, double vkick): Element(name)
    { 
      if (l >=0 )
	{
	  TYPE=string("VKICKER");
	  GROUP=string("");
	  L=l;
          VKICK=vkick;
	}
      else
	{
	  cout<<"Error: VKICKER length should be non-negative."<<endl;
          exit(1); 
	}
    }

void VKICKER::SetP(const char *name, double value)
    { 
      if (strcmp( name, "VKICK" ) == 0) 
	{
	  VKICK=value;
	}
      else 
	{
	  cout<<"VKICKER does not have a parameter of "<<name<<endl; 
          exit(0); 
	} 
    } 

double VKICKER::GetP(const char *name)
    {
      if (strcmp( name, "VKICK" ) == 0) 
	{
	  return VKICK;
	}
      else 
	{
	  cout<<"KICKER does not have a parameter of "<<name<<endl; 
          exit(0); 
	} 
    }  

void VKICKER::Pass(double x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  KICK_Pass(x,L, 0., VKICK);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void VKICKER::DAPass(tps x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  KICK_Pass(x,L, 0., VKICK);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void VKICKER::sPass(double x[9]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  KICK_sPass(x,L, 0., VKICK);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

//----------------------------HACMULT --( horizontal single order ac multipole )-----------------
HACMULT::HACMULT(string name, double l, int norder, double klmax, int tturns, double phi0 ): Element(name)
    { 
      if (l >=0 )
	{
	  TYPE=string("HACMULT");
	  GROUP=string("");
	  L=l;
          Norder=norder;
          KLMAX=klmax;
          TTURNS=tturns;
          PHI0=phi0;
	}
      else
	{
	  cout<<"Error: HACMULT length should be non-negative."<<endl;
          exit(1); 
	}
    }

void HACMULT::SetP(const char *name, double value) 
    {
      if (strcmp( name, "KLMAX" ) == 0) 
	{
	  KLMAX=value;
	}
      else if (strcmp( name, "Norder" ) == 0) 
	{
	  Norder=int(value);
	}
      else if (strcmp( name, "TTURNS" ) == 0) 
	{
	  TTURNS=int(value);
	}
      else if (strcmp( name, "PHI0" ) == 0) 
	{
	  PHI0=value;
	}
      else 
	{
	  cout<<"HACMULT does not have a parameter of "<<name<<endl; 
          exit(0); 
	} 
    }

double HACMULT::GetP(const char *name) 
  {
    if (strcmp( name, "KLMAX" ) == 0) 
      {
	return KLMAX;
      }
    else if (strcmp( name, "Norder" ) == 0) 
      {
	return Norder;
      }
    else if (strcmp( name, "TTURNS" ) == 0) 
      {
	return TTURNS;
      }
    else if (strcmp( name, "PHI0" ) == 0) 
      {
	return PHI0;
      }
    else 
      {
	cout<<"HACMULT  does not have a parameter of "<<name<<endl; 
	exit(0); 
      } 
  }

void HACMULT::Pass(double x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  ACMULT_Pass(x, L, Norder, KLMAX, 0,TTURNS, PHI0);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void HACMULT::DAPass(tps x[6]){
  DRIFT_Pass(x, L);
}

void HACMULT::sPass(double x[9]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  ACMULT_sPass(x, L, Norder, KLMAX, 0, TTURNS, PHI0);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}


//----------------------------VACMULT --( vertical single order ac multipole )-----------------
VACMULT::VACMULT(string name, double l, int norder, double kslmax, int tturns, double phi0 ): Element(name)
    { 
      if (l >=0 )
	{
	  TYPE=string("VACMULT");
	  GROUP=string("");
	  L=l;
          Norder=norder;
          KSLMAX=kslmax;
          TTURNS=tturns;
          PHI0=phi0;
	}
      else
	{
	  cout<<"Error: VACMULT length should be non-negative."<<endl;
          exit(1); 
	}
    }

void VACMULT::SetP(const char *name, double value) 
    {
      if (strcmp( name, "KSLMAX" ) == 0) 
	{
	  KSLMAX=value;
	}
      else if (strcmp( name, "Norder" ) == 0) 
	{
	  Norder=int(value);
	}
      else if (strcmp( name, "TTURNS" ) == 0) 
	{
	  TTURNS=int(value);
	}
      else if (strcmp( name, "PHI0" ) == 0) 
	{
	  PHI0=value;
	}
      else 
	{
	  cout<<"VACMULT does not have a parameter of "<<name<<endl; 
          exit(0); 
	} 
    }

double VACMULT::GetP(const char *name) 
  {
    if (strcmp( name, "KSLMAX" ) == 0) 
      {
	return KSLMAX;
      }
    else if (strcmp( name, "Norder" ) == 0) 
      {
	return Norder;
      }
    else if (strcmp( name, "TTURNS" ) == 0) 
      {
	return TTURNS;
      }
    else if (strcmp( name, "PHI0" ) == 0) 
      {
	return PHI0;
      }
    else 
      {
	cout<<"VACMULT  does not have a parameter of "<<name<<endl; 
	exit(0); 
      } 
  }

void VACMULT::Pass(double x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  ACMULT_Pass(x, L, Norder, 0, KSLMAX,TTURNS, PHI0);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void VACMULT::DAPass(tps x[6]){
  DRIFT_Pass(x, L);
}

void VACMULT::sPass(double x[9]){   //  need updating
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  ACMULT_sPass(x, L, Norder, 0, KSLMAX, TTURNS, PHI0);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

//------------------------------HACDIP ( adiabetically ramping up horizontal kick)--------------------------------

HACDIP::HACDIP(string name, double l, double hkickmax, double nud, double phid, int turns, int turne ): Element(name)
    { 
      if (l >=0 )
	{
	  TYPE=string("HACDIP");
	  GROUP=string("");
	  L=l;
	  HKICKMAX=hkickmax;
          NUD=nud;
	  PHID=phid;
	  TURNS=turns;
          TURNE=turne;
	}
      else
	{
	  cout<<"Error: HACDIP length should be non-negative. "<<endl;
          exit(1); 
	}
    }

void HACDIP::SetP(const char *name, double value)
    { 
      if (strcmp( name, "HKICKMAX" ) == 0) 
	{
	  HKICKMAX=value;
	}
      else if (strcmp( name, "NUD" ) == 0) 
	{
	  NUD=value;
	}
      else if (strcmp( name, "PHID" ) == 0) 
	{
	  PHID=value;
	}
      else if (strcmp( name, "TURNS" ) == 0) 
	{
	  TURNS=int(value);
	}
      else if (strcmp( name, "TURNE" ) == 0) 
	{
	  TURNE=int(value);
	}

      else 
	{
	  cout<<"HACDIP does not have a parameter of "<<name<<endl; 
          exit(0); 
	} 
    } 

double HACDIP::GetP(const char *name)
    { 
      if (strcmp( name, "HKICKMAX" ) == 0) 
	{
	  return HKICKMAX;
	}
      else if (strcmp( name, "TURNS" ) == 0) 
	{
	  return TURNS;
	}
      else if (strcmp( name, "TURNE" ) == 0) 
	{
	  return TURNE;
	}
      else if (strcmp( name, "NUD" ) == 0) 
	{
	  return NUD;
	}
      else if (strcmp( name, "PHID" ) == 0) 
	{
	  return PHID;
	}
      else 
	{
	  cout<<"HACDIP does not have a parameter of "<<name<<endl; 
          exit(0); 
	} 
    } 

void HACDIP::Pass(double x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  ACDIP_Pass(x, L, HKICKMAX, 0., NUD, TURNS, TURNE, PHID);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void HACDIP::DAPass(tps x[6]){
  DRIFT_Pass(x,L);
}

void HACDIP::sPass(double x[9]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  ACDIP_sPass(x, L, HKICKMAX, 0., NUD, TURNS, TURNE, PHID);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

//------------------------------VACDIP------------------------------------------
VACDIP::VACDIP(string name, double l, double vkickmax, double nud, double phid, int turns, int turne ): Element(name)
    { 
      if (l >=0 )
	{
	  TYPE=string("VACDIP");
	  GROUP=string("");
	  L=l;
	  VKICKMAX=vkickmax;
          NUD=nud;
	  PHID=phid;
	  TURNS=turns;
          TURNE=turne;
	  
	}
      else
	{
	  cout<<"Error: VACDIP length should be no-negative."<<endl;
          exit(1); 
	}
    }

void VACDIP::SetP(const char *name, double value)
    { 
      if (strcmp( name, "VKICKMAX" ) == 0) 
	{
	  VKICKMAX=value;
	}
      else if (strcmp( name, "NUD" ) == 0) 
	{
	  NUD=value;
	}
      else if (strcmp( name, "PHID" ) == 0) 
	{
	  PHID=value;
	}
      else if (strcmp( name, "TURNS" ) == 0) 
	{
	  TURNS=int(value);
	}
      else if (strcmp( name, "TURNE" ) == 0) 
	{
	  TURNE=int(value);
	}
      else 
	{
	  cout<<"VACDIP does not have a parameter of "<<name<<endl; 
          exit(0); 
	} 
    } 

double VACDIP::GetP(const char *name)
    { 
      if (strcmp( name, "VKICKMAX" ) == 0) 
	{
	  return VKICKMAX;
	}
      else if (strcmp( name, "TURNS" ) == 0) 
	{
	  return TURNS;
	}
      else if (strcmp( name, "TURNE" ) == 0) 
	{
	  return TURNE;
	}
      else if (strcmp( name, "NUD" ) == 0) 
	{
	  return NUD;
	}
      else if (strcmp( name, "PHID" ) == 0) 
	{
	  return PHID;
	}
      else 
	{
	  cout<<"VACDIP does not have a parameter of "<<name<<endl; 
          exit(0); 
	} 
    } 

void VACDIP::Pass(double x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  ACDIP_Pass(x, L, 0.,VKICKMAX,NUD, TURNS, TURNE, PHID);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void VACDIP::DAPass(tps x[6]){
  DRIFT_Pass(x,L);
}

void VACDIP::sPass(double x[9]){     //  need updating
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  ACDIP_sPass(x, L, 0.,VKICKMAX,NUD, TURNS, TURNE, PHID);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

//---------------------------------BPM-------------------------
BPM::BPM(string name, double l): Element(name)
{ 
    if (l >=0 )
      {
	TYPE=string("BPM");
	GROUP=string("");
	L=l;
      }
    else
      {
	cout<<"Error: BPM length should be non-negative."<<endl;
	exit(1); 
      }
  }

void BPM::SetP(const char *name, double value) 
  {
    cout<<"No parameter to be set for BPM."<<endl;
    exit(1);
  }

double BPM::GetP(const char *name) 
  {
    cout<<"No parameter to be returned for BPM."<<endl;
    exit(1);
  }

void BPM::Pass(double x[6]){
  DRIFT_Pass(x, L);
}

void BPM::DAPass(tps x[6]){
  DRIFT_Pass(x, L);
}

void BPM::sPass(double x[9]){
  DRIFT_Pass(x, L);
}

//----------------------------HBPM-----------------------------------
HBPM::HBPM(string name, double l): Element(name)
    { 
      if (l >=0 )
	{
	  TYPE=string("HBPM");
	  GROUP=string("");
	  L=l;
	}
      else
	{
	  cout<<"Error: HBPM length shoulf be non-negative."<<endl;
          exit(1); 
	}
    }

void  HBPM::SetP(const char *name, double value) 
    {
      cout<<"No parameter to be set for HBPM."<<endl;
      exit(1);
    }

double  HBPM::GetP(const char *name) 
    {
      cout<<"No parameter to be returned for HBPM."<<endl;
      exit(1);
    }

void HBPM::Pass(double x[6]){
  DRIFT_Pass(x, L);
}

void HBPM::DAPass(tps x[6]){
  DRIFT_Pass(x, L);
}

void HBPM::sPass(double x[9]){
  DRIFT_Pass(x, L);
}

//---------------------------------VBPM------------------------------
VBPM::VBPM(string name, double l): Element(name)
    { 
      if (l >=0 )
	{
	  TYPE=string("VBPM");
	  GROUP=string("");
	  L=l;
	}
      else
	{
	  cout<<"Error: VBPM length should be non-negative."<<endl;
          exit(1); 
	}
    }

void VBPM::SetP(const char *name, double value) 
    {
      cout<<"No parameter to be set for VBPM."<<endl;
      exit(1);
    }

double VBPM::GetP(const char *name) 
    {
      cout<<"No parameter to be returned for VBPM."<<endl;
      exit(0);
    }

void VBPM::Pass(double x[6]){
  DRIFT_Pass(x, L);
}

void VBPM::DAPass(tps x[6]){
  DRIFT_Pass(x, L);
}

void VBPM::sPass(double x[9]){
  DRIFT_Pass(x, L);
}

//-----------------------------MARKER--------------------------
MARKER::MARKER(string name, double l): Element(name)
    { 
      if (l >=0 )
	{
	  TYPE=string("MARKER");
	  GROUP=string("");
	  L=l;
	}
      else
	{
	  cout<<"Error: MARKER length should be non-negative."<<endl;
          exit(1); 
	}
    }

void MARKER::SetP(const char *name, double value) 
    {
      cout<<"No parameter to be set for MARKER."<<endl;
      exit(1);
    }

double MARKER::GetP(const char *name) 
    {
      cout<<"No parameter to be returned for MARKER."<<endl;
      exit(1);
    }

void MARKER::Pass(double x[6]){
  DRIFT_Pass(x, L);   
}

void MARKER::DAPass(tps x[6]){
  DRIFT_Pass(x, L);   
}

void MARKER::sPass(double x[9]){
  DRIFT_Pass(x, L);   
}

//----------------------------RFCAV--------------------------------
RFCAV::RFCAV(string name, double l, double vrf, double frf, double phase0): Element(name)
    { 
      if (l >=0 )
	{
	  TYPE=string("RFCAV");
	  GROUP=string("");
	  VRF= vrf; 
	  FRF= frf;
	  PHASE0=phase0;
	  L=l;
	}
      else{
	cout<<"Error: RFCAV length must be non-negative."<<endl;
	exit(1); 
      }
    }

void RFCAV::SetP(const char *name, double value) 
    {
     if (strcmp( name, "VRF" ) == 0) 
	{
	  VRF=value;
	}
     else if (strcmp( name, "FRF" ) == 0) 
	{
	  FRF=value;
	}
     else if (strcmp( name, "PHASE0" ) == 0) 
	{
	  PHASE0=value;
	}
     else 
       {
	  cout<<"RFCAV does not have a parameter of "<<name<<endl; 
          exit(0); 
	} 
    }

double RFCAV::GetP(const char *name) 
    {
     if (strcmp( name, "VRF" ) == 0) 
	{
	  return VRF;
	}
     else if (strcmp( name, "FRF" ) == 0) 
	{
	  return FRF;
	}
     else if (strcmp( name, "PHASE0" ) == 0) 
	{
	  return PHASE0;
	}
      else 
	{
	  cout<<"RFCAV does not have a parameter of "<<name<<endl; 
          exit(0); 
	} 
    }

void RFCAV::Pass(double x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  RFCAV_Pass(x, L, VRF,FRF,PHASE0);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void RFCAV::DAPass(tps x[6]){
  if( VRF != 0.) {
    GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI); 
    DRIFT_Pass(x, L/2.);
    if(GP.gamma > GP.gammat ){
      x[pt_] =x[pt_] + abs(VRF*GP.Q)*sin(2.0*PI*FRF*x[z_]/2.99792458e8)/GP.beta/GP.energy;
    }
    else{
      x[pt_] =x[pt_] + abs(VRF*GP.Q)*sin(-2.0*PI*FRF*x[z_]/2.99792458e8)/GP.beta/GP.energy;
    }
    DRIFT_Pass(x, L/2.);
    LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI); 
  }
  else{
    DRIFT_Pass(x, L);
  }
}
void RFCAV::sPass(double x[9]){    //  need  updating spin tracking
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  RFCAV_Pass(x, L, VRF,FRF,PHASE0);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI); 
}

//----------------------------CRAB CAVITY--------------------------------
CRABRF::CRABRF(string name, double l, double vrf, double frf, double phase0): Element(name)
    { 
      if (l >=0 )
	{
	  TYPE=string("CRABRF");
	  GROUP=string("");
	  VRF= vrf; 
	  FRF= frf;
	  PHASE0=phase0;
	  L=l;
	}
      else{
	cout<<"Error: CRABRF length must be non-negative."<<endl;
	exit(1); 
      }
    }

void CRABRF::SetP(const char *name, double value) 
    {
     if (strcmp( name, "VRF" ) == 0) 
	{
	  VRF=value;
	}
     else if (strcmp( name, "FRF" ) == 0) 
	{
	  FRF=value;
	}
     else if (strcmp( name, "PHASE0" ) == 0) 
	{
	  PHASE0=value;
	}
     else 
       {
	  cout<<"CRABRF does not have a parameter of "<<name<<endl; 
          exit(0); 
	} 
    }

double CRABRF::GetP(const char *name) 
    {
     if (strcmp( name, "VRF" ) == 0) 
	{
	  return VRF;
	}
     else if (strcmp( name, "FRF" ) == 0) 
	{
	  return FRF;
	}
     else if (strcmp( name, "PHASE0" ) == 0) 
	{
	  return PHASE0;
	}
      else 
	{
	  cout<<"CRABRF does not have a parameter of "<<name<<endl; 
          exit(0); 
	} 
    }

void CRABRF::Pass(double x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  CRABRF_Pass(x, L, VRF,FRF,PHASE0);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void CRABRF::DAPass(tps x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  tps dz;
  dz= x[z_] - x[4][0];  // make  sure  there  is  no additional kick to px for map calculation
  
  if( VRF != 0.) {
    GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
    DRIFT_Pass(x, L/2.);
    
    x[px_]=x[px_]- (VRF*GP.Q)*sin(2.0*PI*FRF*dz/2.99792458e8 + PHASE0)/GP.beta/GP.energy;
    x[pt_]=x[pt_]- (VRF*GP.Q)*cos(2.0*PI*FRF*dz/2.99792458e8 + PHASE0)*x[0]*(2.0*PI*FRF/2.99792458e8)/GP.beta/GP.energy;

    //x[px_]=x[px_]- (VRF*GP.Q)*(2.0*PI*FRF*x[z_]/2.99792458e8)/GP.beta/GP.energy;
    //x[pt_]=x[pt_]- (VRF*GP.Q)*(2.0*PI*FRF*x[x_]/2.99792458e8)/GP.beta/GP.energy;
    
    DRIFT_Pass(x, L/2.);
    LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  }
  else{
    DRIFT_Pass(x, L);
  }
   LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void CRABRF::sPass(double x[9]){   //  need updating spin tracking
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  DRIFT_Pass(x, L);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

//----------------------------CRAB CAVITY  MULTI (dipole kick not included)--------------------------------
CCMULT::CCMULT(string name, double l, double vrf, double frf, double phase0, double b1, double a1, double b2, double a2, double b3, double a3, double b4, double a4): Element(name)
    { 
      if (l >=0 )
	{
	  TYPE=string("CCMULT");
	  GROUP=string("");
	  L=l;
	  VRF= vrf; 
	  FRF= frf;
	  PHASE0=phase0;

	  B1=b1; A1=a1;
	  B2=b2; A2=a2;
	  B3=b3; A3=a3;
	  B4=b4; A4=a4;	  
	}
      else{
	cout<<"Error: CCMULT length must be non-negative."<<endl;
	exit(1); 
      }
    }

void CCMULT::SetP(const char *name, double value) 
    {
     if (strcmp( name, "VRF" ) == 0) 
	{
	  VRF=value;
	}
     else if (strcmp( name, "FRF" ) == 0) 
	{
	  FRF=value;
	}
     else if (strcmp( name, "PHASE0" ) == 0) 
	{
	  PHASE0=value;
	}
     else if (strcmp( name, "B1" ) == 0) 
	{
	  B1=value;
	}
     else if (strcmp( name, "A1" ) == 0) 
	{
	  A1=value;
	}
     else if (strcmp( name, "B2" ) == 0) 
	{
	  B2=value;
	}
     else if (strcmp( name, "A2" ) == 0) 
	{
	  A2=value;
	}
     else if (strcmp( name, "B3" ) == 0) 
	{
	  B3=value;
	}
     else if (strcmp( name, "A3" ) == 0) 
       {
	 A3=value;
       }
     else if (strcmp( name, "B4" ) == 0) 
       {
	 B4=value;
       }
     else if (strcmp( name, "A4" ) == 0) 
       {
	 A4=value;
       }
     else 
       {
	  cout<<"CCMULT does not have a parameter of "<<name<<endl; 
          exit(0); 
	} 
    }

double CCMULT::GetP(const char *name) 
    {
     if (strcmp( name, "VRF" ) == 0) 
	{
	  return VRF;
	}
     else if (strcmp( name, "FRF" ) == 0) 
	{
	  return FRF;
	}
     else if (strcmp( name, "PHASE0" ) == 0) 
	{
	  return PHASE0;
	}
     else if (strcmp( name, "B1" ) == 0) 
	{
	  return B1;
	}
     else if (strcmp( name, "A1" ) == 0) 
	{
	  return A1;
	}
     else if (strcmp( name, "B2" ) == 0) 
	{
	  return B2;
	}
     else if (strcmp( name, "A2" ) == 0) 
	{
	  return A2;
	}
     else if (strcmp( name, "B3" ) == 0) 
	{
	  return B3;
	}
     else if (strcmp( name, "A3" ) == 0) 
	{
	  return A3;
	}
     else if (strcmp( name, "B4" ) == 0) 
	{
	  return B4;
	}
     else if (strcmp( name, "A4" ) == 0) 
	{
	  return A4;
	}     
      else 
	{
	  cout<<"CCMULT does not have a parameter of "<<name<<endl; 
          exit(0); 
	} 
    }

void CCMULT::Pass(double x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  CCMULT_Pass(x, L, VRF,FRF,PHASE0,B1,A1,B2,A2,B3,A3,B4,A4);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void CCMULT::DAPass(tps x[6]){  //  Here multipole errors not included
  return;
}

void CCMULT::sPass(double x[9]){   //  here CC effect on spin not  included
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  CCMULT_Pass(x, L, VRF,FRF,PHASE0,B1,A1,B2,A2,B3,A3,B4,A4);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

//----------------------------LBT: Lorentz Boost Transfer--------------------------------
LBT::LBT(string name, double angle): Element(name)
    { 
      TYPE=string("LBT");
      GROUP=string("");
      THETA= angle;
      L=0.;
    }

void LBT::SetP(const char *name, double value) 
    {
     if (strcmp( name, "THETA" ) == 0) 
	{
	  THETA=value;
	}
     else 
       {
	  cout<<"LBT does not have a parameter of "<<name<<endl; 
          exit(0); 
	} 
    }

double LBT::GetP(const char *name) 
    {
     if (strcmp( name, "THETA" ) == 0) 
	{
	  return THETA;
	}
      else 
	{
	  cout<<"LBT does not have a parameter of "<<name<<endl; 
          exit(0); 
	} 
    }

void LBT::Pass(double x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  LBT_Pass(x, THETA);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void LBT::DAPass(tps x[6]){
  return ;
}

void LBT::sPass(double x[9]){  
  LBT_Pass(x, THETA);
}

//----------------------------ILBT: Inverse Lorentz Boost Transfer----------------------------
ILBT::ILBT(string name, double angle): Element(name)
    { 
      TYPE=string("ILBT");
      GROUP=string("");
      THETA= angle;
      L=0;
    }

void ILBT::SetP(const char *name, double value) 
    {
     if (strcmp( name, "THETA" ) == 0) 
	{
	  THETA=value;
	}
     else 
       {
	  cout<<"ILBT does not have a parameter of "<<name<<endl; 
          exit(0); 
	} 
    }

double ILBT::GetP(const char *name) 
    {
     if (strcmp( name, "THETA" ) == 0) 
	{
	  return THETA;
	}
      else 
	{
	  cout<<"ILBT does not have a parameter of "<<name<<endl; 
          exit(0); 
	} 
    }

void ILBT::Pass(double x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  ILBT_Pass(x, THETA);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void ILBT::DAPass(tps x[6]){
  return ; 
}

void ILBT::sPass(double x[9]){   
  ILBT_Pass(x, THETA);
}

//----------------------head-on BB (BEAMBEAM)------------------------------
BEAMBEAM::BEAMBEAM(string name, int treatment, double np, double bbscale, double sigmal, int nslice, double emitx,  
	  double betax, double alfax, double emity, double betay, double alfay): Element(name)
 //     emitx, emity:  un-normalized rms emittance,  sigma=SQRT[ emitx * betax ]  
  { 
    TREATMENT=treatment;
    NP=np;
    BBSCALE=bbscale;
    SIGMAL=sigmal;
    NSLICE=nslice;
    EMITX=emitx;
    EMITY=emity;
    BETAX=betax;
    ALFAX=alfax;
    BETAY=betay;
    ALFAY=alfay;
    TYPE=string("BEAMBEAM");
    GROUP=string("");
    L=0.;
  }

void BEAMBEAM::SetP(const char *name, double value) 
  {
    if (strcmp( name, "TREATMENT" ) == 0) 
      { 
	int temp;
	temp= int(value);
	if(temp == 4 || temp == 6 ) {
	  TREATMENT=temp;}
	else{
	  cout<<"Error: beambeam treatment only can be 4 or 6"<<endl;
	  exit(1);
	}
      }
    else if (strcmp( name, "NP" ) == 0) 
      {
	NP=value;
      }
    else if (strcmp( name, "BBSCALE" ) == 0) 
      {
	BBSCALE=value;
      }
    else if (strcmp( name, "SIGMAL" ) == 0) 
      {
	SIGMAL=value;
      }
    else if (strcmp( name, "NSLICE" ) == 0) 
      {
	NSLICE=int(value);
      }
    else if (strcmp( name, "EMITX" ) == 0) 
      {
	EMITX=value;
      }
    else if (strcmp( name, "EMITY" ) == 0) 
      {
	EMITY=value;
      }
    else if (strcmp( name, "BETAX" ) == 0) 
      {
	BETAX=value;
      }
    else if (strcmp( name, "ALFAX" ) == 0) 
      {
	ALFAX=value;
      }
    else if (strcmp( name, "BETAY" ) == 0) 
      {
	BETAY=value;
      }
     else if (strcmp( name, "ALFAY" ) == 0) 
       {
	 ALFAY=value;
       }
     else 
       {
	 cout<<"BEAMBEAM does not have  parameter of "<<name<<endl; 
	 exit(0); 
       }    
  }

double BEAMBEAM::GetP(const char *name) 
  {
    if (strcmp( name, "TREATMENT" ) == 0) 
      {
	return TREATMENT;
      }
    else if (strcmp( name, "NP" ) == 0) 
      {
	return NP;
      }
    else if (strcmp( name, "BBSCALE" ) == 0) 
      {
	return BBSCALE;
      }
    else if (strcmp( name, "SIGMAL" ) == 0) 
      {
	return SIGMAL;
      }
    else if (strcmp( name, "NSLICE" ) == 0) 
      {
	return NSLICE;
      }
    else if (strcmp( name, "EMITX" ) == 0) 
      {
	return EMITX;
      }
    else if (strcmp( name, "EMITY" ) == 0) 
      {
	return EMITY;
      }
    else if (strcmp( name, "BETAX" ) == 0) 
      {
	return BETAX;
      }
    else if (strcmp( name, "ALFAX" ) == 0) 
      {
	return ALFAX;
      }
    else if (strcmp( name, "BETAY" ) == 0) 
      {
	return BETAY;
      }
    else if (strcmp( name, "ALFAY" ) == 0) 
      {
	return ALFAY;
      }
    else 
      {
	cout<<"BEAMBEAM does not have  parameter of "<<name<<endl; 
	exit(0); 
      }    
  }

void BEAMBEAM::Pass(double x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  if(NP !=0.) BEAMBEAM_Pass(x, TREATMENT,NP,BBSCALE,SIGMAL,NSLICE,EMITX,BETAX,ALFAX,EMITY,BETAY,ALFAY);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void BEAMBEAM::DAPass(tps x[6]){
  double sigmax, sigmay;
  double ksi_x, ksi_y;
  double rp= 1.534698e-18;
  sigmax=  sqrt( EMITX * BETAX );
  sigmay=  sqrt( EMITY * BETAY );
  //---Beam-beam parameter per IP is
  //  xi_z  = N * rp  * beta_z * / ( 2 * PI * gamma * sigma_z * ( sigma_x + sigma_y)
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  ksi_x=   2.0*NP*rp/sigmax/GP.gamma/(sigmax+sigmay);  
  ksi_y=   2.0*NP*rp/sigmay/GP.gamma/(sigmax+sigmay);  
  x[1]=x[1]+ksi_x * x[0] * BBSCALE; 
  x[3]=x[3]+ksi_y * x[2] * BBSCALE;
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void BEAMBEAM::sPass(double x[9]){   // need updating spin tracking
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  if(NP !=0.) BEAMBEAM_Pass(x, TREATMENT,NP,BBSCALE,SIGMAL,NSLICE,EMITX,BETAX,ALFAX,EMITY,BETAY,ALFAY);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

//---------------------LRBB ( 4-d long-range BEAMBEAM  )-----------------------------
LRBB::LRBB(string name, double np, double bbscale, double sepx, double sepy, double sigmax, double sigmay): Element(name)
  { 
    NP=np;
    BBSCALE=bbscale;
    SEPX = sepx;
    SEPY = sepy;
    SIGMAX=sigmax; 
    SIGMAY=sigmay;
    TYPE=string("LRBB");
    GROUP=string("");
    L=0.;
  }

void LRBB::SetP(const char *name, double value) 
  {
    if (strcmp( name, "NP" ) == 0) 
      {
	NP=value;
      }
    else if (strcmp( name, "BBSCALE" ) == 0) 
      {
	BBSCALE=value;
      }
    else if (strcmp( name, "SEPX" ) == 0) 
      {
	SEPX=value;
      }
    else if (strcmp( name, "SEPY" ) == 0) 
      {
	SEPY=value;
      }
    else if (strcmp( name, "SIGMAX" ) == 0) 
      {
	SIGMAX=value;
      }
    else if (strcmp( name, "SIGMAY" ) == 0) 
      {
	SIGMAY=value;
      }
    else 
       {
	 cout<<"LRBB does not have  parameter of "<<name<<endl; 
	 exit(0); 
       }    
  }

double LRBB::GetP(const char *name) 
  {
    if (strcmp( name, "NP" ) == 0) 
      {
	return NP;
      }
    else if (strcmp( name, "BBSCALE" ) == 0) 
      {
	return BBSCALE;
      }
    else if (strcmp( name, "SEPX" ) == 0) 
      {
	return SEPX;
      }
    else if (strcmp( name, "SEPY" ) == 0) 
      {
	return SEPY;
      }
    else if (strcmp( name, "SIGMAX" ) == 0) 
      {
	return SIGMAX;
      }
    else if (strcmp( name, "SIGMAY" ) == 0) 
      {
	return SIGMAY;
      }
    else 
      {
	cout<<"LRBB does not have a parameter of "<<name<<endl; 
	exit(0); 
      }    
  }

void LRBB::Pass(double x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  LRBB_Pass(x, NP, BBSCALE, SEPX, SEPY, SIGMAX, SIGMAY);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void LRBB::DAPass(tps x[6]){
  return; 
}

void LRBB::sPass(double x[9]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  LRBB_Pass(x, NP, BBSCALE, SEPX, SEPY, SIGMAX, SIGMAY);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI); 
}

//---------------------------CRBB ( 6-d BB with horizontal crossing angle )-------------------------
CRBB::CRBB(string name, double np, double theta, double bbscale, double sigmal, int nslice, 
           double emitx, double betax, double alfax, double emity, double betay, double alfay): Element(name)
      // emitx, emity:  un-normalized rms emittance,  sigma=SQRT[ emitx * betax ]  
  { 
    NP=np;
    THETA=theta;
    BBSCALE=bbscale;
    SIGMAL=sigmal;
    NSLICE=nslice;
    EMITX=emitx;
    EMITY=emity;
    BETAX=betax;
    ALFAX=alfax;
    BETAY=betay;
    ALFAY=alfay;
    TYPE=string("CRBB");
    GROUP=string("");
    L=0.;
  }

void CRBB::SetP(const char *name, double value) 
  {
    if (strcmp( name, "NP" ) == 0) 
      {
	NP=value;
      }
    else if (strcmp( name, "THETA" ) == 0) 
      {
	THETA=value;
      }
    else if (strcmp( name, "BBSCALE" ) == 0) 
      {
	BBSCALE=value;
      }
    else if (strcmp( name, "SIGMAL" ) == 0) 
      {
	SIGMAL=value;
      }
    else if (strcmp( name, "NSLICE" ) == 0) 
      {
	NSLICE=int(value);
      }
    else if (strcmp( name, "EMITX" ) == 0) 
      {
	EMITX=value;
      }
    else if (strcmp( name, "EMITY" ) == 0) 
      {
	EMITY=value;
      }
    else if (strcmp( name, "BETAX" ) == 0) 
      {
	BETAX=value;
      }
    else if (strcmp( name, "ALFAX" ) == 0) 
      {
	ALFAX=value;
      }
    else if (strcmp( name, "BETAY" ) == 0) 
      {
	BETAY=value;
      }
     else if (strcmp( name, "ALFAY" ) == 0) 
       {
	 ALFAY=value;
       }
     else 
       {
	 cout<<"CRBB does not have a parameter of "<<name<<endl; 
	 exit(0); 
       }    
  }

double CRBB::GetP(const char *name) 
  {
    if (strcmp( name, "NP" ) == 0) 
      {
	return NP;
      }
    else if (strcmp( name, "THETA" ) == 0) 
      {
	return THETA;
      }
    else if (strcmp( name, "BBSCALE" ) == 0) 
      {
	return BBSCALE;
      }
    else if (strcmp( name, "SIGMAL" ) == 0) 
      {
	return SIGMAL;
      }
    else if (strcmp( name, "NSLICE" ) == 0) 
      {
	return NSLICE;
      }
    else if (strcmp( name, "EMITX" ) == 0) 
      {
	return EMITX;
      }
    else if (strcmp( name, "EMITY" ) == 0) 
      {
	return EMITY;
      }
    else if (strcmp( name, "BETAX" ) == 0) 
      {
	return BETAX;
      }
    else if (strcmp( name, "ALFAX" ) == 0) 
      {
	return ALFAX;
      }
    else if (strcmp( name, "BETAY" ) == 0) 
      {
	return BETAY;
      }
    else if (strcmp( name, "ALFAY" ) == 0) 
      {
	return ALFAY;
      }
    else 
      {
	cout<<"CRBB does not have a parameter of "<<name<<endl; 
	exit(0); 
      }    
  }

void CRBB::Pass(double x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  CRBB_Pass(x, NP,THETA, BBSCALE,SIGMAL,NSLICE,EMITX,BETAX,ALFAX,EMITY,BETAY,ALFAY);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void CRBB::DAPass(tps x[6]){
  double sigmax, sigmay;
  double ksi_x, ksi_y;
  double rp= 1.534698e-18;
  sigmax= sqrt( EMITX * BETAX +  SIGMAL*SIGMAL*pow(tan(THETA),2.0) );
  //sigmax= sigmax * ( 1+ (SIGMAL*tan(THETA)/sigmax)*(SIGMAL*tan(THETA)/sigmax) );
  sigmay=  sqrt( EMITY * BETAY );
  //---Beam-beam parameter per IP is
  //  xi_z  = N * rp  * beta_z * / ( 2 * PI * gamma * sigma_z * ( sigma_x + sigma_y)
  //  this is true only for bunch core
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  ksi_x=   2.0*NP*rp/sigmax/GP.gamma/(sigmax+sigmay);  
  ksi_y=   2.0*NP*rp/sigmay/GP.gamma/(sigmax+sigmay);  
  x[1]=x[1]+ksi_x * x[0] * BBSCALE; 
  x[3]=x[3]+ksi_y * x[2] * BBSCALE;
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void CRBB::sPass(double x[9]){  // need updating spin tracking
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  CRBB_Pass(x, NP,THETA, BBSCALE,SIGMAL,NSLICE,EMITX,BETAX,ALFAX,EMITY,BETAY,ALFAY);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

//-----------------------------ERHICBB----------------------------------------------------
ERHICBB::ERHICBB(string name, double ne, double bbscale, double sigmax, double sigmay): Element(name)
 //  sigmax, sigmay for electron slice at IP:  only for tune-shift calculation purpose, not for tracking.
 //  electron beam only has one slice, however, its beam size  changes around IP  
  { 
    NE=ne;
    BBSCALE=bbscale;
    SIGMAX=sigmax;
    SIGMAY=sigmay;
    TYPE=string("ERHICBB");
    GROUP=string("");
    L=0.;
  }

void ERHICBB::SetP(const char *name, double value) 
  {
    if (strcmp( name, "NE" ) == 0) 
      { 
	NE=value;
      }
    else if (strcmp( name, "BBSCALE" ) == 0) 
      { 
	BBSCALE=value;
      }
    else if (strcmp( name, "SIGMAX" ) == 0) 
      {
	SIGMAX=value;
      }
    else if (strcmp( name, "SIGMAY" ) == 0) 
      {
	SIGMAY=value;
      }
    else 
      {
	cout<<"ERHICBB does not have parameter of "<<name<<endl; 
	exit(0); 
      }    
  }

double ERHICBB::GetP(const char *name) 
  {
    if (strcmp( name, "NE" ) == 0) 
      {
	return NE;
      }
    else if (strcmp( name, "BBSCALE" ) == 0) 
      {
	return BBSCALE;
      }
    else if (strcmp( name, "SIGMAX" ) == 0) 
      {
	return SIGMAX;
      }
    else if (strcmp( name, "SIGMAY" ) == 0) 
      {
	return SIGMAY;
      }
    else 
      {
	cout<<"ERHICBB does not have  parameter of "<<name<<endl; 
	exit(0); 
      }    
  }

void ERHICBB::Pass(double x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  ERHICBB_Pass(x, GP.gamma, NE, BBSCALE);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI); 
}

void ERHICBB::DAPass(tps x[6]){
  double ksi_x, ksi_y;
  double rp= 1.534698e-18;
  //---sigmax and sigmay are only used for tune shift calculation not for tracking
  //  xi_z  = N * rp  * beta_z * / ( 2 * PI * gamma * sigma_z * ( sigma_x + sigma_y)
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  ksi_x= 2.0*NE*rp/SIGMAX/GP.gamma/(SIGMAX + SIGMAY);
  ksi_y= 2.0*NE*rp/SIGMAY/GP.gamma/(SIGMAX + SIGMAY);  
  x[1]=x[1]+ksi_x * x[0] * BBSCALE ; 
  x[3]=x[3]+ksi_y * x[2] * BBSCALE ;
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI); 
}

void ERHICBB::sPass(double x[9]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  ERHICBB_Pass(x, GP.gamma, NE, BBSCALE);    // need updating spin tracking
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI); 
}

//----------------------E-LENS------------------------------
ELENS::ELENS(string name, double l, double ne, double bbscale, int nslice, double betae, double sigmax, double sigmay): Element(name)
    { 
      if (l >=0 )
	{
	  TYPE=string("ELENS");
	  GROUP=string("");
	  L=l;
	  NE=ne;
          BBSCALE=bbscale;
	  NSLICE=nslice;
	  BETAE=betae;
	  SIGMAX=sigmax;
	  SIGMAY=sigmay;
	}
      else
	{
	  cout<<"Error: check signs of L ELENS "<<endl;
          exit(1); 
	}
    }

void ELENS::SetP(const char *name, double value) 
    {
      if (strcmp( name, "NE" ) == 0) 
	{
	  NE=value;
	}
      else if (strcmp( name, "BBSCALE" ) == 0) 
	{
	  BBSCALE=value;
	}
      else if (strcmp( name, "NSLICE" ) == 0) 
	{
	  NSLICE=int(value);
	}
      else if (strcmp( name, "BETAE" ) == 0) 
	{
	  BETAE=value;
	}
      else if (strcmp( name, "SIGMAX" ) == 0) 
	{
	  SIGMAX=value;
	}
      else if (strcmp( name, "SIGMAY" ) == 0) 
	{
	  SIGMAY=value;
	}
      else 
	{
	  cout<<"ELENS does not have  parameter of "<<name<<endl; 
          exit(0); 
	}    
    }

double ELENS::GetP(const char *name) 
    {
      if (strcmp( name, "NE" ) == 0) 
	{
	  return NE;
	}
      else if (strcmp( name, "BBSCALE" ) == 0) 
	{
	  return BBSCALE;
	}
     else if (strcmp( name, "NSLICE" ) == 0) 
	{
	  return NSLICE;
	}
     else if (strcmp( name, "BETAE" ) == 0) 
	{
	  return BETAE;
	}
     else if (strcmp( name, "SIGMAX" ) == 0) 
	{
	  return SIGMAX;
	}
     else if (strcmp( name, "SIGMAY" ) == 0) 
	{
	  return SIGMAY;
	}
      else 
	{
	  cout<<"ELENS does not have  parameter of "<<name<<endl; 
          exit(0); 
	}    
    }

void ELENS::Pass(double x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  ELENS_Pass(x, L, NE, BBSCALE,BETAE, NSLICE, SIGMAX, SIGMAY);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void ELENS::DAPass(tps x[6]){
  int i;
  double ksi_x, ksi_y;
  double rp= 1.534698e-18;
  ksi_x= -2.0*(1.0*NE/NSLICE)*rp*(1+BETAE)/SIGMAX/GP.gamma/(SIGMAX+SIGMAY);  
  ksi_y= -2.0*(1.0*NE/NSLICE)*rp*(1+BETAE)/SIGMAY/GP.gamma/(SIGMAX+SIGMAY); 
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  for(i=0;i<NSLICE;i++){
    DRIFT_Pass(x, L*1.0/NSLICE/2.0);
    x[1]=x[1]+ksi_x * x[0] * BBSCALE;
    x[3]=x[3]+ksi_y * x[2] * BBSCALE;
    DRIFT_Pass(x, L*1.0/NSLICE/2.0);
  }
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void ELENS::sPass(double x[9]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  ELENS_Pass(x, L, NE, BBSCALE,BETAE, NSLICE, SIGMAX, SIGMAY);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

//----------------------HOLLOW E-LENS------------------------------
HELENS::HELENS(string name, double l, double ne, double bbscale, int nslice, double betae, double rout, double rin): Element(name)
    { 
      if (l >=0 )
	{
	  TYPE=string("HELENS");
	  GROUP=string("");
	  L=l;
	  NE=ne;
          BBSCALE=bbscale;
	  NSLICE=nslice;
	  BETAE=betae;
	  ROUT=rout;
	  RIN=rin;
	}
      else
	{
	  cout<<"Error: check signs of L. "<<endl;
          exit(1); 
	}
    }

void HELENS::SetP(const char *name, double value) 
    {
      if (strcmp( name, "NE" ) == 0) 
	{
	  NE=value;
	}
      else if (strcmp( name, "BBSCALE" ) == 0) 
	{
	  BBSCALE=value;
	}
      else if (strcmp( name, "NSLICE" ) == 0) 
	{
	  NSLICE=int(value);
	}
      else if (strcmp( name, "BETAE" ) == 0) 
	{
	  BETAE=value;
	}
      else if (strcmp( name, "ROUT" ) == 0) 
	{
	  ROUT=value;
	}
      else if (strcmp( name, "RIN" ) == 0) 
	{
	  RIN=value;
	}
      else 
	{
	  cout<<"HELENS does not have  parameter of "<<name<<endl; 
          exit(0); 
	}    
    }

double HELENS::GetP(const char *name) 
    {
      if (strcmp( name, "NE" ) == 0) 
	{
	  return NE;
	}
      else if (strcmp( name, "BBSCALE" ) == 0) 
	{
	  return BBSCALE;
	}
     else if (strcmp( name, "NSLICE" ) == 0) 
	{
	  return NSLICE;
	}
     else if (strcmp( name, "BETAE" ) == 0) 
	{
	  return BETAE;
	}
     else if (strcmp( name, "ROUT" ) == 0) 
	{
	  return ROUT;
	}
     else if (strcmp( name, "RIN" ) == 0) 
	{
	  return RIN;
	}
      else 
	{
	  cout<<"HELENS does not have  parameter of "<<name<<endl; 
          exit(0); 
	}    
    }

void HELENS::Pass(double x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  HELENS_Pass(x, L, NE, BBSCALE, BETAE, NSLICE, ROUT, RIN);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void HELENS::DAPass(tps x[6]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  DRIFT_Pass(X,L);
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

void HELENS::sPass(double x[9]){
  GtoL(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
  HELENS_Pass(x, L, NE, BBSCALE,BETAE, NSLICE, ROUT, RIN);  //  need updating spin tracking
  LtoG(x,L,0.,DX,DY,DS,DTHETA,DPHI,DPSI);
}

//---------------ROTAT-----------------------------------
ROTAT::ROTAT(string name, double l, double n[3], double angle): Element(name)
{ 
  if (l >=0 )
    {
      TYPE=string("ROTAT");
      GROUP=string("");
      L=l;
      N[0]=n[0];
      N[1]=n[1];
      N[2]=n[2];
      ANGLE=angle;
    }
  else
    {
      cout<<"Error: ROTATOR length should be non-negative."<<endl;
      exit(1); 
    }
}

void ROTAT::SetP(const char *name, double value) 
{
  if (strcmp( name, "NX" ) == 0) 
    {
      N[0]=value;
    }
  else if (strcmp( name, "NY" ) == 0) 
    {
      N[1]=value;
    }
  else if (strcmp( name, "NS" ) == 0) 
    {
      N[2]=value;
    }
  else if (strcmp( name, "ANGLE" ) == 0) 
    {
      ANGLE=value;
    }
  else{
    cout<<"No parameter to be set for ROTAT."<<endl;
    exit(0);
  }
}

double ROTAT::GetP(const char *name) 
{
  if (strcmp( name, "NX" ) == 0) 
    {
      return N[0];
    }
  else if (strcmp( name, "NY" ) == 0) 
    {
      return N[1];
    }
  else if (strcmp( name, "NS" ) == 0) 
    {
      return N[2];
    }
  else if (strcmp( name, "ANGLE" ) == 0) 
    {
      return ANGLE;
    }
  else{
    cout<<"No parameter to be set for ROTAT."<<endl;
    exit(0);
  }
}

void ROTAT::Pass(double x[6]){
  DRIFT_Pass(x, L);
}

void ROTAT::DAPass(tps x[6]){
  DRIFT_Pass(x, L); 
}

void ROTAT::sPass(double x[9]){
  ROTAT_sPass(x,L,N,ANGLE);
}

//---------------SNAKE-----------------------------------
SNAKE::SNAKE(string name, double l, double n[3], double angle): Element(name)
{ 
  if (l >=0 )
    {
      TYPE=string("SNAKE");
      GROUP=string("");
      L=l;
      N[0]=n[0];
      N[1]=n[1];
      N[2]=n[2];
      ANGLE=angle;
    }
  else
    {
      cout<<"Error: SHAKE length should be non-negative."<<endl;
      exit(1); 
    }
}

void SNAKE::SetP(const char *name, double value) 
{
  if (strcmp( name, "NX" ) == 0) 
    {
      N[0]=value;
    }
  else if (strcmp( name, "NY" ) == 0) 
    {
      N[1]=value;
    }
  else if (strcmp( name, "NS" ) == 0) 
    {
      N[2]=value;
    }
  else if (strcmp( name, "ANGLE" ) == 0) 
    {
      ANGLE=value;
    }
  else{
    cout<<"No parameter to be set for ROTAT."<<endl;
    exit(0);
  }
}

double SNAKE::GetP(const char *name) 
{
  if (strcmp( name, "NX" ) == 0) 
    {
      return N[0];
    }
  else if (strcmp( name, "NY" ) == 0) 
    {
      return N[1];
    }
  else if (strcmp( name, "NS" ) == 0) 
    {
      return N[2];
    }
  else if (strcmp( name, "ANGLE" ) == 0) 
    {
      return ANGLE;
    }
  else{
    cout<<"No parameter to be set for ROTAT."<<endl;
    exit(0);
  }
}

void SNAKE::Pass(double x[6]){
  DRIFT_Pass(x, L);
}

void SNAKE::DAPass(tps x[6]){
  DRIFT_Pass(x, L); 
}

void SNAKE::sPass(double x[9]){
  SNAKE_sPass(x,L,N,ANGLE);
}

//---------------------------DIFFUSE--------------------------------
DIFFUSE::DIFFUSE(string name, double diff_x, double diff_y, double diff_delta ): Element(name)
{ 
  TYPE=string("DIFFUSE");
  GROUP=string("");
  L=0;
  DIFF_X = diff_x;  
  DIFF_Y = diff_y; 
  DIFF_DELTA = diff_delta;
}

void DIFFUSE::SetP(const char *name, double value) 
{
  if (strcmp( name, "DIFF_X" ) == 0) 
    {
      DIFF_X = value;
    }
  else if (strcmp( name, "DIFF_Y" ) == 0) 
    {
      DIFF_Y = value;
    }
  else if (strcmp( name, "DIFF_delta" ) == 0) 
    {
      DIFF_DELTA = value;
    }
  else 
    {
      cout<<"DIFFUSE does not have a parameter of "<<name<<endl; 
      exit(0); 
    } 
}

double DIFFUSE::GetP(const char *name) 
{
  if (strcmp( name, "DIFF_X" ) == 0) 
    {
      return DIFF_X;
    }
  else if (strcmp( name, "DIFF_Y" ) == 0) 
    {
      return DIFF_Y;
    }
  else if (strcmp( name, "DIFF_DELTA" ) == 0) 
    {
      return DIFF_DELTA;
    }
  else 
    {
      cout<<"DIFFUSE does not have a parameter of "<<name<<endl; 
      exit(0); 
    } 
}

void DIFFUSE::Pass(double x[6]){
  DIFFUSE_Pass(x, DIFF_X, DIFF_Y, DIFF_DELTA);
}

void DIFFUSE::DAPass(tps x[6]){ 
  DRIFT_Pass(x, L);
}

void DIFFUSE::sPass(double x[9]){
  DIFFUSE_Pass(x, DIFF_X, DIFF_Y, DIFF_DELTA); 
}

//---------------------------COOLING--------------------------------
COOLING::COOLING(string name, double alpha ): Element(name)
{ 
  TYPE=string("COOLING");
  GROUP=string("");
  L=0;
  ALPHA=alpha;
}

void COOLING::SetP(const char *name, double value) 
{
  if (strcmp( name, "ALPHA" ) == 0) 
    {
      ALPHA=value;
    }
  else 
    {
      cout<<"COOLING does not have a parameter of "<<name<<endl; 
      exit(0); 
    } 
}

double COOLING::GetP(const char *name) 
{
  if (strcmp( name, "ALPHA" ) == 0) 
    {
      return ALPHA;
    }
  else 
    {
      cout<<"COOLING does not have a parameter of "<<name<<endl; 
      exit(0); 
    } 
}

void COOLING::Pass(double x[6]){
  COOLING_Pass(x,ALPHA);  
}

void COOLING::DAPass(tps x[6]){ 
  DRIFT_Pass(x, L);
}

void COOLING::sPass(double x[9]){
  COOLING_Pass(x,ALPHA);  
}

//--------------------TRANS (for coordinate system change)-------------------------------
TRANS::TRANS(string name, double dx, double dy, double ds): Element(name)
{ 
  TYPE=string("TRANS");
  GROUP=string("");
  L=   ds;     
  DX=  dx;
  DY=  dy;
  DS=  ds;
}

void TRANS::SetP(const char *name, double value) 
    {
      if (strcmp( name, "DX" ) == 0) 
	{
	  DX=value;
	}
      else if (strcmp( name, "DY" ) == 0) 
	{
	  DY=value;
	}
      else if (strcmp( name, "DS" ) == 0) 
	{
	  DS=value;
	}
      else 
	{
	  cout<<"TRANS does not have a parameter of "<<name<<endl; 
          exit(0); 
	}        
    }

double TRANS::GetP(const char *name) 
    {
      if (strcmp( name, "DX" ) == 0) 
	{
	  return DX;
	}
      else if (strcmp( name, "DY" ) == 0) 
	{
	  return DY;
	}
      else if (strcmp( name, "DS" ) == 0) 
	{
	  return DS;
	}
      else 
	{
	  cout<<"TRANS does not have a parameter of "<<name<<endl; 
          exit(0); 
	}        
    }

void TRANS::Pass(double x[6]){
  TRANS_Pass(x,DX,DY,DS);
}

void TRANS::DAPass(tps x[6]){
  TRANS_Pass(x,DX,DY,DS);
}

void TRANS::sPass(double x[9]){
  TRANS_Pass(x,DX,DY,DS); 
}

//--------------------SROTATION (for coordinate system change)-------------------------------
SROTAT::SROTAT(string name, double psi): Element(name)
{ 
  TYPE=string("SROTAT");
  GROUP=string("");
  L=0;
  PSI = psi;
}

void SROTAT::SetP(const char *name, double value) 
    {
      if (strcmp( name, "PSI" ) == 0) 
	{
	  PSI=value;
	}
      else 
	{
	  cout<<"SROTAT does not have a parameter of "<<name<<endl; 
          exit(0); 
	}        
    }

double SROTAT::GetP(const char *name) 
    {
      if (strcmp( name, "PSI" ) == 0) 
	{
	  return PSI;
	}
      else 
	{
	  cout<<"SROTAT does not have a parameter of "<<name<<endl; 
          exit(0); 
	}        
    }

void SROTAT::Pass(double x[6]){
  SROTAT_Pass(x, PSI);
}

void SROTAT::DAPass(tps x[6]){
  SROTAT_Pass(x, PSI); 
}

void SROTAT::sPass(double x[9]){
  SROTAT_Pass(x, PSI); 
}

//--------------------YROTATION (for coordinate system change)-------------------------------
YROTAT::YROTAT(string name, double theta): Element(name)
{ 
  TYPE=string("YROTAT");
  GROUP=string("");  
  L=0;
  THETA = theta;
}

void YROTAT::SetP(const char *name, double value) 
    {
      if (strcmp( name, "THETA" ) == 0) 
	{
	  THETA=value;
	}
      else 
	{
	  cout<<"YROTAT does not have a parameter of "<<name<<endl; 
          exit(0); 
	}        
    }

double YROTAT::GetP(const char *name) 
{
  if (strcmp( name, "THETA" ) == 0) 
    {
      return THETA;
    }
  else 
    {
      cout<<"YROTAT does not have a parameter of "<<name<<endl; 
      exit(0); 
    }        
}

void YROTAT::Pass(double x[6]){
  YROTAT_Pass(x, THETA);
}

void YROTAT::DAPass(tps x[6]){
  YROTAT_Pass(x, THETA);
}

void YROTAT::sPass(double x[9]){
  YROTAT_Pass(x, THETA);
}

//--------------------XROTATION (for coordinate system change)-------------------------------
XROTAT::XROTAT(string name, double phi): Element(name)
{ 
  TYPE=string("XROTAT");
  GROUP=string("");
  L=0;
  PHI = phi;
}

void XROTAT::SetP(const char *name, double value) 
    {
      if (strcmp( name, "PHI" ) == 0) 
	{
	  PHI=value;
	}
      else 
	{
	  cout<<"XROTAT does not have a parameter of "<<name<<endl; 
          exit(0); 
	}        
    }

double XROTAT::GetP(const char *name) 
{
  if (strcmp( name, "PHI" ) == 0) 
    {
      return PHI;
    }
  else 
    {
      cout<<"XROTAT does not have a parameter of "<<name<<endl; 
      exit(0); 
    }        
}

void XROTAT::Pass(double x[6]){
  XROTAT_Pass(x, PHI);
}

void XROTAT::DAPass(tps x[6]){
  XROTAT_Pass(x, PHI);
}

void XROTAT::sPass(double x[9]){
  XROTAT_Pass(x, PHI);
}

//--------------------PATCH (for coordinate system change)-------------------------------
PATCH::PATCH(string name, double dx, double dy, double ds, double theta): Element(name)
{ 
  TYPE=string("PATCH");
  GROUP=string("");
  L=   ds;     
  DX=  dx;
  DY=  dy;
  DS=  ds;
  THETA = theta;
}

void PATCH::SetP(const char *name, double value) 
    {
      if (strcmp( name, "DX" ) == 0) 
	{
	  DX=value;
	}
      else if (strcmp( name, "DY" ) == 0) 
	{
	  DY=value;
	}
      else if (strcmp( name, "DS" ) == 0) 
	{
	  DS=value;
	}
      else if (strcmp( name, "THETA" ) == 0) 
	{
	  THETA=value;
	}      
      else 
	{
	  cout<<"PATCH does not have a parameter of "<<name<<endl; 
          exit(0); 
	}        
    }

double PATCH::GetP(const char *name) 
    {
      if (strcmp( name, "DX" ) == 0) 
	{
	  return DX;
	}
      else if (strcmp( name, "DY" ) == 0) 
	{
	  return DY;
	}
      else if (strcmp( name, "DS" ) == 0) 
	{
	  return DS;
	}
      else if (strcmp( name, "THETA" ) == 0) 
	{
	  return THETA;
	}      
      else 
	{
	  cout<<"PATCH does not have a parameter of "<<name<<endl; 
          exit(0); 
	}        
    }

void PATCH::Pass(double x[6]){
  PATCH_Pass(x,DX,DY,DS,THETA);
}

void PATCH::DAPass(tps x[6]){
  PATCH_Pass(x,DX,DY,DS,THETA);
}

void PATCH::sPass(double x[9]){
  PATCH_Pass(x,DX,DY,DS,THETA); 
}

//=========================================
//         
//           Line definition
//
//=========================================

Line::Line()
{
  Ncell=0;
  Length=0;
  Tune1=0.;
  Tune2=0.;
  Tune3=0.;
  SpinTune=0.;
  Chromx1=0.;
  Chromy1=0.;
  Chromx2=0.;
  Chromy2=0.;
  Chromx3=0.;
  Chromy3=0.;
}

Line::~Line()
{
  Empty();
}

void Line::Update()
{
  int i;
  double spointer=0.;
  for(i=0;i<Cell.size();i++){
    spointer=Cell[i]->L+spointer;
    Cell[i]->S=spointer;
  }
  Ncell=Cell.size();  
  Length=spointer;
  frev0 = ( speed_light * GP.beta ) / Length ;
}

void Line::Append(Element * x)
  {
    Cell.push_back( x );
    Update();
  }

void Line::Delete(int i)
  {
    if (i < 0 || i >= Cell.size()) return;
 
    delete Cell[i];
    Cell.erase(Cell.begin() + i);  
    
    Update();
  }

void Line::Insert( int i, Element * temp)
  {
    if (i < 0 || i >= Cell.size()) return;

    Cell.insert(Cell.begin()+i, temp);
    
    Update();
  }

void Line::Replace( int i, Element * temp)
  {
    if (i < 0 || i >= Cell.size()) return;
    
    delete Cell[i];        
    Cell[i] = temp;
    
    Update();
  }

void Line::Rewind(int i)
{
    if (i < 0 || i >= Cell.size()) return;

    std::rotate(Cell.begin(), Cell.begin() + i, Cell.end());

    Update();
}

void Line::Invert()
{
    std::reverse(Cell.begin(), Cell.end());
}

void Line::Empty()
{
  for (auto ptr : Cell) {
    delete ptr;
  }
  Cell.clear();
  Ncell = 0;
}

/*----
void Line_Append_Element(Line & linename, Element * new_element)
{
  linename.Append(new_element);
}

void Line_Append_Elements(Line & linename1, Line & linename2)
{
  int i;
  Element * new_element;
  for(i=0;i<linename2.Ncell;i++){
    new_element=linename2.Cell[i];
    Line_Append_Element(linename1, new_element);
  }
}

void Line_Delete_Element(Line & linename, int k)
{
  int i;
  Element * new_element;
  Line  temp_line;

  if( k > linename.Ncell or k <0 ) {
    cout<<"Element i out of the range of line. exit. "<<endl;
    exit(0);
  }

  for(i=0;i<k;i++){
    new_element = linename.Cell[i];
    temp_line.Append( new_element );
  }
  for(i=k+1;i<linename.Ncell;i++){
    new_element = linename.Cell[i];
    temp_line.Append( new_element );
  }
  linename = temp_line;  //  this is dangerous!!!
}

void Line_Delete_Elements(Line & linename, int i1, int i2)
{
  int i;
  Element * new_element;
  Line  temp_line;

  if( i2 > linename.Ncell or i1 <0  or i1 >= i2) {
    cout<<"Range issue when calling Line_Delete_Elements(). exit. "<<endl;
    exit(0);
  }

  for(i=0;i<i1;i++){
    new_element = linename.Cell[i];
    temp_line.Append( new_element );
  }
  for(i=i2+1;i<linename.Ncell;i++){
    new_element = linename.Cell[i];
    temp_line.Append( new_element );
  }
  linename = temp_line;  //  this is dangerous!!!
}

void Line_Insert_Element(Line & linename, int k, Element * new_element1)
{
  int i;
  Element * new_element;
  Line  temp_line;
  
  if( k > linename.Ncell or k <0 ) {
    cout<<"Range issue when calling Line_Insert_Element(). exit. "<<endl;
    exit(0);
  }
  
  for(i=0;i<k;i++){
    new_element = linename.Cell[i];
    temp_line.Append( new_element );
  }
  
  temp_line.Append( new_element1 );
  
  for(i=k;i<linename.Ncell;i++){
    new_element = linename.Cell[i];
    temp_line.Append( new_element );
  }
  linename = temp_line;    //  this is dangerous!!!
}

void Line_Insert_Elements(Line & linename1, Line & linename2, int i1)
{
  int i;
  Element * new_element;
  
  for(i=0;i<linename2.Ncell;i++) {
    new_element = linename2.Cell[i];
    Line_Insert_Element(linename1, i1+i+1, new_element);
  }  
}

void Line_Replace_Element(Line & linename, int k, Element * new_element1)
{
  int i;
  Element * new_element;
  Line  temp_line;
  
  if( k > linename.Ncell or k <0 ) {
    cout<<"Range issue when calling Line_Replace_Element(). exit. "<<endl;
    exit(0);
  }
  
  for(i=0;i<k;i++){
    new_element = linename.Cell[i];
    temp_line.Append( new_element );
  }
  
  temp_line.Append( new_element1 );
  
  for(i=k+1;i<linename.Ncell;i++){
    new_element = linename.Cell[i];
    temp_line.Append( new_element );
  }
  linename = temp_line;
  linename.Update();
}

void Line_Rewind(Line & linename,  int k )
{
  int i;
  Element * new_element;
  Line      temp_line;
  
  for(i=k;i<linename.Ncell;i++) {
    new_element = linename.Cell[i];
    temp_line.Append( new_element );
  }
  for(i=0; i<k;i++) {
    new_element = linename.Cell[i];
    temp_line.Append( new_element );
  }

  linename=temp_line;
}

void Line_Invert(Line & linename)
{
  int i;
  Element * new_element;
  Line  temp_line;
  
  for(i=0;i<linename.Ncell;i++) {
    new_element = linename.Cell[linename.Ncell-1 - i];
    temp_line.Append( new_element );
  }
  linename=temp_line;
}

void Line_Repeat(Line & linename1, Line & linename2, int n)
{
  int i,j;
  Element * new_element;
  
  for(j=0;j<n;j++){
    for(i=0;i<linename2.Ncell;i++) {
      new_element = linename2.Cell[i];
      linename1.Append( new_element );
    }
  }
}

void Line_Connect(Line & linename1, Line & linename2, Line & linename3)
{
  int i;
  Element * new_element;

  for(i=0;i<linename2.Ncell;i++) {
    new_element = linename2.Cell[i];
    linename1.Append( new_element );
  }
  for(i=0;i<linename3.Ncell;i++) {
    new_element = linename3.Cell[i];
    linename1.Append( new_element );
  }
}

void Line_PickOut_Segment(Line & linename1, Line & linename2, int i1, int i2)
{
  int i;
  Element * new_element;

  for(i=i1;i<i2;i++) {
    new_element = linename1.Cell[i];
    linename2.Append( new_element );
  }  
}

---*/

int Count_One_Name_Element(Line & linename, const char * name)
{
  int i;
  int count=0;
  
  for(i=0;i<linename.Ncell;i++){
    if(linename.Cell[i]->NAME==string(name)){
      count++;
    }
  }

  return count;
}

int Count_One_Type_Elements(Line & linename, const char * name)
{
  int i;
  int count=0;
  for(i=0;i<linename.Ncell;i++){
    if(linename.Cell[i]->TYPE==string(name)){
      count++;
    }
  }
  return count;
}

int Get_Index(Line & linename, const char * name, int k)
{
  int i;
  int count=0;
  for(i=0;i<linename.Ncell;i++){
    if(linename.Cell[i]->NAME==string(name)){
      count++;
      if(count == k ){
	return i;
      }
    }
  }
  if(count == 0) {
    cout<<"Element "<<string(name)<<" not found. Exit."<<endl;
    return -1;
  }
  return 0;
}

void Print_One_Type_Elements(Line & linename, const char * name)
{
  int  i,count=0;

  cout<<"...Element Type :"<< name<<endl;
  for(i=0;i<linename.Ncell;i++){
    if( linename.Cell[i]->TYPE == string(name)  ){
      cout<<linename.Cell[i]->NAME<<endl;
      count++;
    }
  }
  cout<<"    there  are totally "<<count<< "   "<<name<<endl;
  
}

void Print_All_Element_Types(Line & linename)
{
  int  i,j, flag;
  vector<string> types;
  
  for(i=0;i<linename.Ncell;i++){
    
    flag=0;
    for(j=0;j<types.size();j++){
      if( linename.Cell[i]->TYPE == string(types[j]) ){
	flag =1; break;
      }
    }
    
    if( flag==0 ){
      types.push_back( linename.Cell[i]->TYPE );
    }
    
  }
  
  cout<<"....types  in the line: "<<endl;
  for(i=0;i<types.size();i++){
    cout<<"    "<<types[i]<<endl;
  }
  
}

void Print_Sbend_Parameters(Line & linename, const char* name, int  ith)
{
  int  i;
  i=Get_Index(linename, name, ith);
  if( linename.Cell[i]->TYPE != string("SBEND") ){
    cout<<" Print_Sbend_Parameters:  not  a sbend. exit. "<<endl;
  }
  else{
    cout<<linename.Cell[i]->NAME<<":  L="<<linename.Cell[i]->L<<",  S="<<linename.Cell[i]->S<<"  "<<endl;
    cout<<"Angle = "<<linename.Cell[i]->GetP("ANGLE")<<"  ,  "<<" E1 = "<<linename.Cell[i]->GetP("E1")<<"  "<<" E2 = "<<linename.Cell[i]->GetP("E2")<<endl;
  }
}

void Print_Quad_Parameters(Line & linename, const char* name, int  ith)
{
  int  i;
  i=Get_Index(linename, name, ith);
  if( linename.Cell[i]->TYPE != string("QUAD") ){
    cout<<" Print_Quad_Parameters:  not  a quadrupole. exit. "<<endl;
  }
  else{
    cout<<linename.Cell[i]->NAME<<":  L="<<linename.Cell[i]->L<<",  S="<<linename.Cell[i]->S<<"  "<<endl;
    cout<<"K1L = "<<linename.Cell[i]->GetP("K1L")<<"  ,  "<<" K1SL = "<<linename.Cell[i]->GetP("K1SL")<<endl;
  }
}

void Print_Sext_Parameters(Line & linename, const char* name, int  ith)
{
  int  i;
  i=Get_Index(linename, name, ith);
  if( linename.Cell[i]->TYPE != string("SEXT") ){
    cout<<" Print_Sext_Parameters:  not  a sextupole. exit. "<<endl;
  }
  else{
    cout<<linename.Cell[i]->NAME<<":  L="<<linename.Cell[i]->L<<",  S="<<linename.Cell[i]->S<<"  "<<endl;
    cout<<"K2L = "<<linename.Cell[i]->GetP("K2L")<<"  ,  "<<" K2SL = "<<linename.Cell[i]->GetP("K2SL")<<endl;
  }
}

void Print_Mult_Parameters(Line & linename, const char* name, int  ith)
{
  int  i;
  i=Get_Index(linename, name, ith);
  if( linename.Cell[i]->TYPE != string("MULT") ){
    cout<<" Print_Quad_Parameters:  not  a mult. exit. "<<endl;
  }
  else{
    cout<<linename.Cell[i]->NAME<<": L="<<linename.Cell[i]->L<<",  S="<<linename.Cell[i]->S<<"  "<<endl;
    cout<<" knl:={ "<<linename.Cell[i]->GetP("K0L")<<","<<linename.Cell[i]->GetP("K1L")<<","<<linename.Cell[i]->GetP("K2L")<<",";
    cout<<linename.Cell[i]->GetP("K3L")<<","<<linename.Cell[i]->GetP("K4L")<<","<<linename.Cell[i]->GetP("K5L")<<",";  
    cout<<linename.Cell[i]->GetP("K6L")<<","<<linename.Cell[i]->GetP("K7L")<<","<<linename.Cell[i]->GetP("K8L")<<",";
    cout<<linename.Cell[i]->GetP("K9L")<<","<<linename.Cell[i]->GetP("K10L")<<"},"<<endl;
    cout<<" ksl:={ "<<linename.Cell[i]->GetP("K0SL")<<","<<linename.Cell[i]->GetP("K1SL")<<","<<linename.Cell[i]->GetP("K2SL")<<",";
    cout<<linename.Cell[i]->GetP("K3SL")<<","<<linename.Cell[i]->GetP("K4SL")<<","<<linename.Cell[i]->GetP("K5SL")<<",";  
    cout<<linename.Cell[i]->GetP("K6SL")<<","<<linename.Cell[i]->GetP("K7SL")<<","<<linename.Cell[i]->GetP("K8SL")<<",";
    cout<<linename.Cell[i]->GetP("K9SL")<<","<<linename.Cell[i]->GetP("K10SL")<<"};"<<endl;
  }
}

void Print_All_Sbend_Parameters(Line & linename)
{
  int  i;
  for(i=0;i<linename.Ncell;i++){
    if( linename.Cell[i]->TYPE == string("SBEND") ){
      cout<<linename.Cell[i]->NAME<<":  L="<<linename.Cell[i]->L<<",  S="<<linename.Cell[i]->S<<"  "<<endl;
      cout<<"Angle = "<<linename.Cell[i]->GetP("ANGLE")<<"  ,  "<<" E1 = "<<linename.Cell[i]->GetP("E1")<<"  "<<" E2 = "<<linename.Cell[i]->GetP("E2")<<endl;
    }
  }
}

void Print_All_Quad_Parameters(Line & linename)
{
  int  i;
  for(i=0;i<linename.Ncell;i++){
    if( linename.Cell[i]->TYPE == string("QUAD") ){
      cout<<" Print_Quad_Parameters:  not  a quadrupole. exit. "<<endl;
      cout<<linename.Cell[i]->NAME<<":  L="<<linename.Cell[i]->L<<",  S="<<linename.Cell[i]->S<<"  "<<endl;
      cout<<"K1L = "<<linename.Cell[i]->GetP("K1L")<<"  ,  "<<" K1SL = "<<linename.Cell[i]->GetP("K1SL")<<endl;
    }
  }
}

void Print_All_Sext_Parameters(Line & linename)
{
  int  i;
  for(i=0;i<linename.Ncell;i++){
    if( linename.Cell[i]->TYPE == string("SEXT") ){
      cout<<linename.Cell[i]->NAME<<":  L="<<linename.Cell[i]->L<<",  S="<<linename.Cell[i]->S<<"  "<<endl;
      cout<<"K2L = "<<linename.Cell[i]->GetP("K2L")<<"  ,  "<<" K2SL = "<<linename.Cell[i]->GetP("K2SL")<<endl;
    }
  }
}

void Print_All_Mult_Parameters(Line & linename)
{
  int  i;
  for(i=0;i<linename.Ncell;i++){
    for(i=0;i<linename.Ncell;i++){
      if( linename.Cell[i]->TYPE ==string("MULT") ){
	cout<<linename.Cell[i]->NAME<<":   L="<<linename.Cell[i]->L<<", S= "<<linename.Cell[i]->S<<"  "<<endl;
	cout<<" knl:={ "<<linename.Cell[i]->GetP("K0L")<<","<<linename.Cell[i]->GetP("K1L")<<","<<linename.Cell[i]->GetP("K2L")<<",";
	cout<<linename.Cell[i]->GetP("K3L")<<","<<linename.Cell[i]->GetP("K4L")<<","<<linename.Cell[i]->GetP("K5L")<<",";  
	cout<<linename.Cell[i]->GetP("K6L")<<","<<linename.Cell[i]->GetP("K7L")<<","<<linename.Cell[i]->GetP("K8L")<<",";
	cout<<linename.Cell[i]->GetP("K9L")<<","<<linename.Cell[i]->GetP("K10L")<<"},"<<endl;
	cout<<" ksl:={ "<<linename.Cell[i]->GetP("K0SL")<<","<<linename.Cell[i]->GetP("K1SL")<<","<<linename.Cell[i]->GetP("K2SL")<<",";
	cout<<linename.Cell[i]->GetP("K3SL")<<","<<linename.Cell[i]->GetP("K4SL")<<","<<linename.Cell[i]->GetP("K5SL")<<",";  
      cout<<linename.Cell[i]->GetP("K6SL")<<","<<linename.Cell[i]->GetP("K7SL")<<","<<linename.Cell[i]->GetP("K8SL")<<",";
      cout<<linename.Cell[i]->GetP("K9SL")<<","<<linename.Cell[i]->GetP("K10SL")<<"};"<<endl;
      }
      
    }
    
  }
}

void Get_Mult_Parameters( Line & linename, int i, double knl[11], double knsl[11] )
{
 
  if(linename.Cell[i]->TYPE != string("MULT") ){
    cout<<"  Get_MULT_KNl_KNSL(), not a MULT. exit "<<endl;
    exit(0);
  }
  else{
      knl[0] =  linename.Cell[i]->GetP("K0L");
      knl[1] =  linename.Cell[i]->GetP("K1L");
      knl[2] =  linename.Cell[i]->GetP("K2L");
      knl[3] =  linename.Cell[i]->GetP("K3L");
      knl[4] =  linename.Cell[i]->GetP("K4L");
      knl[5] =  linename.Cell[i]->GetP("K5L");
      knl[6] =  linename.Cell[i]->GetP("K6L");
      knl[7] =  linename.Cell[i]->GetP("K7L");
      knl[8] =  linename.Cell[i]->GetP("K8L");
      knl[9] =  linename.Cell[i]->GetP("K9L");
      knl[10]=  linename.Cell[i]->GetP("K10L");
	
      knsl[0] =  linename.Cell[i]->GetP("K0SL");
      knsl[1] =  linename.Cell[i]->GetP("K1SL");
      knsl[2] =  linename.Cell[i]->GetP("K2SL");
      knsl[3] =  linename.Cell[i]->GetP("K3SL");
      knsl[4] =  linename.Cell[i]->GetP("K4SL");
      knsl[5] =  linename.Cell[i]->GetP("K5SL");
      knsl[6] =  linename.Cell[i]->GetP("K6SL");
      knsl[7] =  linename.Cell[i]->GetP("K7SL");
      knsl[8] =  linename.Cell[i]->GetP("K8SL");
      knsl[9] =  linename.Cell[i]->GetP("K9SL");
      knsl[10] = linename.Cell[i]->GetP("K10SL");
  }

}


void Get_Element_Parameter( Line & linename, const char* name, int ith, const char* para ){
  cout<<linename.Cell[ Get_Index( linename, name, ith ) ]->GetP(para) <<endl;
}

void Print_Line_Elements(Line & linename, const char* filename)
{
  
  int i;
  fstream fout;
  fout.open(filename, ios::out);
  for(i=0;i<linename.Ncell;i++){
    fout <<setw(15) <<linename.Cell[i]->NAME<<setw(15) <<linename.Cell[i]->TYPE<<setw(15) <<linename.Cell[i]->L<<setw(15)<<linename.Cell[i]->S<<endl;
  }
  fout.close();
  
}

void Set_Integration_Steps(double bl, double ql, double gl )
{
  BLslice = bl;
  QLslice = ql;
  GLslice = gl;
}

void Change_Integration_Steps(Line & linename, double bl, double ql, double gl )
{
  int i;
  int nint;
  
  for(i=0;i<linename.Ncell;i++){
    if(linename.Cell[i]->TYPE== string("SBEND")  || linename.Cell[i]->TYPE== string("RBEND")  ||
       linename.Cell[i]->TYPE== string("GSBEND") || linename.Cell[i]->TYPE== string("GMULT")  ||
       linename.Cell[i]->TYPE== string("SBENDMULT") || linename.Cell[i]->TYPE== string("GSBENDMULT") ){
      nint=int(linename.Cell[i]->L / bl );
      linename.Cell[i]->SetP("Nint",nint);
    }
    
    if(linename.Cell[i]->TYPE== string("QUAD")  || linename.Cell[i]->TYPE== string("SEXT")  ||
       linename.Cell[i]->TYPE== string("OCT") || linename.Cell[i]->TYPE== string("MULT") ) { 
      nint=int(linename.Cell[i]->L / ql );
      linename.Cell[i]->SetP("Nint",nint);
    }    
  }
  
}

void Change_Nint_Bend(Line & linename, double ds)
{
  int i;

  for(i=0;i<linename.Ncell;i++){
    if( linename.Cell[i]->TYPE == string("SBEND") ){
      linename.Cell[i]->SetP("Nint", int( linename.Cell[i]->L/ds) +1 );
    }
  }
}

void  Change_Nint_Mult(Line & linename, double ds)
{
  int i;

  for(i=0;i<linename.Ncell;i++){
    if( linename.Cell[i]->TYPE == string("QUAD") || linename.Cell[i]->TYPE == string("SEXT")  ||
	linename.Cell[i]->TYPE == string("OCT")  || linename.Cell[i]->TYPE == string("MULT")  ||
        linename.Cell[i]->TYPE == string("SBENDMULT") ||  linename.Cell[i]->TYPE == string("GMULT") ) {
      linename.Cell[i]->SetP("Nint", int( linename.Cell[i]->L/ds) +1 );
    }
  }
}

void Split_Drift(Line & linename, int i, int m )
{
  int j;
  string old_name;
  double length;
  Element * new_element;

  if(m<1) {cout<<" Number of split slices should be => 1 !"<<endl; exit(1);}
  
  if(linename.Cell[i]->TYPE==string("DRIFT")){
    old_name=linename.Cell[i]->NAME;
    length=linename.Cell[i]->L;
    linename.Delete(i);    
    for(j=0;j<m;j++) {
      new_element= new DRIFT(old_name,length/m);
      linename.Insert(i,new_element);
    }
  }
  else{
    cout<<" The i-th element is not DRIFT."<<endl;
    exit(1);
  }
}

void Split_SBend(Line & linename, int i, int m )
{
  int j;
  string old_name;
  double length,angle, e1, e2;
  Element * new_element;

  if(m<1) {cout<<" Number of split slices should be => 1 !"<<endl; exit(1);}
  
  if(linename.Cell[i]->TYPE==string("SBEND")){
    old_name=linename.Cell[i]->NAME;
    length=linename.Cell[i]->L;
    angle=  linename.Cell[i]->GetP("ANGLE");
    e1   =  linename.Cell[i]->GetP("E1");
    e2   =  linename.Cell[i]->GetP("E2");
    linename.Delete(i);    
    for(j=0;j<m;j++) {
      new_element= new SBEND(old_name,length/m, angle/m, 0,0);
      linename.Insert(i,new_element);
    }
    linename.Cell[i]->SetP("E1",e1);
    linename.Cell[i+m-1]->SetP("E2",e2);
  }
  else{
    cout<<" The i-th element is not SBEND ."<<endl;
    exit(1);
  }
}

void Split_Quad(Line & linename, int i, int m )
{
  int j;
  string old_name;
  double length,k1l,k1sl;
  Element * new_element;
  
  if(m<1) {cout<<" Number of split slices should be => 1 !"<<endl; exit(1);}

  if(linename.Cell[i]->TYPE==string("QUAD")){
    old_name=linename.Cell[i]->NAME;
    length=linename.Cell[i]->L;
    k1l =  linename.Cell[i]->GetP("K1L");
    k1sl=  linename.Cell[i]->GetP("K1SL");
    linename.Delete(i);
    for(j=0;j<m;j++) {
      new_element= new QUAD(old_name,length/m, k1l/m, k1sl/m);
      linename.Insert(i,new_element);
    }
  }
  else{
    cout<<" The i-th element is not QUAD." <<endl;
    exit(1);
  }
}

void Split_Sext(Line & linename, int i, int m )
{
  int j;
  string old_name;
  double length,k2l, k2sl;
  Element * new_element;
  string new_name;
  
  if(m<1) {cout<<" Number of split slices should be => 1 !"<<endl; exit(1);}
  
  if(linename.Cell[i]->TYPE==string("SEXT")){
    old_name=linename.Cell[i]->NAME;
    length=linename.Cell[i]->L;
    k2l=  linename.Cell[i]->GetP("K2L");
    k2sl=  linename.Cell[i]->GetP("K2SL");
    linename.Delete(i);    
    for(j=0;j<m;j++) {
      new_element= new SEXT(old_name,length/m, k2l/m, k2sl/m);
      linename.Insert(i,new_element);
    }
  }
  else{
    cout<<" The i-th element is not SEXT." <<endl;
    exit(1);
  }
}

void Split_Mult(Line & linename, int i, int m )
{
  int j;
  string old_name;
  double length;
  double knl[11],knsl[11];
  Element * new_element;

  if(m<1) {cout<<" Number of split slices should be => 1 !"<<endl; exit(1);}
  
  if(linename.Cell[i]->TYPE==string("MULT")){
    old_name=linename.Cell[i]->NAME;
    length=linename.Cell[i]->L;
    
    for(j=0;j<11;j++){
      char name1[125], name2[125];
      sprintf(name1, "K%dL",j);
      sprintf(name2, "K%dSL",j);
      knl[j]= linename.Cell[i]->GetP(name1)/m;
      knsl[j]=linename.Cell[i]->GetP(name2)/m;
    }

    linename.Delete(i);    
    for(j=0;j<m;j++) {
      new_element= new MULT(old_name,length/m,knl, knsl);
      linename.Insert(i,new_element);
    }
  }
  else{
    cout<<" The i-th element is not MULT."<<endl;
    exit(1);
  }
}

void Split_GMULT(Line & linename, int i, int m )
{
  int j;
  string old_name;
  double length;
  double angle,  e1, e2, knl[11],knsl[11];
  Element * new_element;

  if(m<1) {cout<<" Number of split slices should be => 1 !"<<endl; exit(1);}
  
  if(linename.Cell[i]->TYPE==string("GMULT")){
    old_name=linename.Cell[i]->NAME;
    length=linename.Cell[i]->L;
    angle= linename.Cell[i]->GetP("ANGLE");  
    e1   =  linename.Cell[i]->GetP("E1");
    e2   =  linename.Cell[i]->GetP("E2");
    
    for(j=0;j<11;j++){
      char name1[125], name2[125];
      sprintf(name1, "K%dL",j);
      sprintf(name2, "K%dSL",j);
      knl[j]= linename.Cell[i]->GetP(name1)/m;
      knsl[j]=linename.Cell[i]->GetP(name2)/m;
    }
    
    linename.Delete(i);    
    for(j=0;j<m;j++) {
      new_element= new GMULT(old_name, length/m, angle/m, 0, 0, knl, knsl);
      linename.Insert(i,new_element);
    }
    linename.Cell[i]->SetP("E1",e1);
    linename.Cell[i+m-1]->SetP("E2",e2);
  }
  else{
    cout<<" The i-th element is not GMULT."<<endl;
    exit(1);
  }
}

void Split_SBENDMULT(Line & linename, int i, int m )
{
  int j;
  string old_name;
  double length;
  double angle,  e1, e2, knl[11],knsl[11];
  Element * new_element;

  if(m<1) {cout<<" Number of split slices should be => 1 !"<<endl; exit(1);}
  
  if(linename.Cell[i]->TYPE==string("SBENDMULT")){
    old_name=linename.Cell[i]->NAME;
    length=linename.Cell[i]->L;
    angle= linename.Cell[i]->GetP("ANGLE");  
    e1   =  linename.Cell[i]->GetP("E1");
    e2   =  linename.Cell[i]->GetP("E2");
    
    for(j=0;j<11;j++){
      char name1[125], name2[125];
      sprintf(name1, "K%dL",j);
      sprintf(name2, "K%dSL",j);
      knl[j]= linename.Cell[i]->GetP(name1)/m;
      knsl[j]=linename.Cell[i]->GetP(name2)/m;
    }
    
    linename.Delete(i);    
    for(j=0;j<m;j++) {
      new_element= new SBENDMULT(old_name, length/m, angle/m, 0, 0, knl, knsl);
      linename.Insert(i,new_element);
    }
    linename.Cell[i]->SetP("E1",e1);
    linename.Cell[i+m-1]->SetP("E2",e2);
  }
  else{
    cout<<" The i-th element is not SBENDMULT."<<endl;
    exit(1);
  }
}

void Split_SMult(Line & linename, int i, int m )
{
  int j;
  string old_name;
  double length;
  double angle,  e1, e2, knl[11],knsl[11];
  Element * new_element;

  if(m<1) {cout<<" Number of split slices should be => 1 !"<<endl; exit(1);}
  
  if(linename.Cell[i]->TYPE==string("SMULT")){
    old_name=linename.Cell[i]->NAME;
    length=linename.Cell[i]->L;
    angle= linename.Cell[i]->GetP("ANGLE");  
    e1   =  linename.Cell[i]->GetP("E1");
    e2   =  linename.Cell[i]->GetP("E2");
    
    for(j=0;j<11;j++){
      char name1[125], name2[125];
      sprintf(name1, "K%dL",j);
      sprintf(name2, "K%dSL",j);
      knl[j]= linename.Cell[i]->GetP(name1)/m;
      knsl[j]=linename.Cell[i]->GetP(name2)/m;
    }
    
    linename.Delete(i);    
    for(j=0;j<m;j++) {
      new_element= new SMULT(old_name, length/m, angle/m, 0, 0, knl, knsl);
      linename.Insert(i,new_element);
    }
    linename.Cell[i]->SetP("E1",e1);
    linename.Cell[i+m-1]->SetP("E2",e2);
  }
  else{
    cout<<" The i-th element is not SMULT."<<endl;
    exit(1);
  }
}

void Split_Drift_All(Line & linename,int m )
{
  int i;

  i=0;
  do{
    if(linename.Cell[i]->TYPE==string("DRIFT") and linename.Cell[i]->L !=0. ){
      Split_Drift(linename, i, m );
      i=i+m;
    }
    else{
      i++;
    }
  } while( i<linename.Ncell );

}

void Split_Sbend_All(Line & linename,int m )
{
  int i;

  i=0;
  do{
    if(linename.Cell[i]->TYPE==string("SBEND") and linename.Cell[i]->L !=0. ){
      Split_SBend(linename, i, m );
      i=i+m;
    }
    else{
      i++;
    }
  } while( i<linename.Ncell );

}

void Split_Quad_All(Line & linename,int m )
{
  int i;

  i=0;
  do{
    if(linename.Cell[i]->TYPE==string("QUAD") and linename.Cell[i]->L !=0. ){
      Split_Quad(linename, i, m );
      i=i+m;
    }
    else{
      i++;
    }
  } while( i<linename.Ncell );

}

void Split_Sext_All(Line & linename,int m )
{
  int i;

  i=0;
  do{
    if(linename.Cell[i]->TYPE==string("SEXT") and linename.Cell[i]->L !=0. ){
      Split_Sext(linename, i, m );
      i=i+m;
    }
    else{
      i++;
    }
  } while( i<linename.Ncell );

}

void Split_Mult_All(Line & linename,int m )
{
  int i;

  i=0;
  do{
    if(linename.Cell[i]->TYPE==string("MULT") and linename.Cell[i]->L !=0. ){
      Split_Mult(linename, i, m );
      i=i+m;
    }
    else{
      i++;
    }
  } while( i<linename.Ncell );

}

void Split_GMULT_All(Line & linename,int m )
{
  int i;

  i=0;
  do{
    if(linename.Cell[i]->TYPE==string("GMULT") and linename.Cell[i]->L !=0. ){
      Split_GMULT(linename, i, m );
      i=i+m;
    }
    else{
      i++;
    }
  } while( i<linename.Ncell );

}

void Split_SBENDMULT_All(Line & linename,int m )
{
  int i;

  i=0;
  do{
    if(linename.Cell[i]->TYPE==string("SBENDMULT") and linename.Cell[i]->L !=0. ){
      Split_SBENDMULT(linename, i, m );
      i=i+m;
    }
    else{
      i++;
    }
  } while( i<linename.Ncell );

}

void Split_SMult_All(Line & linename,int m )
{
  int i;

  i=0;
  do{
    if(linename.Cell[i]->TYPE==string("SMULT") and linename.Cell[i]->L !=0. ){
      Split_SMult(linename, i, m );
      i=i+m;
    }
    else{
      i++;
    }
  } while( i<linename.Ncell );

}

//----function: concate adjacient DRIFT, doesn't change Twiss
void Concat_Drift(Line & linename)
{
  int i;
  double length;
  
  i=1;
  while(i < linename.Ncell) {
    if( linename.Cell[i-1]->TYPE == string("DRIFT")  and  linename.Cell[i]->TYPE == string("DRIFT") ){
      length=linename.Cell[i]->L;
      linename.Delete(i);  
      linename.Cell[i-1]->L=linename.Cell[i-1]->L + length;
      linename.Update();
    }
    else{
      i++;
    }
  }
}

//----get rid of zero strength element, doesn't change Twiss
void Clean_Up(Line & linename)
{ 
  int i;  
  double length;
  Element *temp_element;

  for(i=0;i<linename.Ncell;i++) {
    if( linename.Cell[i]->TYPE==string("SBEND") and linename.Cell[i]->GetP("ANGLE") ==0 ) {
      length= linename.Cell[i]->L;
      linename.Delete(i);
      temp_element= new DRIFT("TEMPD",length);
      linename.Insert(i,temp_element);
    }
  }

  for(i=0;i<linename.Ncell;i++) {
    if( linename.Cell[i]->TYPE==string("QUAD") and linename.Cell[i]->GetP("K1L") ==0 and linename.Cell[i]->GetP("K1SL") ==0  ) {
      length= linename.Cell[i]->L;
      linename.Delete(i);
      temp_element= new DRIFT("TEMPD",length);
      linename.Insert(i,temp_element);
    }
  }

  for(i=0;i<linename.Ncell;i++) {
    if( linename.Cell[i]->TYPE==string("SKEWQ") and linename.Cell[i]->GetP("K1SL") ==0  ) {
      length= linename.Cell[i]->L;
      linename.Delete(i);
      temp_element= new DRIFT("TEMPD",length);
      linename.Insert(i,temp_element);
    }
  }

  for(i=0;i<linename.Ncell;i++) {
    if( linename.Cell[i]->TYPE==string("SEXT") and linename.Cell[i]->GetP("K2L") ==0 and linename.Cell[i]->GetP("K2SL")  ) {
      length= linename.Cell[i]->L;
      linename.Delete(i);
      temp_element= new DRIFT("TEMPD",length);
      linename.Insert(i,temp_element);
    }
  }

  for(i=0;i<linename.Ncell;i++) {
    if( linename.Cell[i]->TYPE==string("OCT") and linename.Cell[i]->GetP("K3L") ==0 and linename.Cell[i]->GetP("K3SL")  ) {
      length= linename.Cell[i]->L;
      linename.Delete(i);
      temp_element= new DRIFT("TEMPD",length);
      linename.Insert(i,temp_element);
    }
  }

  for(i=0;i<linename.Ncell;i++) {
    if( linename.Cell[i]->TYPE==string("SOLEN") and linename.Cell[i]->GetP("KS") ==0  ) {
      length= linename.Cell[i]->L;
      linename.Delete(i);
      temp_element= new DRIFT("TEMPD",length);
      linename.Insert(i,temp_element);
    }
  }

  for(i=0;i<linename.Ncell;i++) {
    if( linename.Cell[i]->TYPE==string("RFCAV") and linename.Cell[i]->GetP("VRF") ==0  ) {
      length= linename.Cell[i]->L;
      linename.Delete(i);
      temp_element= new DRIFT("TEMPD",length);
      linename.Insert(i,temp_element);
    }
  }
  
  for(i=0;i<linename.Ncell;i++)
    if( linename.Cell[i]->TYPE==string("KICKER") and linename.Cell[i]->GetP("HKICK") ==0.  and  linename.Cell[i]->GetP("VKICK") ==0.  ) {
      length= linename.Cell[i]->L;
      linename.Delete(i);  
      if (length != 0. ) {
	temp_element=new DRIFT("TEMPD", length);
	linename.Insert(i, temp_element);
      }
    }
  
  for(i=0;i<linename.Ncell;i++)
    if( linename.Cell[i]->TYPE==string("HKICKER") and linename.Cell[i]->GetP("HKICK") ==0. ) {
      length= linename.Cell[i]->L;
      linename.Delete(i);  
      if (length != 0. ) {
	temp_element=new DRIFT("TEMPD", length);
	linename.Insert(i, temp_element);
      }
    }
  
  for(i=0;i<linename.Ncell;i++)
    if( linename.Cell[i]->TYPE==string("VKICKER") and linename.Cell[i]->GetP("VKICK") ==0. ) {
      length= linename.Cell[i]->L;
      linename.Delete(i);  
      if (length != 0. ) {
	temp_element=new DRIFT("TEMPD", length);
	linename.Insert(i, temp_element);
      }
    }

  for(i=0;i<linename.Ncell;i++)
    if(    linename.Cell[i]->TYPE==string("MARKER") || linename.Cell[i]->TYPE==string("BPM")  ||
	   linename.Cell[i]->TYPE==string("HBPM")   || linename.Cell[i]->TYPE==string("VBPM")   ){
      length= linename.Cell[i]->L;
      linename.Delete(i);  
      temp_element=new DRIFT("TEMPD", length);
      linename.Insert(i, temp_element);
    }
  
  for(i=0;i<linename.Ncell;i++)
    if( linename.Cell[i]->TYPE==string("MULT") ) {
      length= linename.Cell[i]->L;
      int m;
      double kl=0.;
      double  knl[11], knsl[11];
      for(m=0;m<11;m++){
        char name1[125], name2[125];
        sprintf(name1, "K%dL",m);
        sprintf(name2, "K%dSL",m);
        knl[m] =linename.Cell[i]->GetP(name1);
        knsl[m]=linename.Cell[i]->GetP(name2);
      }
      for(m=0;m<11;m++)  kl=kl+abs(knl[m])+abs(knsl[m]);
      
      if(kl==0.){
	linename.Delete(i);  
	temp_element=new DRIFT("TEMPD",length);
	linename.Insert(i, temp_element);
      }
    }
  Concat_Drift(linename);
}

//--- make thin of some elements  with drift-kick-drift model,  may change Twiss
void Make_Thin(Line & linename) 
{
  int i;
  double length;
  Element * temp_element;

  i=0;
  while(i < linename.Ncell-1 ) {
    if( linename.Cell[i]->TYPE ==string("SEXT")    || linename.Cell[i]->TYPE ==string("OCT")     ||
        linename.Cell[i]->TYPE ==string("MULT")    || linename.Cell[i]->TYPE ==string("RFCAV")   ||
	linename.Cell[i]->TYPE ==string("HKICKER") || linename.Cell[i]->TYPE ==string("VKICKER") ||
        linename.Cell[i]->TYPE ==string("HACMULT") || linename.Cell[i]->TYPE ==string("VACMULT")  )
      {
	length=linename.Cell[i]->L;
	if(length != 0. ) {
	  temp_element=new DRIFT("TEMPD", length/2.0);
	  linename.Insert(i, temp_element);
	  linename.Cell[i+1]->L =0.;
	  temp_element=new DRIFT("TEMPD", length/2.0);
	  linename.Insert(i+2, temp_element);
	}
      }
    i++;  
  }

  Concat_Drift(linename);
}

double Get_KL(Line & linename, const char * name, const char *  kl)
{
  int i;
  for(i=0;i<linename.Ncell;i++) {
    if(linename.Cell[i]->NAME==string(name) ) {
      return linename.Cell[i]->GetP(kl);
    }
  }  
  return 0.;
}

void Set_KL(Line & linename, const char * name, const char *  kl, double strength)
{
  int i;
  for(i=0;i<linename.Ncell;i++) {
    if(linename.Cell[i]->NAME==string(name) )  linename.Cell[i]->SetP(kl, strength);
  }  
}

void Set_KL(Line & linename, int  index,  const char *  kl, double strength)
{
 linename.Cell[index]->SetP(kl, strength);
}


void Set_KL(Line & linename, vector<int> index_elements,  const char *  kl, double strength)
{
  int  i;

  for(i=0;i<index_elements.size();i++) {
    linename.Cell[index_elements[i]]->SetP(kl, strength);
  }
  
}

void Set_dKL(Line & linename, const char * name, const char * kl, double dstrength)
{
  int i;
  for(i=0;i<linename.Ncell;i++) {
    if(linename.Cell[i]->NAME==string(name) ) {
      linename.Cell[i]->SetP(kl, linename.Cell[i]->GetP(kl) + dstrength);
    }
  } 
}

void Set_dKL(Line & linename, int index, const char * kl, double dstrength)
{
  linename.Cell[index]->SetP(kl, linename.Cell[index]->GetP(kl) + dstrength);
}

void Set_dKL(Line & linename, vector<int> index_elements,  const char *  kl, double dstrength)
{
  int  i;

  for(i=0;i<index_elements.size();i++) {
    linename.Cell[ index_elements[i] ]->SetP(kl, linename.Cell[ index_elements[i] ]->GetP(kl) + dstrength);
  }
  
}

void  Cal_CC_Voltage(double theta, double freq_cc, double beta_star, double beta_cc, double energy, double & vrf)
//   theta half cross angle, energy in unit of MeV, vrf in unit of MV, for ideal 90 degrees
{
  double  c =  2.99792458e8;  // speed of light;
  energy    = energy * 1e6;  // to unit of eV 
  
  vrf = c  * energy * ( 2 * theta ) / sqrt( beta_star  * beta_cc )/ 4 /PI/ freq_cc / 1e6 ;
}

double   Cal_CC_Voltage_To_Kick(double Vcc, double freq_cc, double energy)
//  Vcc  in MV,  energy  in MeV, freq in Hz
{
  return  Vcc * 2*PI *freq_cc / energy / 2.99792458e8;
}

double   Cal_CC_Kick_To_Voltage(double kick, double freq_cc, double energy)
//  Vcc  in MV,  energy  in MeV, freq in Hz
{
  return  kick* energy  * 2.99792458e8 /  (  2*PI*freq_cc);
}

void  Read_MADX_Lattice(Line & linename, const char * filename)
{
  int i;
  fstream f1;
  string line;
  string temp_name;
  string temp_type;
  double temp_s;
  double temp_l;
  double temp_angle;
  double temp_e1;
  double temp_e2;
  double temp_tilt;
  double knl[11], knsl[11];
  double temp_sk;
  double hkick=0;
  double vkick=0;
  Element * temp_element;
  int nskip;

  f1.open(filename,ios::in);
  if(!f1)
    {
      cout<<"error in opening the file: "<<filename<<endl;
      exit(0);
    }

  for(i=0;i<100;i++) {
    getline(f1,line,'\n');
    if( line.find( string("KEYWORD") ) != string::npos ){
      break;
    }
  }
  nskip=i+2;
  
  f1.clear();
  f1.seekg(0);
  
  for(i=0;i<nskip;i++) getline(f1,line,'\n');
  
  double s1=0;
  while( getline(f1,line,'\n')){
    istringstream ss(line);
    ss>>temp_name>>temp_type>>temp_s>>temp_l>>temp_angle>>temp_e1>>temp_e2>>temp_tilt
      >>knl[0]>>knsl[0]
      >>knl[1]>>knsl[1]
      >>knl[2]>>knsl[2]
      >>knl[3]>>knsl[3]
      >>knl[4]>>knsl[4]
      >>knl[5]>>knsl[5]
      >>knl[6]>>knsl[6]
      >>knl[7]>>knsl[7]
      >>knl[8]>>knsl[8]
      >>knl[9]>>knsl[9]
      >>knl[10]>>knsl[10]
      >>temp_sk;  //>>hkick>>vkick;
    temp_name=string(temp_name,1,temp_name.size()-2 );     
    temp_type=string(temp_type,1,temp_type.size()-2 );
    
    s1=s1+temp_l;
    if(temp_type==string("DRIFT")){
      temp_element=new DRIFT(temp_name,temp_l);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("SBEND"))    {
      if(knl[0] != 0. and knl[2] == 0.  ) {         //  for sbend with By  !=  ideal By0
	temp_element= new GSBEND(temp_name,temp_l,temp_angle, knl[0], temp_e1, temp_e2);
	temp_element->DPSI = temp_tilt;
	linename.Cell.push_back(temp_element);
      }
      else if (knl[1] == 0. and  knl[2] == 0.  ) {  // for pure sbend, k0l, k1l, k2l all zeros
	temp_element= new SBEND(temp_name,temp_l,temp_angle, temp_e1, temp_e2);
        temp_element->DPSI = temp_tilt;
	linename.Cell.push_back(temp_element);
      }
      else {                                       //  for sbend with multipole errors
	temp_element=new GMULT(temp_name,temp_l, temp_angle, temp_e1, temp_e2, knl, knsl);
	temp_element->DPSI = temp_tilt;
	linename.Cell.push_back(temp_element);
      }
    }
    else if (temp_type==string("RBEND"))    {
      if(knl[1] == 0. and knl[2]==0.) {
	//temp_l = abs( 0.5* temp_l / sin(0.5*temp_angle) * temp_angle ); // madx already ouput arc length
	temp_element= new SBEND(temp_name, temp_l, temp_angle, temp_angle/2, temp_angle/2);
	linename.Cell.push_back(temp_element); }
      else {
	temp_element=new GMULT(temp_name,temp_l, temp_angle, temp_e1, temp_e2, knl, knsl);
	linename.Cell.push_back(temp_element);
      }
    }
    else if (temp_type==string("QUADRUPOLE")){  
      temp_element= new QUAD(temp_name,temp_l, knl[1], knsl[1]);
      temp_element->DPSI = temp_tilt;
      linename.Cell.push_back(temp_element); 
    }
    else if (temp_type==string("SEXTUPOLE")){   
      temp_element= new SEXT(temp_name,temp_l, knl[2], knsl[2]);
      temp_element->DPSI = temp_tilt;
      linename.Cell.push_back(temp_element); 
    }
    else if (temp_type==string("OCTUPOLE")){    
      temp_element= new OCT(temp_name,temp_l,knl[3],knsl[3]);
      temp_element->DPSI = temp_tilt;
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("HKICKER")){
      temp_element= new HKICKER(temp_name,temp_l, hkick );
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("VKICKER")){
      temp_element= new VKICKER(temp_name,temp_l, vkick);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("TKICKER")){
      temp_element= new KICKER(temp_name,temp_l, hkick, vkick);
      linename.Cell.push_back(temp_element);
    }
    
    /*-----
    else if (temp_type==string("MULTIPOLE")){   // length can be non-zero. For madx:  k0l-->angle, hkick, vkick --> correctors
      if(knl[0] == 0 ){
      temp_element=new MULT(temp_name,temp_l, knl, knsl);
      linename.Cell.push_back(temp_element);
      }
      else{
	temp_angle = knl[0];  	knl[0]=0;   temp_l=1.0e-6;
	temp_element=new GMULT(temp_name, temp_l, temp_angle, 0, 0, knl, knsl);
	linename.Cell.push_back(temp_element);
      }	
    }
    ----*/

    /*----
    else if (temp_type==string("MULTIPOLE")){   // In madx:  k0l-->angle, hkick, vkick --> correctors
      if(knl[0] == 0 ){
      temp_element=new MULT(temp_name,temp_l, knl, knsl);
      temp_element->DPSI = temp_tilt;
      linename.Cell.push_back(temp_element);
      }
      else{
	temp_angle = knl[0];  	knl[0]=0;   temp_l=1.0e-6;
	//temp_element=new SBEND(temp_name+"_tb",temp_l,temp_angle, 0, 0);
	//linename.Cell.push_back(temp_element);
	temp_l=0;
	temp_element=new MULT(temp_name,temp_l, knl, knsl); 
	temp_element->DPSI = temp_tilt;
	linename.Cell.push_back(temp_element);
      }	
    }
    ----*/
    
    else if (temp_type==string("MULTIPOLE")){   // In madx:  thin multipoles
      temp_element=new MULT(temp_name,temp_l, knl, knsl);
      temp_element->DPSI = temp_tilt;
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("VMONITOR")){
      temp_element= new VBPM(temp_name,temp_l);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("MONITOR")){
      temp_element= new BPM(temp_name,temp_l);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("MARKER")){
      temp_element= new MARKER(temp_name,temp_l);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("SOLENOID")){  
      if(temp_sk !=0.0  and temp_l !=0. ){    // meaningful input for solenoid 
	temp_element= new SOLEN(temp_name,temp_l, temp_sk);
	linename.Cell.push_back(temp_element);
      }
      else{
	temp_element=new DRIFT(temp_name,temp_l);
	linename.Cell.push_back(temp_element);
      }
    }
    else if (temp_type==string("RFCAVITY")){ 
      temp_element= new RFCAV(temp_name,temp_l, 0.0, 0.0, 0.);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("BEAMBEAM")){
      if(temp_l !=0.){ 
	cout<<"Beam-Beam element should have zero length!"<<endl;
        exit(1);  
      }
      temp_element= new BEAMBEAM(temp_name,6, 0., 1, 0.4545, 11, 2.5e-06, 2.5e-06, 0.53, 0.0, 0.53, 0.0);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("ELENS")){
      temp_element= new ELENS(temp_name, temp_l, 0.0, 1, 0.14, 4, .31e-03, .31e-03);
      linename.Cell.push_back(temp_element);
    }
    else
      {
	temp_element=new DRIFT(temp_name,temp_l);
	linename.Cell.push_back(temp_element);
      }
  }
  linename.Ncell=linename.Cell.size();
  linename.Update();
  GP.circumference = linename.Cell[linename.Ncell-1]->S;
  f1.close();
}

void  Print_MADX_Lattice(Line & linename, const char * filename)
{
  int i, j, flag;
  fstream f2;
  char str[125];
  vector <string> namelist;
  
  f2.open(filename,ios::out);
  if(!f2)
    {
      cout<<"error in opening the file: "<<filename<<endl;
      exit(0);
    }
  
  for(i=0;i<linename.Ncell; i++){
    flag=0;
    for(j=0;j<namelist.size();j++) {
      if(linename.Cell[i]->NAME  == namelist [j] ) {
	flag=1;break;
      }
    }
    if(flag==0) {
      namelist.push_back(linename.Cell[i]->NAME);
      if(linename.Cell[i]->TYPE==string("DRIFT") ){
      	f2<<linename.Cell[i]->NAME<<" : DRIFT, L = "<<setprecision(16)<<linename.Cell[i]->L<<";"<<endl; }
      else if(linename.Cell[i]->TYPE==string("QUAD") )  {
	f2<<linename.Cell[i]->NAME<<" : QUADRUPOLE, L = "<<setprecision(16)<<linename.Cell[i]->L<<", K1 =  ";
        f2<<setprecision(10)<<linename.Cell[i]->GetP("K1L")/ linename.Cell[i]->L<<" ;"<<endl;}
      else if(linename.Cell[i]->TYPE==string("SEXT") ) {
	f2<<linename.Cell[i]->NAME<<" : SEXTUPOLE,  L = "<<setprecision(16)<<linename.Cell[i]->L<<", K2 =  ";
        f2<<setprecision(10)<<linename.Cell[i]->GetP("K2L")/ linename.Cell[i]->L<<" ;"<<endl; }
      else if(linename.Cell[i]->TYPE==string("OCT") ) {
	f2<<linename.Cell[i]->NAME<<" : OCTUPOLE,  L = "<<setprecision(16)<<linename.Cell[i]->L<<", K3 =  ";
        f2<<setprecision(10)<<linename.Cell[i]->GetP("KL")/ linename.Cell[i]->L<<" ;"<<endl; }
      else if(linename.Cell[i]->TYPE==string("SBEND") ) { 
	f2<<linename.Cell[i]->NAME<<" : SBEND,      L = "<<setprecision(16)<<linename.Cell[i]->L<<", ANGLE =  ";
        f2<<setprecision(10)<<linename.Cell[i]->GetP("ANGLE")<<" ; "<< endl;
      }
      else if(linename.Cell[i]->TYPE==string("RFCAV") ) {  
	f2<<linename.Cell[i]->NAME<<" : RFCAVITY,   L = "<<setprecision(16)<<linename.Cell[i]->L<<" ;"<<endl; }
      else if( linename.Cell[i]->TYPE==string("MULT")  ) {
	if( abs(linename.Cell[i]->L) > 1.0e-6 ){
	  cout<<"Warning: "<< " Multipole "<<linename.Cell[i]->NAME<<" has non-zero length. "<<endl;
	}
	f2<<linename.Cell[i]->NAME<<" : MULTIPOLE,  "<<endl;
	f2<<" knl:={ "<<linename.Cell[i]->GetP("K0L")<<","<<linename.Cell[i]->GetP("K1L")<<","<<linename.Cell[i]->GetP("K2L")<<",";
        f2<<linename.Cell[i]->GetP("K3L")<<","<<linename.Cell[i]->GetP("K4L")<<","<<linename.Cell[i]->GetP("K5L")<<",";  
	f2<<linename.Cell[i]->GetP("K6L")<<","<<linename.Cell[i]->GetP("K7L")<<","<<linename.Cell[i]->GetP("K8L")<<",";
        f2<<linename.Cell[i]->GetP("K9L")<<","<<linename.Cell[i]->GetP("K10L")<<"},"<<endl;
	f2<<" ksl:={ "<<linename.Cell[i]->GetP("K0SL")<<","<<linename.Cell[i]->GetP("K1SL")<<","<<linename.Cell[i]->GetP("K2SL")<<",";
        f2<<linename.Cell[i]->GetP("K3SL")<<","<<linename.Cell[i]->GetP("K4SL")<<","<<linename.Cell[i]->GetP("K5SL")<<",";  
	f2<<linename.Cell[i]->GetP("K6SL")<<","<<linename.Cell[i]->GetP("K7SL")<<","<<linename.Cell[i]->GetP("K8SL")<<",";
        f2<<linename.Cell[i]->GetP("K9SL")<<","<<linename.Cell[i]->GetP("K10SL")<<"};"<<endl; }
      else if( linename.Cell[i]->TYPE==string("GMULT")  ) {   //   DIPOLE  with k1l and k2l
	if( abs(linename.Cell[i]->L) > 1.0e-6 ){
	  cout<<"Warning: "<< " GMult "<<linename.Cell[i]->NAME<<" has non-zero length. "<<endl;
	}
	f2<<linename.Cell[i]->NAME<<" : SBEND,      L = "<<setprecision(16)<<linename.Cell[i]->L<<", ANGLE =  "<<setprecision(10)<<linename.Cell[i]->GetP("ANGLE")<<" ,";
	f2<<" K1 = "<<linename.Cell[i]->GetP("K1L")/linename.Cell[i]->L<<","<<" K2 = "<<linename.Cell[i]->GetP("K2L")/linename.Cell[i]->L <<";"<<endl;
      }
      else if(linename.Cell[i]->TYPE==string("MARKER")  and linename.Cell[i]->L==0.   ) {  
	f2<<linename.Cell[i]->NAME<<" : MARKER;"<<endl; }
      else if(linename.Cell[i]->TYPE==string("BPM") ) {  
	f2<<linename.Cell[i]->NAME<<" : MONITOR, L="<<setprecision(16)<<linename.Cell[i]->L<<";"<<endl; }
      else if(linename.Cell[i]->TYPE==string("HBPM") ) {  
	f2<<linename.Cell[i]->NAME<<" : HMONITOR, L="<<setprecision(16)<<linename.Cell[i]->L<<";"<<endl; }
      else if(linename.Cell[i]->TYPE==string("VBPM") ) {  
	f2<<linename.Cell[i]->NAME<<" : VMONITOR, L="<<setprecision(16)<<linename.Cell[i]->L<<";"<<endl; }
      else if(linename.Cell[i]->TYPE==string("KICKER") ) {  
	f2<<linename.Cell[i]->NAME<<" : KICKER, L="<<setprecision(16)<<linename.Cell[i]->L<<","<<"HKICK="<<linename.Cell[i]->GetP("HKICK")<<",";
        f2<<"VKICK="<<linename.Cell[i]->GetP("VKICK")<<";"<<endl; }
      else if(linename.Cell[i]->TYPE==string("HKICKER") ) {  
	f2<<linename.Cell[i]->NAME<<" : HKICKER, L="<<setprecision(16)<<linename.Cell[i]->L<<","<<"KICK="<<linename.Cell[i]->GetP("HKICK")<<";"<<endl; }
      else if(linename.Cell[i]->TYPE==string("VKICKER") ) {  
	f2<<linename.Cell[i]->NAME<<" : VKICKER, L="<<setprecision(16)<<linename.Cell[i]->L<<","<<"KICK="<<linename.Cell[i]->GetP("VKICK")<<";"<<endl; }
      else if(linename.Cell[i]->TYPE==string("BEAMBEAM") ) {  
	f2<<linename.Cell[i]->NAME<<" : BEAMBEAM;"<<endl; }
      else {
	f2<<linename.Cell[i]->NAME<<" : DRIFT, L = "<<setprecision(16)<<linename.Cell[i]->L<<";"<<endl;
        cout<<"Warning: "<< linename.Cell[i]->NAME<<"  "<<linename.Cell[i]->TYPE<<"   "<<setprecision(10)<<linename.Cell[i]->L;
        cout<<"  transferred to DRIFT ."<<endl;
      }
    }
  }
  
  int istart=0;
  int inumber10=0; 
  int nline=linename.Ncell / 1000;

  for(i=0;i<nline+1;i++) {
    sprintf(str,"LIN%d",i);
    f2<<str<<" : LINE=("<<endl;
    inumber10=0;

    do{
      f2<<linename.Cell[istart]->NAME<<",";
      istart++;
      inumber10++;
      if (inumber10 ==10){ f2<< endl; inumber10=0;}
    } while ( istart < 1000*(i+1)  &&  istart < linename.Ncell-1  );
    
    if(i < nline ){
      f2<<linename.Cell[istart]->NAME<<" ); "<<endl;
      istart++;
    }
    else{
      f2<<linename.Cell[linename.Ncell-1]->NAME<<" ); "<<endl;
    }
  }

  f2<<"rhic:  LINE = ( ";
  for(i=0;i<nline;i++) {
    sprintf(str,"LIN%d",i);
    f2<<str<<",";
  }
  sprintf(str,"LIN%d",i);
  f2<<str<<" );"<<endl;

  f2<<"beam, mass:=0.93827, charge:=1, gamma:=268.2, exn:=20.0e-06, eyn:=20.0e-06, sige:=0.001;"<<endl;
  f2<<"use, period=rhic;"<<endl;
  f2<<"select,flag=twiss, clear;"<<endl;  
  f2<<"select, flag=twiss,  column=NAME, KEYWORD,S,L,ANGLE,E1,E2,tilt, K0L,K0SL,K1L,K1SL,K2L,K2SL,K3L,K3SL,K4L,K4SL,K5L,K5SL,K6L,K6SL,K7L,K7SL, K8L,K8SL,K9L,K9SL,K10L,K10SL, KS;"<<endl;
  f2<<"twiss,table=twiss,file=parameters_input;"<<endl;
  f2<<"stop;"<<endl;
  f2.close();
}

void  Print_MADX_Lattice_v2(Line & linename, const char * filename)
{
  int i, j, flag;
  fstream f2;
  char str[125];
  vector <string> namelist;
  
  f2.open(filename,ios::out);
  if(!f2)
    {
      cout<<"error in opening the file: "<<filename<<endl;
      exit(0);
    }
  
  for(i=0;i<linename.Ncell; i++){
    flag=0;
    for(j=0;j<namelist.size();j++) {
      if(linename.Cell[i]->NAME  == namelist [j] ) {
	flag=1;break;
      }
    }
    if(flag==0) {
      namelist.push_back(linename.Cell[i]->NAME);
      if(linename.Cell[i]->TYPE==string("DRIFT") ){
      	f2<<linename.Cell[i]->NAME<<" : DRIFT, L = "<<setprecision(16)<<linename.Cell[i]->L<<";"<<endl; }
      else if(linename.Cell[i]->TYPE==string("QUAD") )  {
	f2<<linename.Cell[i]->NAME<<" : QUADRUPOLE, L = "<<setprecision(16)<<linename.Cell[i]->L<<", K1 =  ";
        f2<<setprecision(10)<<linename.Cell[i]->GetP("K1L")/ linename.Cell[i]->L<<" ;"<<endl;}
      else if(linename.Cell[i]->TYPE==string("SEXT") ) {
	f2<<linename.Cell[i]->NAME<<" : SEXTUPOLE, L = "<<setprecision(16)<<linename.Cell[i]->L<<", K2 =  ";
        f2<<setprecision(10)<<linename.Cell[i]->GetP("K2L")/ linename.Cell[i]->L<<" ;"<<endl; }
      else if(linename.Cell[i]->TYPE==string("OCT") ) {
	f2<<linename.Cell[i]->NAME<<" : OCTUPOLE, L = "<<setprecision(16)<<linename.Cell[i]->L<<", K3 =  ";
        f2<<setprecision(10)<<linename.Cell[i]->GetP("K3L")/ linename.Cell[i]->L<<" ;"<<endl; }
      else if(linename.Cell[i]->TYPE==string("SBEND") ) { 
	f2<<linename.Cell[i]->NAME<<" : SBEND, L = "<<setprecision(16)<<linename.Cell[i]->L
	  <<", ANGLE =  "<<setprecision(10)<<linename.Cell[i]->GetP("ANGLE")
	  <<", E1 =  "   <<setprecision(10)<<linename.Cell[i]->GetP("E1")
	  <<", E2 =  "   <<setprecision(10)<<linename.Cell[i]->GetP("E2")  <<" ; "<< endl;
      }
      else if(linename.Cell[i]->TYPE==string("RFCAV") ) {  
	f2<<linename.Cell[i]->NAME<<" : RFCAVITY,   L = "<<setprecision(16)<<linename.Cell[i]->L<<" ;"<<endl; }
      else if( linename.Cell[i]->TYPE==string("MULT")  ) {   //  NOTE: MADX only takes  zero length multipoles
	if( abs(linename.Cell[i]->L) > 1.0e-6 ){
	  cout<<"Warning: "<< " Multipole "<<linename.Cell[i]->NAME<<" has non-zero length. "<<endl;
	}
	f2<<linename.Cell[i]->NAME<<" : MULTIPOLE,  "<<endl;
	f2<<" knl:={ "<<linename.Cell[i]->GetP("K0L")<<","<<linename.Cell[i]->GetP("K1L")<<","<<linename.Cell[i]->GetP("K2L")<<",";
        f2<<linename.Cell[i]->GetP("K3L")<<","<<linename.Cell[i]->GetP("K4L")<<","<<linename.Cell[i]->GetP("K5L")<<",";  
	f2<<linename.Cell[i]->GetP("K6L")<<","<<linename.Cell[i]->GetP("K7L")<<","<<linename.Cell[i]->GetP("K8L")<<",";
        f2<<linename.Cell[i]->GetP("K9L")<<","<<linename.Cell[i]->GetP("K10L")<<"},"<<endl;
	f2<<" ksl:={ "<<linename.Cell[i]->GetP("K0SL")<<","<<linename.Cell[i]->GetP("K1SL")<<","<<linename.Cell[i]->GetP("K2SL")<<",";
        f2<<linename.Cell[i]->GetP("K3SL")<<","<<linename.Cell[i]->GetP("K4SL")<<","<<linename.Cell[i]->GetP("K5SL")<<",";  
	f2<<linename.Cell[i]->GetP("K6SL")<<","<<linename.Cell[i]->GetP("K7SL")<<","<<linename.Cell[i]->GetP("K8SL")<<",";
        f2<<linename.Cell[i]->GetP("K9SL")<<","<<linename.Cell[i]->GetP("K10SL")<<"};"<<endl; }
      else if( linename.Cell[i]->TYPE==string("GMULT")  ) {   //   NOTE:   only pass  k1, k2  for sbend so that MADX can read in
	if( abs(linename.Cell[i]->L) > 1.0e-6 ){
	  cout<<"Warning: "<< " GMult "<<linename.Cell[i]->NAME<<" transferred to SBEND with K1  annd K2. "<<endl;
	}
	f2<<linename.Cell[i]->NAME<<" : SBEND,      L = "<<setprecision(16)<<linename.Cell[i]->L<<", ANGLE =  "<<setprecision(10)<<linename.Cell[i]->GetP("ANGLE")<<" ,";
	f2<<" K1 = "<<linename.Cell[i]->GetP("K1L")/linename.Cell[i]->L<<","<<" K2 = "<<linename.Cell[i]->GetP("K2L")/linename.Cell[i]->L <<";"<<endl;
      }
      else if(linename.Cell[i]->TYPE==string("MARKER") and linename.Cell[i]->L == 0.  ) {  
	f2<<linename.Cell[i]->NAME<<" : MARKER;"<<endl; }
      else if(linename.Cell[i]->TYPE==string("BPM") and linename.Cell[i]->L == 0.   ) {  
	f2<<linename.Cell[i]->NAME<<" : MONITOR, L="<<setprecision(16)<<linename.Cell[i]->L<<";"<<endl; }
      else if(linename.Cell[i]->TYPE==string("HBPM") and  linename.Cell[i]->L == 0.  ) {  
	f2<<linename.Cell[i]->NAME<<" : HMONITOR, L="<<setprecision(16)<<linename.Cell[i]->L<<";"<<endl; }
      else if(linename.Cell[i]->TYPE==string("VBPM") and  linename.Cell[i]->L == 0.  ) {  
	f2<<linename.Cell[i]->NAME<<" : VMONITOR, L="<<setprecision(16)<<linename.Cell[i]->L<<";"<<endl; }
      else if(linename.Cell[i]->TYPE==string("KICKER") ) {  
	f2<<linename.Cell[i]->NAME<<" : TKICKER, L="<<setprecision(16)<<linename.Cell[i]->L<<","<<"HKICK="<<linename.Cell[i]->GetP("HKICK")<<",";
        f2<<"VKICK="<<linename.Cell[i]->GetP("VKICK")<<";"<<endl; }
      else if(linename.Cell[i]->TYPE==string("HKICKER")  ) {  
	f2<<linename.Cell[i]->NAME<<" : HKICKER, L="<<setprecision(16)<<linename.Cell[i]->L<<","<<"KICK="<<linename.Cell[i]->GetP("HKICK")<<";"<<endl; }
      else if(linename.Cell[i]->TYPE==string("VKICKER")  ) {  
	f2<<linename.Cell[i]->NAME<<" : VKICKER, L="<<setprecision(16)<<linename.Cell[i]->L<<","<<"KICK="<<linename.Cell[i]->GetP("VKICK")<<";"<<endl; }
      else if(linename.Cell[i]->TYPE==string("BEAMBEAM") ) {  
	f2<<linename.Cell[i]->NAME<<" : BEAMBEAM;"<<endl; }
      else {
	f2<<linename.Cell[i]->NAME<<" : DRIFT, L = "<<setprecision(16)<<linename.Cell[i]->L<<";"<<endl;
        cout<<"Warning: "<< linename.Cell[i]->NAME<<"  "<<linename.Cell[i]->TYPE<<"   "<<setprecision(10)<<linename.Cell[i]->L;
        cout<<"  transferred to DRIFT ."<<endl;
      }
    }
  }
  
  int istart=0;
  int inumber10=0; 
  int nline=linename.Ncell / 1000;

  for(i=0;i<nline+1;i++) {
    sprintf(str,"LIN%d",i);
    f2<<str<<" : LINE=("<<endl;
    inumber10=0;

    do{
      f2<<linename.Cell[istart]->NAME<<",";
      istart++;
      inumber10++;
      if (inumber10 ==10){ f2<< endl; inumber10=0;}
    } while ( istart < 1000*(i+1)  &&  istart < linename.Ncell-1  );
    
    if(i < nline ){
      f2<<linename.Cell[istart]->NAME<<" ); "<<endl;
      istart++;
    }
    else{
      f2<<linename.Cell[linename.Ncell-1]->NAME<<" ); "<<endl;
    }
  }

  f2<<"rhic:  LINE = ( ";
  for(i=0;i<nline;i++) {
    sprintf(str,"LIN%d",i);
    f2<<str<<",";
  }
  sprintf(str,"LIN%d",i);
  f2<<str<<" );"<<endl;

  f2<<"beam, mass:=0.93827, charge:=1, gamma:=268.2, exn:=20.0e-06, eyn:=20.0e-06, sige:=0.001;"<<endl;
  f2<<"use, period=rhic;"<<endl;
  f2<<"select,flag=twiss, clear;"<<endl;  
  f2<<"select, flag=twiss,  column=NAME, KEYWORD,S,L,ANGLE,E1,E2,tilt, K0L,K0SL,K1L,K1SL,K2L,K2SL,K3L,K3SL,K4L,K4SL,K5L,K5SL,K6L,K6SL,K7L,K7SL, K8L,K8SL,K9L,K9SL,K10L,K10SL, KS;"<<endl;
  f2<<"twiss,table=twiss,file=parameters_input;"<<endl;
  f2<<"stop;"<<endl;
  f2.close();
}

void  Read_BMAD_Lattice(Line & linename, const char * filename)
{

  int i,m;
  fstream f1;
  string line;
  string temp_name;
  string temp_type;
  double temp_s;
  double temp_l;
  double temp_angle;
  double temp_angle1;
  double temp_e1;
  double temp_e2;
  double temp_tilt;
  double temp_xoff;
  double temp_zoff;  
  double temp_xpitch;
  double knl[11], knsl[11],k1l_mult, k2l_mult;
  double temp_sk;
  double hkick=0;
  double vkick=0;
  Element * temp_element;
  double kl;

  f1.open(filename,ios::in);
  if(!f1)
    {
      cout<<"error in opening the file: "<<filename<<endl;
      exit(0);
    }
  for(i=0;i<1;i++) getline(f1,line,'\n');
  
  double s1=0;
  while( getline(f1,line,'\n')){
    istringstream ss(line);
    ss>>temp_name>>temp_type>>temp_s>>temp_l>>temp_angle>>temp_angle1>>temp_e1>>temp_e2>>temp_tilt>>temp_xoff>>temp_zoff>>temp_xpitch
      >>knl[0]>>knsl[0]
      >>knl[1]>>knsl[1]      
      >>knl[2]>>knsl[2]      
      >>knl[3]>>knsl[3]
      >>knl[4]>>knsl[4]
      >>knl[5]>>knsl[5]
      >>knl[6]>>knsl[6]
      >>knl[7]>>knsl[7]
      >>knl[8]>>knsl[8]
      >>knl[9]>>knsl[9]
      >>knl[10]>>knsl[10]
      >>temp_sk>>hkick>>vkick;
    
    s1=s1+temp_l;

    if(temp_type==string("Drift")){
      temp_element=new DRIFT(temp_name,temp_l);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("Sbend"))    {

      if( abs(temp_angle) < 1.0e-9 ) {             // treated as mult
	temp_element=new MULT(temp_name,temp_l, knl, knsl);
	linename.Cell.push_back(temp_element);
      }
      else if( abs(temp_angle -  knl[0] )  >  1.0e-6 ) {   //  treated as gsbend
	temp_element= new GSBEND(temp_name, temp_l, temp_angle, knl[0], temp_e1, temp_e2);
	linename.Cell.push_back(temp_element);
      }
      else {
	knl[0]=0.;
	if( knl[2] == 0. and  knl[1] == 0. ){  // treated  as  sbend
	  temp_element= new SBEND(temp_name,temp_l,temp_angle, temp_e1, temp_e2);
	  linename.Cell.push_back(temp_element);
	}
	else{                                  // treated  as gmult
	  temp_element=new SBENDMULT(temp_name,temp_l, temp_angle, temp_e1, temp_e2, knl, knsl);
	  linename.Cell.push_back(temp_element);
	}
      }
    }
    else if (temp_type==string("RBEND"))    {
      if(knl[1] == 0. and knl[2]==0.) {
	temp_element= new SBEND(temp_name,temp_l,temp_angle, temp_angle/2., temp_angle/2.);
	linename.Cell.push_back(temp_element); }
      else {
	temp_element=new SBENDMULT(temp_name,temp_l, temp_angle, temp_angle/2., temp_angle/2., knl, knsl);
	linename.Cell.push_back(temp_element);
      }
    }
    else if (temp_type==string("Quadrupole")){  
      temp_element= new QUAD(temp_name,temp_l, knl[1], knsl[1]);
      linename.Cell.push_back(temp_element); 
    }
    else if (temp_type==string("Sextupole")){
      temp_element= new SEXT(temp_name,temp_l, knl[2], knsl[2]);
      linename.Cell.push_back(temp_element); 
    }
    else if (temp_type==string("OCTUPOLE")){    
      temp_element= new OCT(temp_name,temp_l,knl[3],knsl[3]);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("HKICKER")){
      temp_element= new HKICKER(temp_name,temp_l, hkick );
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("VKICKER")){
      temp_element= new VKICKER(temp_name,temp_l, vkick);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("Kicker")){
      temp_element= new KICKER(temp_name,temp_l, hkick, vkick);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("Multipole")){   // In madx:  k0l-->angle, hkick, vkick --> correctors
      //kl = 0.;
      //for(m=0;m<11;m++)  kl=kl+abs(knl[m])+abs(knsl[m]);
      if(knl[0] == 0  ){
	temp_element=new MULT(temp_name,temp_l, knl, knsl);
	linename.Cell.push_back(temp_element);
      }
      else{
	cout<<"   found  one  multipole "<<temp_name<<"  with angle, changed to  SBEND. "<<endl;
	temp_angle = knl[0];  	knl[0]=0;   temp_l=1.0e-6;
	temp_element=new SBEND(temp_name+"_tb",temp_l,temp_angle, 0, 0);
	linename.Cell.push_back(temp_element);
      }	
    }
    else if (temp_type==string("HMONITOR")){
      temp_element= new HBPM(temp_name,temp_l);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("VMONITOR")){
      temp_element= new VBPM(temp_name,temp_l);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("Monitor")){
      temp_element= new BPM(temp_name,temp_l);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("MARKER")){
      temp_element= new MARKER(temp_name,temp_l);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("Solenoid")){
      //if(temp_sk !=0.0  and temp_l !=0. ){    // meaningful input for solenoid
      if(temp_l !=0. ) {
	temp_element= new SOLEN(temp_name,temp_l, temp_sk);
	linename.Cell.push_back(temp_element);
      }
      else{
	temp_element=new DRIFT(temp_name,temp_l);
	linename.Cell.push_back(temp_element);
      }
    }
    else if (temp_type==string("RFcavity")){ 
      temp_element= new RFCAV(temp_name,temp_l, 0.0, 0.0, 0.);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("BEAMBEAM")){
      if(temp_l !=0.){ 
	cout<<"Beam-Beam element should have zero length!"<<endl;
        exit(1);  
      }
      temp_element= new BEAMBEAM(temp_name,6, 0., 1, 0.4545, 11, 2.5e-06, 2.5e-06, 0.53, 0.0, 0.53, 0.0);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("ELENS")){
      temp_element= new ELENS(temp_name, temp_l, 0.0, 1, 0.14, 4, .31e-03, .31e-03);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("Patch")){
      temp_element= new PATCH(temp_name, temp_xoff, 0.0, temp_zoff, -temp_xpitch );
      linename.Cell.push_back(temp_element);   
    }
    else
      {
	temp_element=new DRIFT(temp_name,temp_l);
	linename.Cell.push_back(temp_element);
      }
  }
  linename.Ncell=linename.Cell.size();
  linename.Update();
  GP.circumference = linename.Cell[linename.Ncell-1]->S;
  f1.close();
}


void  Print_Parameters_Input( Line & linename, const char * filename)
{
  int  i;
  fstream f2;
  
  f2.open(filename,ios::out);
  if(!f2)
    {
      cout<<"error in opening the file: "<<filename<<endl;
      exit(0);
    }

  for(i=0;i<47;i++)
    f2<<"..."<<endl;

  for(i=0;i<linename.Ncell;i++){
    if(linename.Cell[i]->TYPE==string("DRIFT") ){ 
      f2<<"\""<<linename.Cell[i]->NAME<<"\" "<<"\"DRIFT\" "<<scientific<<linename.Cell[i]->S<<" "<<linename.Cell[i]->L
	<<" 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "<<endl;
    } 
    else if(linename.Cell[i]->TYPE==string("SBEND") ){
      f2<<"\""<<linename.Cell[i]->NAME<<"\" "<<"\"SBEND\" "<<scientific<<linename.Cell[i]->S<<" "<<linename.Cell[i]->L
	<<" "<<linename.Cell[i]->GetP("ANGLE")<<" "<<linename.Cell[i]->GetP("E1")<<" "<<linename.Cell[i]->GetP("E2")<<" "
        <<" 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"<<endl;
    }
   else if(linename.Cell[i]->TYPE==string("QUAD") ){
      f2<<"\""<<linename.Cell[i]->NAME<<"\" "<<"\"QUADRUPOLE\" "<<scientific<<linename.Cell[i]->S<<" "<<linename.Cell[i]->L
	<<" 0 0 0 0 0 0 "<<linename.Cell[i]->GetP("K1L")<<" "<<linename.Cell[i]->GetP("K1SL")
        <<" 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "<<endl;
   }
   else if(linename.Cell[i]->TYPE==string("SKEWQ") ){
      f2<<"\""<<linename.Cell[i]->NAME<<"\" "<<"\"QUADRUPOLE\" "<<scientific<<linename.Cell[i]->S<<" "<<linename.Cell[i]->L
	<<" 0 0 0 0 0 0 "<<linename.Cell[i]->GetP("K1L")<<" "<<linename.Cell[i]->GetP("K1SL")
        <<" 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "<<endl;
   }
   else if(linename.Cell[i]->TYPE==string("SEXT") ){
      f2<<"\""<<linename.Cell[i]->NAME<<"\" "<<"\"SEXTUPOLE\" "<<scientific<<linename.Cell[i]->S<<" "<<linename.Cell[i]->L
	<<" 0 0 0 0 0 0 0 0 "<<linename.Cell[i]->GetP("K2L")<<" "<<linename.Cell[i]->GetP("K2SL")
        <<" 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "<<endl;
   }   
   else if(linename.Cell[i]->TYPE==string("OCT") ){
      f2<<"\""<<linename.Cell[i]->NAME<<"\" "<<"\"OCTUPOLE\" "<<scientific<<linename.Cell[i]->S<<" "<<linename.Cell[i]->L
	<<" 0 0 0 0 0 0 0 0 0 0 "<<linename.Cell[i]->GetP("K3L")<<" "<<linename.Cell[i]->GetP("K3SL")
        <<" 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "<<endl;
   } 
   else if(linename.Cell[i]->TYPE==string("MULT") ){
      f2<<"\""<<linename.Cell[i]->NAME<<"\" "<<"\"MULTIPOLE\" "<<scientific<<linename.Cell[i]->S<<" "<<linename.Cell[i]->L
	<<" 0 0 0 0 "
        <<linename.Cell[i]->GetP("K0L") <<" "<<linename.Cell[i]->GetP("K0SL")<<" "
        <<linename.Cell[i]->GetP("K1L") <<" "<<linename.Cell[i]->GetP("K1SL")<<" "
        <<linename.Cell[i]->GetP("K2L") <<" "<<linename.Cell[i]->GetP("K2SL")<<" "
        <<linename.Cell[i]->GetP("K3L") <<" "<<linename.Cell[i]->GetP("K3SL")<<" "
        <<linename.Cell[i]->GetP("K4L") <<" "<<linename.Cell[i]->GetP("K4SL")<<" "
        <<linename.Cell[i]->GetP("K5L") <<" "<<linename.Cell[i]->GetP("K5SL")<<" "
        <<linename.Cell[i]->GetP("K6L") <<" "<<linename.Cell[i]->GetP("K6SL")<<" "
        <<linename.Cell[i]->GetP("K7L") <<" "<<linename.Cell[i]->GetP("K7SL")<<" "
        <<linename.Cell[i]->GetP("K8L") <<" "<<linename.Cell[i]->GetP("K8SL")<<" "
        <<linename.Cell[i]->GetP("K9L") <<" "<<linename.Cell[i]->GetP("K9SL")<<" "
        <<linename.Cell[i]->GetP("K10L")<<" "<<linename.Cell[i]->GetP("K10SL")<<" 0 "<<endl;
   } 
   else if(linename.Cell[i]->TYPE==string("BPM") ){
      f2<<"\""<<linename.Cell[i]->NAME<<"\" "<<"\"MONITOR\" "<<scientific<<linename.Cell[i]->S<<" "<<linename.Cell[i]->L
       <<" 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "<<endl;
   } 
   else if(linename.Cell[i]->TYPE==string("HBPM") ){
      f2<<"\""<<linename.Cell[i]->NAME<<"\" "<<"\"HMONITOR\" "<<scientific<<linename.Cell[i]->S<<" "<<linename.Cell[i]->L
       <<" 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "<<endl;
   } 
   else if(linename.Cell[i]->TYPE==string("VBPM") ){
     f2<<"\""<<linename.Cell[i]->NAME<<"\" "<<"\"VMONITOR\" "<<scientific<<linename.Cell[i]->S<<" "<<linename.Cell[i]->L
       <<" 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "<<endl;
   }
   else if(linename.Cell[i]->TYPE==string("MARKER") ){
      f2<<"\""<<linename.Cell[i]->NAME<<"\" "<<"\"MARKER\" "<<scientific<<linename.Cell[i]->S<<" "<<linename.Cell[i]->L
       <<" 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "<<endl;
   } 
   else if(linename.Cell[i]->TYPE==string("KICKER") ){
      f2<<"\""<<linename.Cell[i]->NAME<<"\" "<<"\"KICKER\" "<<scientific<<linename.Cell[i]->S<<" "<<linename.Cell[i]->L
       <<" 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "<<endl;
   } 
   else if(linename.Cell[i]->TYPE==string("HKICKER") ){
      f2<<"\""<<linename.Cell[i]->NAME<<"\" "<<"\"HKICKER\" "<<scientific<<linename.Cell[i]->S<<" "<<linename.Cell[i]->L
       <<" 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "<<endl;
   } 
   else if(linename.Cell[i]->TYPE==string("VKICKER") ){
      f2<<"\""<<linename.Cell[i]->NAME<<"\" "<<"\"VKICKER\" "<<scientific<<linename.Cell[i]->S<<" "<<linename.Cell[i]->L
       <<" 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "<<endl;
   }
   else if(linename.Cell[i]->TYPE==string("RFCAV") ){
      f2<<"\""<<linename.Cell[i]->NAME<<"\" "<<"\"RFCAVITY\" "<<scientific<<linename.Cell[i]->S<<" "<<linename.Cell[i]->L
       <<" 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "<<endl;
   }
   else if(linename.Cell[i]->TYPE==string("BEAMBEAM") ){
      f2<<"\""<<linename.Cell[i]->NAME<<"\" "<<"\"BEAMBEAM\" "<<scientific<<linename.Cell[i]->S<<" "<<linename.Cell[i]->L
       <<" 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "<<endl;
   }
   else{
     cout<<"Warning: "<<linename.Cell[i]->NAME<<" at "<<linename.Cell[i]->S<<"  is transferred to DRIFT "<<endl;
     f2<<"\""<<linename.Cell[i]->NAME<<"\" "<<"\"DRIFT\" "<<scientific<<linename.Cell[i]->S<<" "<<linename.Cell[i]->L
       <<" 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "<<endl;
   }
  }
  
  f2.close();
}

void Tokenize(const string& str, vector<string>& tokens)
{
    string  delimiters=" ";
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    string::size_type pos     = str.find_first_of(delimiters, lastPos);
   

    while (string::npos != pos || string::npos != lastPos)
    {
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = str.find_first_not_of(delimiters, pos);
        pos = str.find_first_of(delimiters, lastPos);
    }
}

void Read_LatticeSequence(Line & linename, const char * filename)
{
  int i,j;
  fstream f1;
  string input0;
  vector<string> tokens, input1;
  double l,angle, e1,e2, k1l, k1sl, k2l, k2sl, k3l, k3sl, ks;
  double knl[11], knsl[11];
  double xin[6], xout[6], m66[36];
  double hkick, vkick;
  double klmax, kslmax,phi0;
  double hkickmax, vkickmax, nud, phid;  
  double dx, dpx, dy, dpy, tilt;
  double vrf, frf, phase0;
  double np, bbscale, sigmal,emitx, emity,betax, alfax, betay, alfay;
  double ne, betae,sigmax, sigmay,sepx,sepy;
  double n[3];
  double diff_x, diff_y,diff_delta;
  double alfa;

  int norder, tturns;
  int turns, turne;
  int treatment, nslice;


  Element * temp_element;

  f1.open(filename,ios::in);
  if(!f1)
    {
      cout<<"error in opening the file: "<<filename<<endl;
      exit(0);
    }

  getline(f1, input0,'\n');

  input1.clear(); 
  while(getline(f1, input0,'\n') ) {
    
    tokens.clear();  Tokenize(input0, tokens);
    if( tokens[ tokens.size()-1 ] == "," ) {
      for(i=0;i< tokens.size()-1;i++) input1.push_back( tokens[i]);
      }
    else{
      for(i=0;i< tokens.size();i++) input1.push_back( tokens[i]);
      
      if(input1[1]==string("DRIFT")){
	l= atof(input1[2].c_str());
	temp_element=new DRIFT(input1[0],l);
	linename.Cell.push_back(temp_element);
      }
      else if(input1[1]==string("SBEND")){
	l= atof(input1[2].c_str());     angle=atof(input1[3].c_str()); 
        e1=atof(input1[4].c_str());     e2=atof(input1[5].c_str()); 
	temp_element= new SBEND( input1[0],l,angle,e1,e2);
	linename.Cell.push_back(temp_element);
      }
      else if(input1[1]==string("QUAD")){
	l= atof(input1[2].c_str());          
        k1l=atof(input1[3].c_str());    k1sl=atof(input1[4].c_str()); 
	temp_element= new QUAD( input1[0],l,k1l,k1sl);
	linename.Cell.push_back(temp_element);
      }
      else if(input1[1]==string("SKEWQ")){
	l= atof(input1[2].c_str());          
        k1sl=atof(input1[3].c_str()); 
	temp_element= new SKEWQ( input1[0], l,k1sl);
	linename.Cell.push_back(temp_element);
      }
      else if(input1[1]==string("SEXT")){
	l= atof(input1[2].c_str());          
        k2l=atof(input1[3].c_str());    k2sl=atof(input1[4].c_str()); 
	temp_element= new SEXT( input1[0],l,k2l,k2sl);
	linename.Cell.push_back(temp_element);
      }
      else if(input1[1]==string("OCT")){
	l= atof(input1[2].c_str());          
        k3l=atof(input1[3].c_str());    k3sl=atof(input1[4].c_str()); 
	temp_element= new OCT( input1[0],l,k3l,k3sl);
	linename.Cell.push_back(temp_element);
      }
      else if(input1[1]==string("MULT")){
	l= atof(input1[2].c_str());          
        for(j=0;j<11;j++) knl[j]= atof(input1[3+j].c_str());
	for(j=0;j<11;j++) knsl[j]= atof(input1[14+j].c_str());
	temp_element=new MULT(input1[0],l,knl,knsl);
	linename.Cell.push_back(temp_element);
      } 
      else if(input1[1]==string("GMULT")){
	l= atof(input1[2].c_str());          angle=atof(input1[3].c_str()); 
        e1=atof(input1[4].c_str());          e2=atof(input1[5].c_str()); 
        for(j=0;j<11;j++) knl[j]= atof(input1[6+j].c_str());
	for(j=0;j<11;j++) knsl[j]= atof(input1[17+j].c_str());
	temp_element=new GMULT(input1[0],l,angle,e1,e2,knl, knsl);
	linename.Cell.push_back(temp_element);
      }
      else if(input1[1]==string("SMULT")){
	l= atof(input1[2].c_str());          angle=atof(input1[3].c_str()); 
        e1=atof(input1[4].c_str());          e2=atof(input1[5].c_str()); 
        for(j=0;j<11;j++) knl[j]= atof(input1[6+j].c_str());
	for(j=0;j<11;j++) knsl[j]= atof(input1[17+j].c_str());
	temp_element=new SMULT(input1[0],l,angle,e1,e2,knl, knsl);
	linename.Cell.push_back(temp_element);
      }      
      else if(input1[1]==string("SOLEN")){
	l= atof(input1[2].c_str());    ks=atof(input1[3].c_str()); 
        if(ks !=0. ){
	  temp_element=new SOLEN(input1[0],l, ks);
	  linename.Cell.push_back(temp_element);
	}
	else{
          cout<<"transfer a zero-streng SOLEN to DRIFT. "<<endl;
	  temp_element=new DRIFT(input1[0],l);
	  linename.Cell.push_back(temp_element);
	}
      }      
      else if(input1[1]==string("MATRIX")){
	l= atof(input1[2].c_str());   
        for(j=0;j<6;j++)  xin[j]=atof(input1[3+j].c_str());
        for(j=0;j<6;j++)  xout[j]=atof(input1[9+j].c_str());
        for(j=0;j<36;j++) m66[j]=atof(input1[15+j].c_str());
	temp_element=new MATRIX(input1[0],l, xin,xout,m66);
	linename.Cell.push_back(temp_element);
      }
      else if(input1[1]==string("KICKER")){
	l= atof(input1[2].c_str());   
        hkick=atof(input1[3].c_str());  vkick=atof(input1[4].c_str()); 
	temp_element=new KICKER(input1[0],l,hkick,vkick);
	linename.Cell.push_back(temp_element);
      }
      else if(input1[1]==string("HKICKER")){
	l= atof(input1[2].c_str());   
        hkick=atof(input1[3].c_str()); 
	temp_element=new HKICKER(input1[0],l,hkick);
	linename.Cell.push_back(temp_element);
      }
      else if(input1[1]==string("VKICKER")){
	l= atof(input1[2].c_str());   
        vkick=atof(input1[3].c_str()); 
	temp_element=new VKICKER(input1[0],l,vkick);
	linename.Cell.push_back(temp_element);
      }
      else if(input1[1]==string("HACMULT")){
	l= atof(input1[2].c_str());  norder = atoi(input1[3].c_str());       
        klmax=atof(input1[4].c_str());   tturns=atoi(input1[5].c_str());   phi0=atof(input1[6].c_str());
	temp_element=new HACMULT(input1[0],l,norder, klmax, tturns, phi0);
	linename.Cell.push_back(temp_element);
      }
      else if(input1[1]==string("VACMULT")){
	l= atof(input1[2].c_str());  norder = atoi(input1[3].c_str());       
        kslmax=atof(input1[4].c_str());   tturns=atoi(input1[5].c_str());    phi0=atof(input1[6].c_str());
	temp_element=new VACMULT(input1[0],l,norder, kslmax, tturns, phi0);
	linename.Cell.push_back(temp_element);
      }
      else if(input1[1]==string("HACDIP")){
	l= atof(input1[2].c_str()); 
        hkickmax=atof(input1[3].c_str());   nud=atof(input1[4].c_str());    phid=atof(input1[5].c_str());
        turns=atoi(input1[6].c_str());      turne=atoi(input1[7].c_str()); 
	temp_element=new HACDIP(input1[0],l,hkickmax,nud,phid,turns,turne);
	linename.Cell.push_back(temp_element);
      }
      else if(input1[1]==string("VACDIP")){
        vkickmax=atof(input1[3].c_str());   nud=atof(input1[4].c_str());    phid=atof(input1[5].c_str());
        turns=atoi(input1[6].c_str());      turne=atoi(input1[7].c_str()); 
	temp_element=new VACDIP(input1[0],l,vkickmax,nud,phid,turns,turne);
	linename.Cell.push_back(temp_element);
      }
      else if(input1[1]==string("BPM")){
	l= atof(input1[2].c_str()); 
	temp_element=new BPM(input1[0],l);
	linename.Cell.push_back(temp_element);
      }
      else if(input1[1]==string("HBPM")){
	l= atof(input1[2].c_str()); 
	temp_element=new HBPM(input1[0],l);
	linename.Cell.push_back(temp_element);
      }
      else if(input1[1]==string("VBPM")){
	l= atof(input1[2].c_str()); 
	temp_element=new VBPM(input1[0],l);
	linename.Cell.push_back(temp_element);
      }
      else if(input1[1]==string("MARKER")){
	l= atof(input1[2].c_str()); 
	temp_element=new MARKER(input1[0],l);
	linename.Cell.push_back(temp_element);
      }
      else if(input1[1]==string("RFCAV")){
	l= atof(input1[2].c_str());   vrf=atof(input1[3].c_str());
        frf=atof(input1[4].c_str());  phase0=atof(input1[5].c_str()); 
	temp_element=new RFCAV(input1[0],l, vrf, frf, phase0);
	linename.Cell.push_back(temp_element);
      }
      else if(input1[1]==string("BEAMBEAM")){
	treatment=atoi(input1[2].c_str()); np=atof(input1[3].c_str()); bbscale = atof(input1[4].c_str());
        sigmal=atof(input1[5].c_str()); nslice=atoi(input1[6].c_str());
        emitx=atof(input1[7].c_str());  emity=atof(input1[8].c_str());
        betax=atof(input1[9].c_str());  alfax=atof(input1[10].c_str());
        betay=atof(input1[11].c_str()); alfay=atof(input1[12].c_str());
	temp_element=new BEAMBEAM(input1[0], treatment,np,bbscale,sigmal,nslice,emitx,emity,betax,alfax,betay,alfay);
	linename.Cell.push_back(temp_element);
      }
      else if(input1[1]==string("LRBB")){
	np=atof(input1[2].c_str());      bbscale= atof(input1[3].c_str());
        sepx=atof(input1[4].c_str());    sepy=atof(input1[5].c_str());
        sigmax=atof(input1[6].c_str());  sigmay=atof(input1[7].c_str());
	temp_element=new LRBB(input1[0], np,bbscale, sepx, sepy, sigmax, sigmay);
	linename.Cell.push_back(temp_element);
      }
      else if(input1[1]==string("ELENS")){
	l= atof(input1[2].c_str());       ne=atof(input1[3].c_str());     
        bbscale=atof(input1[4].c_str());  nslice=atof(input1[5].c_str());
        betae=atof(input1[6].c_str());  sigmax=atof(input1[7].c_str()); sigmay=atof(input1[8].c_str()); 
	temp_element=new ELENS(input1[0],l,ne, bbscale, nslice,betae,sigmax,sigmay);
	linename.Cell.push_back(temp_element);
      }
      else if(input1[1]==string("ERHICBB")){
	ne=atof(input1[2].c_str());  bbscale=atof(input1[3].c_str()); 
        sigmax=atof(input1[4].c_str()); sigmay=atof(input1[5].c_str()); 
	temp_element=new ERHICBB(input1[0],ne,bbscale,sigmax,sigmay);
	linename.Cell.push_back(temp_element);
      }
      else if(input1[1]==string("ROTAT")){
	l= atof(input1[2].c_str());    
        for(j=0;j<3;j++) n[j]=atof(input1[3+j].c_str());
	angle=atof(input1[6].c_str());
	temp_element=new ROTAT(input1[0],l,n,angle);
	linename.Cell.push_back(temp_element);
      }
      else if(input1[1]==string("SNAKE")){
	l= atof(input1[2].c_str());  
        for(j=0;j<3;j++) n[j]=atof(input1[3+j].c_str());
	angle=atof(input1[6].c_str());
	temp_element=new SNAKE(input1[0],l,n,angle);
	linename.Cell.push_back(temp_element);
      }
      else if(input1[1]==string("DIFFUSE")){
	diff_x= atof(input1[2].c_str()); diff_y= atof(input1[3].c_str());  diff_delta= atof(input1[4].c_str());
	temp_element=new DIFFUSE(input1[0],diff_x,diff_y,diff_delta);
	linename.Cell.push_back(temp_element);
      }
      else if(input1[1]==string("COOLING")){
	alfa= atof(input1[2].c_str());
	temp_element=new COOLING(input1[0],alfa);
	linename.Cell.push_back(temp_element);
      }
      else{
        cout<<"can't recognize ELEMENT type "<<input1[1] << endl;
        cout<<"transfer it to be a DRIFT. "<<endl;
	l= atof(input1[2].c_str());
	temp_element=new DRIFT(input1[0],l);
	linename.Cell.push_back(temp_element);
      }

      input1.clear();  
    }
    
  }
  
  linename.Ncell=linename.Cell.size();
  linename.Update();

  f1.close();
}

void Print_LatticeSequence(Line & linename, const char * filename)
{
  int i;
  fstream f2;
  
  f2.open(filename,ios::out);
  if(!f2)
    {
      cout<<"error in opening the file: "<<filename<<endl;
      exit(0);
    }
  
  f2<<"Information Line: rhic lattice "<<endl;

  for(i=0;i<linename.Ncell;i++){
    if(linename.Cell[i]->TYPE==string("DRIFT") ){
      f2<<linename.Cell[i]->NAME<<"   "<<linename.Cell[i]->TYPE<<" "<<scientific<<linename.Cell[i]->L<<endl;
    }
    else if(linename.Cell[i]->TYPE==string("SBEND")){
      f2<<linename.Cell[i]->NAME<<"   "<<linename.Cell[i]->TYPE<<" "<<scientific<<linename.Cell[i]->L<<" "
        <<linename.Cell[i]->GetP("ANGLE")<<" "<<linename.Cell[i]->GetP("E1")<<" "<<linename.Cell[i]->GetP("E2")<<endl;
    }
    else if(linename.Cell[i]->TYPE==string("QUAD")){
      f2<<linename.Cell[i]->NAME<<"   "<<linename.Cell[i]->TYPE<<" "<<scientific<<linename.Cell[i]->L<<" "
        <<linename.Cell[i]->GetP("K1L")<<" "<<linename.Cell[i]->GetP("K1SL")<<endl;
    }
    else if(linename.Cell[i]->TYPE==string("SKEWQ")){
      f2<<linename.Cell[i]->NAME<<"   "<<linename.Cell[i]->TYPE<<" "<<scientific<<linename.Cell[i]->L<<" "
        <<linename.Cell[i]->GetP("K1SL")<<endl;
    }
    else if(linename.Cell[i]->TYPE==string("SEXT")){
      f2<<linename.Cell[i]->NAME<<"   "<<linename.Cell[i]->TYPE<<" "<<scientific<<linename.Cell[i]->L<<" "
        <<linename.Cell[i]->GetP("K2L")<<" "<<linename.Cell[i]->GetP("K2SL")<<endl;
    }
    else if(linename.Cell[i]->TYPE==string("OCT")){
      f2<<linename.Cell[i]->NAME<<"   "<<linename.Cell[i]->TYPE<<" "<<scientific<<linename.Cell[i]->L<<" "
        <<linename.Cell[i]->GetP("K3L")<<" "<<linename.Cell[i]->GetP("K3SL")<<endl;
    }
    else if(linename.Cell[i]->TYPE==string("MULT")){
      f2<<linename.Cell[i]->NAME<<"   "<<linename.Cell[i]->TYPE<<" "<<scientific<<linename.Cell[i]->L<<" ,  "<<endl;
      f2<<"           "
        <<linename.Cell[i]->GetP("K0L")<<" "<<linename.Cell[i]->GetP("K1L")<<" "<<linename.Cell[i]->GetP("K2L")<<" "
	<<linename.Cell[i]->GetP("K3L")<<" "<<linename.Cell[i]->GetP("K4L")<<" "<<linename.Cell[i]->GetP("K5L")<<" "
        <<linename.Cell[i]->GetP("K6L")<<" "<<linename.Cell[i]->GetP("K7L")<<" "<<linename.Cell[i]->GetP("K8L")<<" "
        <<linename.Cell[i]->GetP("K9L")<<" "<<linename.Cell[i]->GetP("K10L")<<" , "<<endl;
      f2<<"           "
        <<linename.Cell[i]->GetP("K0SL")<<" "<<linename.Cell[i]->GetP("K1SL")<<" "<<linename.Cell[i]->GetP("K2SL")<<" "
	<<linename.Cell[i]->GetP("K3SL")<<" "<<linename.Cell[i]->GetP("K4SL")<<" "<<linename.Cell[i]->GetP("K5SL")<<" "
        <<linename.Cell[i]->GetP("K6SL")<<" "<<linename.Cell[i]->GetP("K7SL")<<" "<<linename.Cell[i]->GetP("K8SL")<<" "
        <<linename.Cell[i]->GetP("K9SL")<<" "<<linename.Cell[i]->GetP("K10SL")<<endl;
    }
    else if(linename.Cell[i]->TYPE==string("GMULT")){
      f2<<linename.Cell[i]->NAME<<"   "<<linename.Cell[i]->TYPE<<" "<<scientific<<linename.Cell[i]->L<<" "
        <<linename.Cell[i]->GetP("ANGLE")<<" "<<linename.Cell[i]->GetP("E1")<<" "<<linename.Cell[i]->GetP("E2")<<" , "<<endl;
      f2<<"                        "
        <<linename.Cell[i]->GetP("K0L")<<" "<<linename.Cell[i]->GetP("K1L")<<" "<<linename.Cell[i]->GetP("K2L")<<" "
	<<linename.Cell[i]->GetP("K3L")<<" "<<linename.Cell[i]->GetP("K4L")<<" "<<linename.Cell[i]->GetP("K5L")<<" "
        <<linename.Cell[i]->GetP("K6L")<<" "<<linename.Cell[i]->GetP("K7L")<<" "<<linename.Cell[i]->GetP("K8L")<<" "
        <<linename.Cell[i]->GetP("K9L")<<" "<<linename.Cell[i]->GetP("K10L")<<" , "<<endl;
      f2<<"                        "
        <<linename.Cell[i]->GetP("K0SL")<<" "<<linename.Cell[i]->GetP("K1SL")<<" "<<linename.Cell[i]->GetP("K2SL")<<" "
	<<linename.Cell[i]->GetP("K3SL")<<" "<<linename.Cell[i]->GetP("K4SL")<<" "<<linename.Cell[i]->GetP("K5SL")<<" "
        <<linename.Cell[i]->GetP("K6SL")<<" "<<linename.Cell[i]->GetP("K7SL")<<" "<<linename.Cell[i]->GetP("K8SL")<<" "
        <<linename.Cell[i]->GetP("K9SL")<<" "<<linename.Cell[i]->GetP("K10SL")<<endl;
    }
    else if(linename.Cell[i]->TYPE==string("SMULT")){
      f2<<linename.Cell[i]->NAME<<"   "<<linename.Cell[i]->TYPE<<" "<<scientific<<linename.Cell[i]->L<<" "
        <<linename.Cell[i]->GetP("ANGLE")<<" "<<linename.Cell[i]->GetP("E1")<<" "<<linename.Cell[i]->GetP("E2")<<" , "<<endl;
      f2<<"                        "
        <<linename.Cell[i]->GetP("K0L")<<" "<<linename.Cell[i]->GetP("K1L")<<" "<<linename.Cell[i]->GetP("K2L")<<" "
	<<linename.Cell[i]->GetP("K3L")<<" "<<linename.Cell[i]->GetP("K4L")<<" "<<linename.Cell[i]->GetP("K5L")<<" "
        <<linename.Cell[i]->GetP("K6L")<<" "<<linename.Cell[i]->GetP("K7L")<<" "<<linename.Cell[i]->GetP("K8L")<<" "
        <<linename.Cell[i]->GetP("K9L")<<" "<<linename.Cell[i]->GetP("K10L")<<" , "<<endl;
      f2<<"                        "
        <<linename.Cell[i]->GetP("K0SL")<<" "<<linename.Cell[i]->GetP("K1SL")<<" "<<linename.Cell[i]->GetP("K2SL")<<" "
	<<linename.Cell[i]->GetP("K3SL")<<" "<<linename.Cell[i]->GetP("K4SL")<<" "<<linename.Cell[i]->GetP("K5SL")<<" "
        <<linename.Cell[i]->GetP("K6SL")<<" "<<linename.Cell[i]->GetP("K7SL")<<" "<<linename.Cell[i]->GetP("K8SL")<<" "
        <<linename.Cell[i]->GetP("K9SL")<<" "<<linename.Cell[i]->GetP("K10SL")<<endl;
    }    
    else if(linename.Cell[i]->TYPE==string("SOLEN")){
      f2<<linename.Cell[i]->NAME<<"   "<<linename.Cell[i]->TYPE<<" "<<scientific<<linename.Cell[i]->L<<" "
        <<linename.Cell[i]->GetP("KS")<<endl;
    }
    else if(linename.Cell[i]->TYPE==string("MATRIX")){
      f2<<linename.Cell[i]->NAME<<"   "<<linename.Cell[i]->TYPE<<" "<<scientific<<linename.Cell[i]->L<<" , "<<endl;
      f2<<"           "
        <<linename.Cell[i]->GetP("XCO_IN_X")<<" "<<linename.Cell[i]->GetP("XCO_IN_PX")<<" "
        <<linename.Cell[i]->GetP("XCO_IN_Y")<<" "<<linename.Cell[i]->GetP("XCO_IN_PY")<<" "
        <<linename.Cell[i]->GetP("XCO_IN_Z")<<" "<<linename.Cell[i]->GetP("XCO_IN_DELTA")<<" , "<<endl;
      f2<<"           "
        <<linename.Cell[i]->GetP("XCO_OUT_X")<<" "<<linename.Cell[i]->GetP("XCO_OUT_PX")<<" "
        <<linename.Cell[i]->GetP("XCO_OUT_Y")<<" "<<linename.Cell[i]->GetP("XCO_OUT_PY")<<" "
        <<linename.Cell[i]->GetP("XCO_OUT_Z")<<" "<<linename.Cell[i]->GetP("XCO_OUT_DELTA")<<" , "<<endl;
     f2<<"            "
        <<linename.Cell[i]->GetP("M11")<<" "<<linename.Cell[i]->GetP("M12")<<" "
	<<linename.Cell[i]->GetP("M13")<<" "<<linename.Cell[i]->GetP("M14")<<" "
       <<linename.Cell[i]->GetP("M15")<<" "<<linename.Cell[i]->GetP("M16")<<" ,  "<<endl;
     f2<<"           "
        <<linename.Cell[i]->GetP("M21")<<" "<<linename.Cell[i]->GetP("M22")<<" "
	<<linename.Cell[i]->GetP("M23")<<" "<<linename.Cell[i]->GetP("M24")<<" "
       <<linename.Cell[i]->GetP("M25")<<" "<<linename.Cell[i]->GetP("M26")<<" , "<<endl;
     f2<<"           "
        <<linename.Cell[i]->GetP("M31")<<" "<<linename.Cell[i]->GetP("M32")<<" "
	<<linename.Cell[i]->GetP("M33")<<" "<<linename.Cell[i]->GetP("M34")<<" "
       <<linename.Cell[i]->GetP("M35")<<" "<<linename.Cell[i]->GetP("M36")<<" , "<<endl;
     f2<<"           "
        <<linename.Cell[i]->GetP("M41")<<" "<<linename.Cell[i]->GetP("M42")<<" "
	<<linename.Cell[i]->GetP("M43")<<" "<<linename.Cell[i]->GetP("M44")<<" "
       <<linename.Cell[i]->GetP("M45")<<" "<<linename.Cell[i]->GetP("M46")<<" , "<<endl;
     f2<<"           "
        <<linename.Cell[i]->GetP("M51")<<" "<<linename.Cell[i]->GetP("M52")<<" "
	<<linename.Cell[i]->GetP("M53")<<" "<<linename.Cell[i]->GetP("M54")<<" "
       <<linename.Cell[i]->GetP("M55")<<" "<<linename.Cell[i]->GetP("M56")<<" , "<<endl;
     f2<<"           "
        <<linename.Cell[i]->GetP("M61")<<" "<<linename.Cell[i]->GetP("M62")<<" "
	<<linename.Cell[i]->GetP("M63")<<" "<<linename.Cell[i]->GetP("M64")<<" "
	<<linename.Cell[i]->GetP("M65")<<" "<<linename.Cell[i]->GetP("M66")<<endl;
    }
    else if(linename.Cell[i]->TYPE==string("KICKER")){
      f2<<linename.Cell[i]->NAME<<"   "<<linename.Cell[i]->TYPE<<" "<<scientific<<linename.Cell[i]->L<<" "
        <<linename.Cell[i]->GetP("HKICK")<<" "<<linename.Cell[i]->GetP("VKICK")<<endl;
    }
    else if(linename.Cell[i]->TYPE==string("HKICKER")){
      f2<<linename.Cell[i]->NAME<<"   "<<linename.Cell[i]->TYPE<<" "<<scientific<<linename.Cell[i]->L<<" "
        <<linename.Cell[i]->GetP("HKICK")<<endl;
    }
    else if(linename.Cell[i]->TYPE==string("VKICKER")){
      f2<<linename.Cell[i]->NAME<<"   "<<linename.Cell[i]->TYPE<<" "<<scientific<<linename.Cell[i]->L<<" "
        <<linename.Cell[i]->GetP("VKICK")<<endl;
    }
    else if(linename.Cell[i]->TYPE==string("HACMULT")){
      f2<<linename.Cell[i]->NAME<<"   "<<linename.Cell[i]->TYPE<<" "<<scientific<<linename.Cell[i]->L<<" "<<linename.Cell[i]->GetP("Norder")<<" "
        <<linename.Cell[i]->GetP("KLMAX")<<" "<<linename.Cell[i]->GetP("TTURNS")<<" " <<linename.Cell[i]->GetP("PHI0")<<endl;
    }
    else if(linename.Cell[i]->TYPE==string("VACMULT")){
      f2<<linename.Cell[i]->NAME<<"   "<<linename.Cell[i]->TYPE<<" "<<scientific<<linename.Cell[i]->L<<" "<<linename.Cell[i]->GetP("Norder")<<" "
        <<linename.Cell[i]->GetP("KSLMAX")<<" "<<linename.Cell[i]->GetP("TTURNS")<<" "<<linename.Cell[i]->GetP("PHI0")<<endl;
    }
    else if(linename.Cell[i]->TYPE==string("HACDIP")){
      f2<<linename.Cell[i]->NAME<<"   "<<linename.Cell[i]->TYPE<<" "<<scientific<<linename.Cell[i]->L<<" "
        <<linename.Cell[i]->GetP("HKICKMAX")<<" "<<linename.Cell[i]->GetP("NUD")<<" "
        <<linename.Cell[i]->GetP("PHID")<<" "<<linename.Cell[i]->GetP("TURNS")<<" "
        <<linename.Cell[i]->GetP("TURNE")<<endl;
    }
    else if(linename.Cell[i]->TYPE==string("VACDIP")){
      f2<<linename.Cell[i]->NAME<<"   "<<linename.Cell[i]->TYPE<<" "<<scientific<<linename.Cell[i]->L<<" "
        <<linename.Cell[i]->GetP("VKICKMAX")<<" "<<linename.Cell[i]->GetP("NUD")<<" "
        <<linename.Cell[i]->GetP("PHID")<<" "<<linename.Cell[i]->GetP("TURNS")<<" "
        <<linename.Cell[i]->GetP("TURNE")<<endl;
    }
    else if(linename.Cell[i]->TYPE==string("BPM")){
      f2<<linename.Cell[i]->NAME<<"   "<<linename.Cell[i]->TYPE<<" "<<scientific<<linename.Cell[i]->L<<endl;
    }
    else if(linename.Cell[i]->TYPE==string("HBPM")){
      f2<<linename.Cell[i]->NAME<<"   "<<linename.Cell[i]->TYPE<<" "<<scientific<<linename.Cell[i]->L<<endl;
    }
    else if(linename.Cell[i]->TYPE==string("VBPM")){
      f2<<linename.Cell[i]->NAME<<"   "<<linename.Cell[i]->TYPE<<" "<<scientific<<linename.Cell[i]->L<<endl;
    }
    else if(linename.Cell[i]->TYPE==string("MARKER")){
      f2<<linename.Cell[i]->NAME<<"   "<<linename.Cell[i]->TYPE<<" "<<scientific<<linename.Cell[i]->L<<endl;
    }
    else if(linename.Cell[i]->TYPE==string("RFCAV")){
      f2<<linename.Cell[i]->NAME<<"   "<<linename.Cell[i]->TYPE<<" "<<scientific<<linename.Cell[i]->L<<" "
        <<linename.Cell[i]->GetP("VRF")<<" "<<linename.Cell[i]->GetP("FRF")<<" "
        <<linename.Cell[i]->GetP("PHASE0")<<endl;
    }
    else if(linename.Cell[i]->TYPE==string("BEAMBEAM")){
      f2<<linename.Cell[i]->NAME<<"   "<<linename.Cell[i]->TYPE<<" "
        <<linename.Cell[i]->GetP("TREATMENT")<<" "<<linename.Cell[i]->GetP("NP")<<" "<<linename.Cell[i]->GetP("BBSCALE")<<" "
        <<linename.Cell[i]->GetP("SIGMAL")<<" "<<int(linename.Cell[i]->GetP("NSLICE"))<<" , "<<endl;
      f2<<"           "
        <<linename.Cell[i]->GetP("EMITX")<<" "<<linename.Cell[i]->GetP("EMITY")<<" "
        <<linename.Cell[i]->GetP("BETAX")<<" "<<linename.Cell[i]->GetP("ALFAX")<<" "
        <<linename.Cell[i]->GetP("BETAY")<<" "<<linename.Cell[i]->GetP("ALFAY")<<endl;
    }
    else if(linename.Cell[i]->TYPE==string("LRBB")){
      f2<<linename.Cell[i]->NAME<<"   "<<linename.Cell[i]->TYPE<<"  "
        <<linename.Cell[i]->GetP("TREATMENT")<<" "<<linename.Cell[i]->GetP("NP")<<" "<<linename.Cell[i]->GetP("BBSCALE")<<" "
        <<linename.Cell[i]->GetP("SIGMAL")<<" "<<int(linename.Cell[i]->GetP("NSLICE"))<<" ,  "<<endl;
      f2<<"           "
        <<linename.Cell[i]->GetP("SEPX")<<" "<<linename.Cell[i]->GetP("SEPY")<<" "
        <<linename.Cell[i]->GetP("EMITX")<<" "<<linename.Cell[i]->GetP("EMITY")<<" "
        <<linename.Cell[i]->GetP("BETAX")<<" "<<linename.Cell[i]->GetP("ALFAX")<<" "
        <<linename.Cell[i]->GetP("BETAY")<<" "<<linename.Cell[i]->GetP("ALFAY")<<endl;
    }
    else if(linename.Cell[i]->TYPE==string("ELENS")){
      f2<<linename.Cell[i]->NAME<<"   "<<linename.Cell[i]->TYPE<<" "<<scientific<<linename.Cell[i]->L<<" "
        <<linename.Cell[i]->GetP("NE")<<" "<<linename.Cell[i]->GetP("BBSCALE")<<" "<<int(linename.Cell[i]->GetP("NSLICE"))<<" "
        <<linename.Cell[i]->GetP("BETAE")<<" "<<linename.Cell[i]->GetP("SIGMAX")<<" "
        <<linename.Cell[i]->GetP("SIGMAY")<<endl;
    }
    else if(linename.Cell[i]->TYPE==string("ERHICBB")){
      f2<<linename.Cell[i]->NAME<<"   "<<linename.Cell[i]->TYPE<<" "
        <<linename.Cell[i]->GetP("NE")<<" "<<linename.Cell[i]->GetP("BBSCALE")<<" "
        <<linename.Cell[i]->GetP("SIGMAX")<<" "<<linename.Cell[i]->GetP("SIGMAY")<<endl;
    }
    else if(linename.Cell[i]->TYPE==string("ROTAT")){
      f2<<linename.Cell[i]->NAME<<"   "<<linename.Cell[i]->TYPE<<" "<<scientific<<linename.Cell[i]->L<<" "
        <<linename.Cell[i]->GetP("NX")<<" "<<linename.Cell[i]->GetP("NY")<<" "
        <<linename.Cell[i]->GetP("NS")<<" "<<linename.Cell[i]->GetP("ANGLE")<<endl;
    }
    else if(linename.Cell[i]->TYPE==string("SNAKE")){
      f2<<linename.Cell[i]->NAME<<"   "<<linename.Cell[i]->TYPE<<" "<<scientific<<linename.Cell[i]->L<<" "
        <<linename.Cell[i]->GetP("NX")<<" "<<linename.Cell[i]->GetP("NY")<<" "
        <<linename.Cell[i]->GetP("NS")<<" "<<linename.Cell[i]->GetP("ANGLE")<<endl;
    }
    else if(linename.Cell[i]->TYPE==string("DIFFUSE")){
      f2<<linename.Cell[i]->NAME<<"   "<<linename.Cell[i]->TYPE<<"  "
        <<linename.Cell[i]->GetP("DIFF_X")<<" "<<linename.Cell[i]->GetP("DIFF_Y")<<" "
        <<linename.Cell[i]->GetP("DIFF_DELTA")<<endl;
    }
    else if(linename.Cell[i]->TYPE==string("COOLING")){
      f2<<linename.Cell[i]->NAME<<"   "<<linename.Cell[i]->TYPE<<" "
        <<linename.Cell[i]->GetP("ALPHA")<<endl;
    }
    else{
      cout<<"Can't save this element: "<<linename.Cell[i]->NAME<<"   "<<linename.Cell[i]->TYPE<<endl;
      exit(0);
    }
  }
  f2.close();
}

//================================
//
//     Beam Dynmaics 
//   
//================================

//---Survey only in horizontal plane, not able to do vertical survey.
void Survey(Line & linename, double X0, double Y0, double Theta0, int dir, const char* filename)
{
  int i;
  double x0, y0, x1,y1,Ox, Oy;
  double dtheta,L,B,R;
  double X1, Y1, Theta1;
  fstream fout;
  
  if(dir !=1 and dir != -1) {
    cout<<"Error: dir in Survey() should be +1 or -1 !"<<endl;
    exit(0);
  }

  fout.open(filename, ios::out);
  fout<<setw(15) <<"NAME"<<setw(15) <<"TYPE"<<setw(15) <<"S"
      <<setw(15) <<"L   "<<setw(15) <<"X    "<<setw(15)<<"Y    "<<setw(15)<<"THETA"<<endl;
  fout<<setw(15) <<"START$POINT"<<setw(15) <<"MARKER"<<setw(15) <<0
	 <<setw(15) <<0<<setw(15) <<X0<<setw(15) <<Y0<<setw(15) <<Theta0<<endl;

  for(i=0;i<linename.Ncell;i++){
    if( (linename.Cell[i]->TYPE==string("SBEND") and linename.Cell[i]->GetP("ANGLE")!= 0. ) or
        (linename.Cell[i]->TYPE==string("GMULT") and linename.Cell[i]->GetP("ANGLE")!= 0. ) or
	(linename.Cell[i]->TYPE==string("SMULT") and linename.Cell[i]->GetP("ANGLE")!= 0. ) ){
      L=linename.Cell[i]->L;
      dtheta=linename.Cell[i]->GetP("ANGLE");
      B=dtheta * GP.brho /L;
      R=abs(GP.brho/B);
      if( dir == 1) {
	if(dtheta < 0) {
	  x0=cos(Theta0); y0=sin(Theta0);
	  rotation(x1, y1, x0, y0, -PI/2);
	  Ox=  X0 - R * x1 ;
	  Oy=  Y0 - R * y1 ;
	  x0= R * x1;
	  y0= R * y1 ;
	  rotation(x1, y1, x0, y0, -dtheta); 
	  X1= Ox + x1;  
	  Y1= Oy + y1;   
	  Theta1= Theta0 - dtheta;}
	else{
	  x0=cos(Theta0); y0=sin(Theta0);
	  rotation(x1, y1, x0, y0, PI/2);
	  Ox=  X0 - R * x1 ;
	  Oy=  Y0 - R * y1 ;
	  x0= R * x1;
	  y0= R * y1 ;
	  rotation(x1, y1, x0, y0, -dtheta); 
	  X1= Ox + x1;  
	  Y1= Oy + y1;
	  Theta1= Theta0 - dtheta;	
	}
      }
      else{
	if(dtheta < 0) {
	  x0=cos(Theta0); y0=sin(Theta0);
	  rotation(x1, y1, x0, y0, PI/2);
	  Ox=  X0 - R * x1 ;
	  Oy=  Y0 - R * y1 ;
	  x0= R * x1;
	  y0= R * y1 ;
	  rotation(x1, y1, x0, y0, dtheta); 
	  X1= Ox + x1;  
	  Y1= Oy + y1;   
	  Theta1= Theta0 + dtheta;}
	else{
	  x0=cos(Theta0); y0=sin(Theta0);
	  rotation(x1, y1, x0, y0, -PI/2);
	  Ox=  X0 - R * x1 ;
	  Oy=  Y0 - R * y1 ;
	  x0= R * x1;
	  y0= R * y1 ;
	  rotation(x1, y1, x0, y0, dtheta); 
	  X1= Ox + x1;  
	  Y1= Oy + y1;
	  Theta1= Theta0 + dtheta;	
	}
      }
    }
    else{
      L=linename.Cell[i]->L;
      X1=  X0 + L * cos(Theta0);
      Y1=  Y0 + L * sin(Theta0);
      Theta1= Theta0;
    }
    fout <<setw(15) <<linename.Cell[i]->NAME<<setw(15) <<linename.Cell[i]->TYPE<<setw(15) <<linename.Cell[i]->S
	 <<setw(15) <<linename.Cell[i]->L<<setw(15) <<X1<<setw(15) <<Y1<<setw(15) <<Theta1<<endl;
    X0=X1; Y0=Y1; Theta0=Theta1;
  }
  fout.close();
}

void Cal_Orbit_Num(Line & linename, double deltap)
{
  int i,j,k,iter=0;
  double x0[6],x01[6], x1[6], dx[4];
  double d=1.0e-09;
  double mat[4][4];
  double chi;
  double codeps=1e-12;
  int Max_iter=20;
  int flag;

  vector<double> voltage;
  for(i=0;i<linename.Ncell;i++){
    if(linename.Cell[i]->TYPE==string("RFCAV") || linename.Cell[i]->TYPE==string("CRABRF")
                                               || linename.Cell[i]->TYPE==string("CCMULT") ){
	voltage.push_back( linename.Cell[i]->GetP("VRF") ) ;
        linename.Cell[i]->SetP("VRF", 0.0);
      }
  }

  for (i=0;i<6;i++) x0[i]=0.;
  do{
    
    iter++;
    x0[4]=0;
    x0[5]=deltap;
    
    for(j=0;j<6;j++) x01[j]=x0[j];
    for(j=0;j<linename.Ncell;j++)  linename.Cell[j]->Pass(x01);
    
    
    for(k=0;k<4;k++){
      for(j=0;j<6;j++) x1[j]=x0[j];
      x1[k]=x1[k]+d;
      for(j=0;j<linename.Ncell;j++) linename.Cell[j]->Pass(x1);
      for(j=0;j<4;j++) mat[j][k]=(x1[j]-x01[j])/d;
    }
    
    for(j=0;j<6;j++) x1[j]=x0[j];
    for(j=0;j<linename.Ncell;j++)  linename.Cell[j]->Pass(x1);
    for(i=0;i<4;i++) {
      dx[i]=x0[i]-x1[i];
    }
    chi=0;
    for(i=0;i<4;i++) chi+=dx[i]*dx[i];
    chi=sqrt(chi/4.);
    
    if( chi > codeps ){
      for(i=0;i<4;i++) mat[i][i]=mat[i][i]-1.0000001;
      flag=mat_inv(&mat[0][0],4);
      if( flag == 0  ){
	cout<<" Failed in closed orbit searching. Exit. "<<endl;
	exit(1);
      }
      for(i=0;i<4;i++){
	for(j=0;j<4;j++) x0[i]+=mat[i][j]*dx[j];  
      } 
    } 
    
  }while (chi>codeps && iter <Max_iter );
  
  if(iter == Max_iter-1 ) 
    {
      cout<<"Failed to find COD."<<endl;
      exit(1);
    }
  else
    {
      x0[4]=0.000;
      x0[5]=deltap;
      for(j=0;j<linename.Ncell; j++) {
	linename.Cell[j]->Pass(x0);
	for(i=0;i<6;i++) linename.Cell[j]->X[i]=x0[i];
      } 
    }

  j=0;
  for(i=0;i<linename.Ncell;i++){
    if(linename.Cell[i]->TYPE==string("RFCAV") || linename.Cell[i]->TYPE==string("CRABRF")
                                               || linename.Cell[i]->TYPE==string("CCMULT")  ) {
      linename.Cell[i]->SetP("VRF", voltage[j]);
      j++;
    }
  }

}

//----4-d  optics calculation, with fixed deltap

void Cal_Orbit_Num(Line & linename, double z, double deltap)   // here  z  is   w.r.t. on-momentum particle
{
  int i,j,k,iter=0;
  double x0[6],x01[6], x1[6], dx[4];
  double d=1.0e-09;
  double mat[4][4];
  double chi;
  double codeps=1e-12;
  int Max_iter=20;
  int flag;
  vector<double> voltage;

  for(i=0;i<linename.Ncell;i++){   //  here we  only block  RFCAV  to keep dp/p0 constant
    if(linename.Cell[i]->TYPE==string("RFCAV") ){
      voltage.push_back( linename.Cell[i]->GetP("VRF") ) ;
      linename.Cell[i]->SetP("VRF", 0.0);
    }
  }

  for (i=0;i<6;i++) x0[i]=0.;
  do{
    
    iter++;
    x0[4]=z;
    x0[5]=deltap;
    
    for(j=0;j<6;j++) x01[j]=x0[j];
    for(j=0;j<linename.Ncell;j++)  linename.Cell[j]->Pass(x01);
    
    
    for(k=0;k<4;k++){
      for(j=0;j<6;j++) x1[j]=x0[j];
      x1[k]=x1[k]+d;
      for(j=0;j<linename.Ncell;j++) linename.Cell[j]->Pass(x1);
      for(j=0;j<4;j++) mat[j][k]=(x1[j]-x01[j])/d;
    }
    
    for(j=0;j<6;j++) x1[j]=x0[j];
    for(j=0;j<linename.Ncell;j++)  linename.Cell[j]->Pass(x1);
    for(i=0;i<4;i++) {
      dx[i]=x0[i]-x1[i];
    }
    chi=0;
    for(i=0;i<4;i++) chi+=dx[i]*dx[i];
    chi=sqrt(chi/4.);
    
    if( chi > codeps ){
      for(i=0;i<4;i++) mat[i][i]=mat[i][i]-1.0000001;
      flag=mat_inv(&mat[0][0],4); 
      if( flag == 0  ){
	cout<<" Failed in closed orbit searching. Exit. "<<endl;
	exit(1);
      }

      for(i=0;i<4;i++){
	for(j=0;j<4;j++) x0[i]+=mat[i][j]*dx[j];  
      } 
    } 
    
  }while (chi>codeps && iter <Max_iter );
  
  if(iter == Max_iter-1 ) 
    {
      cout<<"Failed to find COD."<<endl;
      exit(1);
    }
  else
    {
      x0[4]=z;
      x0[5]=deltap;
      for(j=0;j<linename.Ncell; j++) {
	linename.Cell[j]->Pass(x0);
	for(i=0;i<6;i++) linename.Cell[j]->X[i]=x0[i];
      } 
    }

  j=0;
  for(i=0;i<linename.Ncell;i++){
    if(linename.Cell[i]->TYPE==string("RFCAV") ){
      linename.Cell[i]->SetP("VRF", voltage[j]);
      j++;
    }
  }

}

void Cal_OneTurnMap(Line & linename, double deltap)
{
  int i,j;
  double x[6];
  tps tps1[6];
  linmap m1;
  double m66[36], b66[36];
  int flag;
  double u[6], v[6];

  vector<double> voltage;
  for(i=0;i<linename.Ncell;i++){
    if(linename.Cell[i]->TYPE==string("RFCAV") ){
	voltage.push_back( linename.Cell[i]->GetP("VRF") ) ;
        linename.Cell[i]->SetP("VRF", 0.0);
      }
  }

  for(i=0;i<6;i++) x[i]=linename.Cell[linename.Ncell-1]->X[i];
  x[4]=0.0;   
  x[5]=deltap; 
  m1.identity();
  m1=m1+x;

  for(i=0;i<6;i++) tps1[i]=m1[i];
  for(j=0;j<linename.Ncell; j++) linename.Cell[j]->DAPass(tps1);
  for(i=0;i<6;i++) m1[i]=tps1[i];
  Getmat(m1, m66); 
  for(i=0;i<36;i++) linename.Cell[linename.Ncell-1]->M[i]=m66[i];
  
  if( false ) {
    cout<<"One turn map M:"<<endl;
    for (i=0;i<6;i++) {
      for(j=0;j<6;j++) cout<<setw(12)<<m66[i*6+j]<<"  ";
      cout<<endl;
    }
    for(i=0;i<36;i++) b66[i]=m66[i];
    cout<<"Det of M = "<<mat_det(b66,6)<<endl;
  }

  mat_change_hessenberg(m66, 6);
  flag=mat_root_hessenberg(m66,6,u,v,1.0e-10,60);
  if(flag < 0){
    cout<<"One Turn Map: unstable.  exit. \n"<<endl;
    exit(1);
  }


  j=0;
  for(i=0;i<linename.Ncell;i++){
    if(linename.Cell[i]->TYPE==string("RFCAV")  ) {
      linename.Cell[i]->SetP("VRF", voltage[j]);
      j++;
    }
  }

    
}

void Cal_OneElementMap(Element * temp_element, double x[], double t66[])
{
  int i,j,k;
  tps   tps1[6];
  linmap m1;
  
  m1.identity();
  m1=m1+x;
  
  for(i=0;i<6;i++) tps1[i]=m1[i];
  temp_element->DAPass(tps1);
  for(i=0;i<6;i++) m1[i]=tps1[i];
  Getmat(m1, t66);

  if(false){
    cout<<"Transfer map:"<<endl;
    for (i=0;i<6;i++) {
      for(j=0;j<6;j++) cout<<setw(12)<<t66[i*6+j]<<"  ";
      cout<<endl;
    }
  }
  
}

void Cal_ElementMap(Line & linename, double deltap)
{
  int i,j,k, istart;
  double x[6];
  tps tps1[6];
  linmap m1;
  double t66[36];

  vector<double> voltage;
  for(i=0;i<linename.Ncell;i++){
    if(linename.Cell[i]->TYPE==string("RFCAV") ){
	voltage.push_back( linename.Cell[i]->GetP("VRF") ) ;
        linename.Cell[i]->SetP("VRF", 0.0);
      }
  }

  for(k=0; k<linename.Ncell; k++) 
    {
      istart=k-1;
      if (istart <0) istart=linename.Ncell-1;

      for(i=0;i<6;i++) x[i]=linename.Cell[istart]->X[i];
      x[4]=0.000;
      x[5]=deltap; 
      m1.identity();
      m1=m1+x;
      
      for(i=0;i<6;i++) tps1[i]=m1[i];
      linename.Cell[k]->DAPass(tps1);
      for(i=0;i<6;i++) m1[i]=tps1[i];
      Getmat(m1, t66);
      
      for(j=0;j<36;j++) linename.Cell[k]->T[j]=t66[j]; 
    }

  j=0;
  for(i=0;i<linename.Ncell;i++){
    if(linename.Cell[i]->TYPE==string("RFCAV")  ) {
      linename.Cell[i]->SetP("VRF", voltage[j]);
      j++;
    }
  }
  
}

void  Get_ElementMap(Line & linename, int i1, double m66[])
{
  int  i;

  for(i=0;i<36;i++) m66[i]=linename.Cell[i1]->T[i];

}

void Cal_SectionMap(Line & linename, int i1, int i2, double deltap, double t66[36] )
// not include  the first element i1, but include last element i2
{
  int i;
  double x[6];
  tps tps1[6];
  linmap m1;
  
  if (i1== 0 ) { 
    for(i=0;i<6;i++) x[i] = linename.Cell[linename.Ncell-1]->X[i]; }
  else {
    for(i=0;i<6;i++) x[i] = linename.Cell[i1-1]->X[i]; }
  x[4]=0.000;
  x[5]=deltap; 
  m1.identity();
  m1=m1+x;
  
  for(i=0;i<6;i++) tps1[i]=m1[i];
  for(i=i1+1;i<=i2;i++)  linename.Cell[i]->DAPass(tps1);
  for(i=0;i<6;i++) m1[i]=tps1[i];
  Getmat(m1, t66);
}

void Cal_SectionMap(Line & linename, int i1, int i2, double x0[], double t66[36] )
// not include  the first element i1, but include last element i2
{
  int i;
  double x[6];
  tps tps1[6];
  linmap m1;
  
  for(i=0;i<6;i++) x[i] = x0[i];
  m1.identity();
  m1=m1+x;
  
  for(i=0;i<6;i++) tps1[i]=m1[i];
  for(i=i1+1;i<=i2;i++)  linename.Cell[i]->DAPass(tps1);
  for(i=0;i<6;i++) m1[i]=tps1[i];
  Getmat(m1, t66);
}

void Cal_A(Line & linename, double deltap)
{
  
  int i,j;
  double m44[4][4];

  double wr[4], wi[4], vr[4][4], vi[4][4], temp_wr, temp_wi, temp_vr[4], temp_vi[4];  
  double temp1, temp2;
  double theta1, theta2, tempx, tempy, x1,y1;
  double a44[4][4], b44[4][4];
  double det, scale;
  
  //---solve eigen problem to construct A[4*4]
  for(i=0;i<4;i++)
   for(j=0;j<4;j++) m44[i][j]=linename.Cell[linename.Ncell-1]->M[i*6+j];
  EigenSolver(m44, wr, wi, vr, vi);

  //---link eigen vetcors to  H/V planes 
  temp1=abs(vr[0][0] ) + abs( vi[0][0] ) + abs( vr[0][1] ) + abs( vi[0][1] );
  temp2=abs(vr[2][0] ) + abs( vi[2][0] ) + abs( vr[2][1] ) + abs( vi[2][1] );
  if (temp2 > temp1 ) {
      temp_wr    =  wr[0] ;
      temp_wi    =  wi[0] ;
      wr[0]      =  wr[2] ;
      wi[0]      =  wi[2] ;
      wr[2]      =  temp_wr;
      wi[2]      =  temp_wi;
      wr[1]      =  wr[0];
      wi[1]      = -wi[0]; 
      wr[3]      =  wr[2];
      wi[3]      = -wi[2]; 
    for(i=0;i<4;i++){
      temp_vr[i]    = vr[0][i];
      temp_vi[i]    = vi[0][i];
      vr[0][i]      = vr[2][i];
      vi[0][i]      = vi[2][i];
      vr[2][i]      = temp_vr[i];
      vi[2][i]      = temp_vi[i];
      vr[1][i]      = vr[0][i];
      vi[1][i]      =-vi[0][i]; 
      vr[3][i]      = vr[2][i];
      vi[3][i]      =-vi[2][i]; 
    }
  }

  //----rotating to make A12=0, A34=0
  theta1=atan2(vi[0][0], vr[0][0]);
  x1  =cos(-PI/2-theta1);
  y1  =sin(-PI/2-theta1);
  for(i=0;i<4;i++){
    tempx=vr[0][i]*x1-vi[0][i]*y1;
    tempy=vi[0][i]*x1+vr[0][i]*y1;
    vr[0][i]=tempx;
    vi[0][i]=tempy; }
  x1 =cos(PI/2+theta1);
  y1 =sin(PI/2+theta1);
  for(i=0;i<4;i++){
    tempx=vr[1][i]*x1-vi[1][i]*y1;
    tempy=vi[1][i]*x1+vr[1][i]*y1;
    vr[1][i]=tempx;
    vi[1][i]=tempy;}

  theta2=atan2(vi[2][2], vr[2][2]);
  x1 =cos(-PI/2-theta2);
  y1 =sin(-PI/2-theta2);
  for(i=0;i<4;i++){
    tempx=vr[2][i]*x1-vi[2][i]*y1;
    tempy=vi[2][i]*x1+vr[2][i]*y1;
    vr[2][i]=tempx;
    vi[2][i]=tempy; }
  x1  = cos(PI/2+theta2);
  y1=sin(PI/2+theta2);
  for(i=0;i<4;i++){
    tempx=vr[3][i]*x1-vi[3][i]*y1;
    tempy=vi[3][i]*x1+vr[3][i]*y1;
    vr[3][i]=tempx;
    vi[3][i]=tempy;}

  //----produce A
  for(i=0;i<4;i++) a44[0][i]=-(vi[0][i]-vi[1][i]);
  for(i=0;i<4;i++) a44[1][i]=  vr[0][i]+vr[1][i];
  for(i=0;i<4;i++) a44[2][i]=-(vi[2][i]-vi[3][i]);
  for(i=0;i<4;i++) a44[3][i]=  vr[2][i]+vr[3][i];

  /*----
  if(true){
    if(a44[1][1] < 0. ) {
      for(i=0;i<4;i++) a44[i][1]=- a44[i][1];
      //cout<<"Mode I : negative beta reverted."<<endl; 
   }
    if(a44[3][3] < 0. ) {
      for(i=0;i<4;i++) a44[i][3]=- a44[i][3];
      //cout<<"Mode II : negative beta reverted."<<endl; 
    }
  }
  ------*/
  
  if(true){
    if(a44[1][1] < 0. ) {
      for(i=0;i<4;i++) a44[1][i]=- a44[1][i];
      //cout<<"Mode I : negative beta reverted."<<endl; 
   }
    if(a44[3][3] < 0. ) {
      for(i=0;i<4;i++) a44[3][i]=- a44[3][i];
      //cout<<"Mode II : negative beta reverted."<<endl; 
    }
  }
  
  
  if(false) {
    cout<<"A matrix (before normalizing):"<<endl;
    for (i=0;i<4;i++) {
      for(j=0;j<4;j++) cout<<setw(12)<<a44[i][j]<<"  ";
      cout<<endl;
    }

    for(i=0;i<4;i++) {
      for(j=0;j<4;j++) b44[i][j]=a44[i][j];
    }
    det= mat_det( &b44[0][0], 4);
    cout<<" Det of A44  =  "<< det <<endl;
  }
  
  //---normalizing A
   
  double  s44[4][4];
  for(i=0;i<4;i++)
    for(j=0;j<4;j++) s44[i][j]=0.;
  s44[0][1]= 1.0;
  s44[1][0]=-1.0;
  s44[2][3]= 1.0;
  s44[3][2]=-1.0;  

  if(false){
    cout<<"S matrix :"<<endl;
    for (i=0;i<4;i++) {
      for(j=0;j<4;j++) cout<<setw(12)<<s44[i][j]<<"  ";
      cout<<endl;
    }
  }

  temp1=(-a44[0][1]*a44[1][0] + a44[0][0]*a44[1][1] -a44[0][3]*a44[1][2] +a44[0][2]*a44[1][3] );
  temp2=(-a44[2][1]*a44[3][0] + a44[2][0]*a44[3][1] -a44[2][3]*a44[3][2] +a44[2][2]*a44[3][3] );  

  //cout<<"...normalize new v1 and v2, v3 adn v4:"<<endl;
  //cout<<temp1<<" "<<temp2<<endl;
  
  for(i=0;i<4;i++){
    a44[0][i]=a44[0][i]/sqrt(temp1);
    a44[1][i]=a44[1][i]/sqrt(temp1);
    a44[2][i]=a44[2][i]/sqrt(temp2);
    a44[3][i]=a44[3][i]/sqrt(temp2);
  }
  
  temp1=(-a44[0][1]*a44[1][0] + a44[0][0]*a44[1][1] -a44[0][3]*a44[1][2] +a44[0][2]*a44[1][3] ); 
  temp2=(-a44[2][1]*a44[3][0] + a44[2][0]*a44[3][1] -a44[2][3]*a44[3][2] +a44[2][2]*a44[3][3] );

  //cout<<"...normalize new v1 and v2, v3 adn v4:"<<endl;
  //cout<<temp1<<" "<<temp2<<endl; 

  if(false) {
    cout<<"A matrix :"<<endl;
    for (i=0;i<4;i++) {
      for(j=0;j<4;j++) cout<<setw(12)<<a44[j][i]<<"  ";
      cout<<endl;
    }
    
    for(i=0;i<4;i++) {
      for(j=0;j<4;j++) b44[i][j]=a44[j][i];
    }
    det= mat_det( &b44[0][0], 4);
    cout<<" Det of A  =  "<< det <<endl;
  }

  //----save
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      linename.Cell[linename.Ncell-1]->A[i*6+j]=a44[j][i];
  
}

void Trace_A(Line & linename, double deltap)
{
  int i,j,k;
  double T[4][4], A[4][4], G[4][4];
  double dphi1, dphi2;
  double scale1, scale2;
  double mu1, mu2;
 
  for(i=0;i<4;i++)
    for(j=0;j<4;j++) A[i][j]= linename.Cell[linename.Ncell-1]->A[i*6+j];

  mu1=0.;
  mu2=0.;
  for(k=0;k<linename.Ncell;k++){
    
    for (i=0;i<4;i++){
      for(j=0;j<4;j++) 	T[i][j]=linename.Cell[k]->T[i*6+j];
    }
    
    mat_mult( &T[0][0], &A[0][0], &G[0][0], 4,4,4 );
    
    dphi1=atan2(G[0][1], G[0][0]);
    dphi2=atan2(G[2][3], G[2][2]);
    scale1=sqrt( G[0][0]*G[0][0] + G[0][1]*G[0][1]);
    scale2=sqrt( G[2][2]*G[2][2] + G[2][3]*G[2][3]);

    A[0][0]= ( G[0][0]*G[0][0]+G[0][1]*G[0][1] ) / scale1;
    A[0][1]= (-G[0][0]*G[0][1]+G[0][1]*G[0][0] ) / scale1;
    A[1][0]= ( G[1][0]*G[0][0]+G[1][1]*G[0][1] ) / scale1;
    A[1][1]= (-G[1][0]*G[0][1]+G[1][1]*G[0][0] ) / scale1;
                                                                          
    A[2][0]= ( G[2][0]*G[0][0]+G[2][1]*G[0][1] ) / scale1;
    A[2][1]= (-G[2][0]*G[0][1]+G[2][1]*G[0][0] ) / scale1;
    A[3][0]= ( G[3][0]*G[0][0]+G[3][1]*G[0][1] ) / scale1;
    A[3][1]= (-G[3][0]*G[0][1]+G[3][1]*G[0][0] ) / scale1; 
                                                                          
    A[0][2]= ( G[0][2]*G[2][2]+G[0][3]*G[2][3] ) / scale2;
    A[0][3]= (-G[0][2]*G[2][3]+G[0][3]*G[2][2] ) / scale2;
    A[1][2]= ( G[1][2]*G[2][2]+G[1][3]*G[2][3] ) / scale2;
    A[1][3]= (-G[1][2]*G[2][3]+G[1][3]*G[2][2] ) / scale2;
                                                                          
    A[2][2]= ( G[2][2]*G[2][2]+G[2][3]*G[2][3] ) / scale2;
    A[2][3]= (-G[2][2]*G[2][3]+G[2][3]*G[2][2] ) / scale2;
    A[3][2]= ( G[3][2]*G[2][2]+G[3][3]*G[2][3] ) / scale2;
    A[3][3]= (-G[3][2]*G[2][3]+G[3][3]*G[2][2] ) / scale2; 

    for(i=0;i<4;i++){
      for(j=0;j<4;j++) linename.Cell[k]->A[i*6+j]=A[i][j];
    }
    
    mu1=mu1+dphi1/ 2./PI;
    mu2=mu2+dphi2/ 2./PI;
    linename.Cell[k]->Mu1= mu1;
    linename.Cell[k]->Mu2= mu2;
  }

  linename.Tune1=  mu1;
  linename.Tune2=  mu2;
}

void Cal_Twiss(Line & linename, double deltap)
//  here  input "deltap" is always dp/p0 no matter GP.H_expand, 4-d  treatment
{
  int i,j,k;
  double pt, temp1, temp2,r;
  double A11[4], A12[4], A22[4], A22I[4], C[4];
  int    flag;

  if(GP.H_expand == true ){
    Cal_Orbit_Num(linename,deltap);
    Cal_OneTurnMap(linename,deltap);
    Cal_A(linename,deltap);
    Cal_ElementMap(linename,deltap);
    Trace_A(linename,deltap);
  }
  else{
    pt=DeltaToPt(deltap);
    Cal_Orbit_Num(linename,pt);
    Cal_OneTurnMap(linename,pt);
    Cal_A(linename,pt);
    Cal_ElementMap(linename,pt);
    Trace_A(linename,pt);
  }

  for(k=0; k<linename.Ncell;k++){
    linename.Cell[k]->Beta1= linename.Cell[k]->A[0*6+0]/linename.Cell[k]->A[1*6+1];
    linename.Cell[k]->Alfa1=-linename.Cell[k]->A[1*6+0]/linename.Cell[k]->A[1*6+1];
    linename.Cell[k]->Beta2= linename.Cell[k]->A[2*6+2]/linename.Cell[k]->A[3*6+3];
    linename.Cell[k]->Alfa2=-linename.Cell[k]->A[3*6+2]/linename.Cell[k]->A[3*6+3];
    
    for(i=0;i<2;i++){  //  updated March 26, 2022
      for(j=0;j<2;j++) {
	A11[i*2+j]= linename.Cell[k]->A[6*i+j];
	A12[i*2+j]= linename.Cell[k]->A[6*i+(j+2)];
	A22[i*2+j]= linename.Cell[k]->A[6*(i+2)+(j+2)];
      }    
    }
    
    temp1=mat_det(A11,2);
    temp2=mat_det(A22,2);
    if( abs(temp1 -temp2) > 1.0e-10) {
      //----switched off
      //cout<<"Warning:  during calculating r Cal_Twiss(). "<<endl;
    }
    r=sqrt( mat_det( A22,2 ) );
    
    A22I[0]  =  A22[3]/( A22[0] * A22[3] );
    A22I[1]  = -A22[1]/( A22[0] * A22[3] );
    A22I[2]  = -A22[2]/( A22[0] * A22[3] );
    A22I[3]  =  A22[0]/( A22[0] * A22[3] );

    mat_mult(A12, A22I, C,2,2,2);
    for(i=0;i<4;i++) C[i]=C[i]*r;
    
    linename.Cell[k]->r=r;
    linename.Cell[k]->c11=C[0];
    linename.Cell[k]->c12=C[1];
    linename.Cell[k]->c21=C[2];
    linename.Cell[k]->c22=C[3];
  }
  
}

//-----for  Twiss  calcualtion  with z  constant

void Cal_OneTurnMap_Z(Line & linename, double z0)
{
  int i,j;
  double x[6];
  tps tps1[6];
  linmap m1;
  double m66[36], b66[36];
  int flag;
  double u[6], v[6];

  vector<double> voltage;
  for(i=0;i<linename.Ncell;i++){
    if(linename.Cell[i]->TYPE==string("RFCAV") ){
	voltage.push_back( linename.Cell[i]->GetP("VRF") ) ;
        linename.Cell[i]->SetP("VRF", 0.0);
      }
  }

  for(i=0;i<6;i++) x[i]=linename.Cell[linename.Ncell-1]->X[i];
  x[4]=z0;   
  x[5]=0.0; 
  m1.identity();
  m1=m1+x;

  for(i=0;i<6;i++) tps1[i]=m1[i];
  for(j=0;j<linename.Ncell; j++) linename.Cell[j]->DAPass(tps1);
  for(i=0;i<6;i++) m1[i]=tps1[i];
  Getmat(m1, m66); 
  for(i=0;i<36;i++) linename.Cell[linename.Ncell-1]->M[i]=m66[i];
  
  if( false ) {
    cout<<"One turn map M:"<<endl;
    for (i=0;i<6;i++) {
      for(j=0;j<6;j++) cout<<setw(12)<<m66[i*6+j]<<"  ";
      cout<<endl;
    }
    for(i=0;i<36;i++) b66[i]=m66[i];
    cout<<"Det of M = "<<mat_det(b66,6)<<endl;
  }

  mat_change_hessenberg(m66, 6);
  flag=mat_root_hessenberg(m66,6,u,v,1.0e-10,60);
  if(flag < 0){
    cout<<"One Turn Map: unstable.  exit. \n"<<endl;
    exit(1);
  }


  j=0;
  for(i=0;i<linename.Ncell;i++){
    if(linename.Cell[i]->TYPE==string("RFCAV")  ) {
      linename.Cell[i]->SetP("VRF", voltage[j]);
      j++;
    }
  }

    
}

void Cal_ElementMap_Z(Line & linename, double z0)
{
  int i,j,k, istart;
  double x[6];
  tps tps1[6];
  linmap m1;
  double t66[36];

  vector<double> voltage;
  for(i=0;i<linename.Ncell;i++){
    if(linename.Cell[i]->TYPE==string("RFCAV") ){
	voltage.push_back( linename.Cell[i]->GetP("VRF") ) ;
        linename.Cell[i]->SetP("VRF", 0.0);
      }
  }

  for(k=0; k<linename.Ncell; k++) 
    {
      istart=k-1;
      if (istart <0) istart=linename.Ncell-1;

      for(i=0;i<6;i++) x[i]=linename.Cell[istart]->X[i];
      x[4]=z0;
      x[5]=0.; 
      m1.identity();
      m1=m1+x;
      
      for(i=0;i<6;i++) tps1[i]=m1[i];
      linename.Cell[k]->DAPass(tps1);
      for(i=0;i<6;i++) m1[i]=tps1[i];
      Getmat(m1, t66);
      
      for(j=0;j<36;j++) linename.Cell[k]->T[j]=t66[j]; 
    }

  j=0;
  for(i=0;i<linename.Ncell;i++){
    if(linename.Cell[i]->TYPE==string("RFCAV")  ) {
      linename.Cell[i]->SetP("VRF", voltage[j]);
      j++;
    }
  }
  
}

void Cal_Twiss_Z(Line & linename, double z0)
//  here  deltp is always dp/p0 no matter GP.H_expand, 4-d  treatment
{
  int i,j,k;
  double pt, temp1, temp2,r;
  double A11[4], A12[4], A22[4], A22I[4], C[4];
  int    flag;
  
  Cal_Orbit_Num(linename,z0,0);
  Cal_OneTurnMap_Z(linename,z0);
  Cal_A(linename,0);
  Cal_ElementMap_Z(linename,z0);
  Trace_A(linename,0);

  for(k=0; k<linename.Ncell;k++){
    linename.Cell[k]->Beta1= linename.Cell[k]->A[0*6+0]/linename.Cell[k]->A[1*6+1];
    linename.Cell[k]->Alfa1=-linename.Cell[k]->A[1*6+0]/linename.Cell[k]->A[1*6+1];
    linename.Cell[k]->Beta2= linename.Cell[k]->A[2*6+2]/linename.Cell[k]->A[3*6+3];
    linename.Cell[k]->Alfa2=-linename.Cell[k]->A[3*6+2]/linename.Cell[k]->A[3*6+3];
    
    for(i=0;i<2;i++){  //  updated March 26, 2022
      for(j=0;j<2;j++) {
	A11[i*2+j]= linename.Cell[k]->A[6*i+j];
	A12[i*2+j]= linename.Cell[k]->A[6*i+(j+2)];
	A22[i*2+j]= linename.Cell[k]->A[6*(i+2)+(j+2)];
      }    
    }
    
    temp1=mat_det(A11,2);
    temp2=mat_det(A22,2);
    if( abs(temp1 -temp2) > 1.0e-10) {
      //----switched, quite annoying.
      //cout<<"Warning:  during calculating r Cal_Twiss(). "<<endl;
    }
    r=sqrt( mat_det( A22,2 ) );
    
    A22I[0]  =  A22[3]/( A22[0] * A22[3] );
    A22I[1]  = -A22[1]/( A22[0] * A22[3] );
    A22I[2]  = -A22[2]/( A22[0] * A22[3] );
    A22I[3]  =  A22[0]/( A22[0] * A22[3] );

    mat_mult(A12, A22I, C,2,2,2);
    for(i=0;i<4;i++) C[i]=C[i]*r;
    
    linename.Cell[k]->r=r;
    linename.Cell[k]->c11=C[0];
    linename.Cell[k]->c12=C[1];
    linename.Cell[k]->c21=C[2];
    linename.Cell[k]->c22=C[3];
  }
  
}

void Cal_Tune_vs_Z(Line & linename, const char *filename)
{
  int i;
  double qx[21],qy[21];
  double z0;
  fstream fout;

  for(i=0;i<21;i++){
    z0 =  GP.step_deltaz * (i-10);
    Cal_Twiss_Z(linename, z0);
    qx[i]=linename.Tune1;
    qy[i]=linename.Tune2;
  }

  fout.open(filename, ios::out);
  for (i = 0; i < 21; i++)
    fout << setw(10) <<GP.step_deltaz * (i-10)
  	 << scientific << setw(15) << qx[i]
  	 << scientific << setw(15) << qy[i]<<endl;
   fout.close();
}

void Cal_Beta_Star_vs_Z(Line & linename, const char *filename)
{
  int i;
  double betx[21],bety[21];
  double z0;
  fstream fout;

  for(i=0;i<21;i++){
    z0 = GP.step_deltaz * (i-10);
    Cal_Twiss_Z(linename, z0);
    betx[i]=linename.Cell[linename.Ncell-1]->Beta1;
    bety[i]=linename.Cell[linename.Ncell-1]->Beta2;
  }

  fout.open(filename, ios::out);
  for (i = 0; i < 21; i++)
    fout << setw(10) <<GP.step_deltaz * (i-10)
  	 << scientific << setw(15) << betx[i]
  	 << scientific << setw(15) << bety[i]<<endl;
   fout.close();
}

//-----trace the orbit through the ring to the end 
void Trace_Orbit(Line & linename, int istart, double x[6])
{
  int i,j;
  for(i=istart;i<linename.Ncell;i++){
    linename.Cell[i]->Pass(x);
    for(j=0;j<6;j++) linename.Cell[i]->X[j]=x[j];
  }
}

void Trace_Twiss(Line & linename, double deltap, double x[4], double Beta1, double Beta2, double Alfa1, double Alfa2, double c11, double c12, double c21, double c22)
{
  double r;
  double x0[6];
  int i,j,k;
  double pt,temp;
  double A12[4], A22[4], C[4];
  int  flag;

  double T[4][4], A[4][4], G[4][4];
  double dphi1, dphi2;
  double mu1, mu2;
  
  pt=DeltaToPt(deltap);
  
  for(i=0;i<4;i++)  x0[i]=x[i];
  x0[4]=0.000;
  x0[5]=pt;
  for(j=0;j<linename.Ncell; j++) {
    linename.Cell[j]->Pass(x0);
    for(i=0;i<6;i++) linename.Cell[j]->X[i]=x0[i];
  } 
  
  r= sqrt(1 - ( c11*c22-c12*c21));
  linename.Cell[linename.Ncell-1]->A[0*6+0]= r * sqrt(Beta1);
  linename.Cell[linename.Ncell-1]->A[0*6+1]= 0; 
  linename.Cell[linename.Ncell-1]->A[0*6+2]= c11*sqrt(Beta2)- c12*Alfa2/sqrt(Beta2);
  linename.Cell[linename.Ncell-1]->A[0*6+3]= c12 / sqrt(Beta2);
  linename.Cell[linename.Ncell-1]->A[1*6+0]= -Alfa1*r / sqrt(Beta1);
  linename.Cell[linename.Ncell-1]->A[1*6+1]= r/sqrt(Beta1); 
  linename.Cell[linename.Ncell-1]->A[1*6+2]= c21*sqrt(Beta2)-c22*Alfa2/sqrt(Beta2);
  linename.Cell[linename.Ncell-1]->A[1*6+3]= c22 / sqrt(Beta2);
  linename.Cell[linename.Ncell-1]->A[2*6+0]=  -c12*Alfa1 / sqrt(Beta1) - c22*sqrt(Beta1);
  linename.Cell[linename.Ncell-1]->A[2*6+1]=  c12/sqrt(Beta1); 
  linename.Cell[linename.Ncell-1]->A[2*6+2]= r*sqrt(Beta2);
  linename.Cell[linename.Ncell-1]->A[2*6+3]= 0.;
  linename.Cell[linename.Ncell-1]->A[3*6+0]= c11*Alfa1 / sqrt(Beta1) + c21*sqrt(Beta1);
  linename.Cell[linename.Ncell-1]->A[3*6+1]= -c11/sqrt(Beta1); 
  linename.Cell[linename.Ncell-1]->A[3*6+2]= -Alfa2*r/sqrt(Beta2);
  linename.Cell[linename.Ncell-1]->A[3*6+3]= r/sqrt(Beta2);
  
  Cal_ElementMap(linename, pt);
  Trace_A(linename, pt);

  for(k=0; k<linename.Ncell;k++){
    linename.Cell[k]->Beta1= linename.Cell[k]->A[0*6+0]/linename.Cell[k]->A[1*6+1];
    linename.Cell[k]->Alfa1=-linename.Cell[k]->A[1*6+0]/linename.Cell[k]->A[1*6+1];
    linename.Cell[k]->Beta2= linename.Cell[k]->A[2*6+2]/linename.Cell[k]->A[3*6+3];
    linename.Cell[k]->Alfa2=-linename.Cell[k]->A[3*6+2]/linename.Cell[k]->A[3*6+3];
    
    temp=linename.Cell[k]->A[0*6+0]*linename.Cell[k]->A[1*6+1];
    temp=temp*linename.Cell[k]->A[2*6+2]*linename.Cell[k]->A[3*6+3];
    temp=sqrt(sqrt(temp));
    linename.Cell[k]->r=temp;
    
    for(i=0;i<2;i++){
      for(j=0;j<2;j++) {
	A12[i*2+j]= linename.Cell[k]->A[6*i+(j+2)];
	A22[i*2+j]= linename.Cell[k]->A[6*(i+2)+(j+2)];
      }    
    }
    
    //mat_inv(A22,2);
    A22[0] = 1.0/A22[0];
    A22[1] = 0.;
    A22[2] =-A22[2];
    A22[3] = 1.0/A22[0];
    
    mat_mult(A12, A22,C,2,2,2);
    linename.Cell[k]->c11=temp*C[0];
    linename.Cell[k]->c12=temp*C[1];
    linename.Cell[k]->c21=temp*C[2];
    linename.Cell[k]->c22=temp*C[3];
  }
  
}

//-----6-d  optics calcualtion

void Cal_Orbit_Num_6D(Line & linename)
{
  int i,j,k,iter=0;
  double x0[6], x01[6], x1[6], dx[6];
  double d=1.0e-09;
  double mat[6][6];
  double chi;
  double codeps=1e-12;
  int Max_iter=20;

  vector<double> voltage;
  for(i=0;i<linename.Ncell;i++){
    if( linename.Cell[i]->TYPE==string("CRABRF") || linename.Cell[i]->TYPE==string("CCMULT") ){
	voltage.push_back( linename.Cell[i]->GetP("VRF") ) ;
        linename.Cell[i]->SetP("VRF", 0.0);
      }
  }
  
  for (i=0;i<6;i++) x0[i]=0.;

  do{
    iter++;
    
    for(j=0;j<6;j++) x01[j]=x0[j];
    for(j=0;j<linename.Ncell;j++)  linename.Cell[j]->Pass(x01);
    
    for(k=0;k<6;k++){
      for(j=0;j<6;j++) x1[j]=x0[j];
      x1[k]=x1[k]+d;
      for(j=0;j<linename.Ncell;j++) linename.Cell[j]->Pass(x1);
      for(j=0;j<6;j++) mat[j][k]=(x1[j]-x01[j])/d;
    }
    
    for(j=0;j<6;j++) x1[j]=x0[j];
    for(j=0;j<linename.Ncell;j++)  linename.Cell[j]->Pass(x1);
    for(i=0;i<6;i++) {
      dx[i]=x0[i]-x1[i];
    }
    chi=0;
    for(i=0;i<6;i++) chi+=dx[i]*dx[i];
    chi=sqrt(chi/6.);
    
    if( chi > codeps ){
      for(i=0;i<6;i++) mat[i][i]=mat[i][i]-1.0000001;
      mat_inv(&mat[0][0],6); 
      for(i=0;i<6;i++){
	for(j=0;j<6;j++) x0[i]+=mat[i][j]*dx[j];  
      } 
    } 
    
  }while (chi>codeps && iter <Max_iter );
  
  if(iter == Max_iter-1 ) 
    {
      cout<<"Failed to find COD."<<endl;
      exit(1);
    }
  else
    {
      for(j=0;j<linename.Ncell; j++) {
	linename.Cell[j]->Pass(x0);
	for(i=0;i<6;i++) linename.Cell[j]->X[i]=x0[i];
      } 
    }

  j=0;
  for(i=0;i<linename.Ncell;i++){
    if( linename.Cell[i]->TYPE==string("CRABRF") || linename.Cell[i]->TYPE==string("CCMULT")  ) {
      linename.Cell[i]->SetP("VRF", voltage[j]);
      j++;
    }
  }

  
}

void Cal_OneTurnMap_6D(Line & linename)
{
  int i,j;
  double x[6];
  tps tps1[6];
  linmap m1;
  double m66[36], b66[36];
  int flag;
  double u[6], v[6];

  for(i=0;i<6;i++) x[i]=linename.Cell[linename.Ncell-1]->X[i];
  m1.identity();
  m1=m1+x;
  
  for(i=0;i<6;i++) tps1[i]=m1[i];
  for(j=0;j<linename.Ncell; j++) linename.Cell[j]->DAPass(tps1);
  for(i=0;i<6;i++) m1[i]=tps1[i];
  Getmat(m1, m66); 

  for(i=0;i<36;i++) linename.Cell[linename.Ncell-1]->M[i]=m66[i];
  
  if(true) {
    cout<<"One turn map M:"<<endl;
    for (i=0;i<6;i++) {
      for(j=0;j<6;j++) cout<<setw(12)<<m66[i*6+j]<<"  ";
      cout<<endl;
    }
    for(i=0;i<36;i++) b66[i]=m66[i];
    cout<<"Det of M = "<<mat_det(b66,6)<<endl;
  }

  mat_change_hessenberg(m66, 6);
  flag=mat_root_hessenberg(m66,6,u,v,1.0e-10,60);
  if(flag < 0){
    cout<<"One Turn Map: unstable.  exit. \n"<<endl;
    exit(1);
  }
}

void Cal_OneTurnMap_Element_6D(Line & linename, int k)
{
  int i,j;
  double x[6];
  tps tps1[6];
  linmap m1;
  double m66[36], b66[36];
  int flag;
  double u[6], v[6];

  if(k<0 or k>linename.Ncell-1 ){
    cout<<" Index out of ring range. exit. "<<endl;
    return;
  }

  for(i=0;i<6;i++) x[i]=linename.Cell[k]->X[i];
  m1.identity();
  m1=m1+x;
  
  for(i=0;i<6;i++) tps1[i]=m1[i];
  for(j=k+1;j<linename.Ncell; j++) linename.Cell[j]->DAPass(tps1);
  for(j=0;j<k+1; j++) linename.Cell[j]->DAPass(tps1);  
  for(i=0;i<6;i++) m1[i]=tps1[i];
  Getmat(m1, m66); 
  
  for(i=0;i<36;i++) linename.Cell[k]->M[i]=m66[i];

}

void Cal_ElementMap_6D(Line & linename)
{
  int i,j,k, istart;
  double x[6];
  tps tps1[6];
  linmap m1;
  double t66[36];

  for(k=0; k<linename.Ncell; k++) 
    {
      istart=k-1;
      if (istart <0) istart=linename.Ncell-1;

      for(i=0;i<6;i++) x[i]=linename.Cell[istart]->X[i];
      m1.identity();
      m1=m1+x;
      
      for(i=0;i<6;i++) tps1[i]=m1[i];
      linename.Cell[k]->DAPass(tps1);
      for(i=0;i<6;i++) m1[i]=tps1[i];
      Getmat(m1, t66);
      
      for(j=0;j<36;j++) linename.Cell[k]->T[j]=t66[j]; 
    }
}

void Cal_A_6D(Line & linename)
{
  
  int i,j;
  double m66[6][6];
  double wr[6], wi[6], vr[6][6], vi[6][6], temp_wr, temp_wi, temp_vr[6], temp_vi[6];  
  double temp1, temp2, temp3;
  double theta1, theta2, theta3, tempx, tempy, x1,y1;
  double a66[6][6], b66[6][6], c44[4][4];
  double det, scale;
  
  //---solve eigen problem to construct A[6*6]
  for(i=0;i<6;i++)
   for(j=0;j<6;j++) m66[i][j]=linename.Cell[linename.Ncell-1]->M[i*6+j];
  EigenSolver_6D(m66, wr, wi, vr, vi);

  if(false){
    cout<<"Eigen values and eigen vectors: "<<endl;
    for(i=0;i<6;i++){
      cout<<"    "<<i<<"  :  "<<wr[i]<< "  +j  "<<wi[i]<<endl;
      for(j=0;j<6;j++){
        cout<<sqrt(vr[i][j]*vr[i][j] +vi[i][j]*vi[i][j]  ) <<"  "<<atan2(vi[i][j],vr[i][j])<<endl; 
	//cout<<vr[i][j]<<" +j "<<vi[i][j]<<endl;
      }
    }
  }

 //---link eigen vetcors S planes 
  temp1=abs(atan2(wi[0], wr[0]))/2/PI;
  temp2=abs(atan2(wi[2], wr[2]))/2/PI; 
  temp3=abs(atan2(wi[4], wr[4]))/2/PI;
  cout<<"Fractional Tunes: "<<temp1<<"  "<<temp2<<"   "<<temp3<<endl;
  if( temp1 < temp2 and temp1 < temp3) {  // assuming Qs < Qx and Qs < Qy, which is  true normally
      temp_wr    =  wr[4] ;
      temp_wi    =  wi[4] ;
      wr[4]      =  wr[0] ;
      wi[4]      =  wi[0] ;
      wr[0]      =  temp_wr;
      wi[0]      =  temp_wi;
      wr[1]      =  wr[0];
      wi[1]      = -wi[0]; 
      wr[5]      =  wr[4];
      wi[5]      = -wi[4]; 
    for(i=0;i<6;i++){
      temp_vr[i]    = vr[4][i];
      temp_vi[i]    = vi[4][i];
      vr[4][i]      = vr[0][i];
      vi[4][i]      = vi[0][i];
      vr[0][i]      = temp_vr[i];
      vi[0][i]      = temp_vi[i];
      vr[1][i]      = vr[0][i];
      vi[1][i]      =-vi[0][i]; 
      vr[5][i]      = vr[4][i];
      vi[5][i]      =-vi[4][i]; 
    }
  }
  
  if( temp2 < temp1 and temp2 < temp3) {
      temp_wr    =  wr[4] ;
      temp_wi    =  wi[4] ;
      wr[4]      =  wr[2] ;
      wi[4]      =  wi[2] ;
      wr[2]      =  temp_wr;
      wi[2]      =  temp_wi;
      wr[3]      = wr[2];
      wi[3]      =-wi[2]; 
      wr[5]      = wr[4];
      wi[5]      =-wi[4]; 
    for(i=0;i<6;i++){
      temp_vr[i]    = vr[4][i];
      temp_vi[i]    = vi[4][i];
      vr[4][i]      = vr[2][i];
      vi[4][i]      = vi[2][i];
      vr[2][i]      = temp_vr[i];
      vi[2][i]      = temp_vi[i];
      vr[3][i]      = vr[2][i];
      vi[3][i]      =-vi[2][i]; 
      vr[5][i]      = vr[4][i];
      vi[5][i]      =-vi[4][i]; 
    }
  }

  if(false){
    cout<<"Eigen values and eigen vectors after s mode fixed: "<<endl;
    for(i=0;i<6;i++){
      cout<<"    "<<i<<"  :  "<<wr[i]<< "  +j  "<<wi[i]<<endl;
      for(j=0;j<6;j++) cout<<vr[i][j]<<" +j "<<vi[i][j]<<endl;
    }
  }

 //---link eigen vetcors H/V planes 
  //temp1=abs(vr[0][0] ) + abs( vi[0][0] ) + abs( vr[0][1] ) + abs( vi[0][1] );  // horizontal 'projection'
  //temp2=abs(vr[2][0] ) + abs( vi[2][0] ) + abs( vr[2][1] ) + abs( vi[2][1] );  // vertical 'projection'
  temp1=abs(vr[0][0] ) + abs( vi[0][0] ) ;  // mode I  to horizontal 'projection'
  temp2=abs(vr[2][0] ) + abs( vi[2][0] ) ;  // mode II to horizontal 'projection'
  
  if (temp2 > temp1 ) {     
      temp_wr    =  wr[0] ;
      temp_wi    =  wi[0] ;
      wr[0]      =  wr[2] ;
      wi[0]      =  wi[2] ;
      wr[2]      =  temp_wr;
      wi[2]      =  temp_wi;
      wr[1]      =  wr[0];
      wi[1]      = -wi[0]; 
      wr[3]      =  wr[2];
      wi[3]      = -wi[2]; 
    for(i=0;i<6;i++){
      temp_vr[i]    = vr[0][i];
      temp_vi[i]    = vi[0][i];
      vr[0][i]      = vr[2][i];
      vi[0][i]      = vi[2][i];
      vr[2][i]      = temp_vr[i];
      vi[2][i]      = temp_vi[i];
      vr[1][i]      = vr[0][i];
      vi[1][i]      =-vi[0][i]; 
      vr[3][i]      = vr[2][i];
      vi[3][i]      =-vi[2][i]; 
    }
  }

  if(false){
    cout<<"Eigen values and eigen vectors after H/V/S modes fixed: "<<endl;
    for(i=0;i<6;i++){
      cout<<"    "<<i<<"  :  "<<wr[i]<< "  +j  "<<wi[i]<<endl;
      for(j=0;j<6;j++){
	cout<<sqrt(vr[i][j]*vr[i][j] +vi[i][j]*vi[i][j]  ) <<"  "<<atan2(vi[i][j],vr[i][j])<<endl;
	//cout<<vr[i][j]<<" +j "<<vi[i][j]<<endl;
      }
    }
  }

  //----rotating to make A12=0, A34=0, A56=0
  theta1=atan2(vi[0][0], vr[0][0]);
  x1  =cos(-PI/2-theta1);
  y1  =sin(-PI/2-theta1);
  for(i=0;i<6;i++){
    tempx=vr[0][i]*x1-vi[0][i]*y1;
    tempy=vi[0][i]*x1+vr[0][i]*y1;
    vr[0][i]=tempx;
    vi[0][i]=tempy; }
  x1 =cos(PI/2+theta1);
  y1 =sin(PI/2+theta1);
  for(i=0;i<6;i++){
    tempx=vr[1][i]*x1-vi[1][i]*y1;
    tempy=vi[1][i]*x1+vr[1][i]*y1;
    vr[1][i]=tempx;
    vi[1][i]=tempy;}

  theta2=atan2(vi[2][2], vr[2][2]);
  x1 =cos(-PI/2-theta2);
  y1 =sin(-PI/2-theta2);
  for(i=0;i<6;i++){
    tempx=vr[2][i]*x1-vi[2][i]*y1;
    tempy=vi[2][i]*x1+vr[2][i]*y1;
    vr[2][i]=tempx;
    vi[2][i]=tempy; }
  x1=cos(PI/2+theta2);
  y1=sin(PI/2+theta2);
  for(i=0;i<6;i++){
    tempx=vr[3][i]*x1-vi[3][i]*y1;
    tempy=vi[3][i]*x1+vr[3][i]*y1;
    vr[3][i]=tempx;
    vi[3][i]=tempy;}

  theta3=atan2(vi[4][4], vr[4][4]);
  x1 =cos(-PI/2-theta3);
  y1 =sin(-PI/2-theta3);
  for(i=0;i<6;i++){
    tempx=vr[4][i]*x1-vi[4][i]*y1;
    tempy=vi[4][i]*x1+vr[4][i]*y1;
    vr[4][i]=tempx;
    vi[4][i]=tempy; }
  x1=cos(PI/2+theta3);
  y1=sin(PI/2+theta3);
  for(i=0;i<6;i++){
    tempx=vr[5][i]*x1-vi[5][i]*y1;
    tempy=vi[5][i]*x1+vr[5][i]*y1;
    vr[5][i]=tempx;
    vi[5][i]=tempy;}

  if(false){
    cout<<"Eigen values and eigen vectors after phase fixed: "<<endl;
    for(i=0;i<6;i++){
      cout<<"    "<<i<<"  :  "<<wr[i]<< "  +j  "<<wi[i]<<endl;
      for(j=0;j<6;j++){
	cout<<sqrt(vr[i][j]*vr[i][j] +vi[i][j]*vi[i][j]  ) <<"  "<<atan2(vi[i][j],vr[i][j])<<endl;
	//cout<<vr[i][j]<<" +j "<<vi[i][j]<<endl;
      }
    }
  }

  //----produce A and normalize
  for(i=0;i<6;i++) a66[0][i]=-(vi[0][i]-vi[1][i]);
  for(i=0;i<6;i++) a66[1][i]=  vr[0][i]+vr[1][i];
  for(i=0;i<6;i++) a66[2][i]=-(vi[2][i]-vi[3][i]);
  for(i=0;i<6;i++) a66[3][i]=  vr[2][i]+vr[3][i];
  for(i=0;i<6;i++) a66[4][i]=-(vi[4][i]-vi[5][i]);
  for(i=0;i<6;i++) a66[5][i]=  vr[4][i]+vr[5][i];
  
  if(true){
    if(a66[1][1] < 0. ) {
      for(i=0;i<6;i++) a66[i][1]= -a66[i][1];
      cout<<"Mode I : negative beta reverted."<<endl; 
   }
    if(a66[3][3] < 0. ) {
      for(i=0;i<6;i++) a66[i][3]= -a66[i][3];
      cout<<"Mode II : negative beta reverted."<<endl; 
    }
    if(a66[5][5] < 0. ) {
      for(i=0;i<6;i++) a66[i][5]= -a66[i][5];
      cout<<"Mode III: negative beta reverted."<<endl; 
    }
  }

  if(false) {
    cout<<"A matrix :"<<endl;
    for (i=0;i<6;i++) {
      for(j=0;j<6;j++) cout<<setw(12)<<a66[i][j]<<"  ";
      cout<<endl;
    }

    for(i=0;i<6;i++) {
      for(j=0;j<6;j++) b66[i][j]=a66[i][j];
    }
    det= mat_det( &b66[0][0], 6);
    cout<<" Det of A66  =  "<< det <<endl;
  }
  
  //---normalizing A
   
  double  s66[6][6];
  for(i=0;i<6;i++)
    for(j=0;j<6;j++) s66[i][j]=0.;
  s66[0][1]= 1.0;
  s66[1][0]=-1.0;
  s66[2][3]= 1.0;
  s66[3][2]=-1.0;  
  s66[4][5]= 1.0;
  s66[5][4]=-1.0;  

  if(false){
    cout<<"S matrix :"<<endl;
    for (i=0;i<6;i++) {
      for(j=0;j<6;j++) cout<<setw(12)<<s66[i][j]<<"  ";
      cout<<endl;
    }
  }

  temp1=(-a66[0][1]*a66[1][0] + a66[0][0]*a66[1][1] -a66[0][3]*a66[1][2] +a66[0][2]*a66[1][3] - a66[0][5]*a66[1][4] +a66[0][4]*a66[1][5]);
  temp2=(-a66[2][1]*a66[3][0] + a66[2][0]*a66[3][1] -a66[2][3]*a66[3][2] +a66[2][2]*a66[3][3] - a66[2][5]*a66[3][4] +a66[2][4]*a66[3][5]);  
  temp3=(-a66[4][1]*a66[5][0] + a66[4][0]*a66[5][1] -a66[4][3]*a66[5][2] +a66[4][2]*a66[5][3] - a66[4][5]*a66[5][4] +a66[4][4]*a66[5][5]);

  //cout<<"...normalize new v1 and v2, v3 adn v4, v5 and v6 :"<<endl;
  //cout<<temp1<<" "<<temp2<<" "<<temp3<<endl;
  //temp3=-temp3;
  
  for(i=0;i<6;i++){
    a66[0][i]=a66[0][i]/sqrt(temp1);
    a66[1][i]=a66[1][i]/sqrt(temp1);
    a66[2][i]=a66[2][i]/sqrt(temp2);
    a66[3][i]=a66[3][i]/sqrt(temp2);
    a66[4][i]=a66[4][i]/sqrt(temp3);
    a66[5][i]=a66[5][i]/sqrt(temp3);
  }

  if(false){
    temp1=(-a66[0][1]*a66[1][0] + a66[0][0]*a66[1][1] -a66[0][3]*a66[1][2] +a66[0][2]*a66[1][3] - a66[0][5]*a66[1][4] +a66[0][4]*a66[1][5]);
    temp2=(-a66[2][1]*a66[3][0] + a66[2][0]*a66[3][1] -a66[2][3]*a66[3][2] +a66[2][2]*a66[3][3] - a66[2][5]*a66[3][4] +a66[2][4]*a66[3][5]);  
    temp3=(-a66[4][1]*a66[5][0] + a66[4][0]*a66[5][1] -a66[4][3]*a66[5][2] +a66[4][2]*a66[5][3] - a66[4][5]*a66[5][4] +a66[4][4]*a66[5][5]);

    cout<<"...normalize new v1 and v2, v3 adn v4, v5 and v6 :"<<endl;
    cout<<temp1<<" "<<temp2<<" "<<temp3<<endl;
  }

  if(true) {
    cout<<"A matrix :"<<endl;
    for (i=0;i<6;i++) {
      for(j=0;j<6;j++) cout<<setw(12)<<a66[j][i]<<"  ";
      cout<<endl;
    }
    
    for(i=0;i<6;i++) {
      for(j=0;j<6;j++) b66[i][j]=a66[j][i];
    }
    det= mat_det( &b66[0][0], 6);
    cout<<" Det of A  =  "<< det <<endl;
  }

  //----save 
  for(i=0;i<6;i++)
    for(j=0;j<6;j++)
     linename.Cell[linename.Ncell-1]->A[i*6+j]=a66[j][i];
  
}

void Trace_A_6D(Line & linename)
{
  int i,j,k;
  double T[6][6], A[6][6], G[6][6], RI[6][6];
  double dphi1, dphi2, dphi3;
  double mu1, mu2, mu3;
 
  for(i=0;i<6;i++)
    for(j=0;j<6;j++) A[i][j]= linename.Cell[linename.Ncell-1]->A[i*6+j];

  mu1=0.;
  mu2=0.;
  mu3=0.;
  for(k=0;k<linename.Ncell;k++){   
    for (i=0;i<6;i++){
      for(j=0;j<6;j++) 	T[i][j]=linename.Cell[k]->T[i*6+j];
    }
  
    mat_mult( &T[0][0], &A[0][0], &G[0][0], 6,6,6 );
    dphi1=atan2(G[0][1], G[0][0]);
    dphi2=atan2(G[2][3], G[2][2]);
    dphi3=atan2(G[4][5], G[4][4]);

    //cout<<k<<"  "<<linename.Cell[k]->NAME<<" "<<linename.Cell[k]->TYPE<<"  "<<dphi1<<" "<<dphi2<<"  "<<dphi3<<endl;

    for(i=0;i<6;i++)
      for(j=0;j<6;j++) RI[i][j]=0;

    RI[0][0]  =  cos(dphi1);   
    RI[0][1]  = -sin(dphi1);  
    RI[1][0]  =  sin(dphi1);  
    RI[1][1]  =  cos(dphi1);

    RI[2][2]  =  cos(dphi2);   
    RI[2][3]  = -sin(dphi2);  
    RI[3][2]  =  sin(dphi2);  
    RI[3][3]  =  cos(dphi2);

    RI[4][4]  =  cos(dphi3);   
    RI[4][5]  = -sin(dphi3);  
    RI[5][4]  =  sin(dphi3);  
    RI[5][5]  =  cos(dphi3);

   mat_mult( &G[0][0], &RI[0][0], &A[0][0], 6,6,6 );
   for(i=0;i<6;i++){
     for(j=0;j<6;j++) linename.Cell[k]->A[i*6+j]=A[i][j];
   }
    
    mu1=mu1+dphi1/ 2./PI;
    mu2=mu2+dphi2/ 2./PI;
    mu3=mu3+dphi3/ 2./PI;
    linename.Cell[k]->Mu1= mu1;
    linename.Cell[k]->Mu2= mu2;
    linename.Cell[k]->Mu3= mu3;
  }

  linename.Tune1=  mu1;
  linename.Tune2=  mu2;
  linename.Tune3=  mu3;
}

void Cal_Twiss_6D(Line & linename)
{
  int i,j,k;
  double temp;
  double A12[4], A22[4], C[4];

  Cal_Orbit_Num_6D(linename);
  Cal_OneTurnMap_6D(linename);
  Cal_A_6D(linename);
  Cal_ElementMap_6D(linename);
  Trace_A_6D(linename);
  
  for(k=0; k<linename.Ncell;k++){
    linename.Cell[k]->Beta1= linename.Cell[k]->A[0*6+0]/linename.Cell[k]->A[1*6+1];
    linename.Cell[k]->Alfa1=-linename.Cell[k]->A[1*6+0]/linename.Cell[k]->A[1*6+1];
    linename.Cell[k]->Beta2= linename.Cell[k]->A[2*6+2]/linename.Cell[k]->A[3*6+3];
    linename.Cell[k]->Alfa2=-linename.Cell[k]->A[3*6+2]/linename.Cell[k]->A[3*6+3];
    linename.Cell[k]->Beta3= linename.Cell[k]->A[4*6+4]/linename.Cell[k]->A[5*6+5];
    linename.Cell[k]->Alfa3=-linename.Cell[k]->A[5*6+4]/linename.Cell[k]->A[5*6+5];
    
    temp=linename.Cell[k]->A[0*6+0]*linename.Cell[k]->A[1*6+1];
    temp=temp*linename.Cell[k]->A[2*6+2]*linename.Cell[k]->A[3*6+3];
    temp=sqrt(sqrt(temp));
    linename.Cell[k]->r=temp;
    
    for(i=0;i<2;i++){
      for(j=0;j<2;j++) {
	A12[i*2+j]= linename.Cell[k]->A[6*i+(j+2)];
	A22[i*2+j]= linename.Cell[k]->A[6*(i+2)+(j+2)];
      }    
    }
    
    //mat_inv(A22,2);
    
    //C[1]=-A22[1];
    //C[2]=-A22[2];
    //C[3]= A22[0];
    //det = (A22[0]*A22[3]-A22[1]*A22[2]);  
    //for(i=0;i<3;i++) A22[i]= C[0] / det ;
    //if(det==0){
    //  cout<<" Can't invert A22 in Cal_Twiss(). "<<endl;
    //  exit(1);
    //}
    
    A22[0] = 1.0/A22[0];
    A22[1] = 0.;
    A22[2] =-A22[2];
    A22[3] = 1.0/A22[0];

    mat_mult(A12, A22,C,2,2,2);
    linename.Cell[k]->c11=temp*C[0];
    linename.Cell[k]->c12=temp*C[1];
    linename.Cell[k]->c21=temp*C[2];
    linename.Cell[k]->c22=temp*C[3];
  }
}

void Cal_A_4D(double m[], double a[])
{
  int i,j;
  double m44[4][4];
  double wr[4], wi[4], vr[4][4], vi[4][4], temp_wr, temp_wi, temp_vr[4], temp_vi[4];  
  double temp1, temp2;
  double theta1, theta2, tempx, tempy, x1,y1;
  double a44[4][4], b44[4][4];
  double det, scale;
  
  //---solve eigen problem to construct A[4*4]
  for(i=0;i<4;i++)
    for(j=0;j<4;j++) m44[i][j]=m[i*4+j];


    cout<<"M matrix read in :"<<endl;
    for (i=0;i<4;i++) {
      for(j=0;j<4;j++) cout<<setw(12)<<m44[i][j]<<"  ";
      cout<<endl;
    }
  
  EigenSolver(m44, wr, wi, vr, vi);

  //---link eigen vetcors to  H/V planes 
  temp1=abs(vr[0][0] ) + abs( vi[0][0] ) + abs( vr[0][1] ) + abs( vi[0][1] );
  temp2=abs(vr[2][0] ) + abs( vi[2][0] ) + abs( vr[2][1] ) + abs( vi[2][1] );
  if (temp2 > temp1 ) {
      temp_wr    =  wr[0] ;
      temp_wi    =  wi[0] ;
      wr[0]      =  wr[2] ;
      wi[0]      =  wi[2] ;
      wr[2]      =  temp_wr;
      wi[2]      =  temp_wi;
      wr[1]      =  wr[0];
      wi[1]      = -wi[0]; 
      wr[3]      =  wr[2];
      wi[3]      = -wi[2]; 
    for(i=0;i<4;i++){
      temp_vr[i]    = vr[0][i];
      temp_vi[i]    = vi[0][i];
      vr[0][i]      = vr[2][i];
      vi[0][i]      = vi[2][i];
      vr[2][i]      = temp_vr[i];
      vi[2][i]      = temp_vi[i];
      vr[1][i]      = vr[0][i];
      vi[1][i]      =-vi[0][i]; 
      vr[3][i]      = vr[2][i];
      vi[3][i]      =-vi[2][i]; 
    }
  }

  //------check eigentunes
  if(true){
    temp1=abs(atan2(wi[0], wr[0]))/2/PI;
    temp2=abs(atan2(wi[2], wr[2]))/2/PI; 
    cout<<"Fractional Tunes: "<<temp1<<"  "<<temp2<<endl;
    cout<<"-------------------------------"<<endl;
  }
  
  //----rotating to make A12=0, A34=0
  theta1=atan2(vi[0][0], vr[0][0]);
  x1  = cos(-PI/2-theta1);
  y1  =sin(-PI/2-theta1);
  for(i=0;i<4;i++){
    tempx=vr[0][i]*x1-vi[0][i]*y1;
    tempy=vi[0][i]*x1+vr[0][i]*y1;
    vr[0][i]=tempx;
    vi[0][i]=tempy; }
  x1  = cos(PI/2+theta1);
  y1 =sin(PI/2+theta1);
  for(i=0;i<4;i++){
    tempx=vr[1][i]*x1-vi[1][i]*y1;
    tempy=vi[1][i]*x1+vr[1][i]*y1;
    vr[1][i]=tempx;
    vi[1][i]=tempy;}

  theta2=atan2(vi[2][2], vr[2][2]);
  x1  = cos(-PI/2-theta2);
  y1 =sin(-PI/2-theta2);
  for(i=0;i<4;i++){
    tempx=vr[2][i]*x1-vi[2][i]*y1;
    tempy=vi[2][i]*x1+vr[2][i]*y1;
    vr[2][i]=tempx;
    vi[2][i]=tempy; }
  x1  = cos(PI/2+theta2);
  y1=sin(PI/2+theta2);
  for(i=0;i<4;i++){
    tempx=vr[3][i]*x1-vi[3][i]*y1;
    tempy=vi[3][i]*x1+vr[3][i]*y1;
    vr[3][i]=tempx;
    vi[3][i]=tempy;}

  //----produce A
  for(i=0;i<4;i++) a44[0][i]=-(vi[0][i]-vi[1][i]);
  for(i=0;i<4;i++) a44[1][i]=  vr[0][i]+vr[1][i];
  for(i=0;i<4;i++) a44[2][i]=-(vi[2][i]-vi[3][i]);
  for(i=0;i<4;i++) a44[3][i]=  vr[2][i]+vr[3][i];

  
  if(false){
    if(a44[1][1] < 0. ) {
      for(i=0;i<4;i++) a44[i][1]=- a44[i][1];
      cout<<"Mode I : negative beta reverted."<<endl; 
   }
    if(a44[3][3] < 0. ) {
      for(i=0;i<4;i++) a44[i][3]=- a44[i][3];
      cout<<"Mode II : negative beta reverted."<<endl; 
    }
  }
  
  if(false) {
    cout<<"A matrix (before normalizing):"<<endl;
    for (i=0;i<4;i++) {
      for(j=0;j<4;j++) cout<<setw(12)<<a44[i][j]<<"  ";
      cout<<endl;
    }

    for(i=0;i<4;i++) {
      for(j=0;j<4;j++) b44[i][j]=a44[i][j];
    }
    det= mat_det( &b44[0][0], 4);
    cout<<" Det of A44  =  "<< det <<endl;
  }
  
  //---normalizing A
   
  double  s44[4][4];
  for(i=0;i<4;i++)
    for(j=0;j<4;j++) s44[i][j]=0.;
  s44[0][1]= 1.0;
  s44[1][0]=-1.0;
  s44[2][3]= 1.0;
  s44[3][2]=-1.0;  

  if(false){
    cout<<"S matrix :"<<endl;
    for (i=0;i<4;i++) {
      for(j=0;j<4;j++) cout<<setw(12)<<s44[i][j]<<"  ";
      cout<<endl;
    }
  }

  temp1=(-a44[0][1]*a44[1][0] + a44[0][0]*a44[1][1] -a44[0][3]*a44[1][2] +a44[0][2]*a44[1][3] );
  temp2=(-a44[2][1]*a44[3][0] + a44[2][0]*a44[3][1] -a44[2][3]*a44[3][2] +a44[2][2]*a44[3][3] );  

  //cout<<"...normalize new v1 and v2, v3 adn v4:"<<endl;
  //cout<<temp1<<" "<<temp2<<endl;
  
  for(i=0;i<4;i++){
    a44[0][i]=a44[0][i]/sqrt(temp1);
    a44[1][i]=a44[1][i]/sqrt(temp1);
    a44[2][i]=a44[2][i]/sqrt(temp2);
    a44[3][i]=a44[3][i]/sqrt(temp2);
  }
  
  temp1=(-a44[0][1]*a44[1][0] + a44[0][0]*a44[1][1] -a44[0][3]*a44[1][2] +a44[0][2]*a44[1][3] ); 
  temp2=(-a44[2][1]*a44[3][0] + a44[2][0]*a44[3][1] -a44[2][3]*a44[3][2] +a44[2][2]*a44[3][3] );

  //cout<<"...normalize new v1 and v2, v3 adn v4:"<<endl;
  //cout<<temp1<<" "<<temp2<<endl; 

  if(true) {
    cout<<"A matrix (after normalizing) :"<<endl;
    for (i=0;i<4;i++) {
      for(j=0;j<4;j++) cout<<setw(12)<<a44[j][i]<<"  ";
      cout<<endl;
    }
    
    for(i=0;i<4;i++) {
      for(j=0;j<4;j++) b44[i][j]=a44[j][i];
    }
    det= mat_det( &b44[0][0], 4);
    cout<<" Det of A  =  "<< det <<endl;
  }

  //----return
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      a[i*4+j]=  a44[j][i]; 
}

void Cal_A_6D(double m[], double a[])
{
  int i,j;
  double m66[6][6];
  double wr[6], wi[6], vr[6][6], vi[6][6], temp_wr, temp_wi, temp_vr[6], temp_vi[6];  
  double temp1, temp2, temp3;
  double theta1, theta2, theta3, tempx, tempy, x1,y1;
  double a66[6][6], b66[6][6], c44[4][4];
  double det, scale;


    cout<<"M3 matrix :"<<endl;
    for (i=0;i<6;i++) {
      for(j=0;j<6;j++) cout<<setw(12)<<m[i*6+j]<<"  ";
      cout<<endl;
    }

  //---solve eigen problem to construct A[6*6]
  for(i=0;i<6;i++)
    for(j=0;j<6;j++)
      m66[i][j]=m[i*6+j];
  
  EigenSolver_6D(m66, wr, wi, vr, vi);
  if(false){
    cout<<"Eigen values and eigen vectors: "<<endl;
    for(i=0;i<6;i++){
      cout<<"    "<<i<<"  :  "<<wr[i]<< "  +j  "<<wi[i]<<endl;
      for(j=0;j<6;j++){
        cout<<sqrt(vr[i][j]*vr[i][j] +vi[i][j]*vi[i][j]  ) <<"  "<<atan2(vi[i][j],vr[i][j])<<endl; 
	//cout<<vr[i][j]<<" +j "<<vi[i][j]<<endl;
      }
    }
  }

 //---link eigen vetcors S planes 
  temp1=abs(atan2(wi[0], wr[0]))/2/PI;
  temp2=abs(atan2(wi[2], wr[2]))/2/PI; 
  temp3=abs(atan2(wi[4], wr[4]))/2/PI;
  cout<<"Fractional Tunes: "<<temp1<<"  "<<temp2<<"   "<<temp3<<endl;
  if( temp1 < temp2 and temp1 < temp3) {  // assuming Qs < Qx and Qs < Qy, which is  true normally
      temp_wr    =  wr[4] ;
      temp_wi    =  wi[4] ;
      wr[4]      =  wr[0] ;
      wi[4]      =  wi[0] ;
      wr[0]      =  temp_wr;
      wi[0]      =  temp_wi;
      wr[1]      =  wr[0];
      wi[1]      = -wi[0]; 
      wr[5]      =  wr[4];
      wi[5]      = -wi[4]; 
    for(i=0;i<6;i++){
      temp_vr[i]    = vr[4][i];
      temp_vi[i]    = vi[4][i];
      vr[4][i]      = vr[0][i];
      vi[4][i]      = vi[0][i];
      vr[0][i]      = temp_vr[i];
      vi[0][i]      = temp_vi[i];
      vr[1][i]      = vr[0][i];
      vi[1][i]      =-vi[0][i]; 
      vr[5][i]      = vr[4][i];
      vi[5][i]      =-vi[4][i]; 
    }
  }
  
  if( temp2 < temp1 and temp2 < temp3) {
      temp_wr    =  wr[4] ;
      temp_wi    =  wi[4] ;
      wr[4]      =  wr[2] ;
      wi[4]      =  wi[2] ;
      wr[2]      =  temp_wr;
      wi[2]      =  temp_wi;
      wr[3]      = wr[2];
      wi[3]      =-wi[2]; 
      wr[5]      = wr[4];
      wi[5]      =-wi[4]; 
    for(i=0;i<6;i++){
      temp_vr[i]    = vr[4][i];
      temp_vi[i]    = vi[4][i];
      vr[4][i]      = vr[2][i];
      vi[4][i]      = vi[2][i];
      vr[2][i]      = temp_vr[i];
      vi[2][i]      = temp_vi[i];
      vr[3][i]      = vr[2][i];
      vi[3][i]      =-vi[2][i]; 
      vr[5][i]      = vr[4][i];
      vi[5][i]      =-vi[4][i]; 
    }
  }

  if(false){
    cout<<"Eigen values and eigen vectors after s mode fixed: "<<endl;
    for(i=0;i<6;i++){
      cout<<"    "<<i<<"  :  "<<wr[i]<< "  +j  "<<wi[i]<<endl;
      for(j=0;j<6;j++) cout<<vr[i][j]<<" +j "<<vi[i][j]<<endl;
    }
  }

 //---link eigen vetcors H/V planes 
  //temp1=abs(vr[0][0] ) + abs( vi[0][0] ) + abs( vr[0][1] ) + abs( vi[0][1] );  // horizontal 'projection'
  //temp2=abs(vr[2][0] ) + abs( vi[2][0] ) + abs( vr[2][1] ) + abs( vi[2][1] );  // vertical 'projection'
  temp1=abs(vr[0][0] ) + abs( vi[0][0] ) ;  // mode I  to horizontal 'projection'
  temp2=abs(vr[2][0] ) + abs( vi[2][0] ) ;  // mode II to horizontal 'projection'
  
  if (temp2 > temp1 ) {     
      temp_wr    =  wr[0] ;
      temp_wi    =  wi[0] ;
      wr[0]      =  wr[2] ;
      wi[0]      =  wi[2] ;
      wr[2]      =  temp_wr;
      wi[2]      =  temp_wi;
      wr[1]      =  wr[0];
      wi[1]      = -wi[0]; 
      wr[3]      =  wr[2];
      wi[3]      = -wi[2]; 
    for(i=0;i<6;i++){
      temp_vr[i]    = vr[0][i];
      temp_vi[i]    = vi[0][i];
      vr[0][i]      = vr[2][i];
      vi[0][i]      = vi[2][i];
      vr[2][i]      = temp_vr[i];
      vi[2][i]      = temp_vi[i];
      vr[1][i]      = vr[0][i];
      vi[1][i]      =-vi[0][i]; 
      vr[3][i]      = vr[2][i];
      vi[3][i]      =-vi[2][i]; 
    }
  }

  if(false){
    cout<<"Eigen values and eigen vectors after H/V/S modes fixed: "<<endl;
    for(i=0;i<6;i++){
      cout<<"    "<<i<<"  :  "<<wr[i]<< "  +j  "<<wi[i]<<endl;
      for(j=0;j<6;j++){
	cout<<sqrt(vr[i][j]*vr[i][j] +vi[i][j]*vi[i][j]  ) <<"  "<<atan2(vi[i][j],vr[i][j])<<endl;
	//cout<<vr[i][j]<<" +j "<<vi[i][j]<<endl;
      }
    }
  }

  //----rotating to make A12=0, A34=0, A56=0
  theta1=atan2(vi[0][0], vr[0][0]);
  x1  =cos(-PI/2-theta1);
  y1  =sin(-PI/2-theta1);
  for(i=0;i<6;i++){
    tempx=vr[0][i]*x1-vi[0][i]*y1;
    tempy=vi[0][i]*x1+vr[0][i]*y1;
    vr[0][i]=tempx;
    vi[0][i]=tempy; }
  x1 =cos(PI/2+theta1);
  y1 =sin(PI/2+theta1);
  for(i=0;i<6;i++){
    tempx=vr[1][i]*x1-vi[1][i]*y1;
    tempy=vi[1][i]*x1+vr[1][i]*y1;
    vr[1][i]=tempx;
    vi[1][i]=tempy;}

  theta2=atan2(vi[2][2], vr[2][2]);
  x1 =cos(-PI/2-theta2);
  y1 =sin(-PI/2-theta2);
  for(i=0;i<6;i++){
    tempx=vr[2][i]*x1-vi[2][i]*y1;
    tempy=vi[2][i]*x1+vr[2][i]*y1;
    vr[2][i]=tempx;
    vi[2][i]=tempy; }
  x1=cos(PI/2+theta2);
  y1=sin(PI/2+theta2);
  for(i=0;i<6;i++){
    tempx=vr[3][i]*x1-vi[3][i]*y1;
    tempy=vi[3][i]*x1+vr[3][i]*y1;
    vr[3][i]=tempx;
    vi[3][i]=tempy;}

  theta3=atan2(vi[4][4], vr[4][4]);
  x1 =cos(-PI/2-theta3);
  y1 =sin(-PI/2-theta3);
  for(i=0;i<6;i++){
    tempx=vr[4][i]*x1-vi[4][i]*y1;
    tempy=vi[4][i]*x1+vr[4][i]*y1;
    vr[4][i]=tempx;
    vi[4][i]=tempy; }
  x1=cos(PI/2+theta3);
  y1=sin(PI/2+theta3);
  for(i=0;i<6;i++){
    tempx=vr[5][i]*x1-vi[5][i]*y1;
    tempy=vi[5][i]*x1+vr[5][i]*y1;
    vr[5][i]=tempx;
    vi[5][i]=tempy;}

  if(false){
    cout<<"Eigen values and eigen vectors after[ phase fixed: "<<endl;
    for(i=0;i<6;i++){
      cout<<"    "<<i<<"  :  "<<wr[i]<< "  +j  "<<wi[i]<<endl;
      for(j=0;j<6;j++){
	cout<<sqrt(vr[i][j]*vr[i][j] +vi[i][j]*vi[i][j]  ) <<"  "<<atan2(vi[i][j],vr[i][j])<<endl;
	//cout<<vr[i][j]<<" +j "<<vi[i][j]<<endl;
      }
    }
  }

  //----produce A and normalize
  for(i=0;i<6;i++) a66[0][i]=-(vi[0][i]-vi[1][i]);
  for(i=0;i<6;i++) a66[1][i]=  vr[0][i]+vr[1][i];
  for(i=0;i<6;i++) a66[2][i]=-(vi[2][i]-vi[3][i]);
  for(i=0;i<6;i++) a66[3][i]=  vr[2][i]+vr[3][i];
  for(i=0;i<6;i++) a66[4][i]=-(vi[4][i]-vi[5][i]);
  for(i=0;i<6;i++) a66[5][i]=  vr[4][i]+vr[5][i];
  
  if(true){
    if(a66[1][1] < 0. ) {
      for(i=0;i<6;i++) a66[i][1]=- a66[i][1];
      //cout<<"Mode I : negative beta reverted."<<endl; 
   }
    if(a66[3][3] < 0. ) {
      for(i=0;i<6;i++) a66[i][3]=- a66[i][3];
      //cout<<"Mode II : negative beta reverted."<<endl; 
    }
    if(a66[5][5] < 0. ) {
      for(i=0;i<6;i++) a66[i][5]=- a66[i][5];
      //cout<<"Mode III: negative beta reverted."<<endl; 
    }
  }
  
  if(false) {
    cout<<"Intermediate A matrix :"<<endl;
    for (i=0;i<6;i++) {
      for(j=0;j<6;j++) cout<<setw(12)<<a66[i][j]<<"  ";
      cout<<endl;
    }

    for(i=0;i<6;i++) {
      for(j=0;j<6;j++) b66[i][j]=a66[i][j];
    }
    det= mat_det( &b66[0][0], 6);
    cout<<" Det of A66  =  "<< det <<endl;
  }
  
  //---normalizing A: method I
   
  double  s66[6][6];
  for(i=0;i<6;i++)
    for(j=0;j<6;j++) s66[i][j]=0.;
  s66[0][1]= 1.0;
  s66[1][0]=-1.0;
  s66[2][3]= 1.0;
  s66[3][2]=-1.0;  
  s66[4][5]= 1.0;
  s66[5][4]=-1.0;  

  if(false){
    cout<<"S matrix :"<<endl;
    for (i=0;i<6;i++) {
      for(j=0;j<6;j++) cout<<setw(12)<<s66[i][j]<<"  ";
      cout<<endl;
    }
  }

  temp1=(-a66[0][1]*a66[1][0] + a66[0][0]*a66[1][1] -a66[0][3]*a66[1][2] +a66[0][2]*a66[1][3] - a66[0][5]*a66[1][4] +a66[0][4]*a66[1][5]);
  temp2=(-a66[2][1]*a66[3][0] + a66[2][0]*a66[3][1] -a66[2][3]*a66[3][2] +a66[2][2]*a66[3][3] - a66[2][5]*a66[3][4] +a66[2][4]*a66[3][5]);  
  temp3=(-a66[4][1]*a66[5][0] + a66[4][0]*a66[5][1] -a66[4][3]*a66[5][2] +a66[4][2]*a66[5][3] - a66[4][5]*a66[5][4] +a66[4][4]*a66[5][5]);

  //cout<<"...normalize new v1 and v2, v3 adn v4, v5 and v6 :"<<endl;
  //cout<<temp1<<" "<<temp2<<" "<<temp3<<endl;
  
  for(i=0;i<6;i++){
    a66[0][i]=a66[0][i]/sqrt(temp1);
    a66[1][i]=a66[1][i]/sqrt(temp1);
    a66[2][i]=a66[2][i]/sqrt(temp2);
    a66[3][i]=a66[3][i]/sqrt(temp2);
    a66[4][i]=a66[4][i]/sqrt(temp3);
    a66[5][i]=a66[5][i]/sqrt(temp3);
  }

  if(false){
    temp1=(-a66[0][1]*a66[1][0] + a66[0][0]*a66[1][1] -a66[0][3]*a66[1][2] +a66[0][2]*a66[1][3] - a66[0][5]*a66[1][4] +a66[0][4]*a66[1][5]);
    temp2=(-a66[2][1]*a66[3][0] + a66[2][0]*a66[3][1] -a66[2][3]*a66[3][2] +a66[2][2]*a66[3][3] - a66[2][5]*a66[3][4] +a66[2][4]*a66[3][5]);  
    temp3=(-a66[4][1]*a66[5][0] + a66[4][0]*a66[5][1] -a66[4][3]*a66[5][2] +a66[4][2]*a66[5][3] - a66[4][5]*a66[5][4] +a66[4][4]*a66[5][5]);

    cout<<"...normalize new v1 and v2, v3 adn v4, v5 and v6 :"<<endl;
    cout<<temp1<<" "<<temp2<<" "<<temp3<<endl;
  }

  if(true) {
    cout<<"A3 matrix :"<<endl;
    for (i=0;i<6;i++) {
      for(j=0;j<6;j++) cout<<setw(12)<<a66[j][i]<<"  ";
      cout<<endl;
    }
    
    for(i=0;i<6;i++) {
      for(j=0;j<6;j++) b66[i][j]=a66[j][i];
    }
    det= mat_det( &b66[0][0], 6);
    cout<<" Det of A3  =  "<< det <<endl;
  }
  
  //----return
  for(i=0;i<6;i++)
    for(j=0;j<6;j++)
      a[i*6+j]=  a66[j][i]; 
  
}

//-----calculate 4-d  chrom, dispersion, etc.

void Cal_Chrom( Line & linename)
{
  double deltap;
  double qx0, qy0, qxp, qyp, qxm, qym;
  
  if( GP.H_expand == true ){
    deltap = GP.step_deltap;
  }
  else{
    deltap = DeltaToPt( GP.step_deltap );
  }
    
  Cal_Twiss(linename, deltap);
  qxp=linename.Tune1;
  qyp=linename.Tune2;
  
  Cal_Twiss(linename, -deltap);
  qxm=linename.Tune1;
  qym=linename.Tune2;

  Cal_Twiss(linename, 0.);
  qx0=linename.Tune1;
  qy0=linename.Tune2;

  linename.Chromx1=(qxp-qxm)/2./GP.step_deltap;
  linename.Chromy1=(qyp-qym)/2./GP.step_deltap;
  linename.Chromx2=(qxp+qxm-2*qx0)/2./GP.step_deltap/GP.step_deltap;
  linename.Chromy2=(qyp+qym-2*qy0)/2./GP.step_deltap/GP.step_deltap;
}

void Cal_Chrom_Dispersion_Num( Line & linename)
{
  int i,j;
  double deltap;
  double qx0, qy0, qxp, qyp, qxm, qym;
  int    nelement = linename.Ncell;
  double xp[4*nelement], xm[4*nelement];
  
  if( GP.H_expand == true ){
    deltap = GP.step_deltap;
  }
  else{
    deltap = DeltaToPt( GP.step_deltap );
  }

  //cout<<"  I am  here. "<<endl;
  
  Cal_Twiss(linename, deltap);
  qxp=linename.Tune1;
  qyp=linename.Tune2;
  for(i=0;i<linename.Ncell;i++) {
    for(j=0;j<4;j++) {
      xp[i*4+j] = linename.Cell[i]->X[j];
    }
  }

  //cout<<"  I am  here. "<<endl;
				  
  Cal_Twiss(linename, -deltap);
  qxm=linename.Tune1;
  qym=linename.Tune2;
  for(i=0;i<linename.Ncell;i++) {
    for(j=0;j<4;j++) {
      xm[i*4+j] = linename.Cell[i]->X[j];
    }
  }
  
  Cal_Twiss(linename, 0.);
  qx0=linename.Tune1;
  qy0=linename.Tune2;

  linename.Chromx1=(qxp-qxm)/2./GP.step_deltap;
  linename.Chromy1=(qyp-qym)/2./GP.step_deltap;
  linename.Chromx2=(qxp+qxm-2*qx0)/2./GP.step_deltap/GP.step_deltap;
  linename.Chromy2=(qyp+qym-2*qy0)/2./GP.step_deltap/GP.step_deltap;
  
  for(i=0;i<linename.Ncell;i++) {
    linename.Cell[i]->Etax  = ( xp[i*4+0 ]  -  xm[i*4+0]  ) / 2 / GP.step_deltap;
    linename.Cell[i]->Etaxp = ( xp[i*4+1 ]  -  xm[i*4+1]  ) / 2 / GP.step_deltap;
    linename.Cell[i]->Etay  = ( xp[i*4+2 ]  -  xm[i*4+2]  ) / 2 / GP.step_deltap;
    linename.Cell[i]->Etayp = ( xp[i*4+3 ]  -  xm[i*4+3]  ) / 2 / GP.step_deltap;
  }

}

void Cal_Dispersion(Line & linename)
{
  int i,j, k;
  double m66[6][6];
  double r44[4][4];
  double x[6], xtemp[6];  
  double T66[6][6];
  
  for(i=0;i<6;i++) 
    for(j=0;j<6;j++) m66[i][j]= linename.Cell[linename.Ncell-1]->M[i*6+j];
  for(i=0;i<4;i++) 
    for(j=0;j<4;j++) r44[i][j]= linename.Cell[linename.Ncell-1]->M[i*6+j];
  for(i=0;i<4;i++) r44[i][i]  = r44[i][i]-1.0;
  mat_inv(&r44[0][0],4); 
  
  for(i=0;i<4;i++) 
    for(j=0;j<4;j++) r44[i][j]= -1.0 * r44[i][j];
  linename.Cell[linename.Ncell-1]->Etax  = r44[0][0] * m66[0][5] +  r44[0][1] * m66[1][5] +   r44[0][2] * m66[2][5] +  r44[0][3] * m66[3][5];
  linename.Cell[linename.Ncell-1]->Etaxp = r44[1][0] * m66[0][5] +  r44[1][1] * m66[1][5] +   r44[1][2] * m66[2][5] +  r44[1][3] * m66[3][5];
  linename.Cell[linename.Ncell-1]->Etay  = r44[2][0] * m66[0][5] +  r44[2][1] * m66[1][5] +   r44[2][2] * m66[2][5] +  r44[2][3] * m66[3][5];
  linename.Cell[linename.Ncell-1]->Etayp = r44[3][0] * m66[0][5] +  r44[3][1] * m66[1][5] +   r44[3][2] * m66[2][5] +  r44[3][3] * m66[3][5];
  
  x[0]= linename.Cell[linename.Ncell-1]->Etax  ;
  x[1]= linename.Cell[linename.Ncell-1]->Etaxp ;
  x[2]= linename.Cell[linename.Ncell-1]->Etay  ;
  x[3]= linename.Cell[linename.Ncell-1]->Etayp ;
  for(k=0;k<linename.Ncell; k++){
    x[4]= 0.;  x[5]= 1;
    for(i=0;i<6;i++) 
      for(j=0;j<6;j++) T66[i][j]= linename.Cell[k]->T[i*6+j];
    xtemp[0]=  x[0] * T66[0][0] +  x[1] * T66[0][1] +   x[2] * T66[0][2] +  x[3] * T66[0][3]  +   x[5] * T66[0][5]   ;
    xtemp[1]=  x[0] * T66[1][0] +  x[1] * T66[1][1] +   x[2] * T66[1][2] +  x[3] * T66[1][3]  +   x[5] * T66[1][5]   ;
    xtemp[2]=  x[0] * T66[2][0] +  x[1] * T66[2][1] +   x[2] * T66[2][2] +  x[3] * T66[2][3]  +   x[5] * T66[2][5]   ;
    xtemp[3]=  x[0] * T66[3][0] +  x[1] * T66[3][1] +   x[2] * T66[3][2] +  x[3] * T66[3][3]  +   x[5] * T66[3][5]   ;
    linename.Cell[k]->Etax  =  xtemp[0];
    linename.Cell[k]->Etaxp =  xtemp[1];
    linename.Cell[k]->Etay  =  xtemp[2];
    linename.Cell[k]->Etayp =  xtemp[3];
    for(i=0;i<4;i++) x[i]= xtemp[i];
  }

  if(GP.H_expand==false){
    for(k=0;k<linename.Ncell;k++) {
      linename.Cell[k]->Etax  = linename.Cell[k]->Etax  * GP.beta  ;
      linename.Cell[k]->Etaxp = linename.Cell[k]->Etaxp * GP.beta  ;
      linename.Cell[k]->Etay  = linename.Cell[k]->Etay  * GP.beta  ;
      linename.Cell[k]->Etayp = linename.Cell[k]->Etayp * GP.beta  ;
    }
  }
}

void Cal_Momentum_Dispersion(Line & linename)
{
  int i,j, k;
  double m66[6][6];
  double r44[4][4];
  double x[6], xtemp[6];  
  double T66[6][6];

  for(i=0;i<6;i++) 
    for(j=0;j<6;j++) m66[i][j]= linename.Cell[linename.Ncell-1]->M[i*6+j];
  for(i=0;i<4;i++) 
    for(j=0;j<4;j++) r44[i][j]= linename.Cell[linename.Ncell-1]->M[i*6+j];
  for(i=0;i<4;i++) r44[i][i]  = r44[i][i]-1.0;
  mat_inv(&r44[0][0],4); 
  
  for(i=0;i<4;i++) 
    for(j=0;j<4;j++) r44[i][j]= -1.0 * r44[i][j];
  linename.Cell[linename.Ncell-1]->Etax  = r44[0][0] * m66[0][5] +  r44[0][1] * m66[1][5] +   r44[0][2] * m66[2][5] +  r44[0][3] * m66[3][5];
  linename.Cell[linename.Ncell-1]->Etaxp = r44[1][0] * m66[0][5] +  r44[1][1] * m66[1][5] +   r44[1][2] * m66[2][5] +  r44[1][3] * m66[3][5];
  linename.Cell[linename.Ncell-1]->Etay  = r44[2][0] * m66[0][5] +  r44[2][1] * m66[1][5] +   r44[2][2] * m66[2][5] +  r44[2][3] * m66[3][5];
  linename.Cell[linename.Ncell-1]->Etayp = r44[3][0] * m66[0][5] +  r44[3][1] * m66[1][5] +   r44[3][2] * m66[2][5] +  r44[3][3] * m66[3][5];
  
  x[0]= linename.Cell[linename.Ncell-1]->Etax  ;
  x[1]= linename.Cell[linename.Ncell-1]->Etaxp ;
  x[2]= linename.Cell[linename.Ncell-1]->Etay  ;
  x[3]= linename.Cell[linename.Ncell-1]->Etayp ;
  for(k=0;k<linename.Ncell; k++){
    x[4]= 0.;  x[5]= 1;
    for(i=0;i<6;i++) 
      for(j=0;j<6;j++) T66[i][j]= linename.Cell[k]->T[i*6+j];
    xtemp[0]=  x[0] * T66[0][0] +  x[1] * T66[0][1] +   x[2] * T66[0][2] +  x[3] * T66[0][3]  +   x[5] * T66[0][5]   ;
    xtemp[1]=  x[0] * T66[1][0] +  x[1] * T66[1][1] +   x[2] * T66[1][2] +  x[3] * T66[1][3]  +   x[5] * T66[1][5]   ;
    xtemp[2]=  x[0] * T66[2][0] +  x[1] * T66[2][1] +   x[2] * T66[2][2] +  x[3] * T66[2][3]  +   x[5] * T66[2][5]   ;
    xtemp[3]=  x[0] * T66[3][0] +  x[1] * T66[3][1] +   x[2] * T66[3][2] +  x[3] * T66[3][3]  +   x[5] * T66[3][5]   ;
    linename.Cell[k]->Etax  =  xtemp[0];
    linename.Cell[k]->Etaxp =  xtemp[1];
    linename.Cell[k]->Etay  =  xtemp[2];
    linename.Cell[k]->Etayp =  xtemp[3];
    for(i=0;i<4;i++) x[i]= xtemp[i];
  }

  if(GP.H_expand==false){
    for(k=0;k<linename.Ncell;k++) {
      linename.Cell[k]->Etax  = linename.Cell[k]->Etax  * GP.beta  ;
      linename.Cell[k]->Etaxp = linename.Cell[k]->Etaxp * GP.beta  ;
      linename.Cell[k]->Etay  = linename.Cell[k]->Etay  * GP.beta  ;
      linename.Cell[k]->Etayp = linename.Cell[k]->Etayp * GP.beta  ;
    }
  }

}

void Trace_Momentum_Dispersion(Line & linename, int i1, int i2, double Eta1[4])
{
  int i,j, k;
  double m66[6][6];
  double r44[4][4];
  double x[6], xtemp[6];  
  double T66[6][6];
  
  x[0]= Eta1[0] ;
  x[1]= Eta1[1] ;
  x[2]= Eta1[2] ;
  x[3]= Eta1[3] ;
  for(k=i1+1;k<=i2; k++){
    x[4]= 0.;  x[5]= 1;
    for(i=0;i<6;i++) 
      for(j=0;j<6;j++) T66[i][j]= linename.Cell[k]->T[i*6+j];
    xtemp[0]=  x[0] * T66[0][0] +  x[1] * T66[0][1] +   x[2] * T66[0][2] +  x[3] * T66[0][3]  +   x[5] * T66[0][5]   ;
    xtemp[1]=  x[0] * T66[1][0] +  x[1] * T66[1][1] +   x[2] * T66[1][2] +  x[3] * T66[1][3]  +   x[5] * T66[1][5]   ;
    xtemp[2]=  x[0] * T66[2][0] +  x[1] * T66[2][1] +   x[2] * T66[2][2] +  x[3] * T66[2][3]  +   x[5] * T66[2][5]   ;
    xtemp[3]=  x[0] * T66[3][0] +  x[1] * T66[3][1] +   x[2] * T66[3][2] +  x[3] * T66[3][3]  +   x[5] * T66[3][5]   ;
    linename.Cell[k]->Etax  =  xtemp[0];
    linename.Cell[k]->Etaxp =  xtemp[1];
    linename.Cell[k]->Etay  =  xtemp[2];
    linename.Cell[k]->Etayp =  xtemp[3];
    for(i=0;i<4;i++) x[i]= xtemp[i];
  }

}

void Cal_Timeflight_Dispersion(Line & linename)
{
  int i,j, k;
  double m66[6][6];
  double r44[4][4];
  double x[6], xtemp[6];  
  double T66[6][6];

  for(i=0;i<6;i++) 
    for(j=0;j<6;j++) m66[i][j]= linename.Cell[linename.Ncell-1]->M[i*6+j];
  for(i=0;i<4;i++) 
    for(j=0;j<4;j++) r44[i][j]= linename.Cell[linename.Ncell-1]->M[i*6+j];
  for(i=0;i<4;i++) r44[i][i]  = r44[i][i]-1.0;
  mat_inv(&r44[0][0],4); 
  
  for(i=0;i<4;i++) 
    for(j=0;j<4;j++) r44[i][j]= -1.0 * r44[i][j];
  linename.Cell[linename.Ncell-1]->Ksix  = r44[0][0] * m66[0][4] +  r44[0][1] * m66[1][4] +   r44[0][2] * m66[2][4] +  r44[0][3] * m66[3][4];
  linename.Cell[linename.Ncell-1]->Ksixp = r44[1][0] * m66[0][4] +  r44[1][1] * m66[1][4] +   r44[1][2] * m66[2][4] +  r44[1][3] * m66[3][4];
  linename.Cell[linename.Ncell-1]->Ksiy  = r44[2][0] * m66[0][4] +  r44[2][1] * m66[1][4] +   r44[2][2] * m66[2][4] +  r44[2][3] * m66[3][4];
  linename.Cell[linename.Ncell-1]->Ksiyp = r44[3][0] * m66[0][4] +  r44[3][1] * m66[1][4] +   r44[3][2] * m66[2][4] +  r44[3][3] * m66[3][4];
  
  x[0]= linename.Cell[linename.Ncell-1]->Ksix  ;
  x[1]= linename.Cell[linename.Ncell-1]->Ksixp ;
  x[2]= linename.Cell[linename.Ncell-1]->Ksiy  ;
  x[3]= linename.Cell[linename.Ncell-1]->Ksiyp ;

  for(k=0;k<linename.Ncell; k++){
    x[4]= 1.;  x[5]= 0;
    for(i=0;i<6;i++) 
      for(j=0;j<6;j++) T66[i][j]= linename.Cell[k]->T[i*6+j];
    xtemp[0]=  x[0] * T66[0][0] +  x[1] * T66[0][1] +   x[2] * T66[0][2] +  x[3] * T66[0][3]  +   x[4] * T66[0][4]   ;
    xtemp[1]=  x[0] * T66[1][0] +  x[1] * T66[1][1] +   x[2] * T66[1][2] +  x[3] * T66[1][3]  +   x[4] * T66[1][4]   ;
    xtemp[2]=  x[0] * T66[2][0] +  x[1] * T66[2][1] +   x[2] * T66[2][2] +  x[3] * T66[2][3]  +   x[4] * T66[2][4]   ;
    xtemp[3]=  x[0] * T66[3][0] +  x[1] * T66[3][1] +   x[2] * T66[3][2] +  x[3] * T66[3][3]  +   x[4] * T66[3][4]   ;
    linename.Cell[k]->Ksix  =  xtemp[0];
    linename.Cell[k]->Ksixp =  xtemp[1];
    linename.Cell[k]->Ksiy  =  xtemp[2];
    linename.Cell[k]->Ksiyp =  xtemp[3];
    for(i=0;i<4;i++) x[i]= xtemp[i];
  }

  if(GP.H_expand==false){
    for(k=0;k<linename.Ncell;k++) {
      linename.Cell[k]->Etax  = linename.Cell[k]->Etax  * GP.beta  ;
      linename.Cell[k]->Etaxp = linename.Cell[k]->Etaxp * GP.beta  ;
      linename.Cell[k]->Etay  = linename.Cell[k]->Etay  * GP.beta  ;
      linename.Cell[k]->Etayp = linename.Cell[k]->Etayp * GP.beta  ;
    }
  }

}  

void Trace_Timeflight_Dispersion(Line & linename, int i1,  int  i2, double Ksi1[4])
{
  int i,j, k;
  double m66[6][6];
  double r44[4][4];
  double x[6], xtemp[6];  
  double T66[6][6];

  x[0]= Ksi1[0]  ;
  x[1]= Ksi1[1] ;
  x[2]= Ksi1[2]  ;
  x[3]= Ksi1[3] ;

  for(k=i1+1;k<=i2; k++){
    x[4]= 1.;  x[5]= 0;
    for(i=0;i<6;i++) 
      for(j=0;j<6;j++) T66[i][j]= linename.Cell[k]->T[i*6+j];
    xtemp[0]=  x[0] * T66[0][0] +  x[1] * T66[0][1] +   x[2] * T66[0][2] +  x[3] * T66[0][3]  +   x[4] * T66[0][4]   ;
    xtemp[1]=  x[0] * T66[1][0] +  x[1] * T66[1][1] +   x[2] * T66[1][2] +  x[3] * T66[1][3]  +   x[4] * T66[1][4]   ;
    xtemp[2]=  x[0] * T66[2][0] +  x[1] * T66[2][1] +   x[2] * T66[2][2] +  x[3] * T66[2][3]  +   x[4] * T66[2][4]   ;
    xtemp[3]=  x[0] * T66[3][0] +  x[1] * T66[3][1] +   x[2] * T66[3][2] +  x[3] * T66[3][3]  +   x[4] * T66[3][4]   ;
    linename.Cell[k]->Ksix  =  xtemp[0];
    linename.Cell[k]->Ksixp =  xtemp[1];
    linename.Cell[k]->Ksiy  =  xtemp[2];
    linename.Cell[k]->Ksiyp =  xtemp[3];
    for(i=0;i<4;i++) x[i]= xtemp[i];
  }

}  



void Cal_Optics(Line & linename)  // will give chrom, Twiss(dp/p0=0), and dispersions
{
  Cal_Chrom(linename);
  Cal_Dispersion(linename);
}

//-----  fitting tunes and linear chroms

void Fit_Tune(Line & linename, double q1, double q2, const char * qf_name, const char * qd_name)
{
  double tunex0, tuney0, tunex1=q1, tuney1=q2, dtunex, dtuney;
  double qf_k1l_0, qd_k1l_0;
  double dk1l_qf, dtunex_qf,  dtuney_qf, dk1l_qd, dtunex_qd,  dtuney_qd;
  double scale_qf, scale_qd;

  Cal_Twiss(linename,0.0);
  tunex0=linename.Tune1;
  tuney0=linename.Tune2;

  while( (tunex1-tunex0)*(tunex1-tunex0) + (tuney1-tuney0)*(tuney1-tuney0)  > 1.0e-10 ) {
    qf_k1l_0= Get_KL(linename,qf_name, "K1L");
    qd_k1l_0= Get_KL(linename, qd_name,"K1L"); 
    
    dk1l_qf=  qf_k1l_0 * 0.001;
    Set_dKL(linename,qf_name, "K1L", dk1l_qf);
    Cal_Twiss(linename,0.0);
    dtunex_qf=linename.Tune1 - tunex0;
    dtuney_qf=linename.Tune2 - tuney0;
    Set_dKL(linename,qf_name, "K1L",-dk1l_qf);
    
    dk1l_qd=  qd_k1l_0 * 0.001;
    Set_dKL(linename,qd_name, "K1L", dk1l_qd);
    Cal_Twiss(linename,0.0);
    dtunex_qd=linename.Tune1 - tunex0;
    dtuney_qd=linename.Tune2 - tuney0;
    Set_dKL(linename,qd_name, "K1L",-dk1l_qd);
    
    dtunex=tunex1- tunex0;
    dtuney=tuney1- tuney0;
    
    LinearEquations(dtunex_qf, dtunex_qd,dtunex, dtuney_qf, dtuney_qd, dtuney, scale_qf, scale_qd);
    Set_dKL(linename,qf_name, "K1L", dk1l_qf * scale_qf);
    Set_dKL(linename,qd_name, "K1L", dk1l_qd * scale_qd);
    
    Cal_Twiss(linename,0.0);
    tunex0=linename.Tune1;
    tuney0=linename.Tune2;
  }
}

void Fit_Tune_RHICelens(Line & linename, double q1, double q2)
{
  double tunex0, tuney0, tunex1=q1, tuney1=q2, dtunex, dtuney;
  double qf_k1l_0, qd_k1l_0;
  double dk1l_qf, dtunex_qf,  dtuney_qf, dk1l_qd, dtunex_qd,  dtuney_qd;
  double scale_qf, scale_qd;

  Cal_Twiss(linename,0.0);
  tunex0=linename.Tune1;
  tuney0=linename.Tune2;

  while( (tunex1-tunex0)*(tunex1-tunex0) + (tuney1-tuney0)*(tuney1-tuney0)  > 1.0e-10 ) {
    qf_k1l_0= Get_KL(linename,"QF", "K1L");
    qd_k1l_0= Get_KL(linename, "QD","K1L"); 
    
    dk1l_qf=  qf_k1l_0 * 0.001;
    Set_dKL(linename,"QF", "K1L",  dk1l_qf);
    Set_dKL(linename,"QF9", "K1L", dk1l_qf);

    Cal_Twiss(linename,0.0);
    dtunex_qf=linename.Tune1 - tunex0;
    dtuney_qf=linename.Tune2 - tuney0;
    Set_dKL(linename,"QF", "K1L",  -dk1l_qf);
    Set_dKL(linename,"QF9", "K1L", -dk1l_qf);
    
    dk1l_qd=  qd_k1l_0 * 0.001;
    Set_dKL(linename,"QD", "K1L",  dk1l_qd);
    Set_dKL(linename,"QD9", "K1L", dk1l_qd);

    Cal_Twiss(linename,0.0);
    dtunex_qd=linename.Tune1 - tunex0;
    dtuney_qd=linename.Tune2 - tuney0;
    Set_dKL(linename,"QD", "K1L",  -dk1l_qd);
    Set_dKL(linename,"QD9", "K1L", -dk1l_qd);

    dtunex=tunex1- tunex0;
    dtuney=tuney1- tuney0;
    
    LinearEquations(dtunex_qf, dtunex_qd,dtunex, dtuney_qf, dtuney_qd, dtuney, scale_qf, scale_qd);
    Set_dKL(linename,"QF", "K1L",   dk1l_qf * scale_qf);
    Set_dKL(linename,"QF9", "K1L",  dk1l_qf * scale_qf);

    Set_dKL(linename,"QD", "K1L",  dk1l_qd * scale_qd);
    Set_dKL(linename,"QD9", "K1L", dk1l_qd * scale_qd);
  
    Cal_Twiss(linename,0.0);
    tunex0=linename.Tune1;
    tuney0=linename.Tune2;
  }
}


double  Get_RHIC_QF_K1L(Line & linename)
{
  int i;
  double temp;
  
  for(i=0;i<linename.Ncell;i++){
    if( linename.Cell[i]->TYPE==string("QUAD")  and linename.Cell[i]->L == 1.11 ){
      if ( linename.Cell[i]->GetP("K1L") > 0.) {
	temp= linename.Cell[i]->GetP("K1L");
	break;
      }
    }
  }

  return temp;
}


double  Get_RHIC_QD_K1L(Line & linename)
{
  int i;
  double temp;
  
  for(i=0;i<linename.Ncell;i++){
    if( linename.Cell[i]->TYPE==string("QUAD")  and linename.Cell[i]->L == 1.11 ){
      if ( linename.Cell[i]->GetP("K1L") < 0.) {
	temp= linename.Cell[i]->GetP("K1L");
	break;
      }
    }
  }

  return temp;
}


void  Set_RHIC_QF_dK1L(Line & linename, double dk1l)
{
  int i;
  double temp;
  
  for(i=0;i<linename.Ncell;i++){
    if( linename.Cell[i]->TYPE==string("QUAD")  and linename.Cell[i]->L == 1.11 ){
      if ( linename.Cell[i]->GetP("K1L") > 0.) {
	linename.Cell[i]->SetP( "K1L", linename.Cell[i]->GetP("K1L") + dk1l );
      }
    }
  }

}

void  Set_RHIC_QD_dK1L(Line & linename, double dk1l)
{
  int i;
  double temp;
  
  for(i=0;i<linename.Ncell;i++){
    if( linename.Cell[i]->TYPE==string("QUAD")  and linename.Cell[i]->L == 1.11 ){
      if ( linename.Cell[i]->GetP("K1L") < 0.) {
	linename.Cell[i]->SetP( "K1L", linename.Cell[i]->GetP("K1L") + dk1l );
      }
    }
  }

}


void Fit_Tune_RHIC(Line & linename, double q1, double q2 )
{
  double tunex0, tuney0, tunex1=q1, tuney1=q2, dtunex, dtuney;
  double qf_k1l_0, qd_k1l_0;
  double dk1l_qf, dtunex_qf,  dtuney_qf, dk1l_qd, dtunex_qd,  dtuney_qd;
  double scale_qf, scale_qd;

  Cal_Twiss(linename,0.0);
  tunex0=linename.Tune1;
  tuney0=linename.Tune2;

  while( (tunex1-tunex0)*(tunex1-tunex0) + (tuney1-tuney0)*(tuney1-tuney0)  > 1.0e-10 ) {
    qf_k1l_0= Get_RHIC_QF_K1L(linename);
    qd_k1l_0= Get_RHIC_QD_K1L(linename);
    
    dk1l_qf=  qf_k1l_0 * 0.001;
    Set_RHIC_QF_dK1L(linename,  dk1l_qf);
    Cal_Twiss(linename,0.0);
    dtunex_qf=linename.Tune1 - tunex0;
    dtuney_qf=linename.Tune2 - tuney0;
    Set_RHIC_QF_dK1L(linename,  -dk1l_qf);
    
    dk1l_qd=  qd_k1l_0 * 0.001;
    Set_RHIC_QD_dK1L(linename,  dk1l_qd);
    Cal_Twiss(linename,0.0);
    dtunex_qd=linename.Tune1 - tunex0;
    dtuney_qd=linename.Tune2 - tuney0;
    Set_RHIC_QD_dK1L(linename,  -dk1l_qd);
    
    dtunex=tunex1- tunex0;
    dtuney=tuney1- tuney0;
    
    LinearEquations(dtunex_qf, dtunex_qd,dtunex, dtuney_qf, dtuney_qd, dtuney, scale_qf, scale_qd);
    Set_RHIC_QF_dK1L(linename,  dk1l_qf * scale_qf);
    Set_RHIC_QD_dK1L(linename,  dk1l_qd * scale_qd);
    
    Cal_Twiss(linename,0.0);
    tunex0=linename.Tune1;
    tuney0=linename.Tune2;
  }
  
}

void Fit_Chrom(Line & linename, double chrom1x_want, double chrom1y_want, const char * sf_name, const char * sd_name )
{
  double chrom1x0, chrom1y0,  dchrom1x, dchrom1y;
  double dk2l_sf, dk2l_sd, dchrom1x_sf,  dchrom1y_sf, dchrom1x_sd, dchrom1y_sd;
  double scale_sf, scale_sd;

  Cal_Chrom(linename) ; 
  chrom1x0= linename.Chromx1;
  chrom1y0= linename.Chromy1;

  while( (chrom1x_want-chrom1x0)*(chrom1x_want-chrom1x0) + (chrom1y_want-chrom1y0)*(chrom1y_want-chrom1y0)  > 0.0001 ) {
    //sf_k2l_0= Get_KL(linename,sf_name,"K2L");
    //sd_k2l_0= Get_KL(linename,sd_name,"K2L"); 
    
    dk2l_sf= 0.3 * 0.005;
    Set_dKL(linename,sf_name, "K2L", dk2l_sf);
    Cal_Chrom(linename);
    dchrom1x_sf=linename.Chromx1 - chrom1x0;
    dchrom1y_sf=linename.Chromy1 - chrom1y0;
    Set_dKL(linename,sf_name, "K2L",-dk2l_sf);
    
    dk2l_sd= -0.5 * 0.005;
    Set_dKL(linename,sd_name, "K2L", dk2l_sd);
    Cal_Chrom(linename);
    dchrom1x_sd=linename.Chromx1 - chrom1x0;
    dchrom1y_sd=linename.Chromy1 - chrom1y0;
    Set_dKL(linename,sd_name, "K2L",-dk2l_sd);

    dchrom1x=chrom1x_want- chrom1x0;
    dchrom1y=chrom1y_want- chrom1y0;
    
    LinearEquations(dchrom1x_sf, dchrom1x_sd, dchrom1x, dchrom1y_sf, dchrom1y_sd, dchrom1y, scale_sf, scale_sd);
    Set_dKL(linename,sf_name, "K2L", dk2l_sf * scale_sf);
    Set_dKL(linename,sd_name, "K2L", dk2l_sd * scale_sd);
    
    Cal_Chrom(linename);
    chrom1x0=linename.Chromx1;
    chrom1y0=linename.Chromy1;
  }
}


void Fit_Chrom_RHIC8fam(Line & linename, double chrom1x_want, double chrom1y_want )
{
  int i;
  double chrom1x0, chrom1y0,  dchrom1x, dchrom1y;
  double dk2l_sf, dk2l_sd, dchrom1x_sf,  dchrom1y_sf, dchrom1x_sd, dchrom1y_sd;
  double scale_sf, scale_sd;

  Cal_Chrom(linename);
  chrom1x0=linename.Chromx1;
  chrom1y0=linename.Chromy1;

  while( (chrom1x_want-chrom1x0)*(chrom1x_want-chrom1x0) + (chrom1y_want-chrom1y0)*(chrom1y_want-chrom1y0)  > 0.0001 ) {
    
    dk2l_sf =  0.3 * 0.005;
    for(i=0;i<linename.Ncell;i++) 
      if(linename.Cell[i]->NAME.find("SF") != string::npos and linename.Cell[i]->TYPE ==string("SEXT"))
	linename.Cell[i]->SetP("K2L",  linename.Cell[i]->GetP("K2L") +  dk2l_sf );
    Cal_Chrom(linename);
    dchrom1x_sf=linename.Chromx1 - chrom1x0;
    dchrom1y_sf=linename.Chromy1 - chrom1y0;  
    for(i=0;i<linename.Ncell;i++) 
      if(linename.Cell[i]->NAME.find("SF") != string::npos and linename.Cell[i]->TYPE ==string("SEXT")) 
	linename.Cell[i]->SetP( "K2L",linename.Cell[i]->GetP("K2L") -  dk2l_sf );

    dk2l_sd =  -0.5 * 0.005;
    for(i=0;i<linename.Ncell;i++) 
      if(linename.Cell[i]->NAME.find("SD") != string::npos and linename.Cell[i]->TYPE ==string("SEXT")) 
	linename.Cell[i]->SetP( "K2L", linename.Cell[i]->GetP("K2L") +  dk2l_sd );
    Cal_Chrom(linename);
    dchrom1x_sd=linename.Chromx1 - chrom1x0;
    dchrom1y_sd=linename.Chromy1 - chrom1y0;
    for(i=0;i<linename.Ncell;i++) 
      if(linename.Cell[i]->NAME.find("SD") != string::npos and linename.Cell[i]->TYPE ==string("SEXT")) 
	linename.Cell[i]->SetP( "K2L", linename.Cell[i]->GetP("K2L") -  dk2l_sd );

    dchrom1x=chrom1x_want- chrom1x0;
    dchrom1y=chrom1y_want- chrom1y0;
    LinearEquations(dchrom1x_sf, dchrom1x_sd, dchrom1x, dchrom1y_sf, dchrom1y_sd, dchrom1y, scale_sf, scale_sd);
    for(i=0;i<linename.Ncell;i++) 
      if(linename.Cell[i]->NAME.find("SF") != string::npos and linename.Cell[i]->TYPE ==string("SEXT"))
	linename.Cell[i]->SetP( "K2L", linename.Cell[i]->GetP("K2L") + dk2l_sf * scale_sf );
    for(i=0;i<linename.Ncell;i++) 
      if(linename.Cell[i]->NAME.find("SD") != string::npos and linename.Cell[i]->TYPE ==string("SEXT")) 
	linename.Cell[i]->SetP( "K2L", linename.Cell[i]->GetP("K2L") +dk2l_sd * scale_sd );
    
    Cal_Chrom(linename);
    chrom1x0=linename.Chromx1;
    chrom1y0=linename.Chromy1;
  }
}

void Fit_Chrom_RHIC(Line & linename, double chrom1x_want, double chrom1y_want )
{
  int i;
  double chrom1x0, chrom1y0,  dchrom1x, dchrom1y;
  double dk2l_sf, dk2l_sd, dchrom1x_sf,  dchrom1y_sf, dchrom1x_sd, dchrom1y_sd;
  double scale_sf, scale_sd;

  Cal_Chrom(linename);
  chrom1x0=linename.Chromx1;
  chrom1y0=linename.Chromy1;

  while( (chrom1x_want-chrom1x0)*(chrom1x_want-chrom1x0) + (chrom1y_want-chrom1y0)*(chrom1y_want-chrom1y0)  > 0.0001 ) {
    
    dk2l_sf =  0.3 * 0.005;
    for(i=0;i<linename.Ncell;i++) 
      if(linename.Cell[i]->NAME.find("SXF") != string::npos and linename.Cell[i]->TYPE ==string("SEXT"))
	linename.Cell[i]->SetP("K2L",  linename.Cell[i]->GetP("K2L") +  dk2l_sf );
    Cal_Chrom(linename);
    dchrom1x_sf=linename.Chromx1 - chrom1x0;
    dchrom1y_sf=linename.Chromy1 - chrom1y0;  
    for(i=0;i<linename.Ncell;i++) 
      if(linename.Cell[i]->NAME.find("SXF") != string::npos and linename.Cell[i]->TYPE ==string("SEXT")) 
	linename.Cell[i]->SetP( "K2L",linename.Cell[i]->GetP("K2L") -  dk2l_sf );

    dk2l_sd =  -0.5 * 0.005;
    for(i=0;i<linename.Ncell;i++) 
      if(linename.Cell[i]->NAME.find("SXD") != string::npos and linename.Cell[i]->TYPE ==string("SEXT")) 
	linename.Cell[i]->SetP( "K2L", linename.Cell[i]->GetP("K2L") +  dk2l_sd );
    Cal_Chrom(linename);
    dchrom1x_sd=linename.Chromx1 - chrom1x0;
    dchrom1y_sd=linename.Chromy1 - chrom1y0;
    for(i=0;i<linename.Ncell;i++) 
      if(linename.Cell[i]->NAME.find("SXD") != string::npos and linename.Cell[i]->TYPE ==string("SEXT")) 
	linename.Cell[i]->SetP( "K2L", linename.Cell[i]->GetP("K2L") -  dk2l_sd );

    dchrom1x=chrom1x_want- chrom1x0;
    dchrom1y=chrom1y_want- chrom1y0;
    LinearEquations(dchrom1x_sf, dchrom1x_sd, dchrom1x, dchrom1y_sf, dchrom1y_sd, dchrom1y, scale_sf, scale_sd);
    for(i=0;i<linename.Ncell;i++) 
      if(linename.Cell[i]->NAME.find("SXF") != string::npos and linename.Cell[i]->TYPE ==string("SEXT"))
	linename.Cell[i]->SetP( "K2L", linename.Cell[i]->GetP("K2L") + dk2l_sf * scale_sf );
    for(i=0;i<linename.Ncell;i++) 
      if(linename.Cell[i]->NAME.find("SXD") != string::npos and linename.Cell[i]->TYPE ==string("SEXT")) 
	linename.Cell[i]->SetP( "K2L", linename.Cell[i]->GetP("K2L") +dk2l_sd * scale_sd );
    
    Cal_Chrom(linename);
    chrom1x0=linename.Chromx1;
    chrom1y0=linename.Chromy1;
  }
}

void  chrom_fit(double qx[],double qy[],double & chromx1,double & chromy1,double & chromx2,double & chromy2,double & chromx3,double & chromy3 )
{
  int i;
  double temp;
  double xa[21],ya[21];
  double coeff[8];

  temp=qx[10];
  for (i=0; i<21;i++){
    xa[i]=(i-10)*GP.step_deltap;
    ya[i]=qx[i]-temp;
  }
  pfit(xa, ya, 21, coeff, 7);
  chromx1=coeff[1];
  chromx2=coeff[2];
  chromx3=coeff[3];  
  
  temp=qy[10];
  for (i=0; i<21;i++){
    xa[i]=(i-10)*GP.step_deltap;
    ya[i]=qy[i]-temp;
  }
  pfit(xa, ya, 21, coeff,7);
  chromy1=coeff[1];
  chromy2=coeff[2];
  chromy3=coeff[3]; 
}

void Cal_Chrom_Num( Line & linename)
{
  int i;
  double qx[21],qy[21];
  double deltap;
  double chromx1,chromy1,chromx2,chromy2,chromx3,chromy3;  

  for(i=0;i<21;i++){
    deltap=GP.step_deltap*(i-10);
    Cal_Twiss(linename, deltap);
    qx[i]=linename.Tune1;
    qy[i]=linename.Tune2;
  }

  chrom_fit(qx,qy,chromx1,chromy1,chromx2,chromy2,chromx3,chromy3);
  
  linename.Chromx1 =  chromx1;
  linename.Chromy1 =  chromy1;
  linename.Chromx2 =  chromx2;
  linename.Chromy2 =  chromy2;
  linename.Chromx3 =  chromx3;
  linename.Chromy3 =  chromy3;

  Cal_Twiss(linename, 0.0);
}

void Correct_Chrom_Manual( Line & linename)
{
  int i;
  char fam1[16], fam2[16];
  double step;
  
  for(i=0;i<8;i++) {
    cout<<"Input two Sextupole families and change step:"<<endl;
    cin>>fam1>>fam2>>step;
    Set_dKL(linename,fam1,"K2L",step);
    Set_dKL(linename,fam2,"K2L",-step);
    Fit_Chrom_RHIC8fam(linename, 1.0, 1.0);
    Cal_Chrom(linename);  //  Cal_Chrom_Num(linename);
    cout<<linename.Chromx1<<"  "<<linename.Chromy1<<endl;
    cout<<linename.Chromx2<<"  "<<linename.Chromy2<<endl;
    //cout<<linename.Chromx3<<"  "<<linename.Chromy3<<endl;
  }
}

void Cal_Tune_vs_Deltap(Line & linename, const char *filename)
{
  int i;
  double qx[21],qy[21];
  double deltap;
  fstream fout;

  for(i=0;i<21;i++){
    deltap=GP.step_deltap*(i-10);
    Cal_Twiss(linename, deltap);
    qx[i]=linename.Tune1;
    qy[i]=linename.Tune2;
  }

  fout.open(filename, ios::out);
  for (i = 0; i < 21; i++)
    fout << setw(10) <<GP.step_deltap*(i-10)
  	 << scientific << setw(15) << qx[i]
  	 << scientific << setw(15) << qy[i]<<endl;
   fout.close();
}

void Plot_Tune_vs_Deltap(Line & linename, const char* filename)
{
  char command[256];
  fstream fout;

  fout.open("temp333.p", ios::out);
  fout<<"set term post color enhanced 20 solid "<<endl;
  fout<<"set output 'tune_vs_delta.ps' "<<endl;
  fout<<"set xlabel  'dp/p_0 [10^{-3}] '  " <<endl;
  fout<<"set ylabel 'Q_x' " <<endl;
  fout<<"set y2label 'Q_y' "<<endl;
  fout<<"set ytics nomirror"<<endl;
  fout<<"set y2tics"<<endl;
  
  sprintf(command, "plot '%s' u ($1*1000):2 tit 'Q_x' w l lt 1  lw 2,\\", filename );
  fout<<command<<endl;
  sprintf(command, "     '%s' u ($1*1000):3 axes x1y2 tit 'Q_y' w l lt 3  lw 2", filename );
  fout<<command<<endl;
  fout<<"exit"<<endl;
  fout.close();
  //system("gnuplot temp222.p");
  //system("rm temp222.p");
}

void Cal_Beta_Star_vs_Deltap(Line & linename, const char *filename)
{
  int i;
  double betx[21],bety[21];
  double deltap;
  fstream fout;

  for(i=0;i<21;i++){
    deltap=GP.step_deltap*(i-10);
    Cal_Twiss(linename, deltap);
    betx[i]=linename.Cell[linename.Ncell-1]->Beta1;
    bety[i]=linename.Cell[linename.Ncell-1]->Beta2;
  }

  fout.open(filename, ios::out);
  for (i = 0; i < 21; i++)
    fout << setw(10) <<GP.step_deltap*(i-10)
  	 << scientific << setw(15) << betx[i]
  	 << scientific << setw(15) << bety[i]<<endl;
   fout.close();
}

void Plot_Beta_Star_vs_Deltap(Line & linename, const char* filename)
{
  char command[256];
  fstream fout;

  fout.open("temp444.p", ios::out);
  fout<<"set term post color enhanced 20 "<<endl;
  fout<<"set output 'beta_star_vs_delta.ps' "<<endl;
  fout<<"set xlabel  'dp/p_0 [10^{-3}]'  " <<endl;
  fout<<"set ylabel '{/Symbol \142}*_x' " <<endl;
  fout<<"set y2label '{/Symbol \142}*_y' "<<endl;
  fout<<"set ytics nomirror"<<endl;
  fout<<"set y2tics"<<endl;
  
  sprintf(command, "plot '%s' u ($1*1000):2 tit '{/Symbol \142}*_x' w l lt 1  lw 2,\\", filename );
  fout<<command<<endl;
  sprintf(command, "     '%s' u ($1*1000):3 axes x1y2 tit '{/Symbol \142}*_y' w l lt 3  lw 2", filename );
  fout<<command<<endl;
  fout<<"exit"<<endl;
  fout.close();
  //system("gnuplot temp222.p");
  //system("rm temp222.p");
 }

void Cal_Beta_vs_Deltap(Line & linename, const char *filename)
{
  int i,j;
  double delta[21];
  double betax[21][linename.Ncell], betay[21][linename.Ncell];
  double DBXDD[linename.Ncell], DBYDD[linename.Ncell], DBXDD2[linename.Ncell], DBYDD2[linename.Ncell], DBXDD3[linename.Ncell], DBYDD3[linename.Ncell];
  double input1[21],input2[21];
  double term11, term12,term21, term22, term31, term32; 
  fstream fout;
  
  for(i=0;i<21;i++) delta[i]=GP.step_deltap*(i-10);
  
  for(i=0;i<21;i++){
    Cal_Twiss(linename, delta[i]);
    for(j=0;j<linename.Ncell;j++){
      betax[i][j]=linename.Cell[j]->Beta1;
      betay[i][j]=linename.Cell[j]->Beta2;
    }
  }

  for(j=0;j<linename.Ncell;j++) {
    for(i=0;i<21;i++){
      input1[i]=  betax[i][j];
      input2[i]=  betay[i][j];
    }
    chrom_fit(input1,input2,term11, term12,term21, term22, term31, term32);
    DBXDD[j]=term11;   DBXDD2[j]=term21*2;  DBXDD3[j]=term31*6;
    DBYDD[j]=term12;   DBYDD2[j]=term22*2;  DBYDD3[j]=term32*6;
  }
  
  fout.open(filename, ios::out);
  for(i=0;i<linename.Ncell;i++) 
    fout<<linename.Cell[i]->NAME<<setw(20)<<linename.Cell[i]->S
        <<setw(20)<<DBXDD[i]<<setw(20)<<DBXDD2[i]<<setw(20)<<DBXDD3[i]<<setw(20)<<DBYDD[i]<<setw(20)<<DBYDD2[i]<<setw(20)<<DBYDD3[i]<<endl;
  fout.close();
}

void Cal_Dispersion_vs_Deltap(Line & linename, const char *filename)
{
  int i,j;
  double delta[21];
  double xco[21][linename.Ncell], yco[21][linename.Ncell];
  double DX[linename.Ncell], DX2[linename.Ncell], DX3[linename.Ncell];
  double input1[21],input2[21];
  double term11, term12,term21, term22, term31, term32; 
  fstream fout;
 
  for(i=0;i<21;i++) delta[i]=GP.step_deltap*(i-10);
  
  for(i=0;i<21;i++){
    Cal_Twiss(linename, delta[i]);
    for(j=0;j<linename.Ncell;j++){
      xco[i][j]=linename.Cell[j]->X[0];
    }
  }
  
  for(j=0;j<linename.Ncell;j++) {
    for(i=0;i<21;i++){
      input1[i]=  xco[i][j];
      input2[i]=  yco[i][j];
    }
    chrom_fit(input1,input2,term11, term12,term21, term22, term31, term32);
    DX[j]=term11;  DX2[j]=term21*2;  DX3[j]=term31*6;
  }

  fout.open(filename, ios::out);
  for(i=0;i<linename.Ncell;i++) 
    fout<<setw(16)<<linename.Cell[i]->NAME<<setw(20)<<linename.Cell[i]->S<<setw(20)<< DX[i]<<setw(20)<<DX2[i]<<setw(20)<<DX3[i]<<endl;
  fout.close();
}

void Cal_Chromatic_Functions(Line & linename, const char* filename )
{
  int i,j;
  double delta[21];
  double xco[21][linename.Ncell], yco[21][linename.Ncell];
  double betax[21][linename.Ncell], betay[21][linename.Ncell];
  double DX[linename.Ncell], DX2[linename.Ncell], DX3[linename.Ncell];
  double DBXDD[linename.Ncell], DBYDD[linename.Ncell], DBXDD2[linename.Ncell], DBYDD2[linename.Ncell];
  double input1[21],input2[21];
  double term11, term12,term21, term22, term31, term32; 
  fstream fout;  

  fout.open(filename, ios::out);

  for(i=0;i<21;i++) delta[i]=GP.step_deltap*(i-10);
  for(i=0;i<21;i++){
    Cal_Twiss(linename, delta[i]);
    for(j=0;j<linename.Ncell;j++){
      xco[i][j]=linename.Cell[j]->X[0];
      betax[i][j]=linename.Cell[j]->Beta1;
      betay[i][j]=linename.Cell[j]->Beta2;
    }
  }

  for(j=0;j<linename.Ncell;j++) {
    for(i=0;i<21;i++){
      input1[i]=  xco[i][j];
      input2[i]=  yco[i][j];
    }
    chrom_fit(input1,input2,term11, term12,term21, term22, term31, term32);
    DX[j]=term11;  DX2[j]=term21*2;  DX3[j]=term31*6;
  }

  for(j=0;j<linename.Ncell;j++) {
    for(i=0;i<21;i++){
      input1[i]=  betax[i][j];
      input2[i]=  betay[i][j];
    }
    chrom_fit(input1,input2,term11, term12,term21, term22, term31, term32);
    DBXDD[j]=term11;   DBXDD2[j]=term21*2; 
    DBYDD[j]=term12;   DBYDD2[j]=term22*2; 
  }

  for(i=1;i<linename.Ncell;i++){
      fout<<setw(16)<<linename.Cell[i]->NAME<<setw(20)<<linename.Cell[i]->TYPE<<setw(20)<<linename.Cell[i]->S
          <<setw(20)<<DX[i]<<setw(20)<<DX2[i]<<setw(20)<<DX3[i]
	  <<setw(20)<<DBXDD[i]<<setw(20)<<DBXDD2[i]<<setw(20)<<DBYDD[i]<<setw(20)<<DBYDD2[i]<<endl;
  }
  
  fout.close();
}

void Cal_W_Functions(Line & linename, const char *filename)
{
  int i,j;
  double delta[21];
  double BX[linename.Ncell], BY[linename.Ncell],AX[linename.Ncell], AY[linename.Ncell];
  double betax[21][linename.Ncell], betay[21][linename.Ncell], alfax[21][linename.Ncell], alfay[21][linename.Ncell];
  double DBXDD[linename.Ncell], DBYDD[linename.Ncell],  DAXDD[linename.Ncell], DAYDD[linename.Ncell];
  double input1[21],input2[21];
  double term11, term12,term21, term22, term31, term32; 
  fstream fout;
  double b,a, WX, WY;
  
  Cal_Twiss(linename,0.);
  for(j=0;j<linename.Ncell;j++){
    BX[j]=linename.Cell[j]->Beta1;
    BY[j]=linename.Cell[j]->Beta2;
    AX[j]=linename.Cell[j]->Alfa1;
    AY[j]=linename.Cell[j]->Alfa2;
  }

  for(i=0;i<21;i++) delta[i]=GP.step_deltap*(i-10);
  for(i=0;i<21;i++){
    Cal_Twiss(linename, delta[i]);
    for(j=0;j<linename.Ncell;j++){
      betax[i][j]=linename.Cell[j]->Beta1;
      betay[i][j]=linename.Cell[j]->Beta2;
      alfax[i][j]=linename.Cell[j]->Alfa1;
      alfay[i][j]=linename.Cell[j]->Alfa2;
    }
  }

  for(j=0;j<linename.Ncell;j++) {
    for(i=0;i<21;i++){
      input1[i]=  betax[i][j];
      input2[i]=  betay[i][j];
    }
    chrom_fit(input1,input2,term11, term12,term21, term22, term31, term32);
    DBXDD[j]=term11; 
    DBYDD[j]=term12; 
  }
  for(j=0;j<linename.Ncell;j++) {
    for(i=0;i<21;i++){
      input1[i]=  alfax[i][j];
      input2[i]=  alfay[i][j];
    }
    chrom_fit(input1,input2,term11, term12,term21, term22, term31, term32);
    DAXDD[j]=term11; 
    DAYDD[j]=term12; 
  }  
  
  fout.open(filename, ios::out);
  for(i=0;i<linename.Ncell;i++){ 
    b = DBXDD[i]/BX[i] ;
    a = DAXDD[i] - AX[i]*DBXDD[i]/BX[i];
    WX= sqrt(a*a + b*b);
    b = DBYDD[i]/BY[i] ;
    a = DAYDD[i]-AY[i]*DBYDD[i]/BY[i];
    WY= sqrt(a*a + b*b);
    fout<<setw(16)<<linename.Cell[i]->NAME<<setw(20)<<linename.Cell[i]->S
        <<setw(20)<<DBXDD[i]<<setw(20)<<DBYDD[i]<<setw(20)<<WX<<setw(20)<<WY<<endl;
  }
  fout.close();

}



void Cal_Half_Integer_RDT(Line & linename,  const char* filename)
//  horizontal: h20001, vertical: h00021
{
  int i,j;
  fstream fout;
  double h_real, h_imag, v_real, v_imag;
  double k1l, k2l, betx, bety, Dx,dphix, dphiy; 

  Cal_Twiss(linename,0.0);
  Cal_Dispersion(linename);
  fout.open(filename, ios::out);

  for(i=0;i<linename.Ncell;i++){
    h_real=0; h_imag=0;
    v_real=0; v_imag=0;
    
    for(j=0;j<linename.Ncell;j++){
      if( linename.Cell[j]->TYPE == string("QUAD") ){
	k1l=linename.Cell[j]->GetP("K1L");
	betx=linename.Cell[j]->Beta1;
	bety=linename.Cell[j]->Beta2;
	Dx=linename.Cell[j]->Etax;
	dphix= abs(linename.Cell[j]->Mu1 - linename.Cell[i]->Mu1);
	dphiy= abs(linename.Cell[j]->Mu2 - linename.Cell[i]->Mu2);        
	h_real = h_real  -  k1l*betx*cos(  2*dphix * 2 * PI ) ; 
	h_imag = h_imag  -  k1l*betx*sin(  2*dphix * 2 * PI );	
	v_real = v_real  +  k1l*bety*cos(  2*dphiy * 2 * PI ) ; 
	v_imag = v_imag  +  k1l*bety*sin(  2*dphiy * 2 * PI ) ; 
      }
      else if(linename.Cell[j]->TYPE == string("SEXT") ){
	k2l=linename.Cell[j]->GetP("K2L");
	betx=linename.Cell[j]->Beta1;
	bety=linename.Cell[j]->Beta2;
	Dx=linename.Cell[j]->Etax;
	dphix= abs(linename.Cell[j]->Mu1 - linename.Cell[i]->Mu1);
	dphiy= abs(linename.Cell[j]->Mu2 - linename.Cell[i]->Mu2);  
	h_real = h_real +  k2l*betx*Dx*cos(  2*dphix * 2 * PI ) ;  
	h_imag = h_imag +  k2l*betx*Dx*sin(  2*dphix * 2 * PI );	
	v_real = v_real -  k2l*bety*Dx*cos(  2*dphiy * 2 * PI ) ;  
	v_imag = v_imag -  k2l*bety*Dx*sin(  2*dphiy * 2 * PI ) ;  
      }
      else{}
    }
    fout<<setw(16)<<linename.Cell[i]->NAME<<setw(20)<<linename.Cell[i]->S<<setw(20)<<linename.Cell[i]->Beta1<<setw(20)<<linename.Cell[i]->Beta2<<setw(20)<<sqrt(h_real*h_real+h_imag*h_imag)<<"  "<<sqrt(v_real*v_real+v_imag*v_imag)<<endl;
  }

  fout.close();
}

void Cal_Half_Integer_RDT_SextFamily(Line & linename,  const char* sextname)
{
  int j;
  double h_real, h_imag, v_real, v_imag;
  double k2l, betx, bety, Dx,dphix, dphiy; 
  
  Cal_Twiss(linename,0.0);  Cal_Dispersion(linename);
  
  h_real=0; h_imag=0;  v_real=0; v_imag=0;
  for(j=0;j<linename.Ncell;j++){
    if( linename.Cell[j]->TYPE  == string("SEXT")  and linename.Cell[j]->NAME  == string(sextname)  ){
      k2l=linename.Cell[j]->GetP("K2L");
      betx=linename.Cell[j]->Beta1;
      bety=linename.Cell[j]->Beta2;
      Dx=linename.Cell[j]->Etax;
      dphix= linename.Cell[j]->Mu1;
      dphiy= linename.Cell[j]->Mu2;
      h_real = h_real +  k2l*betx*Dx*cos(  2*dphix * 2 * PI ) ;  
      h_imag = h_imag +  k2l*betx*Dx*sin(  2*dphix * 2 * PI );	
      v_real = v_real -  k2l*bety*Dx*cos(  2*dphiy * 2 * PI ) ;  
      v_imag = v_imag -  k2l*bety*Dx*sin(  2*dphiy * 2 * PI ) ;  
    }
  }
  cout<<setw(15)<<sextname<<setw(15)<<h_real<<setw(15)<<h_imag<<setw(15)<<atan2(h_imag, h_real)*180/PI
      <<setw(15)<<v_real<<setw(15)<<v_imag<<setw(15)<<atan2(v_imag, v_real)*180/PI<<endl;
}



void Cal_Q2_Source(Line & linename, const char* filename )
{
  int i,j;
  double delta[21];
  double xco[21][linename.Ncell], yco[21][linename.Ncell];
  double betax[21][linename.Ncell], betay[21][linename.Ncell];
  double DX[linename.Ncell], DX2[linename.Ncell], DX3[linename.Ncell];
  double DBXDD[linename.Ncell], DBYDD[linename.Ncell], DBXDD2[linename.Ncell], DBYDD2[linename.Ncell];
  double input1[21],input2[21];
  double term11, term12,term21, term22, term31, term32; 
  fstream fout;  

  fout.open(filename, ios::out);

  for(i=0;i<21;i++) delta[i]=GP.step_deltap*(i-10);
  for(i=0;i<21;i++){
    Cal_Twiss(linename, delta[i]);
    for(j=0;j<linename.Ncell;j++){
      xco[i][j]=linename.Cell[j]->X[0];
      betax[i][j]=linename.Cell[j]->Beta1;
      betay[i][j]=linename.Cell[j]->Beta2;
    }
  }

  for(j=0;j<linename.Ncell;j++) {
    for(i=0;i<21;i++){
      input1[i]=  xco[i][j];
      input2[i]=  yco[i][j];
    }
    chrom_fit(input1,input2,term11, term12,term21, term22, term31, term32);
    DX[j]=term11;  DX2[j]=term21*2;  DX3[j]=term31*6;
  }

  for(j=0;j<linename.Ncell;j++) {
    for(i=0;i<21;i++){
      input1[i]=  betax[i][j];
      input2[i]=  betay[i][j];
    }
    chrom_fit(input1,input2,term11, term12,term21, term22, term31, term32);
    DBXDD[j]=term11;   DBXDD2[j]=term21*2; 
    DBYDD[j]=term12;   DBYDD2[j]=term22*2; 
  }

  double chromx1,   chromy1, chromx2,  chromy2,  chromx3, chromy3;
  Cal_Twiss(linename,0.);

  for(i=1;i<linename.Ncell;i++){
    
    chromx1=0;  chromy1=0 ; chromx2=0 ;  chromy2=0;  chromx3=0; chromy3=0;
    
    if(linename.Cell[i]->TYPE ==string("SEXT") ) {
      
      chromx1 +=  1.0*( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5*DX[i]*linename.Cell[i]->GetP("K2L")/4/PI;
      chromy1 +=- 1.0*( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5*DX[i]*linename.Cell[i]->GetP("K2L")/4/PI;
      
      chromx2 += -1.0*( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5*DX[i]*linename.Cell[i]->GetP("K2L")/4/PI/2;
      chromy2 += +1.0*( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5*DX[i]*linename.Cell[i]->GetP("K2L")/4/PI/2;
      
      chromx2 +=  linename.Cell[i]->GetP("K2L") * DX[i]  * DBXDD[i]/4/PI/2;
      chromy2 += -linename.Cell[i]->GetP("K2L") * DX[i]  * DBYDD[i]/4/PI/2;
      
      chromx2 +=  linename.Cell[i]->GetP("K2L") * ( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5 * DX2[i]/4/PI/2;
      chromy2 += -linename.Cell[i]->GetP("K2L") * ( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5 * DX2[i]/4/PI/2;
      
      chromx3 += +1.0*( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5*DX[i]*linename.Cell[i]->GetP("K2L")/8/PI/3;
      chromy3 += -1.0*( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5*DX[i]*linename.Cell[i]->GetP("K2L")/8/PI/3;
      
      chromx3 += -linename.Cell[i]->GetP("K2L") * DX[i]  * DBXDD[i]/8/PI/3;
      chromy3 += +linename.Cell[i]->GetP("K2L") * DX[i]  * DBYDD[i]/8/PI/3;
      
      chromx3 += -linename.Cell[i]->GetP("K2L") * ( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5 * DX2[i]/8/PI/3;
      chromy3 += +linename.Cell[i]->GetP("K2L") * ( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5 * DX2[i]/8/PI/3; 
      
      chromx3 +=  linename.Cell[i]->GetP("K2L") *  DX[i] * DBXDD2[i]/8/PI/3;
      chromy3 += -linename.Cell[i]->GetP("K2L") *  DX[i] * DBYDD2[i]/8/PI/3;
      
      chromx3 += -linename.Cell[i]->GetP("K2L") *  DX[i] * DBXDD[i]/8/PI/3;
      chromy3 += +linename.Cell[i]->GetP("K2L") *  DX[i] * DBYDD[i]/8/PI/3;       
      
      
      chromx3 += +linename.Cell[i]->GetP("K2L") *  DX2[i] * DBXDD[i]/8/PI/3;
      chromy3 += -linename.Cell[i]->GetP("K2L") *  DX2[i] * DBYDD[i]/8/PI/3; 
      
      chromx3 += +linename.Cell[i]->GetP("K2L") *  ( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5 * DX3[i]/8/PI/3;
      chromy3 += -linename.Cell[i]->GetP("K2L") *  ( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5 * DX3[i]/8/PI/3;
      
      chromx3 += -linename.Cell[i]->GetP("K2L") *  ( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5 * DX2[i]/8/PI/3;
      chromy3 += +linename.Cell[i]->GetP("K2L") *  ( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5 * DX2[i]/8/PI/3;
      
      chromx3 += +linename.Cell[i]->GetP("K2L") *  DBXDD[i] * DX2[i]/8/PI/3;
      chromy3 += -linename.Cell[i]->GetP("K2L") *  DBYDD[i] * DX2[i]/8/PI/3;        
      
      fout<<setw(16)<<linename.Cell[i]->NAME<<setw(20)<<linename.Cell[i]->TYPE<<setw(20)<<linename.Cell[i]->S
          <<setw(20)<<chromx1<<setw(20)<<chromy1<<setw(20)<<chromx2<<setw(20)<<chromy2<<setw(20)<<chromx3<<setw(20)<<chromy3<<endl;
    }
    
    if(linename.Cell[i]->TYPE ==string("QUAD") ) {
      chromx1 += - 1.0*( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5 * linename.Cell[i]->GetP("K1L")/4/PI;
      chromy1 += + 1.0*( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5 * linename.Cell[i]->GetP("K1L")/4/PI;
      
      chromx2 += +1.0*( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5 * linename.Cell[i]->GetP("K1L")/4/PI/2;
      chromy2 += -1.0*( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5 * linename.Cell[i]->GetP("K1L")/4/PI/2;
      
      chromx2 += -linename.Cell[i]->GetP("K1L") * DBXDD[i]/4/PI/2;
      chromy2 += +linename.Cell[i]->GetP("K1L") * DBYDD[i]/4/PI/2;
      
      chromx3 += -1.0*( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5 * linename.Cell[i]->GetP("K1L")/8/PI/3;
      chromy3 += +1.0*( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5 * linename.Cell[i]->GetP("K1L")/8/PI/3;
      
      chromx3 += +linename.Cell[i]->GetP("K1L") * DBXDD[i]/8/PI/3;
      chromy3 += -linename.Cell[i]->GetP("K1L") * DBYDD[i]/8/PI/3;
      
      chromx3 += -linename.Cell[i]->GetP("K1L") * DBXDD2[i]/8/PI/3;
      chromy3 += +linename.Cell[i]->GetP("K1L") * DBYDD2[i]/8/PI/3;
      
      chromx3 += +linename.Cell[i]->GetP("K1L") * DBXDD[i]/8/PI/3;
      chromy3 += -linename.Cell[i]->GetP("K1L") * DBYDD[i]/8/PI/3;
      
      fout<<setw(16)<<linename.Cell[i]->NAME<<setw(20)<<linename.Cell[i]->TYPE<<setw(20)<<linename.Cell[i]->S
          <<setw(20)<<chromx1<<setw(20)<<chromy1<<setw(20)<<chromx2<<setw(20)<<chromy2<<setw(20)<<chromx3<<setw(20)<<chromy3<<endl;
    }
  }
  fout.close();
}


void Cal_Q2_Source_Section(Line & linename, const char* filename )
{
  int i,j;
  double delta[21];
  double xco[21][linename.Ncell], yco[21][linename.Ncell];
  double betax[21][linename.Ncell], betay[21][linename.Ncell];
  double DX[linename.Ncell], DX2[linename.Ncell], DX3[linename.Ncell];
  double DBXDD[linename.Ncell], DBYDD[linename.Ncell], DBXDD2[linename.Ncell], DBYDD2[linename.Ncell];
  double input1[21],input2[21];
  double term11, term12,term21, term22, term31, term32; 
  
  const char * name[13]={"IR6", "ARC0604","IR4", "ARC0402", "IR2", "ARC0212","IR12", "ARC1210","IR10", "ARC1008","IR8", "ARC0806","SUM"};
  double chrom1x[13], chrom1y[13], chrom2x[13], chrom2y[13], chrom3x[13], chrom3y[13];
  double s[12]={143,484,782,1133,1421,1770,2060,2412,2699,3052,3338,3690};
  double chromx1, chromy1, chromx1_sum=0, chromy1_sum=0;
  double chromx2, chromy2, chromx2_sum=0, chromy2_sum=0;
  double chromx3, chromy3, chromx3_sum=0, chromy3_sum=0;

  fstream fout;  

  fout.open(filename, ios::out);

  for(i=0;i<21;i++) delta[i]=GP.step_deltap*(i-10);
  for(i=0;i<21;i++){
    Cal_Twiss(linename, delta[i]);
    for(j=0;j<linename.Ncell;j++){
      xco[i][j]=linename.Cell[j]->X[0];
      betax[i][j]=linename.Cell[j]->Beta1;
      betay[i][j]=linename.Cell[j]->Beta2;
    }
  }

  for(j=0;j<linename.Ncell;j++) {
    for(i=0;i<21;i++){
      input1[i]=  xco[i][j];
      input2[i]=  yco[i][j];
    }
    chrom_fit(input1,input2,term11, term12,term21, term22, term31, term32);
    DX[j]=term11;  DX2[j]=term21*2;  DX3[j]=term31*6;
  }

  for(j=0;j<linename.Ncell;j++) {
    for(i=0;i<21;i++){
      input1[i]=  betax[i][j];
      input2[i]=  betay[i][j];
    }
    chrom_fit(input1,input2,term11, term12,term21, term22, term31, term32);
    DBXDD[j]=term11;   DBXDD2[j]=term21*2;  
    DBYDD[j]=term12;   DBYDD2[j]=term22*2;  
  }

  Cal_Twiss(linename,0.);
 
  //------Q' x,y calculation
  chromx1 = 0.; chromy1 = 0. ;
  for(i=0;i<linename.Ncell;i++){
    if( linename.Cell[i]->S < s[0] or  linename.Cell[i]->S > s[11] ) {
      if(linename.Cell[i]->TYPE ==string("SEXT") ) {
	chromx1 +=  1.0*linename.Cell[i]->Beta1*linename.Cell[i]->Etax*linename.Cell[i]->GetP("K2L")/4/PI;
	chromy1 +=- 1.0*linename.Cell[i]->Beta2*linename.Cell[i]->Etax*linename.Cell[i]->GetP("K2L")/4/PI;
      }
      
      if(linename.Cell[i]->TYPE ==string("QUAD") ) {
	chromx1 += - 1.0*( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5 * linename.Cell[i]->GetP("K1L")/4/PI;
	chromy1 += + 1.0*( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5 * linename.Cell[i]->GetP("K1L")/4/PI;
      }
    }
  }
  chrom1x[0]= chromx1;  chrom1y[0]= chromy1;
  chromx1_sum += chromx1;  chromy1_sum += chromy1; 
  
  for(j=0;j<11;j++){
    chromx1 = 0.; chromy1 = 0. ;
    for(i=0;i<linename.Ncell;i++){ 
      if( linename.Cell[i]->S > s[j]  and   linename.Cell[i]->S < s[j+1] ) {
	
	if(linename.Cell[i]->TYPE ==string("SEXT") ) {
	  chromx1 += + 1.0*linename.Cell[i]->Beta1*linename.Cell[i]->Etax*linename.Cell[i]->GetP("K2L")/4/PI;
	  chromy1 += - 1.0*linename.Cell[i]->Beta2*linename.Cell[i]->Etax*linename.Cell[i]->GetP("K2L")/4/PI;
	}
	
	if(linename.Cell[i]->TYPE ==string("QUAD") ) {
	  chromx1 += - 1.0*( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5 * linename.Cell[i]->GetP("K1L")/4/PI;
	  chromy1 += + 1.0*( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5 * linename.Cell[i]->GetP("K1L")/4/PI;
	}
      }
    }
    chrom1x[j+1]= chromx1;  chrom1y[j+1]= chromy1;
    chromx1_sum += chromx1;  chromy1_sum += chromy1;    
  }

  chrom1x[12]= chromx1_sum;  chrom1y[12]= chromy1_sum;

  //-----Q''x,y
  chromx2=0; chromy2=0;
  for(i=0;i<linename.Ncell;i++){
    if( linename.Cell[i]->S < s[0] or  linename.Cell[i]->S > s[11] ) {
      if(linename.Cell[i]->TYPE ==string("SEXT") ) {
	chromx2 += -1.0*( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5*DX[i]*linename.Cell[i]->GetP("K2L")/4/PI/2;
	chromy2 += +1.0*( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5*DX[i]*linename.Cell[i]->GetP("K2L")/4/PI/2;
	
	chromx2 +=  linename.Cell[i]->GetP("K2L") * DX[i]  * DBXDD[i]/4/PI/2;
	chromy2 += -linename.Cell[i]->GetP("K2L") * DX[i]  * DBYDD[i]/4/PI/2;
	
	chromx2 +=  linename.Cell[i]->GetP("K2L") * ( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5 * DX2[i]/4/PI/2;
	chromy2 += -linename.Cell[i]->GetP("K2L") * ( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5 * DX2[i]/4/PI/2;      
      }
      if(linename.Cell[i]->TYPE ==string("QUAD") ) {
	chromx2 += +1.0*( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5 * linename.Cell[i]->GetP("K1L")/4/PI/2;
	chromy2 += -1.0*( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5 * linename.Cell[i]->GetP("K1L")/4/PI/2;
	
	chromx2 += -linename.Cell[i]->GetP("K1L") * DBXDD[i]/4/PI/2;
	chromy2 += +linename.Cell[i]->GetP("K1L") * DBYDD[i]/4/PI/2;
      }
    }
  }
  chrom2x[0]= chromx2;  chrom2y[0]= chromy2;
  chromx2_sum += chromx2;  chromy2_sum += chromy2; 
  
  for(j=0;j<11;j++){
    chromx2=0; chromy2=0;
    for(i=0;i<linename.Ncell;i++){ 
      if( linename.Cell[i]->S > s[j]  and   linename.Cell[i]->S < s[j+1] ) {
	if(linename.Cell[i]->TYPE ==string("SEXT") ) {
	  chromx2 += -1.0*( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5*DX[i]*linename.Cell[i]->GetP("K2L")/4/PI/2;
	  chromy2 += +1.0*( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5*DX[i]*linename.Cell[i]->GetP("K2L")/4/PI/2;
	  
	  chromx2 +=  linename.Cell[i]->GetP("K2L") * DX[i]  * DBXDD[i]/4/PI/2;
	  chromy2 += -linename.Cell[i]->GetP("K2L") * DX[i]  * DBYDD[i]/4/PI/2;
	  
	  chromx2 +=  linename.Cell[i]->GetP("K2L") * ( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5 * DX2[i]/4/PI/2;
	  chromy2 += -linename.Cell[i]->GetP("K2L") * ( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5 * DX2[i]/4/PI/2;      
	}
	if(linename.Cell[i]->TYPE ==string("QUAD") ) {
	  chromx2 += +1.0*( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5 * linename.Cell[i]->GetP("K1L")/4/PI/2;
	  chromy2 += -1.0*( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5 * linename.Cell[i]->GetP("K1L")/4/PI/2;
	  
	  chromx2 += -linename.Cell[i]->GetP("K1L") * DBXDD[i]/4/PI/2;
	  chromy2 += +linename.Cell[i]->GetP("K1L") * DBYDD[i]/4/PI/2;
	}
      }
    }
    chrom2x[j+1]= chromx2;  chrom2y[j+1]= chromy2;
    chromx2_sum += chromx2;  chromy2_sum += chromy2;
  }

  chrom2x[12]= chromx2_sum;  chrom2y[12]= chromy2_sum;

 //-----Q'''x,y
  chromx3=0; chromy3=0;
  for(i=0;i<linename.Ncell;i++){
    if( linename.Cell[i]->S < s[0] or  linename.Cell[i]->S > s[11] ) {
      if(linename.Cell[i]->TYPE ==string("SEXT") ) {
	
        //---contribution from second order
	chromx3 += +1.0*( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5*DX[i]*linename.Cell[i]->GetP("K2L")/8/PI/3;
	chromy3 += -1.0*( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5*DX[i]*linename.Cell[i]->GetP("K2L")/8/PI/3;
	
	chromx3 += -linename.Cell[i]->GetP("K2L") * DX[i]  * DBXDD[i]/8/PI/3;
	chromy3 += +linename.Cell[i]->GetP("K2L") * DX[i]  * DBYDD[i]/8/PI/3;
	
	chromx3 += -linename.Cell[i]->GetP("K2L") * ( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5 * DX2[i]/8/PI/3;
	chromy3 += +linename.Cell[i]->GetP("K2L") * ( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5 * DX2[i]/8/PI/3; 
	
	//---contribution from third order
	
        chromx3 +=  linename.Cell[i]->GetP("K2L") *  DX[i] * DBXDD2[i]/8/PI/3;
	chromy3 += -linename.Cell[i]->GetP("K2L") *  DX[i] * DBYDD2[i]/8/PI/3;
	
	chromx3 += -linename.Cell[i]->GetP("K2L") *  DX[i] * DBXDD[i]/8/PI/3;
	chromy3 += +linename.Cell[i]->GetP("K2L") *  DX[i] * DBYDD[i]/8/PI/3;       
	
        
	chromx3 += +linename.Cell[i]->GetP("K2L") *  DX2[i] * DBXDD[i]/8/PI/3;
	chromy3 += -linename.Cell[i]->GetP("K2L") *  DX2[i] * DBYDD[i]/8/PI/3; 
	
	chromx3 += +linename.Cell[i]->GetP("K2L") *  ( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5 * DX3[i]/8/PI/3;
	chromy3 += -linename.Cell[i]->GetP("K2L") *  ( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5 * DX3[i]/8/PI/3;
	
	chromx3 += -linename.Cell[i]->GetP("K2L") *  ( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5 * DX2[i]/8/PI/3;
	chromy3 += +linename.Cell[i]->GetP("K2L") *  ( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5 * DX2[i]/8/PI/3;
	
	chromx3 += +linename.Cell[i]->GetP("K2L") *  DBXDD[i] * DX2[i]/8/PI/3;
	chromy3 += -linename.Cell[i]->GetP("K2L") *  DBYDD[i] * DX2[i]/8/PI/3;        
	
      }

      if(linename.Cell[i]->TYPE ==string("QUAD") ) {
	
	//----contribition from second order 
	chromx3 += -1.0*( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5 * linename.Cell[i]->GetP("K1L")/8/PI/3;
	chromy3 += +1.0*( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5 * linename.Cell[i]->GetP("K1L")/8/PI/3;
	
	chromx3 += +linename.Cell[i]->GetP("K1L") * DBXDD[i]/8/PI/3;
	chromy3 += -linename.Cell[i]->GetP("K1L") * DBYDD[i]/8/PI/3;
	
        //----contrbution from third order
	
        chromx3 += -linename.Cell[i]->GetP("K1L") * DBXDD2[i]/8/PI/3;
	chromy3 += +linename.Cell[i]->GetP("K1L") * DBYDD2[i]/8/PI/3;
	
	chromx3 += +linename.Cell[i]->GetP("K1L") * DBXDD[i]/8/PI/3;
	chromy3 += -linename.Cell[i]->GetP("K1L") * DBYDD[i]/8/PI/3;
	
      }
    }
  }
  chrom3x[0]= chromx3;  chrom3y[0]= chromy3;
  chromx3_sum += chromx3;  chromy3_sum += chromy3; 
  
  for(j=0;j<11;j++){
    chromx3=0; chromy3=0;
    for(i=0;i<linename.Ncell;i++){ 
      if( linename.Cell[i]->S > s[j]  and   linename.Cell[i]->S < s[j+1] ) {
	if(linename.Cell[i]->TYPE ==string("SEXT") ) {
	  
	  //---contribution from second order
	  chromx3 += +1.0*( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5*DX[i]*linename.Cell[i]->GetP("K2L")/8/PI/3;
	  chromy3 += -1.0*( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5*DX[i]*linename.Cell[i]->GetP("K2L")/8/PI/3;
	  
	  chromx3 += -linename.Cell[i]->GetP("K2L") * DX[i]  * DBXDD[i]/8/PI/3;
	  chromy3 += +linename.Cell[i]->GetP("K2L") * DX[i]  * DBYDD[i]/8/PI/3;
	  
	  chromx3 += -linename.Cell[i]->GetP("K2L") * ( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5 * DX2[i]/8/PI/3;
	  chromy3 += +linename.Cell[i]->GetP("K2L") * ( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5 * DX2[i]/8/PI/3; 
	  
	  //---contribution from third order
	  
	  chromx3 +=  linename.Cell[i]->GetP("K2L") *  DX[i] * DBXDD2[i]/8/PI/3;
	  chromy3 += -linename.Cell[i]->GetP("K2L") *  DX[i] * DBYDD2[i]/8/PI/3;
	  
	  chromx3 += -linename.Cell[i]->GetP("K2L") *  DX[i] * DBXDD[i]/8/PI/3;
	  chromy3 += +linename.Cell[i]->GetP("K2L") *  DX[i] * DBYDD[i]/8/PI/3;       
	  
	  
	  chromx3 += +linename.Cell[i]->GetP("K2L") *  DX2[i] * DBXDD[i]/8/PI/3;
	  chromy3 += -linename.Cell[i]->GetP("K2L") *  DX2[i] * DBYDD[i]/8/PI/3; 
	  
	  chromx3 += +linename.Cell[i]->GetP("K2L") *  ( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5 * DX3[i]/8/PI/3;
	  chromy3 += -linename.Cell[i]->GetP("K2L") *  ( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5 * DX3[i]/8/PI/3;
	  
	  chromx3 += -linename.Cell[i]->GetP("K2L") *  ( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5 * DX2[i]/8/PI/3;
	  chromy3 += +linename.Cell[i]->GetP("K2L") *  ( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5 * DX2[i]/8/PI/3;
	  
	  chromx3 += +linename.Cell[i]->GetP("K2L") *  DBXDD[i] * DX2[i]/8/PI/3;
	  chromy3 += -linename.Cell[i]->GetP("K2L") *  DBYDD[i] * DX2[i]/8/PI/3;        
	  
	}

	if(linename.Cell[i]->TYPE ==string("QUAD") ) {
	  
	  //----contribition from second order 
	  chromx3 += -1.0*( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5 * linename.Cell[i]->GetP("K1L")/8/PI/3;
	  chromy3 += +1.0*( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5 * linename.Cell[i]->GetP("K1L")/8/PI/3;
	  
	  chromx3 += +linename.Cell[i]->GetP("K1L") * DBXDD[i]/8/PI/3;
	  chromy3 += -linename.Cell[i]->GetP("K1L") * DBYDD[i]/8/PI/3;
	  
	  //----contrbution from third order
	  
	  chromx3 += -linename.Cell[i]->GetP("K1L") * DBXDD2[i]/8/PI/3;
	  chromy3 += +linename.Cell[i]->GetP("K1L") * DBYDD2[i]/8/PI/3;
	  
	  chromx3 += +linename.Cell[i]->GetP("K1L") * DBXDD[i]/8/PI/3;
	  chromy3 += -linename.Cell[i]->GetP("K1L") * DBYDD[i]/8/PI/3;
	  
	}
      }
    }
    chrom3x[j+1]= chromx3;  chrom3y[j+1]= chromy3;
    chromx3_sum += chromx3;  chromy3_sum += chromy3;
  }

  chrom3x[12]= chromx3_sum;  chrom3y[12]= chromy3_sum;

  //----save results
  for(i=0;i<13;i++){
    chrom1x[i] =  int(chrom1x[i] *100 )/100.;
    chrom2x[i] =  int(chrom2x[i] *100 )/100.;
    chrom3x[i] =  int(chrom3x[i] *100 )/100.;
    chrom1y[i] =  int(chrom1y[i] *100 )/100.;
    chrom2y[i] =  int(chrom2y[i] *100 )/100.;
    chrom3y[i] =  int(chrom3y[i] *100 )/100.;
  }


  cout<<"  "<<endl;
  cout<<"...chrom sucrce:"<<endl;
  cout<<"----------------------------------------------------------------------------------------------------------"<<endl;
  cout<<"                    Q'x          Q'y        Qx''/2       Qy''/2       Qx'''/6       Qy'''/6      "<<endl;
  cout<<"----------------------------------------------------------------------------------------------------------"<<endl;
  cout<<"   actual values:"<<endl;
  for(i=0;i<12;i++)
    cout<<setw(13)<<name[i]<<setw(13)<<chrom1x[i]<<setw(13)<<chrom1y[i]<<setw(13)<<chrom2x[i]<<setw(13)<<chrom2y[i]<<setw(13)<<chrom3x[i]<<setw(13)<<chrom3y[i]<<endl; 
  cout<<setw(13)<<name[i]<<setw(13)<<chrom1x[12]<<setw(13)<<chrom1y[12]<<setw(13)<<chrom2x[12]<<setw(13)<<chrom2y[12]<<setw(13)<<chrom3x[12]<<setw(13)<<chrom3y[12]<<endl;
  cout<<"  percentages :"<<setprecision(5)<<endl;
  for(i=0;i<12;i++)
    cout<<setw(13)<<name[i]<<setw(13)<<chrom1x[i]/abs(chrom1x[12])<<setw(13)<<chrom1y[i]/abs(chrom1y[12])<<setw(13)<<chrom2x[i]/abs(chrom2x[12])<<setw(13)<<chrom2y[i]/abs(chrom2y[12])<<setw(13)<<chrom3x[i]/abs(chrom3x[12])<<setw(13)<<chrom3y[i]/abs(chrom3y[12])<<endl; 
  cout<<setw(13)<<name[i]<<setw(13)<<chrom1x[12]/abs(chrom1x[12])<<setw(13)<<chrom1y[12]/abs(chrom1y[12])<<setw(13)<<chrom2x[12]/abs(chrom2x[12])<<setw(13)<<chrom2y[12]/abs(chrom2y[12])<<setw(13)<<chrom3x[12]/abs(chrom3x[12])<<setw(13)<<chrom3y[12]/abs(chrom3y[12])<<endl;
  cout<<"----------------------------------------------------------------------------------------------------------"<<endl; 
  cout<<"   Summary:"<<endl;
  cout<<"   IR6 and IR8 contribute:  "<<( chrom1x[0] + chrom1x[10] ) / abs(chrom1x[12])<<"  to Qx' "<<endl;
  cout<<"   IR6 and IR8 contribute:  "<<( chrom1y[0] + chrom1y[10] ) / abs(chrom1y[12])<<"  to Qy' "<<endl;
  cout<<"   IR6 and IR8 contribute:  "<<( chrom2x[0] + chrom2x[10] ) / abs(chrom2x[12])<<"  to Qx'' "<<endl;
  cout<<"   IR6 and IR8 contribute:  "<<( chrom2y[0] + chrom2y[10] ) / abs(chrom2y[12])<<"  to Qy'' "<<endl;
  cout<<"   IR6 and IR8 contribute:  "<<( chrom3x[0] + chrom3x[10] ) / abs(chrom3x[12])<<"  to Qx''' "<<endl;
  cout<<"   IR6 and IR8 contribute:  "<<( chrom3y[0] + chrom3y[10] ) / abs(chrom3y[12])<<"  to Qy''' "<<endl; 
  cout<<"  "<<endl;
  cout<<"   Other IRs contribute:  "<<( chrom1x[4] + chrom1x[6] +  chrom1x[8] +  chrom1x[2] )/ abs(chrom1x[12]) <<"  to Qx' "<<endl;
  cout<<"   Other IRs contribute:  "<<( chrom1y[4] + chrom1y[6] +  chrom1y[8] +  chrom1y[2] )/ abs(chrom1y[12]) <<"  to Qy' "<<endl;
  cout<<"   Other IRs contribute:  "<<( chrom2x[4] + chrom2x[6] +  chrom2x[8] +  chrom2x[2] )/ abs(chrom2x[12]) <<"  to Qx'' "<<endl;
  cout<<"   Other IRs contribute:  "<<( chrom2y[4] + chrom2y[6] +  chrom2y[8] +  chrom2y[2] )/ abs(chrom2y[12]) <<"  to Qy'' "<<endl;
  cout<<"   Other IRs contribute:  "<<( chrom3x[4] + chrom3x[6] +  chrom3x[8] +  chrom3x[2] )/ abs(chrom3x[12]) <<"  to Qx''' "<<endl;
  cout<<"   Other IRs contribute:  "<<( chrom3y[4] + chrom3y[6] +  chrom3y[8] +  chrom3y[2] )/ abs(chrom3y[12]) <<"  to Qy''' "<<endl;  
  cout<<"  "<<endl;
  cout<<"   Arcs contribute:  "<< ( chrom1x[1] + chrom1x[3] +  chrom1x[5] +  chrom1x[7]  +  chrom1x[9]  + chrom1x[11] )/ abs(chrom1x[12]) <<"  to Qx' "<<endl;
  cout<<"   Arcs contribute:  " <<( chrom1y[1] + chrom1y[3] +  chrom1y[5] +  chrom1y[7]  +  chrom1y[9]  + chrom1y[11] )/ abs(chrom1y[12])<<"  to Qy' "<<endl;
  cout<<"   Arcs  contribute:  "<<( chrom2x[1] + chrom2x[3] +  chrom2x[5] +  chrom2x[7]  +  chrom2x[9]  + chrom2x[11] )/ abs(chrom2x[12])<<"  to Qx'' "<<endl;
  cout<<"   Arcs  contribute:  "<<( chrom2y[1] + chrom2y[3] +  chrom2y[5] +  chrom2y[7]  +  chrom2y[9]  + chrom2y[11] )/ abs(chrom2y[12])<<"  to Qy'' "<<endl;
  cout<<"   Arcs  contribute:  "<<( chrom3x[1] + chrom3x[3] +  chrom3x[5] +  chrom3x[7]  +  chrom3x[9]  + chrom3x[11] )/ abs(chrom3x[12])<<"  to Qx''' "<<endl;
  cout<<"   Arcs  contribute:  "<<( chrom3y[1] + chrom3y[3] +  chrom3y[5] +  chrom3y[7]  +  chrom3y[9]  + chrom3y[11] )/ abs(chrom3y[12])<<"  to Qy''' "<<endl;  
  cout<<"----------------------------------------------------------------------------------------------------------"<<endl; 
 
  fout.close();
} 

void Cal_Coupling_Coefficient( Line & linename )
// coupling coefficient at the starting point: here solenoid's contribution is excluded
{
  int i;
  double creal=0.0, cimag=0.0;
  double angle, sin_angle, cos_angle;
  
  for(i=0; i< linename.Ncell;i++)
    if( linename.Cell[i]->TYPE == string("SKEWQ") ){
      creal= creal + sqrt(linename.Cell[i]->Beta1 * linename.Cell[i]->Beta2) *
	cos( (linename.Cell[i]->Mu1 - linename.Cell[i]->Mu2)* 2.0 * PI ) *
        linename.Cell[i]->GetP("K1SL") / 2. / PI  ;
      cimag= cimag + sqrt(linename.Cell[i]->Beta1 * linename.Cell[i]->Beta2) *
	sin( (linename.Cell[i]->Mu1 - linename.Cell[i]->Mu2)* 2.0 * PI ) *
        linename.Cell[i]->GetP("K1SL") / 2. / PI  ;
    }

   for(i=0; i< linename.Ncell;i++)
    if( linename.Cell[i]->TYPE == string("QUAD")  and linename.Cell[i]->DPSI !=0 ){
      creal= creal + sqrt(linename.Cell[i]->Beta1 * linename.Cell[i]->Beta2) *
	cos( (linename.Cell[i]->Mu1 - linename.Cell[i]->Mu2)* 2.0 * PI ) *
        linename.Cell[i]->GetP("K1L") * linename.Cell[i]->DPSI*2 / 2. / PI  ;
      cimag= cimag + sqrt(linename.Cell[i]->Beta1 * linename.Cell[i]->Beta2) *
	sin( (linename.Cell[i]->Mu1 - linename.Cell[i]->Mu2)* 2.0 * PI ) *
        linename.Cell[i]->GetP("K1L") * linename.Cell[i]->DPSI*2 / 2. / PI  ;
    }
 
  angle=  atan( abs(cimag/creal) );
  sin_angle=cimag/ sqrt(cimag*cimag+ creal*creal) ;
  cos_angle=creal/ sqrt(cimag*cimag+ creal*creal);
  
  if( sin_angle  >= 0. and  cos_angle >=0 ) angle=angle;
  if( sin_angle >= 0. and cos_angle <=0 )  angle=3.14159265 - angle;
  if( sin_angle < 0.  and cos_angle <=0)   angle=3.14159265 + angle;
  if( sin_angle < 0.  and cos_angle >=0)   angle=2*3.14159265-angle;

  cout<<" C_real  = "<<creal<<endl;  
  cout<<" C_imag  = "<<cimag<<endl;
  cout<<" C_amp   =  "<<sqrt(cimag*cimag+ creal*creal)<<endl;
  cout<<" C_phase =  "<< angle <<endl;
}

void Cal_Coupling_Coefficient_Source( Line & linename, const char * filename)
// single element's coupling contribution
{
  int i;
  double creal=0.0, cimag=0.0;
  double angle, sin_angle, cos_angle;
  fstream  f1;

  f1.open(filename, ios::out);
      f1<<setw(15)<<" NAME "<<setw(15)<<" TYPE "
	<<setw(15)<<" creal "<<setw(15)<<" cimag " <<setw(15)<<" C_amp "<<setw(15)<<" C_phase "<<endl;

  for(i=0; i< linename.Ncell;i++) {

    if( linename.Cell[i]->TYPE == string("SKEWQ") ){
      creal= sqrt(linename.Cell[i]->Beta1 * linename.Cell[i]->Beta2) *
	cos( (linename.Cell[i]->Mu1 - linename.Cell[i]->Mu2)* 2.0 * PI ) / 2. / PI  ;
      cimag= sqrt(linename.Cell[i]->Beta1 * linename.Cell[i]->Beta2) *
	sin( (linename.Cell[i]->Mu1 - linename.Cell[i]->Mu2)* 2.0 * PI ) / 2. / PI  ;
      
      angle=  atan( abs(cimag/creal) );
      sin_angle=cimag/ sqrt(cimag*cimag+ creal*creal) ;
      cos_angle=creal/ sqrt(cimag*cimag+ creal*creal);
      if( sin_angle  >= 0. and  cos_angle >=0 ) angle=angle;
      if( sin_angle >= 0. and cos_angle <=0 )  angle=3.14159265 - angle;
      if( sin_angle < 0.  and cos_angle <=0)   angle=3.14159265 + angle;
      if( sin_angle < 0.  and cos_angle >=0)   angle=2*3.14159265-angle;
      f1<<setw(15)<<linename.Cell[i]->NAME<<setw(15)<<setw(15)<<linename.Cell[i]->TYPE
	  <<setw(15)<<creal<<setw(15)<<cimag<<setw(15)<<sqrt(cimag*cimag+ creal*creal)<<setw(15)<<angle<<endl;
    }
    
    if( linename.Cell[i]->TYPE == string("QUAD")  and linename.Cell[i]->DPSI !=0 ){
      creal= sqrt(linename.Cell[i]->Beta1 * linename.Cell[i]->Beta2) *
	cos( (linename.Cell[i]->Mu1 - linename.Cell[i]->Mu2)* 2.0 * PI ) / 2. / PI  ;
      cimag= sqrt(linename.Cell[i]->Beta1 * linename.Cell[i]->Beta2) *
	sin( (linename.Cell[i]->Mu1 - linename.Cell[i]->Mu2)* 2.0 * PI ) / 2. / PI  ;
      
      angle=  atan( abs(cimag/creal) );
      sin_angle=cimag/ sqrt(cimag*cimag+ creal*creal) ;
      cos_angle=creal/ sqrt(cimag*cimag+ creal*creal);
      if( sin_angle  >= 0. and  cos_angle >=0 ) angle=angle;
      if( sin_angle >= 0. and cos_angle <=0 )  angle=3.14159265 - angle;
      if( sin_angle < 0.  and cos_angle <=0)   angle=3.14159265 + angle;
      if( sin_angle < 0.  and cos_angle >=0)   angle=2*3.14159265-angle;
      f1<<setw(15)<<linename.Cell[i]->NAME<<setw(15)<<setw(15)<<linename.Cell[i]->TYPE
	<<setw(15)<<creal<<setw(15)<<cimag<<setw(15)<<sqrt(cimag*cimag+ creal*creal)<<setw(15)<<angle<<endl;
    }

  }
   f1.close();
}

void  Cal_BeamSize_SigmaMatrx(Line & linename, const char* filename )
{
  int i,l;
  double emit1=2.5e-6/(25/0.9879), emit2=2.5e-6/(25/0.9879), emit3=0* 6e-2 * 6.8e-4;
  double a66[36], a66t[36], sigman[36], sigma[36], b66[36];
  fstream  f1;

  f1.open(filename, ios::out);

  for(i=0;i<36;i++) sigman[i]=0.;
  sigman[0*6+0]= emit1;
  sigman[1*6+1]= emit1;
  sigman[2*6+2]= emit2;
  sigman[3*6+3]= emit2;
  sigman[4*6+4]= emit3;
  sigman[5*6+5]= emit3;

   //---get sigma matrix

  for(l=0;l<linename.Ncell;l++){
    
    for(i=0;i<36;i++) a66[i] = linename.Cell[l]->A[i];
    mat_transpose(a66, a66t, 6,6);
    mat_mult(a66, sigman, b66, 6,6,6);
    mat_mult(b66, a66t, sigma, 6,6,6);

    f1<<linename.Cell[l]->S<<"  "<<sqrt( sigma[0*6+0] ) *1e6<<"    "<<sqrt( sigma[2*6+2] ) *1e6<<endl;

  }

  f1.close();

}


void  Twiss_Propagation_Drift(double  L,  double beta0[], double alfa0[], double beta1[], double alfa1[], double dphase[] )
{
  double gama0[2], gama1[2];

  gama0[0]= (1+ alfa0[0]*alfa0[0])/ beta0[0];
  gama0[1]= (1+ alfa0[1]*alfa0[1])/ beta0[1];  
  
  beta1[0] = beta0[0] -2*L*alfa0[0] + L*L * gama0[0];
  alfa1[0] = alfa0[0] -L * gama0[0];
  gama1[0] = gama0[0];

  beta1[1] = beta0[1] -2*L*alfa0[1] + L*L * gama0[1];
  alfa1[1] = alfa0[1] -L * gama0[1];
  gama1[1] = gama0[1];

  dphase[0]=atan( L / (beta0[0]*1.0 -alfa0[0]*L)  ) / 2. / PI;
  dphase[1]=atan( L / (beta0[1]*1.0 -alfa0[1]*L)  ) / 2. / PI;

}

void  Twiss_Propagation_Matrix(double  T[],  double beta0[], double alfa0[], double beta1[], double alfa1[], double dphase[] )
//  to be  corrected:  not suitable for 4-d case
{
  double gama0[2], gama1[2];

  gama0[0]= (1+ alfa0[0]*alfa0[0])/ beta0[0];
  gama0[1]= (1+ alfa0[1]*alfa0[1])/ beta0[1];  
  
  dphase[0]=atan( T[1,2] / (beta0[0]*T[1,1] -alfa0[0]*T[1,2])  );
  dphase[1]=atan( T[1,2] / (beta0[1]*T[1,1] -alfa0[1]*T[1,2])  );

}

void Coupling_From_Solenoid(double L, double KS, int nstep, double beta0[], double alfa0[], double phasex0, double phasey0, double & creal, double & cimag)
{
  int    i, k;
  double beta1[2], alfa1[2], dphase1[2];
  double betax, betay, alfax, alfay;
  double sqrtbxby, dphi, cosdphi, sindphi,real1, imag1;

  //---calculate  Twiss at  those steps
  
  creal=0;  cimag=0;
  
  for(i=0;i<nstep;i++){
    //cout<<"step: "<<i<<endl;
    Twiss_Propagation_Drift( L/nstep, beta0, alfa0, beta1, alfa1, dphase1);
    betax= (beta0[0] + beta1[0] )/2.0;
    betay= (beta0[1] + beta1[1] )/2.0;
    alfax= (alfa0[0] + alfa1[0] )/2.0;
    alfay= (alfa0[1] + alfa1[1] )/2.0;
    //cout<<beta1[0]<<"  "<<alfa1[0]<<"   "<<beta1[1]<<"  "<<alfa1[1]<<"  "<<dphase1[0]<<"  "<<dphase1[1]<<endl;
    sqrtbxby = sqrt( betax * betay );
    real1=(alfax/betax - alfay/betay);
    imag1=-(1.0/betax + 1./betay);
    dphi = ( phasex0 + dphase1[0]/2 - phasey0- dphase1[1]/2 ) * 2* PI; 
    cosdphi=cos( dphi );
    sindphi=sin( dphi );
    //cout<<sqrtbxby<<endl;
    //cout<<real1<<"  "<<imag1<<endl;
    //cout<<cosdphi<<"  "<<sindphi<<endl;
    creal =creal +  (real1 * cosdphi - imag1*sindphi )*sqrtbxby* (KS/2) * (L/nstep) / 2. / PI ;
    cimag =cimag +  (real1 * sindphi + imag1*cosdphi )*sqrtbxby* (KS/2) * (L/nstep) / 2. / PI ;
    for(k=0;k<2;k++){
      beta0[k] = beta1[k];
      alfa0[k] = alfa1[k];
    }
    phasex0 = phasex0 + dphase1[0];
    phasey0 = phasey0 + dphase1[1];
    //cout<<creal<<"  "<<cimag<<endl;
  }

  //cout<<" C- from solenoid: "<<creal<<"  + j "<<cimag<<endl;
  //cout<<" | C-| =  "<<sqrt( creal*creal +cimag*cimag)<<endl;
}

void Cal_Coupling_Coefficient_Updated( Line & linename, double & creal, double & cimag )
// coupling coefficient at the starting point, before the first element
{
  int i;
  double angle, sin_angle, cos_angle;
  
  double L, KS, beta0[2], alfa0[2], phasex0, phasey0, creal1, cimag1;
  int    nstep=20;
  
  for(i=0; i< linename.Ncell;i++)
    if( linename.Cell[i]->TYPE == string("SKEWQ") ){
      creal= creal + sqrt(linename.Cell[i]->Beta1 * linename.Cell[i]->Beta2) *
	cos( (linename.Cell[i]->Mu1 - linename.Cell[i]->Mu2)* 2.0 * PI ) *
        linename.Cell[i]->GetP("K1SL") / 2. / PI  ;
      cimag= cimag + sqrt(linename.Cell[i]->Beta1 * linename.Cell[i]->Beta2) *
	sin( (linename.Cell[i]->Mu1 - linename.Cell[i]->Mu2)* 2.0 * PI ) *
        linename.Cell[i]->GetP("K1SL") / 2. / PI  ;
    }

   for(i=0; i< linename.Ncell;i++)
    if( linename.Cell[i]->TYPE == string("QUAD")  and linename.Cell[i]->DPSI !=0 ){
      creal= creal + sqrt(linename.Cell[i]->Beta1 * linename.Cell[i]->Beta2) *
	cos( (linename.Cell[i]->Mu1 - linename.Cell[i]->Mu2)* 2.0 * PI ) *
        linename.Cell[i]->GetP("K1L") * (-1.0*linename.Cell[i]->DPSI)*2 / 2. / PI  ;
      cimag= cimag + sqrt(linename.Cell[i]->Beta1 * linename.Cell[i]->Beta2) *
	sin( (linename.Cell[i]->Mu1 - linename.Cell[i]->Mu2)* 2.0 * PI ) *
        linename.Cell[i]->GetP("K1L") * (-1.0*linename.Cell[i]->DPSI)*2 / 2. / PI  ;
    }
   
   for(i=0; i< linename.Ncell;i++)
   if( linename.Cell[i]->TYPE == string("SOLEN") ) {
     //cout<<" here is one solenoid: "<<i<<endl;
     L = linename.Cell[i]->L;
     KS = linename.Cell[i]->GetP("KS");
     //cout<<" ks=="<<L<<" "<< KS << endl;
     if( i == 0 ) {
       beta0[0]= linename.Cell[linename.Ncell-1]->Beta1; beta0[1]= linename.Cell[linename.Ncell-1]->Beta2;
       alfa0[0]= linename.Cell[linename.Ncell-1]->Alfa1; alfa0[1]= linename.Cell[linename.Ncell-1]->Alfa2;
       phasex0 = 0;   phasey0 = 0;
     }
     else{
       beta0[0]= linename.Cell[i-1]->Beta1; beta0[1]= linename.Cell[i-1]->Beta2;
       alfa0[0]= linename.Cell[i-1]->Alfa1; alfa0[1]= linename.Cell[i-1]->Alfa2;
       phasex0 = linename.Cell[i-1]->Mu1;   phasey0 = linename.Cell[i-1]->Mu2;
     }
     //cout<<i<<"   "<<beta0[0]<<"  "<<beta0[1]<<"   "<<alfa0[0]<<"  "<<alfa0[1]<<endl;
     Coupling_From_Solenoid(L, KS, nstep, beta0, alfa0, phasex0, phasey0, creal1, cimag1) ;
     //cout<<i<<"   "<<beta0[0]<<"  "<<beta0[1]<<"   "<<alfa0[0]<<"  "<<alfa0[1]<<endl;
     creal= creal + creal1;
     cimag= cimag + cimag1;
     //cout<<creal1<<"  "<<cimag1<<"  "<<atan2(cimag, creal)<<endl;
   }
   
   //cout<<"...from analytical calculation:"<<endl;
   //cout<<" C^{-}  = "<<creal<<" + J "<<cimag<<endl;
   //cout<<" C^{-}  = "<<sqrt(cimag*cimag+ creal*creal)<<" / "<<atan2(cimag, creal)<<endl;
}

void Cal_Coupling_One_Location(Line linename, int  i, double & creal, double & cimag )
//   calaculate  global coupling at one element with index  i, reference  point is just at entrance of element i
{
  int     j;
  double  dphasex, dphasey, sqrtbxby_ksl, L,KS, beta0[2], alfa0[2], phasex0, phasey0;
  int     nstep=10;
  double  creal1=0, cimag1=0;
  
  creal=0. ;  cimag=0.;

  if( i ==0 ) {     //  for default starting point
    Cal_Coupling_Coefficient_Updated(linename, creal, cimag );
  }
  else{           // for i != 0
  
    //----count elements after observation point
    for( j=i;j< linename.Ncell; j++){
      if( linename.Cell[j]->TYPE == string("SKEWQ") ){
	dphasex = (linename.Cell[j]->Mu1-linename.Cell[i-1]->Mu1)*2*PI;
	dphasey = (linename.Cell[j]->Mu2-linename.Cell[i-1]->Mu2)*2*PI;
	sqrtbxby_ksl =  sqrt(linename.Cell[j]->Beta1 * linename.Cell[j]->Beta2) * linename.Cell[j]->GetP("K1SL");
	creal= creal + sqrtbxby_ksl*cos(dphasex - dphasey) / 2. / PI  ;
	cimag= cimag + sqrtbxby_ksl*sin(dphasex - dphasey) / 2. / PI  ;
      }
      if( linename.Cell[j]->TYPE == string("QUAD") and  linename.Cell[j]->DPSI !=0.  ) {
	dphasex = (linename.Cell[j]->Mu1-linename.Cell[i-1]->Mu1)*2*PI;
	dphasey = (linename.Cell[j]->Mu2-linename.Cell[i-1]->Mu2)*2*PI;
	sqrtbxby_ksl =  sqrt(linename.Cell[j]->Beta1 * linename.Cell[j]->Beta2) * linename.Cell[j]->GetP("K1L") * (-1.0* linename.Cell[j]->DPSI)*2;
	creal= creal + sqrtbxby_ksl*cos( dphasex - dphasey) / 2. / PI  ;
	cimag= cimag + sqrtbxby_ksl*sin( dphasex - dphasey) / 2. / PI  ;
      }
      if( linename.Cell[j]->TYPE == string("SOLEN") ) {
	//cout<<"here is a solenoid: "<<linename.Cell[j]->NAME<<endl;
	L = linename.Cell[j]->L;
	KS = linename.Cell[j]->GetP("KS");
	//cout<<j-1<<"  "<<linename.Cell[j-1]->NAME<<endl;
	//cout<<"HI: "<<linename.Cell[j-1]->Beta1<<"  "<<linename.Cell[j-1]->Beta2<<endl;
	//cout<<"HI: "<<linename.Cell[j-1]->Alfa1<<"  "<<linename.Cell[j-1]->Alfa2<<endl;
	//cout<<"HI: "<<linename.Cell[j-1]->Mu1<<"  "<<linename.Cell[j-1]->Mu2<<endl;
	beta0[0]= linename.Cell[j-1]->Beta1; beta0[1]= linename.Cell[j-1]->Beta2;
	alfa0[0]= linename.Cell[j-1]->Alfa1; alfa0[1]= linename.Cell[j-1]->Alfa2;
	phasex0 = linename.Cell[j-1]->Mu1 - linename.Cell[i-1]->Mu1;
	phasey0 = linename.Cell[j-1]->Mu2 - linename.Cell[i-1]->Mu2;
	Coupling_From_Solenoid(L, KS, nstep, beta0, alfa0, phasex0, phasey0, creal1, cimag1) ;
	creal= creal + creal1;
	cimag= cimag + cimag1;
	//cout<<creal1<<"  "<<cimag1<<"  "<<creal<<"  "<<cimag<<endl;
      }
    
    }
  
    //------count elements before observation point
    for( j=0;j<i;j++){   // here i can't be 0
      if( linename.Cell[j]->TYPE == string("SKEWQ") ){
	//cout<<"...there is an skewq."<<endl;
	dphasex = linename.Tune1 -  abs(linename.Cell[j]->Mu1-linename.Cell[i-1]->Mu1);
	dphasey = linename.Tune2 -  abs(linename.Cell[j]->Mu2-linename.Cell[i-1]->Mu2);
	sqrtbxby_ksl =  sqrt(linename.Cell[j]->Beta1 * linename.Cell[j]->Beta2) * linename.Cell[j]->GetP("K1SL");
	creal= creal + sqrtbxby_ksl*cos( (dphasex - dphasey)* 2.0 * PI ) / 2. / PI  ;
	cimag= cimag + sqrtbxby_ksl*sin( (dphasex - dphasey)* 2.0 * PI ) / 2. / PI  ;
      }
      if( linename.Cell[j]->TYPE == string("QUAD") and  linename.Cell[j]->DPSI !=0.  ) {
	//cout<<"...there is a quad with roll"<<endl;
	dphasex = linename.Tune1 -  abs(linename.Cell[j]->Mu1-linename.Cell[i-1]->Mu1);
	dphasey = linename.Tune2 -  abs(linename.Cell[j]->Mu2-linename.Cell[i-1]->Mu2);
	sqrtbxby_ksl =  sqrt(linename.Cell[j]->Beta1 * linename.Cell[j]->Beta2) * linename.Cell[j]->GetP("K1L") * (-1.0*linename.Cell[j]->DPSI)*2;
	creal= creal + sqrtbxby_ksl*cos( (dphasex - dphasey)* 2.0 * PI ) / 2. / PI  ;
	cimag= cimag + sqrtbxby_ksl*sin( (dphasex - dphasey)* 2.0 * PI ) / 2. / PI  ;
      }
      if( linename.Cell[j]->TYPE == string("SOLEN")  ) {
	//cout<<"here is a solenoid: "<<linename.Cell[j]->NAME<<"  "<<j<<endl;
	L = linename.Cell[j]->L;
	KS = linename.Cell[j]->GetP("KS");
	//cout<<j-1<<"  "<<linename.Cell[j-1]->NAME<<endl;
	//cout<<"HI: "<<linename.Cell[j-1]->Beta1<<"  "<<linename.Cell[j-1]->Beta2<<endl;
	//cout<<"HI: "<<linename.Cell[j-1]->Alfa1<<"  "<<linename.Cell[j-1]->Alfa2<<endl;
	//cout<<"HI: "<<linename.Cell[j-1]->Mu1<<"  "<<linename.Cell[j-1]->Mu2<<endl;
	if(j==0) {
	  beta0[0]= linename.Cell[linename.Ncell-1]->Beta1; beta0[1]= linename.Cell[linename.Ncell-1]->Beta2;
	  alfa0[0]= linename.Cell[linename.Ncell-1]->Alfa1; alfa0[1]= linename.Cell[linename.Ncell-1]->Alfa2;
	  phasex0 = linename.Tune1 - linename.Cell[i-1]->Mu1;
	  phasey0 = linename.Tune2 - linename.Cell[i-1]->Mu2;
	}
	else{
	  beta0[0]= linename.Cell[j-1]->Beta1; beta0[1]= linename.Cell[j-1]->Beta2;
	  alfa0[0]= linename.Cell[j-1]->Alfa1; alfa0[1]= linename.Cell[j-1]->Alfa2;
	  phasex0 = linename.Tune1 -  abs(linename.Cell[i-1]->Mu1-linename.Cell[j-1]->Mu1);
	  phasey0 = linename.Tune2 -  abs(linename.Cell[i-1]->Mu2-linename.Cell[j-1]->Mu2);
	}
	Coupling_From_Solenoid(L, KS, nstep, beta0, alfa0, phasex0, phasey0, creal1, cimag1) ;
	creal= creal + creal1;
	cimag= cimag + cimag1;
	//cout<<creal1<<"  "<<cimag1<<"  "<<creal<<"  "<<cimag<<endl;
      }
    }

  }
  
}

void Cal_Coupling_Along_Ring(Line linename, const char* filename )
{
  int  i,j;
  double creal=0, cimag=0;
  fstream  f1;

  f1.open(filename, ios::out);
  for(i=0;i< linename.Ncell; i++){
    Cal_Coupling_One_Location(linename, i, creal, cimag );
    f1<<linename.Cell[i]->NAME<<"  "<<linename.Cell[i]->S<<"   "<<creal<<"  "<<cimag<<"  "<<sqrt(creal*creal + cimag*cimag)<<"  "<<atan2(cimag, creal)<<endl;
  }
  f1.close();

}

void Cal_Coupling_Along_Ring_Old(Line linename, const char* filename )
{
  int  i,j;
  double dphasex, dphasey, sqrtbxby_ksl;
  double creal=0, cimag=0;

  double L,KS, beta0[2], alfa0[2], phasex0, phasey0, creal1, cimag1;
  int    nstep=1;

  fstream  f1;

  f1.open(filename, ios::out);
  
  for(i=0;i< linename.Ncell; i++){
    
    //if( linename.Cell[i]->TYPE == string("DRIFT") ) {   //  calculate c- at all sextupoles-)

      creal=0. ;  cimag=0.;

      //----count elements after observation point
      for( j=i+1;j< linename.Ncell; j++){
	if( linename.Cell[j]->TYPE == string("SKEWQ") ){
	  cout<<"...there is a skewQ"<<endl;
	  dphasex = (linename.Cell[j]->Mu1-linename.Cell[i]->Mu1)*2*PI;
	  dphasey = (linename.Cell[j]->Mu2-linename.Cell[i]->Mu2)*2*PI;
	  sqrtbxby_ksl =  sqrt(linename.Cell[j]->Beta1 * linename.Cell[j]->Beta2) * linename.Cell[j]->GetP("K1SL");
	  creal= creal + sqrtbxby_ksl*cos(dphasex - dphasey) / 2. / PI  ;
	  cimag= cimag + sqrtbxby_ksl*sin(dphasex - dphasey) / 2. / PI  ;
	}
	if( linename.Cell[j]->TYPE == string("QUAD") and  linename.Cell[j]->DPSI !=0.  ) {
	  cout<<"...there is a quad with roll "<<endl;
	  dphasex = (linename.Cell[j]->Mu1-linename.Cell[i]->Mu1)*2*PI;
	  dphasey = (linename.Cell[j]->Mu2-linename.Cell[i]->Mu2)*2*PI;
	  sqrtbxby_ksl =  sqrt(linename.Cell[j]->Beta1 * linename.Cell[j]->Beta2) * linename.Cell[j]->GetP("K1L") * linename.Cell[j]->DPSI*2;
	  creal= creal + sqrtbxby_ksl*cos( dphasex - dphasey) / 2. / PI  ;
	  cimag= cimag + sqrtbxby_ksl*sin( dphasex - dphasey) / 2. / PI  ;
	}
	if( linename.Cell[j]->TYPE == string("SOLEN") ) {
	  cout<<"here is a solenoid: "<<linename.Cell[j]->NAME<<endl;
	  L = linename.Cell[j]->L;
	  KS = linename.Cell[j]->GetP("KS");
	  cout<<j-1<<"  "<<linename.Cell[j-1]->NAME<<endl;
	  cout<<"HI: "<<linename.Cell[j-1]->Beta1<<"  "<<linename.Cell[j-1]->Beta2<<endl;
	  cout<<"HI: "<<linename.Cell[j-1]->Alfa1<<"  "<<linename.Cell[j-1]->Alfa2<<endl;
	  cout<<"HI: "<<linename.Cell[j-1]->Mu1<<"  "<<linename.Cell[j-1]->Mu2<<endl;
          beta0[0]= linename.Cell[j-1]->Beta1; beta0[1]= linename.Cell[j-1]->Beta2;
          alfa0[0]= linename.Cell[j-1]->Alfa1; alfa0[1]= linename.Cell[j-1]->Alfa2;
	  phasex0 = linename.Cell[j-1]->Mu1 - linename.Cell[i]->Mu1;
	  phasey0 = linename.Cell[j-1]->Mu2 - linename.Cell[i]->Mu2;
	  cout<<"HI: "<<phasex0<<"  "<<phasey0<<endl;
	  cout<<L<<"  "<<KS<<"  "<<beta0[0]<<"  "<<beta0[1]<<"  "<<alfa0[0]<<"  "<<alfa0[1]<<"  "<<phasex0<<"  "<<phasey0<<endl;
	  Coupling_From_Solenoid(L, KS, nstep, beta0, alfa0, phasex0, phasey0, creal1, cimag1) ;
	  creal= creal + creal1;
	  cimag= cimag + cimag1;
	  cout<<creal1<<"  "<<cimag1<<"  "<<creal<<"  "<<cimag<<endl;
	}
	
      }

      //------count elements before observation point
      for( j=0;j<i;j++){
	if( linename.Cell[j]->TYPE == string("SKEWQ") ){
	  cout<<"...there is an skewq."<<endl;
	  dphasex = linename.Tune1 -  abs(linename.Cell[j]->Mu1-linename.Cell[i]->Mu1);
	  dphasey = linename.Tune2 -  abs(linename.Cell[j]->Mu2-linename.Cell[i]->Mu2);
	  sqrtbxby_ksl =  sqrt(linename.Cell[j]->Beta1 * linename.Cell[j]->Beta2) * linename.Cell[j]->GetP("K1SL");
	  creal= creal + sqrtbxby_ksl*cos( (dphasex - dphasey)* 2.0 * PI ) / 2. / PI  ;
	  cimag= cimag + sqrtbxby_ksl*sin( (dphasex - dphasey)* 2.0 * PI ) / 2. / PI  ;
	}
	if( linename.Cell[j]->TYPE == string("QUAD") and  linename.Cell[i]->DPSI !=0.  ) {
	  cout<<"...there is a quad with roll"<<endl;
	  dphasex = linename.Tune1 -  abs(linename.Cell[j]->Mu1-linename.Cell[i]->Mu1);
	  dphasey = linename.Tune2 -  abs(linename.Cell[j]->Mu2-linename.Cell[i]->Mu2);
	  sqrtbxby_ksl =  sqrt(linename.Cell[j]->Beta1 * linename.Cell[j]->Beta2) * linename.Cell[j]->GetP("K1L") * linename.Cell[j]->DPSI*2;
	  creal= creal + sqrtbxby_ksl*cos( (dphasex - dphasey)* 2.0 * PI ) / 2. / PI  ;
	  cimag= cimag + sqrtbxby_ksl*sin( (dphasex - dphasey)* 2.0 * PI ) / 2. / PI  ;
	}
	if( linename.Cell[j]->TYPE == string("SOLEN") ) {
	  cout<<"here is a solenoid: "<<linename.Cell[j]->NAME<<endl;
	  L = linename.Cell[j]->L;
	  KS = linename.Cell[j]->GetP("KS");
	  cout<<j-1<<"  "<<linename.Cell[j-1]->NAME<<endl;
	  cout<<"HI: "<<linename.Cell[j-1]->Beta1<<"  "<<linename.Cell[j-1]->Beta2<<endl;
	  cout<<"HI: "<<linename.Cell[j-1]->Alfa1<<"  "<<linename.Cell[j-1]->Alfa2<<endl;
	  cout<<"HI: "<<linename.Cell[j-1]->Mu1<<"  "<<linename.Cell[j-1]->Mu2<<endl;
          beta0[0]= linename.Cell[j-1]->Beta1; beta0[1]= linename.Cell[j-1]->Beta2;
          alfa0[0]= linename.Cell[j-1]->Alfa1; alfa0[1]= linename.Cell[j-1]->Alfa2;
	  phasex0 = linename.Tune1 -  (linename.Cell[i]->Mu1-linename.Cell[j-1]->Mu1);
	  phasey0 = linename.Tune2 -  (linename.Cell[i]->Mu2-linename.Cell[j-1]->Mu2);
	  cout<<"HI: "<<phasex0<<"  "<<phasey0<<endl;
	  cout<<L<<"  "<<KS<<"  "<<beta0[0]<<"  "<<beta0[1]<<"  "<<alfa0[0]<<"  "<<alfa0[1]<<"  "<<phasex0<<"  "<<phasey0<<endl;
	  Coupling_From_Solenoid(L, KS, nstep, beta0, alfa0, phasex0, phasey0, creal1, cimag1) ;
	  creal= creal + creal1;
	  cimag= cimag + cimag1;
	  cout<<creal1<<"  "<<cimag1<<"  "<<creal<<"  "<<cimag<<endl;
	}
      }
       
      f1<<linename.Cell[i]->NAME<<"  "<<linename.Cell[i]->S<<"   "<<creal<<"  "<<cimag<<"  "<<sqrt(creal*creal + cimag*cimag)<<endl;

      
      //    }
       }

  f1.close();

}

void Harmonic_Analysis(double xtemp[], int nturn, double tune, double & peak, double &phase )
{
  int  i;
  double sumr=0, sumi=0;
  
  for(i=0;i<nturn;i++){
    sumr = sumr + xtemp[i] * cos(2*PI*i*tune)  ;
    sumi = sumi + xtemp[i] * sin(2*PI*i*tune)  ;
  } 
  peak  = sqrt(sumr*sumr + sumi*sumi);
  phase = atan2(-sumi, sumr);
  
}

void Meas_Coupling_OneDualPlaneTBT(int nturn, int m, double xtemp[], double ytemp[], double &camp, double & cphase )
{
  int i;
  double  tune1, tune2;
  double  peak1x,  peak1y,   peak2x,  peak2y;
  double  phase1x, phase1y,  phase2x, phase2y;
  double  r1, r2, dphase1, dphase2;
  
  FineTuneFinder(1024, xtemp, 0.69, 0.75, tune1); 
  FineTuneFinder(1024, ytemp, 0.65, 0.69, tune2);

  Harmonic_Analysis(xtemp, nturn, tune1, peak1x, phase1x);
  Harmonic_Analysis(xtemp, nturn, tune2, peak2x, phase2x);
  Harmonic_Analysis(ytemp, nturn, tune1, peak1y, phase1y);
  Harmonic_Analysis(ytemp, nturn, tune2, peak2y, phase2y);
  
  r1 =  peak1y /  peak1x ;
  r2 =  peak2x /  peak2y ;
  dphase1 = phase1y - phase1x ;  
  dphase2 = phase2x - phase2y ;  
  camp =   2 *sqrt( r1*r2 ) * abs(tune1-tune2) / (1 + r1 * r2) ;
  cphase = dphase1 ;
}

void Meas_Coupling_Along_Ring(Line linename, int nbpm, int nturn, int m, const  char* filename)  // based on TBT dual BPM
{
  int    i,j;
  int    ncount ;
  double x[6];
  double *xtbt_holder = new double [nbpm*nturn];
  double *ytbt_holder = new double [nbpm*nturn];
  double xtemp[nturn], ytemp[nturn];
  double camp, cphase;
  fstream f1, f2;

  cout<<"I  am here ! "<<endl;
  
  //----generate the TBT  data from tracking
  
  for(i=0;i<6;i++)  x[i]=0;
  x[1]=1e-4; x[3]=1e-4;

  cout<<"I  am here ! "<<endl;
  
  ncount=0;
  for(GP.turn=0; GP.turn < nturn; GP.turn++) {
    for(j=0;j<linename.Ncell;j++) {
      linename.Cell[j]->Pass(x);
      xtbt_holder[ncount]= x[0];
      ytbt_holder[ncount]= x[2];
      ncount++;
    }
  }

  cout<<" ...I  am here!!"<<endl;
  
  //----save TBT data
  f1.open("./all_tbt_bpm.dat", ios::out);
  for(i=0;i<nbpm*nturn;i++){
    f1<<xtbt_holder[i]<<"   "<<ytbt_holder[i]<<endl;
  }
  f1.close();
  
  //----analyze TBT data at each BPM

  f2.open(filename, ios::out);
  for(i=0;i<nbpm; i++) {
    for(j=0;j<nturn;j++) xtemp[j] = xtbt_holder[j*nbpm+ i];
    for(j=0;j<nturn;j++) ytemp[j] = ytbt_holder[j*nbpm+ i];    
    Meas_Coupling_OneDualPlaneTBT( 1024, 10, xtemp, ytemp, camp,  cphase );
    f2<<linename.Cell[i]->S<<"   "<<camp<<"   "<<cphase<<endl;
  }
  f2.close();
  
  delete [] xtbt_holder;
  delete [] ytbt_holder;   
  
}

void Meas_Coupling_Along_Ring_Old(Line linename , int nbpm, int nturn)  // based on TBT dual BPM
{
  int  i,j, ncount;
  double x[6];
  double xtbt_holder[nbpm*nturn], ytbt_holder[nbpm*nturn];
  double xtemp[nturn], ytemp[nturn];
  double camp, cphase;
  
  //----generate the TBT  data from tracking
  
  for(i=0;i<6;i++)  x[i]=0;
  x[1]=1e-4; x[3]=1e-4;
  
  ncount=0;
  for(GP.turn=0; GP.turn < nturn; GP.turn++) {
    for(j=0;j<linename.Ncell;j++) {
      linename.Cell[j]->Pass(x);
      if( linename.Cell[j]->TYPE == string("BPM") ){
        xtbt_holder[ncount]= x[0];
	ytbt_holder[ncount]= x[2];
	ncount++;
      }
    }
  }
  
  //----analyze TBT data at each BPM
  
  for(i=0;i<nbpm; i++) {
    for(j=0;j<nturn;j++) xtemp[j] = xtbt_holder[j*nbpm+ i];
    for(j=0;j<nturn;j++) ytemp[j] = ytbt_holder[j*nbpm+ i];    
    Meas_Coupling_OneDualPlaneTBT( 1024, 10, xtemp, ytemp, camp,  cphase );
    cout<<linename.Cell[i]->S<<"   "<<camp<<"   "<<cphase<<endl;
  }
  
}
  
void SaveData_To_File(double xtemp[], int nturn, const char*  filename )
{
  int i;
  fstream fout;
  
  fout.open(filename, ios::out);
  for(i=0;i<nturn;i++){
    fout<<i<<"  "<<xtemp[i]<<endl;
  }
  fout.close();
}

void FFT_Save_File(double xtemp[], int nturn, int m,  const char*  filename )
{
  int i;
  int nfft;
  fstream fout;
  
  nfft = 1;
  for (i=0;i<m;i++) 
    nfft *= 2;
  
  double  xreal[nfft], ximag[nfft];
  for(i=0;i<nfft;i++){
    xreal[i]= xtemp[i];
    ximag[i]= 0. ;
  }
  fft(m, xreal, ximag);
  
  fout.open(filename, ios::out);
  for(i=0;i<nfft;i++){
    fout<<i<<"  "<<i*1.0/nfft<<"  "<<xreal[i]<<"  "<<ximag[i]<<"   "<<sqrt(xreal[i]*xreal[i]+ ximag[i]*ximag[i] )<<endl;
  }
  fout.close();
  
}

template <class T> void SOLEN_Pass_Old1(T x[], double L, int Nint, double KS)
{
  int i, j, k;	  
  double g=KS/2, theta=KS*L/2, costheta= cos( theta ), sintheta=  sin( theta ) ;
  double cctheta = costheta * costheta, sstheta=sintheta * sintheta, cstheta=sintheta * costheta;
  double M[36];
  T      xtemp[6];
  
  if( KS == 0. ) {
    DRIFT_Pass(x, L);}
  else{
    //----calculate the matrix
    for(i=0;i<36;i++) M[i]= 0.0;

    
    M[0*6+0] = costheta * costheta;
    M[0*6+1] = sintheta * costheta / g; 
    M[0*6+2] = sintheta * costheta;
    M[0*6+3] = sintheta *  sintheta / g;
    
    M[1*6+0] = -g * sintheta * costheta ;
    M[1*6+1] =      costheta * costheta ;
    M[1*6+2] = -g*  sintheta * sintheta ;
    M[1*6+3] =      sintheta *  costheta;
    
    M[2*6+0] = -sintheta * costheta     ;
    M[2*6+1] = -sintheta * sintheta / g ;
    M[2*6+2] =  costheta * costheta     ;
    M[2*6+3] =  sintheta * costheta / g ; 
    
    M[3*6+0] = g *  sintheta * sintheta ;
    M[3*6+1] =     -sintheta * costheta ;
    M[3*6+2] = -g * sintheta * costheta ;
    M[3*6+3] =      costheta * costheta ;
    
    M[4*6+4] = 1.0;
    M[5*6+5] = 1.0;
    
    /*---
      M[0*6+0] = cctheta;
      M[0*6+1] = cstheta / g; 
      M[0*6+2] = cstheta;
      M[0*6+3] = sstheta / g;
    
      M[1*6+0] = -g * cstheta ;
      M[1*6+1] =      cctheta ;
      M[1*6+2] = -g * sstheta ;
      M[1*6+3] =      cstheta;
    
      M[2*6+0] = -cstheta     ;
      M[2*6+1] = -sstheta / g ;
      M[2*6+2] =  cctheta     ;
      M[2*6+3] =  cstheta / g ; 
    
      M[3*6+0] =  g * sstheta ;
      M[3*6+1] =     -cstheta ;
      M[3*6+2] = -g * cstheta ;
      M[3*6+3] =      cctheta ;
    
      M[4*6+4] = 1.0;
      M[5*6+5] = 1.0;
      ---*/

    //---transfer
    x[1]=x[1]/(1+x[5]);  // x'=px / (1+delta), is it needed ?
    x[3]=x[3]/(1+x[5]);
    for(i=0;i<6;i++) xtemp[i]=x[i];
    
    for(k=0;k<6;k++) {
      x[k]=0.;
      for(j=0;j<6;j++) x[k]= x[k]+ M[k*6+j]*xtemp[j];
    }
    x[1]=x[1]*(1+x[5]); // px=x' * (1+delta), is it needed ?
    x[3]=x[3]*(1+x[5]);
    
  }
  
}


//-----------------supporting  file---------------------------

void Extract_Orbits(Line linename, double holder[] )
{
  int i,j;
  
  for(i=0;i<linename.Ncell;i++){
    for(j=0;j<6;j++){
      holder[i*6+j ] = linename.Cell[i]->X[j];
    }
  }
   
}


void Make_Orbit_Difference(Line linename, double orbit0[], double orbit1[], double diff[] )
{
  int i,j;

  for(i=0;i<linename.Ncell*6;i++){
    diff[i] = orbit1[i]  - orbit0[i];
  }

}

void Save_Orbit_Difference(Line linename, double orbit0[], double orbit1[], const char* filename)
{
  int i,j;
  fstream f1;
  
  f1.open(filename, ios::out);
  for(i=0;i<linename.Ncell;i++){
    f1<<setw(20)<<linename.Cell[i]->NAME<<setw(20)<<linename.Cell[i]->S;
    for(j=0;j<6;j++) f1<<setw(20)<<orbit1[6*i+j] - orbit0[6*i+j];
    f1<<endl;
  }
  f1.close();
  
}

void Cal_Sext_RDTs( Line & linename )
// sextupole linear RDT at the staring point
{
  int i;
  double h2100_c, h2100_s;
  double h3000_c, h3000_s;
  double h1011_c, h1011_s;  
  double h1002_c, h1002_s;
  double h1020_c, h1020_s;
  double amp, phi, k2l, betx, bety, mux, muy;

  h2100_c=0 ; h2100_s =0 ;   h3000_c=0 ; h3000_s =0 ;    h1011_c=0 ; h1011_s =0 ;   
  h1002_c=0 ; h1002_s =0 ;   h1020_c=0 ; h1020_s =0 ;

  for(i=0;i<linename.Ncell;i++) {
    if( linename.Cell[i]->TYPE==string("SEXT") ||  linename.Cell[i]->TYPE==string("MULT") ){
      k2l=linename.Cell[i]->GetP("K2L");
      betx=linename.Cell[i]->Beta1;
      bety=linename.Cell[i]->Beta2;
      mux =linename.Cell[i]->Mu1 * 2.0 * PI;
      muy =linename.Cell[i]->Mu2 * 2.0 * PI;
      
      amp = k2l * pow(betx,1.5);      phi = mux;      h2100_c += amp * cos(phi);      h2100_s += amp * sin(phi);
      amp = k2l * pow(betx,1.5);      phi = 3* mux;   h3000_c += amp * cos(phi);      h3000_s += amp * sin(phi);
      amp = k2l * pow(betx,0.5) * bety;  phi = mux;   h1011_c += amp * cos(phi);      h1011_s += amp * sin(phi);
      amp = k2l * pow(betx,0.5) * bety;  phi = mux - 2* muy;  h1002_c += amp * cos(phi);      h1002_s += amp * sin(phi);
      amp = k2l * pow(betx,0.5) * bety;  phi = mux + 2* muy;  h1020_c += amp * cos(phi);      h1020_s += amp * sin(phi);
     }
  }

  h2100_c=   h2100_c * (-1.0/8);  h2100_s=   h2100_s * (-1.0/8);
  h3000_c=   h3000_c * (-1.0/24); h3000_s=   h3000_s * (-1.0/24);
  h1011_c=   h1011_c * ( 1.0/ 4); h1011_s=   h1011_s * ( 1.0/ 4);
  h1002_c=   h1002_c * ( 1.0/ 8); h1002_s=   h1002_s * ( 1.0/ 8);
  h1020_c=   h1020_c * ( 1.0/ 8); h1020_s=   h1020_s * ( 1.0/ 8);

  cout<<">>>1st order sextupole geometric resonance driving terms: "<<endl;
  cout<<" h2100 :  "<<h2100_c <<" +i  "<<h2100_s<<" , "<<sqrt(h2100_c*h2100_c + h2100_s*h2100_s)<<" / "<<atan2(h2100_s,h2100_c)<<endl;
  cout<<" h3000 :  "<<h3000_c <<" +i  "<<h3000_s<<" , "<<sqrt(h3000_c*h3000_c + h3000_s*h3000_s)<<" / "<<atan2(h3000_s,h3000_c)<<endl;
  cout<<" h1011 :  "<<h1011_c <<" +i  "<<h1011_s<<" , "<<sqrt(h1011_c*h1011_c + h1011_s*h1011_s)<<" / "<<atan2(h1011_s,h1011_c)<<endl;
  cout<<" h1002 :  "<<h1002_c <<" +i  "<<h1002_s<<" , "<<sqrt(h1002_c*h1002_c + h1002_s*h1002_s)<<" / "<<atan2(h1002_s,h1002_c)<<endl;
  cout<<" h1020 :  "<<h1020_c <<" +i  "<<h1020_s<<" , "<<sqrt(h1020_c*h1020_c + h1020_s*h1020_s)<<" / "<<atan2(h1020_s,h1020_c)<<endl;
}


void Cal_3Qx_RDTs_Ring( Line & linename, const char* filename )
// 3Qx RDT along the ring
{
  int i,j;
  double h3000_c, h3000_s;
  double k2l, betx, dmux;
  fstream f1;
  
  f1.open(filename, ios::out);
  
  for(i=0;i<linename.Ncell;i++){
    h3000_c=0;
    h3000_s=0;
    for(j=0;j<linename.Ncell;j++){
      if( linename.Cell[j]->TYPE==string("MULT") || linename.Cell[j]->TYPE==string("SEXT") ){
	k2l=linename.Cell[j]->GetP("K2L");
	betx=linename.Cell[j]->Beta1;
	if(j>i){
	  dmux =linename.Cell[j]->Mu1 * 2.0 * PI  - linename.Cell[i]->Mu1 * 2.0 * PI;  
	  h3000_c +=  k2l * pow(betx,1.5) * cos( dmux * 3 );
	  h3000_s +=  k2l * pow(betx,1.5) * sin( dmux * 3 );
	}
	else{
	  dmux =linename.Tune1 * 2 *PI  - abs(linename.Cell[i]->Mu1 * 2.0 * PI  - linename.Cell[j]->Mu1 * 2.0 * PI);  
	  h3000_c +=  k2l * pow(betx,1.5) * cos( dmux * 3 );
	  h3000_s +=  k2l * pow(betx,1.5) * sin( dmux * 3 );
	}
      }
    }
    h3000_c=   h3000_c * (-1.0/24); h3000_s=   h3000_s * (-1.0/24);
    f1<<linename.Cell[i]->S<<" "<<  h3000_c<< "  "<<  h3000_s<<"  "<<sqrt( h3000_c* h3000_c +h3000_s*h3000_s)<<"  "<<atan2(h3000_s, h3000_c)<<endl;
  }
  f1.close();
}

void Cal_3Qx_RDTs_Source( Line & linename )
// 3Qx RDT contribution along the ring
{
  int j;
  double h3000_c, h3000_s;
  double k2l, betx, dmux;
  
  for(j=0;j<linename.Ncell;j++){
    if( linename.Cell[j]->TYPE==string("MULT") || linename.Cell[j]->TYPE==string("SEXT") ){
      k2l=linename.Cell[j]->GetP("K2L");
      betx=linename.Cell[j]->Beta1;
      dmux =linename.Cell[j]->Mu1 * 2.0 * PI;
      h3000_c =  k2l * pow(betx,1.5) * cos( dmux * 3 );
      h3000_s =  k2l * pow(betx,1.5) * sin( dmux * 3 );
      h3000_c =  h3000_c * (-1.0/24); 
      h3000_s =  h3000_s * (-1.0/24);
      cout<<linename.Cell[j]->S<<" "<<linename.Cell[j]->NAME<<"  "<<linename.Cell[j]->TYPE<<"  "<<h3000_c<< "  "<<  h3000_s<<"  "<<sqrt( h3000_c* h3000_c +h3000_s*h3000_s)<<"  "<<atan2(h3000_s, h3000_c)<<endl;
    }
  }
}

void Cal_SkewSext_RDTs( Line & linename ) 
// skew sextupole linear RDT at the staring point
{
  int i;
  double h0030_c, h0030_s;
  double h0021_c, h0021_s;
  double h1110_c, h1110_s;  
  double h2001_c, h2001_s;
  double h2010_c, h2010_s;
  double amp, phi, k2sl, betx, bety, mux, muy;

  h0030_c=0 ; h0030_s =0 ;   h0021_c=0 ; h0021_s =0 ;    h1110_c=0 ; h1110_s =0 ;   
  h2001_c=0 ; h2001_s =0 ;   h2010_c=0 ; h2010_s =0 ;

  for(i=0;i<linename.Ncell;i++) {
    if( linename.Cell[i]->TYPE==string("SEXT") ||  linename.Cell[i]->TYPE==string("MULT") ){
      k2sl=linename.Cell[i]->GetP("K2SL");
      betx=linename.Cell[i]->Beta1;
      bety=linename.Cell[i]->Beta2;
      mux =linename.Cell[i]->Mu1 * 2.0 * PI;
      muy =linename.Cell[i]->Mu2 * 2.0 * PI;
      
      amp = k2sl * pow(bety,1.5);      phi = 3* muy;      h0030_c += amp * cos(phi);      h0030_s += amp * sin(phi);
      amp = k2sl * pow(bety,1.5);      phi =    muy;      h0021_c += amp * cos(phi);      h0021_s += amp * sin(phi);
      amp = k2sl * pow(bety,0.5) * betx;  phi = muy;      h1110_c += amp * cos(phi);      h1110_s += amp * sin(phi);
      amp = k2sl * pow(bety,0.5) * betx;  phi = 2*mux - muy;  h2001_c += amp * cos(phi);      h2001_s += amp * sin(phi);
      amp = k2sl * pow(bety,0.5) * betx;  phi = 2*mux + muy;  h2010_c += amp * cos(phi);      h2010_s += amp * sin(phi);
     }
  }

  h0030_c=   h0030_c * (-1.0/24);  h0030_s=   h0030_s * (-1.0/24);
  h0021_c=   h0021_c * (1.0/8);    h0021_s=   h0021_s * (1.0/8);
  h1110_c=   h1110_c * ( -1.0/ 4); h1110_s=   h1110_s * ( -1.0/ 4);
  h2001_c=   h2001_c * ( -1.0/ 8); h2001_s=   h2001_s * ( -1.0/ 8);
  h2010_c=   h2010_c * ( -1.0/ 8); h2010_s=   h2010_s * ( -1.0/ 8);

  cout<<">>>1st order skew sextupole geometric resonance driving terms: "<<endl;
  cout<<" h0030 :  "<<h0030_c <<" +i  "<<h0030_s<<" , "<<sqrt(h0030_c*h0030_c + h0030_s*h0030_s)<<" / "<<atan2(h0030_s,h0030_c)<<endl;
  cout<<" h0021 :  "<<h0021_c <<" +i  "<<h0021_s<<" , "<<sqrt(h0021_c*h0021_c + h0021_s*h0021_s)<<" / "<<atan2(h0021_s,h0021_c)<<endl;
  cout<<" h1110 :  "<<h1110_c <<" +i  "<<h1110_s<<" , "<<sqrt(h1110_c*h1110_c + h1110_s*h1110_s)<<" / "<<atan2(h1110_s,h1110_c)<<endl;
  cout<<" h2001 :  "<<h2001_c <<" +i  "<<h2001_s<<" , "<<sqrt(h2001_c*h2001_c + h2001_s*h2001_s)<<" / "<<atan2(h2001_s,h2001_c)<<endl;
  cout<<" h2010 :  "<<h2010_c <<" +i  "<<h2010_s<<" , "<<sqrt(h2010_c*h2010_c + h2010_s*h2010_s)<<" / "<<atan2(h2010_s,h2010_c)<<endl;
}

void Cal_3Qy_RDTs_Ring( Line & linename, const char* filename  )
// 3Qy RDT along the ring
{
  int i,j;
  double h0030_c, h0030_s;
  double k2sl, bety, dmuy;
  fstream f1;
  
  f1.open(filename, ios::out);
  
  for(i=0;i<linename.Ncell;i++){
    h0030_c=0;
    h0030_s=0;
    for(j=0;j<linename.Ncell;j++){
      if( linename.Cell[j]->TYPE==string("MULT") || linename.Cell[j]->TYPE==string("SEXT") ){
	k2sl=linename.Cell[j]->GetP("K2SL");
	bety=linename.Cell[j]->Beta2;
	if(j>i){
	  dmuy =linename.Cell[j]->Mu2 * 2.0 * PI  - linename.Cell[i]->Mu2 * 2.0 * PI;  
	  h0030_c +=  k2sl * pow(bety,1.5) * cos( dmuy * 3 );
	  h0030_s +=  k2sl * pow(bety,1.5) * sin( dmuy * 3 );
	}
	else{
	  dmuy =linename.Tune2 * 2 *PI  - abs(linename.Cell[i]->Mu2 * 2.0 * PI  - linename.Cell[j]->Mu2 * 2.0 * PI);  
	  h0030_c +=  k2sl * pow(bety,1.5) * cos( dmuy * 3 );
	  h0030_s +=  k2sl * pow(bety,1.5) * sin( dmuy * 3 );
	}
      }
    }
    h0030_c=   h0030_c * (-1.0/24);  h0030_s=   h0030_s * (-1.0/24);
    f1<<linename.Cell[i]->S<<" "<<  h0030_c<< "  "<<  h0030_s<<"  "<<sqrt( h0030_c* h0030_c +h0030_s*h0030_s)<<"  "<<atan2(h0030_s, h0030_c)<<endl;
  }
  f1.close();
}

void Cal_3Qy_RDTs_Source( Line & linename )
// 3Qy RDT contribution along the ring
{
  int j;
  double h0030_c, h0030_s;
  double k2sl, bety, dmuy;
  
  for(j=0;j<linename.Ncell;j++){
    if( linename.Cell[j]->TYPE==string("MULT") || linename.Cell[j]->TYPE==string("SEXT") ){
      k2sl=linename.Cell[j]->GetP("K2SL");
      k2sl=1;
      bety=linename.Cell[j]->Beta2;
      dmuy=linename.Cell[j]->Mu2 * 2.0 * PI;
      h0030_c =  k2sl * pow(bety,1.5) * cos( dmuy * 3 );
      h0030_s =  k2sl * pow(bety,1.5) * sin( dmuy * 3 );
      h0030_c=   h0030_c * (-1.0/24);  
      h0030_s=   h0030_s * (-1.0/24);
      cout<<linename.Cell[j]->S<<" "<<linename.Cell[j]->NAME<<"  "<<linename.Cell[j]->TYPE<<"  "<< h0030_c<< "  "<<  h0030_s<<"  "<<sqrt( h0030_c* h0030_c +h0030_s*h0030_s)<<"  "<<atan2(h0030_s, h0030_c)<<endl;
    }
  }
  
}

void Cal_Detuning_Sext( Line & linename )
{
  int i,j;
  double axx=0, axy=0, ayy=0;
  double k2l_1, betx_1, bety_1, mux_1, muy_1;
  double k2l_2, betx_2, bety_2, mux_2, muy_2;
  double dmux, dmuy;
  double Qx0=linename.Tune1, Qy0=linename.Tune2;

  for(i=0;i<linename.Ncell;i++) {
    if(linename.Cell[i]->TYPE==string("SEXT") ||  linename.Cell[i]->TYPE==string("MULT") ) {
      k2l_1=linename.Cell[i]->GetP("K2L");
      betx_1=linename.Cell[i]->Beta1;
      bety_1=linename.Cell[i]->Beta2;
      mux_1 =linename.Cell[i]->Mu1 * 2.0 * PI;
      muy_1 =linename.Cell[i]->Mu2 * 2.0 * PI;
      
      for(j=0;j<linename.Ncell;j++) {
	if(linename.Cell[j]->TYPE==string("SEXT") ||  linename.Cell[j]->TYPE==string("MULT") ) {
	  k2l_2=linename.Cell[j]->GetP("K2L");
	  betx_2=linename.Cell[j]->Beta1;
	  bety_2=linename.Cell[j]->Beta2;
	  mux_2 =linename.Cell[j]->Mu1 * 2.0 * PI;
	  muy_2 =linename.Cell[j]->Mu2 * 2.0 * PI;
	  
	  dmux=abs(mux_1 -  mux_2);   dmuy=abs(muy_1 - muy_2) ;

	  axx += k2l_1 *k2l_2 * pow(betx_1, 1.5) * pow(betx_2, 1.5) 
	    *(      cos(3*dmux -3 *PI*Qx0) /sin(3*PI*Qx0) 
		    + 3* cos(dmux -PI*Qx0 )/sin(PI*Qx0) );

	  axy += k2l_1 *k2l_2 *sqrt(betx_1*betx_2) * bety_1 * bety_2 
	    *(    cos( 2*dmuy + dmux - PI*(2*Qy0+ Qx0))/sin(PI*(2*Qy0+Qx0))
	       +  cos( 2*dmuy - dmux - PI*(2*Qy0- Qx0))/sin(PI*(2*Qy0-Qx0)) )
	       
                -k2l_1 *k2l_2 *sqrt(betx_1*betx_2) * bety_1 *betx_2
	         *2*cos(dmux-PI*Qx0)/sin(PI*Qx0);

	  ayy += k2l_1 *k2l_2 *sqrt(betx_1 *betx_2) * bety_1 * bety_2 
	    *(  cos(2*dmuy+dmux-PI*(2*Qy0+Qx0))/sin(PI*(2*Qy0+Qx0)) 
	      - cos(2*dmuy-dmux-PI*(2*Qy0-Qx0))/sin(PI*(2*Qy0-Qx0)) 
	      + 4*cos(dmux-PI*Qx0) /sin(PI*Qx0)  );
	}
      }
    }
  }
  axx=axx*(-1.0/64/PI); ayy=ayy*(-1.0/64/PI); axy= axy*(-1.0/32/PI);
  cout<<">>>Amplitude dependent tune shifts from sextupoles:" <<endl;
  cout<<"   dQx =  axx * Jx + axy *Jy "<<endl;
  cout<<"   dQy =  axy * Jx + ayy *Jy "<<endl;
  cout<<"  "<<endl;
  cout<<"  a_xx = "<<axx<<endl;
  cout<<"  a_xy = "<<axy<<endl;
  cout<<"  a_yy = "<<ayy<<endl;
}

void Cal_Detuning_Oct( Line & linename )
{
  int i;
  double axx=0, axy=0, ayy=0;
  double k3l, betx, bety;

  for(i=0;i<linename.Ncell;i++) {
    if(linename.Cell[i]->TYPE==string("OCT") ||  linename.Cell[i]->TYPE==string("MULT") ) {
      k3l=linename.Cell[i]->GetP("K3L");
      betx=linename.Cell[i]->Beta1;
      bety=linename.Cell[i]->Beta2;
      axx +=  k3l*betx*betx;
      axy +=  k3l*betx*bety;
      ayy +=  k3l*bety*bety;
    }
  }
  axx=axx*(1.0/16/PI); axy= axy*(-1.0/8/PI); ayy=ayy*(1.0/16/PI);
  cout<<">>>Amplitude dependent tune shifts from octupoles:" <<endl;
  cout<<"   dQx =  axx * Jx +  axy * Jy "<<endl;
  cout<<"   dQy =  axy * Jx +  ayy * Jy "<<endl;
  cout<<"  "<<endl;
  cout<<"  a_xx = "<<axx<<endl;
  cout<<"  a_xy = "<<axy<<endl;
  cout<<"  a_yy = "<<ayy<<endl;
}

void Cal_Q2_Oct( Line & linename )  
{
  int i;
  double k3l, betx, bety, dx;
  double chrom2x=0, chrom2y=0;

  for(i=0;i<linename.Ncell;i++) {
    if(linename.Cell[i]->TYPE==string("OCT") ||  linename.Cell[i]->TYPE==string("MULT") ) {
      k3l=linename.Cell[i]->GetP("K3L");
      betx=linename.Cell[i]->Beta1;
      bety=linename.Cell[i]->Beta2;
      dx  =linename.Cell[i]->Etax;
      chrom2x +=  k3l*betx*dx*dx;
      chrom2y -=  k3l*bety*dx*dx;
    }
  }

  chrom2x += chrom2x / 8 /PI / 2;
  chrom2y -= chrom2y / 8 /PI / 2;
  cout<<">>>Second order chromaticity from octupoles:" <<endl;
  cout<<"   0.5* d^2 Qx / d delta^2  =   "<< chrom2x <<endl;
  cout<<"   0.5* d^2 Qy / d delta^2  =   "<< chrom2y <<endl;
  cout<<"  "<<endl;
}

//------ path length and gamma-t

double  Cal_Pathlength( Line & linename, double deltap)
{
  int i;
  double l0=0;
  double x0[6];

  Cal_Orbit_Num(linename,0);  // 5-d closed orbit with deltap fixed

  for(i=0;i<6;i++) x0[i]=linename.Cell[linename.Ncell-1]->X[i];
  x0[4]=0;
  x0[5]=0;

  cout<<"Start line: "<<endl;
  for(i=0;i<6;i++) cout<<x0[i]<<"  ";
  cout<<endl;

  for(i=0;i<linename.Ncell;i++){
    if(linename.Cell[i]->TYPE != string("RFCAV") and  linename.Cell[i]->TYPE != string("CRABRF") ){
      linename.Cell[i]->Pass(x0);
      l0= l0+ linename.Cell[i]->L;
    }
      
  }

  cout<<"Finish line: "<<endl;
  for(i=0;i<6;i++) cout<<x0[i]<<"  ";
  cout<<endl;
  
  cout<<"...reference frame path-length:"<<endl;
  cout<<l0<<endl;
  
  l0=l0 - x0[4] * GP.beta ;
  
  cout<<"...on-momentum particle's path-length : "<<endl;
  cout<<l0<<endl;

  return l0;

}


void Cal_Gammat(Line & linename)
{
  int i,j;
  double x[6];
  double deltap,pt;
  double input1[21],input2[21];
  double coeff[8];
  double alfa0, alfa1,alfa2;
  fstream fout;

  if(GP.H_expand == true){
    for(i=0;i<21;i++){
      deltap=GP.step_deltap*(i-10);
      Cal_Orbit_Num(linename, deltap);
      x[0]=linename.Cell[linename.Ncell-1]->X[0];
      x[1]=linename.Cell[linename.Ncell-1]->X[1];
      x[2]=linename.Cell[linename.Ncell-1]->X[2];
      x[3]=linename.Cell[linename.Ncell-1]->X[3];
      x[4]=0.000;
      x[5]=deltap;
      for(j=0;j<linename.Ncell; j++) linename.Cell[j]->Pass(x);
      input1[i] = deltap;
      input2[i] = -x[4]*GP.beta/linename.Length;
    } 
  }
  else{
    for(i=0;i<21;i++){
      deltap=GP.step_deltap*(i-10);
      pt=DeltaToPt(deltap);
      Cal_Orbit_Num(linename,pt);
      x[0]=linename.Cell[linename.Ncell-1]->X[0];
      x[1]=linename.Cell[linename.Ncell-1]->X[1];
      x[2]=linename.Cell[linename.Ncell-1]->X[2];
      x[3]=linename.Cell[linename.Ncell-1]->X[3];
      x[4]=0.000;
      x[5]=pt;
      for(j=0;j<linename.Ncell; j++) linename.Cell[j]->Pass(x);
      input1[i] = deltap;
      input2[i] = -x[4]*GP.beta/linename.Length;
    } 
  }

  fout.open("dT_vs_deltap.dat", ios::out);
  for (i = 0; i < 21; i++)
    fout << scientific << setw(15) << input1[i]<< scientific << setw(15) << input2[i]<<endl;
  fout.close();

  pfit(input1, input2, 21, coeff, 7);
  alfa0=coeff[1]; alfa1=coeff[2];   alfa2=coeff[3]; 
  
  linename.Gammat= 1.0/sqrt(alfa0 + 1.0/GP.gamma/GP.gamma ); 
  linename.Slip=alfa0;

  cout<<alfa0<<"   "<<alfa1<<"   "<<alfa2<<endl;
}

void Cal_Orbit_Length(Line & linename, double deltap)
{
  int j;
  double x[6];
  
  Cal_Orbit_Num(linename, deltap);
  x[0]=linename.Cell[linename.Ncell-1]->X[0];
  x[1]=linename.Cell[linename.Ncell-1]->X[1];
  x[2]=linename.Cell[linename.Ncell-1]->X[2];
  x[3]=linename.Cell[linename.Ncell-1]->X[3];
  x[4]=0.000;
  x[5]=deltap;
  for(j=0;j<linename.Ncell; j++) linename.Cell[j]->Pass(x);
  linename.Orbit_Length =  -x[4] * 3.0e8 +  linename.Length;
}

//------ longitudinal optics calculations

void  Cal_Qs(Line & linename)
{
  int i;
  double Vrf_tot=0,  phi_s =0;
  
  for(i=0;i<linename.Ncell;i++)
    if(linename.Cell[i]->TYPE==string("RFCAV") ){
      Vrf_tot +=linename.Cell[i]->GetP("VRF");
    }
  linename.Vrf_tot= Vrf_tot;
  
  Cal_Gammat(linename);
  linename.Qs= sqrt( GP.harm * abs(GP.Q) * Vrf_tot * abs(linename.Slip * cos(phi_s) )/ 2/PI/GP.beta/GP.beta/GP.energy ) ;
}

void Cal_Bucket_Area(Line & linename)
{
  int i;
  double Vrf_tot=0;
  
  for(i=0;i<linename.Ncell;i++)
    if(linename.Cell[i]->TYPE==string("RFCAV") ){
      Vrf_tot +=linename.Cell[i]->GetP("VRF");
    }

  linename.Bucket_length= linename.Length * 1.0e9 / GP.harm / ( GP.beta * speed_light) ;
  linename.Bucket_height= 2*sqrt(abs(GP.Q)*Vrf_tot/2/PI/GP.beta/GP.beta/GP.energy/GP.harm/abs(linename.Slip)); 
  linename.Bucket_area=16*1.0e6*sqrt( GP.beta*GP.beta*GP.energy * abs(GP.Q) * Vrf_tot /2/PI/(2*PI*linename.frev0)/(2*PI*linename.frev0)/GP.harm/GP.harm/GP.harm/linename.Slip ) / GP.A ;
}

double RF_F_function(double phi_s, double phi_right, double phi_left)
{
  int i;
  double dphi,phi;
  double sum=0;

  dphi=(phi_right - phi_left) / 1000;
  for(i=0;i<1000;i++){
    phi= phi_left + dphi * (i+1);
    sum += sqrt( abs(cos(phi_left)-cos(phi) + (phi_left-phi) * sin(phi_s)  ) ) * dphi;
  }
  return sum*sqrt(2.0)/8;
}

void Cal_Bunch_Area(Line & linename, double full_length)
// full length is +/-3sigma_l, that is, 6 sigma_l, in units of ns
{
  double phi_s=0, delta_phi, phi_right, phi_left;

  linename.Bunch_length =  full_length ; 
  delta_phi =  full_length * 2 *PI / linename.Bucket_length ;
  phi_right =   delta_phi /2  ;
  phi_left =   -delta_phi /2  ;
  linename.Bunch_area  =   RF_F_function( 0, phi_right, phi_left) * linename.Bucket_area ;
  linename.Bunch_height =  sqrt(abs(cos(phi_right)  - cos(phi_s) + ( phi_right - phi_s ) * sin(phi_s) ) )
    *  linename.Bucket_area * GP.harm * ( 2*PI*linename.frev0)/ 8/ sqrt(2.0) / (GP.energy*1.0e6)*GP.beta*GP.beta;
}

void Cal_Bunch_Height(Line & linename, double bunch_area)
{
  int i;
  double guess,  bunch_area_0,  bunch_area_1, scale;
  
  Cal_Qs(linename);
  Cal_Bucket_Area(linename);

  i=0;
  guess = linename.Bucket_length*0.5; 
  do{
    Cal_Bunch_Area( linename, guess);
    bunch_area_0  = linename.Bunch_area ;   
    Cal_Bunch_Area( linename, guess+0.05);
    bunch_area_1  = linename.Bunch_area ; 
    scale = 0.05 / ( bunch_area_1 -  bunch_area_0 );
    
    guess = ( bunch_area -  bunch_area_0 ) * scale + guess; 
    Cal_Bunch_Area( linename, guess);
    i++;
  } while (i < 30 && abs(bunch_area - linename.Bunch_area ) > 0.1 );   // there may be a problem for this function.

}

void Print_Longitudinal_Summary( Line & linename)
{
  cout<<"------------------------------------------------------"<<endl;
  cout<<"Beam energy  = "<<GP.energy<<"  MeV "<<endl;
  cout<<"gamma        = "<<GP.gamma<<endl;
  cout<<"beta         = "<<GP.beta<<endl;
  cout<<"circumference= "<<linename.Length<<"  m   "<<endl;
  cout<<"Revolution frequency = "<< linename.frev0 <<"  Hz  "<<endl;
  cout<<"particle's A : "<<GP.A<<endl;
  cout<<"particle's Q  : "<<GP.Q<<endl;
  cout<<"harmnic number = "<<GP.harm<<endl;
  cout<<"RF total voltage = "<< linename.Vrf_tot <<"  MV "<<endl; 
  cout<<"Alpha_p  = "<< linename.Alfa0<<endl;
  cout<<"GammaT  = "<< linename.Gammat<<endl;
  cout<<"Phase slip factor = " << linename.Slip<<endl;
  cout<<"Qs= "<<linename.Qs<<endl;  
  cout<<"bucket length = "<<linename.Bucket_length<<"  ns   "<<endl; 
  cout<<"bucket height  (dp/p0_max) = "<<linename.Bucket_height<<endl; 
  cout<<"bucket area (un-normalized) = "<<linename.Bucket_area << "  eV.s/n  "<<endl;  
  cout<<"bunch length = "<<linename.Bunch_length<<"  ns  "<<endl; 
  cout<<"bunch area (un-normalized) = "<<linename.Bunch_area<<"  eV.s/n   "<<endl; 
  cout<<"bunch height (dp/p0_max) = "<<linename.Bunch_height<<endl; 
  cout<<"------------------------------------------------------"<<endl;
}

//-------- Print and Plot
void Print_Optics(Line & linename,  int  i)
{
  cout<<"-----------------optics  at one point----------------"<<endl;
  cout<<setw(6) <<"NAME"<<setw(6) <<"TYPE"<<setw(6) <<"L"<<setw(6) <<"S"<<" : "
      <<linename.Cell[i]->NAME<<" "<<linename.Cell[i]->TYPE<<" "<<linename.Cell[i]->L<<" "<<linename.Cell[i]->S<<endl;
  cout<<setw(6) <<"BETA1"<<setw(6) <<"ALFA1"<<setw(6) <<"MU1"<<" : "
      <<linename.Cell[i]->Beta1<<" "<<linename.Cell[i]->Alfa1<<" "<<linename.Cell[i]->Mu1<<endl;
  cout<<setw(6) <<"BETA2"<<setw(6) <<"ALFA2"<<setw(6) <<"MU2"<<" : "
      <<linename.Cell[i]->Beta2<<" "<<linename.Cell[i]->Alfa2<<" "<<linename.Cell[i]->Mu2<<endl;
  cout<<setw(6) <<"ETAX"<<setw(6)  <<"ETAXP"<<setw(6) <<"ETAY"<<setw(6)  <<"ETAYP"<<" : "
      <<linename.Cell[i]->Etax<<" "<<linename.Cell[i]->Etaxp<<" "<<linename.Cell[i]->Etay<<" "<<linename.Cell[i]->Etayp<<endl;
  cout<<setw(6)<<"r"<<setw(6)<<"C11"<<setw(6)  <<"C12"<<setw(6) <<"C21"<<setw(6)  <<"C22"<<" : "
      <<linename.Cell[i]->r<<" "<<linename.Cell[i]->c11<<" "<<linename.Cell[i]->c12<<" "<<linename.Cell[i]->c21<<" "<<linename.Cell[i]->c22<<endl;
}

void Print_Orbit(Line & linename, const char* filename)
{
 int i;
 fstream fout;
 
 fout.open(filename, ios::out);
 fout <<setw(15) <<"NAME"<<setw(15) <<"TYPE"<<setw(15) <<"S"
      <<setw(15) <<"X   "<<setw(15)  <<"PX    "<<setw(15)<<"Y    "<<setw(15)<<"PY    "<<setw(15)<<"Z    "<<setw(15)<<"Pt    "<<endl;
 for(i=0;i<linename.Ncell;i++)
   fout <<setw(15) <<linename.Cell[i]->NAME<<setw(15) <<linename.Cell[i]->TYPE<<setw(15) <<linename.Cell[i]->S
	<<setw(15) <<linename.Cell[i]->X[0]<<setw(15) <<linename.Cell[i]->X[1]
        <<setw(15) <<linename.Cell[i]->X[2]<<setw(15) <<linename.Cell[i]->X[3]
        <<setw(15) <<linename.Cell[i]->X[4]<<setw(15) <<linename.Cell[i]->X[5]
        <<endl;
 fout.close();
}

void Print_Twiss(Line & linename, const char* filename)
{
 int i;
 fstream fout;
 
 fout.open(filename, ios::out);
 fout <<setw(15) <<"NAME"<<setw(15) <<"TYPE"<<setw(15) <<"S"
      <<setw(15) <<"BETA1"<<setw(15) <<"ALFA1"<<setw(15) <<"MU1"
      <<setw(15) <<"BETA2"<<setw(15) <<"ALFA2"<<setw(15) <<"MU2"
      <<setw(15) <<"ETAX"<<setw(15)  <<"ETAXP"<<setw(15) <<"ETAY"<<setw(15)  <<"ETAYP"<<endl;
 for(i=0;i<linename.Ncell;i++)
   fout <<setw(15) <<linename.Cell[i]->NAME<<setw(15) <<linename.Cell[i]->TYPE<<setw(15) <<linename.Cell[i]->S
	<<setw(15) <<linename.Cell[i]->Beta1<<setw(15) <<linename.Cell[i]->Alfa1<<setw(15) <<linename.Cell[i]->Mu1
	<<setw(15) <<linename.Cell[i]->Beta2<<setw(15) <<linename.Cell[i]->Alfa2<<setw(15) <<linename.Cell[i]->Mu2
	<<setw(15) <<linename.Cell[i]->Etax<<setw(15) <<linename.Cell[i]->Etaxp<<setw(15) <<linename.Cell[i]->Etay<<setw(15) <<linename.Cell[i]->Etayp<<endl;
 fout.close();
}

void Print_Twiss_6D(Line & linename, const char* filename)
{
 int i;
 fstream fout;
 
 fout.open(filename, ios::out);
 fout <<setw(15) <<"NAME"<<setw(15) <<"TYPE"<<setw(15) <<"S"
      <<setw(15) <<"X   "<<setw(15)  <<"Y    "<<setw(15)<<"Z    "<<setw(15)<<"Pt"
      <<setw(15) <<"BETA1"<<setw(15) <<"ALFA1"<<setw(15) <<"MU1"
      <<setw(15) <<"BETA2"<<setw(15) <<"ALFA2"<<setw(15) <<"MU2"
      <<setw(15) <<"BETA3"<<setw(15) <<"ALFA3"<<setw(15) <<"MU3"
      <<setw(15) <<"ETAX"<<setw(15)  <<"ETAXP"<<setw(15) <<"ETAY"<<setw(15)  <<"ETAYP"<<endl;
 for(i=0;i<linename.Ncell;i++)
   fout <<setw(15) <<linename.Cell[i]->NAME<<setw(15) <<linename.Cell[i]->TYPE<<setw(15) <<linename.Cell[i]->S
	<<setw(15) <<linename.Cell[i]->X[0]<<setw(15) <<linename.Cell[i]->X[2]<<setw(15) <<linename.Cell[i]->X[4]<<setw(15) <<linename.Cell[i]->X[5]
	<<setw(15) <<linename.Cell[i]->Beta1<<setw(15) <<linename.Cell[i]->Alfa1<<setw(15) <<linename.Cell[i]->Mu1
	<<setw(15) <<linename.Cell[i]->Beta2<<setw(15) <<linename.Cell[i]->Alfa2<<setw(15) <<linename.Cell[i]->Mu2
	<<setw(15) <<linename.Cell[i]->Beta3<<setw(15) <<linename.Cell[i]->Alfa3<<setw(15) <<linename.Cell[i]->Mu3
	<<setw(15) <<linename.Cell[i]->Etax<<setw(15) <<linename.Cell[i]->Etaxp<<setw(15) <<linename.Cell[i]->Etay<<setw(15) <<linename.Cell[i]->Etayp<<endl;
 fout.close();
}

void Print_Twiss_Coupling(Line & linename, const char* filename)
{
  int i;
  fstream fout;
 
  fout.open(filename, ios::out);
  fout <<setw(15) <<"NAME"<<setw(15) <<"TYPE"<<setw(15) <<"S"
       <<setw(15) <<"BETA1"<<setw(15) <<"ALFA1"<<setw(15) <<"MU1"
       <<setw(15) <<"BETA2"<<setw(15) <<"ALFA2"<<setw(15) <<"MU2"
       <<setw(15) <<"c11"<<setw(15) <<"c12"<<setw(15) <<"c21"
       <<setw(15) <<"c22"<<setw(15) <<"r"<<setw(15)
       <<setw(15) <<"ETAX"<<setw(15)  <<"ETAXP"<<setw(15) <<"ETAY"<<setw(15)  <<"ETAYP"<<endl;
  for(i=0;i<linename.Ncell;i++)
    fout <<setw(15) <<linename.Cell[i]->NAME<<setw(15) <<linename.Cell[i]->TYPE<<setw(15) <<linename.Cell[i]->S
	 <<setw(15) <<linename.Cell[i]->Beta1<<setw(15) <<linename.Cell[i]->Alfa1<<setw(15) <<linename.Cell[i]->Mu1
	 <<setw(15) <<linename.Cell[i]->Beta2<<setw(15) <<linename.Cell[i]->Alfa2<<setw(15) <<linename.Cell[i]->Mu2
	 <<setw(15) <<linename.Cell[i]->c11<<setw(15) <<linename.Cell[i]->c12<<setw(15) <<linename.Cell[i]->c21
	 <<setw(15) <<linename.Cell[i]->c22<<setw(15) <<linename.Cell[i]->r
	 <<setw(15) <<linename.Cell[i]->Etax<<setw(15) <<linename.Cell[i]->Etaxp<<setw(15) <<linename.Cell[i]->Etay<<setw(15) <<linename.Cell[i]->Etayp<<endl;
  fout.close();
}



void Print_Dispersion(Line & linename, const char* filename)
{
 int i;
 fstream fout;
 
 fout.open(filename, ios::out);
 fout <<setw(15) <<"NAME"<<setw(15) <<"TYPE"<<setw(15) <<"S"
      <<setw(15) <<"ETAX   "<<setw(15)  <<"ETAXP   "<<setw(15)<<"ETAY    "<<setw(15)<<"ETAYP    "
      <<setw(15) <<"KSIX   "<<setw(15)  <<"KSIXP   "<<setw(15)<<"KSIY    "<<setw(15)<<"KSIYP    "<<endl;
 for(i=0;i<linename.Ncell;i++)
   fout <<setw(15) <<linename.Cell[i]->NAME<<setw(15) <<linename.Cell[i]->TYPE<<setw(15) <<linename.Cell[i]->S
	<<setw(15) <<linename.Cell[i]->Etax<<setw(15) <<linename.Cell[i]->Etaxp
        <<setw(15) <<linename.Cell[i]->Etay<<setw(15) <<linename.Cell[i]->Etayp
	<<setw(15) <<linename.Cell[i]->Ksix<<setw(15) <<linename.Cell[i]->Ksixp
        <<setw(15) <<linename.Cell[i]->Ksiy<<setw(15) <<linename.Cell[i]->Ksiyp
        <<endl;
 fout.close();
}

void Print_A_Matrix(Line & linename, const char* filename)
{
  int i,j;
 fstream fout;
 
 fout.open(filename, ios::out);
 for(i=0;i<linename.Ncell;i++) {
   fout<<">>>>"<<setw(15) <<linename.Cell[i]->NAME<<setw(15) <<linename.Cell[i]->TYPE<<setw(15) <<linename.Cell[i]->L<<setw(15) <<linename.Cell[i]->S<<endl;
   for(j=0;j<6;j++)
     fout<<setw(15) <<linename.Cell[i]->A[j*6+0]<<setw(15) <<linename.Cell[i]->A[j*6+1]<<setw(15) <<linename.Cell[i]->A[j*6+2]<<setw(15)
	 <<setw(15) <<linename.Cell[i]->A[j*6+3]<<setw(15) <<linename.Cell[i]->A[j*6+4]<<setw(15) <<linename.Cell[i]->A[j*6+5]<<endl;
 }
 fout.close();
}

void Print_Optics_Summary(Line & linename)
{
  cout<<"-----------------------------------------------------------"<<endl;
  cout<<setw(12)<<" Optics Summary:"<<endl;
  cout<<setw(12)<<" Length:"<<setw(25)<<setprecision(15)<<linename.Length<<endl; 
  cout<<setw(12)<<" x, px   :"<<setw(25)<<setprecision(15)<<linename.Cell[linename.Ncell-1]->X[0]<<setw(25)<<setprecision(15)<<linename.Cell[linename.Ncell-1]->X[1]<<endl;
  cout<<setw(12)<<" y, py   :"<<setw(25)<<setprecision(15)<<linename.Cell[linename.Ncell-1]->X[2]<<setw(25)<<setprecision(15)<<linename.Cell[linename.Ncell-1]->X[3]<<endl;
  cout<<setw(12)<<" Tunes :" <<setw(25)<<setprecision(15)<<linename.Tune1<<setw(25)<<setprecision(15)<<linename.Tune2<<endl;
  cout<<setw(12)<<" Chrom1:" <<setw(25)<<setprecision(15)<<linename.Chromx1<<setw(25)<<setprecision(15)<<linename.Chromy1<<endl;
  cout<<setw(12)<<" Chrom2:" <<setw(25)<<setprecision(15)<<linename.Chromx2<<setw(25)<<setprecision(15)<<linename.Chromy2<<endl;
  cout<<setw(12)<<" Chrom3:" <<setw(25)<<setprecision(15)<<linename.Chromx3<<setw(25)<<setprecision(15)<<linename.Chromy3<<endl;
  cout<<setw(12)<<" Beta0 :" <<setw(25)<<setprecision(15)<<linename.Cell[linename.Ncell-1]->Beta1<<setw(25)<<setprecision(15)<<linename.Cell[linename.Ncell-1]->Beta2<<endl;
  cout<<setw(12)<<" Alfa0 :" <<setw(25)<<setprecision(15)<<linename.Cell[linename.Ncell-1]->Alfa1<<setw(25)<<setprecision(15)<<linename.Cell[linename.Ncell-1]->Alfa2<<endl;
  cout<<setw(12)<<" Etax0 :" <<setw(25)<<setprecision(15)<<linename.Cell[linename.Ncell-1]->Etax<<setw(25)<<setprecision(15)<<linename.Cell[linename.Ncell-1]->Etaxp<<endl;
  cout<<setw(12)<<" Etay0 :" <<setw(25)<<setprecision(15)<<linename.Cell[linename.Ncell-1]->Etay<<setw(25)<<setprecision(15)<<linename.Cell[linename.Ncell-1]->Etayp<<endl;
  cout<<setw(12)<<" Ksix0 :" <<setw(25)<<setprecision(15)<<linename.Cell[linename.Ncell-1]->Ksix<<setw(25)<<setprecision(15)<<linename.Cell[linename.Ncell-1]->Ksixp<<endl;
  cout<<setw(12)<<" Ksiy0 :" <<setw(25)<<setprecision(15)<<linename.Cell[linename.Ncell-1]->Ksiy<<setw(25)<<setprecision(15)<<linename.Cell[linename.Ncell-1]->Ksiyp<<endl;
  cout<<setw(12)<<" ||C|| :" <<setw(25)<<(linename.Cell[linename.Ncell-1]->c11 * linename.Cell[linename.Ncell-1]->c22 - linename.Cell[linename.Ncell-1]->c12 * linename.Cell[linename.Ncell-1]->c21 )<<endl;
  int  i;
  double betax_max=0, betay_max=0;
  for(i=0;i<linename.Ncell;i++){
    if ( linename.Cell[i]->Beta1 > betax_max ) betax_max = linename.Cell[i]->Beta1 ;
    if ( linename.Cell[i]->Beta2 > betay_max ) betay_max = linename.Cell[i]->Beta2 ;
  }
  cout<<setw(12)<<" Beta_max:" <<setw(25)<<setprecision(15)<<betax_max<<setw(25)<<setprecision(15)<<betay_max<<endl;
  cout<<"-----------------------------------------------------------"<<endl;  
}

void Plot_Twiss(Line & linename)
{
  int i;
  fstream fout;
  fout.open("temp_twiss", ios::out);
  for(i=0;i<linename.Ncell;i++)
   fout<<linename.Cell[i]->TYPE<<"  "<<linename.Cell[i]->L<<"  "<<linename.Cell[i]->S<<"    "<<linename.Cell[i]->Beta1<<"    "<<linename.Cell[i]->Beta2<<"  "<<linename.Cell[i]->Etax<<endl;
 fout.close();

 fout.open("temp111.p", ios::out);
 fout<<"set term post color enhanced 20 solid "<<endl;
 fout<<"set output 'twiss.ps' "<<endl;
 fout<<"set xlabel 's [m]'  " <<endl;
 fout<<"set ylabel '{/Symbol \142}_x , {/Symbol \142}_y [m]  " <<endl;
 fout<<"set y2label 'Dx [m]' "<<endl;
 fout<<"set ytics nomirror"<<endl;
 fout<<"set y2tics"<<endl;
 fout<<"plot 'temp_twiss' u 3:4 tit '{/Symbol \142}_x' w l lt 1  lw 2,\\"<<endl;
 fout<<"    'temp_twiss' u 3:5 tit '{/Symbol \142}_y' w l lt 3 lw 2,\\"<<endl;
 fout<<"    'temp_twiss' u 3:6 axes x1y2 tit 'D_x'  w l lt 2 lw 2"<<endl;
 fout<<"exit"<<endl;
 fout.close();

 //system("gnuplot temp1.p");
 //system("rm temp1.p");
 //system("rm temp_twiss");
}

void Plot_Orbit(Line & linename)
{
 int i;
 fstream fout;
 fout.open("temp_orbit", ios::out);
 for(i=0;i<linename.Ncell;i++)
   fout<<linename.Cell[i]->TYPE<<"  "<<linename.Cell[i]->L<<"  "<<linename.Cell[i]->S<<"    "<<linename.Cell[i]->X[0]<<"    "<<linename.Cell[i]->X[1]<<endl;
 fout.close();
 
 fout.open("temp222.p", ios::out);
 fout<<"set term post color enhanced 20 "<<endl;
 fout<<"set output 'orbit.ps' "<<endl;
 fout<<"set xlabel 's [m]'  " <<endl;
 fout<<"set ylabel 'x_{co}, y_{co}  [m]  " <<endl;
 fout<<"plot 'temp_orbit' u 3:4 tit 'x_{co}' w l lw 2,\\"<<endl;
 fout<<"     'temp_orbit' u 3:5 tit 'y_{co}'  w l lt 2 lw 2"<<endl;
 fout<<"exit"<<endl;
 fout.close();

 //system("gnuplot temp1.p");
 //system("rm temp1.p");
 //system("rm temp_orbit");
}

//--------Artifical phase rotator matrix

void Print_PhaseAdvances(Line & linename, int index1, int index2)
{
  cout<<"Phases between "<<linename.Cell[index1]->NAME<<" and "<< linename.Cell[index2]->NAME<<" : ";
  cout<<(linename.Cell[index2]->Mu1- linename.Cell[index1]->Mu1)<<"  ,   "<<(linename.Cell[index2]->Mu2- linename.Cell[index1]->Mu2)<<endl;
}

void Add_Phaser(Line & linename, int loc, const char * name,  double mux, double muy)
{
  int i;
  Element * temp_element;

  double bx, ax, gx,  by, ay, gy, dx, dxp;
  double R11, R12, R21,R22, R33, R34,R43,R44, A,B,C;
  double m[6][6];
  double xco_in[6], xco_out[6];

  for(i=0;i<6;i++) xco_in[i]=0.;
  for(i=0;i<6;i++) xco_out[i]=0.;
  
  bx= linename.Cell[loc-1]->Beta1;
  ax= linename.Cell[loc-1]->Alfa1;
  gx= (1+ax*ax)/bx; 
  
  by= linename.Cell[loc-1]->Beta2;
  ay= linename.Cell[loc-1]->Alfa2;
  gy= (1+ay*ay)/by; 
  
  dx = linename.Cell[loc-1]->Etax;
  dxp= linename.Cell[loc-1]->Etaxp; 
  
  R11=cos(mux) + ax*sin(mux);
  R12=bx*sin(mux);
  R21=-gx*sin(mux);
  R22=cos(mux)-ax*sin(mux);
  
  R33=cos(muy) + ay*sin(muy);
  R34=by*sin(muy);
  R43=-gy*sin(muy);
  R44=cos(muy)-ay*sin(muy);
  
  A=R21*dx-R11*dxp+dxp;
  B=-R12*dxp+R22*dx-dx;
  C=0;
  
  m[0][0]=R11;
  m[0][1]=R12;
  m[0][2]=0;
  m[0][3]=0;
  m[0][4]=0;
  m[0][5]=dx-R11*dx-R12*dxp;
  
  m[1][0]=R21;
  m[1][1]=R22;
  m[1][2]=0;
  m[1][3]=0;
  m[1][4]=0;
  m[1][5]=dxp-R21*dx-R22*dxp;
  
  m[2][0]=0;
  m[2][1]=0;
  m[2][2]=R33;
  m[2][3]=R34;
  m[2][4]=0;
  m[2][5]=0;
  
  m[3][0]=0;
  m[3][1]=0;
  m[3][2]=R43;
  m[3][3]=R44;
  m[3][4]=0;
  m[3][5]=0;
  
  m[4][0]=A;
  m[4][1]=B;
  m[4][2]=0;
  m[4][3]=0;
  m[4][4]=1;
  m[4][5]=C-A*dx-B*dxp;
  
  m[5][0]=0;
  m[5][1]=0;
  m[5][2]=0;
  m[5][3]=0;
  m[5][4]=0;
  m[5][5]=1;
  
  temp_element= new MATRIX(name, 0.0, xco_in, xco_out, &m[0][0]);
  linename.Insert(loc, temp_element);
}

//-------ORBIT correction
void Correct_Orbit_SVD(Line & linename, int m, int n, vector<int> bpm_index, vector<int> kicker_index, int plane)
{
  /*----
  if( n > m ) {
    cout<<" Numer of BPM should be larger than the number of correctors."<<endl;
    exit(0);
  }
  ----*/
  
  int     i,j, k, p;
  double  A[m][n], U[m][m], VT[n][n], AI[n][m];
  double  eps, cut_scale=1000;
  vector<double> temp_vector; 
  vector< vector<double> > V1, U1;
  double read[m], dkick[n];
  
  if( plane ==0 ) {
      for(i=0;i<m;i++) 
	for(j=0;j<n;j++)  A[i][j]=sqrt( linename.Cell[bpm_index[i]]->Beta1 * linename.Cell[kicker_index[j]]->Beta1 ) 
	  *cos( abs(  linename.Cell[bpm_index[i]]->Mu1 - linename.Cell[kicker_index[j]]->Mu1 )*2.0*PI - PI* linename.Tune1 )
	  /2.0/sin( PI * linename.Tune1);
  }
  else{
    for(i=0;i<m;i++) 
      for(j=0;j<n;j++)  A[i][j]=sqrt( linename.Cell[bpm_index[i]]->Beta2 * linename.Cell[kicker_index[j]]->Beta2 ) 
	*cos( abs(  linename.Cell[bpm_index[i]]->Mu2*2*PI    - linename.Cell[kicker_index[j]]->Mu2*2*PI  ) - PI* linename.Tune2 )
	/2.0/sin( PI * linename.Tune2);
  }
  
  eps=1.0e-10;
  i=bmuav(&A[0][0],m,n,&U[0][0], &VT[0][0],eps,m+1);
  if( i<= 0) { cout<<" SVD failed"<<endl;    exit(0);  }
  
  for(i=0;i<n;i++) {
    if( A[i][i] < A[0][0]/cut_scale ) {
      p=i;
      break;
    }
  }
  for(i=0;i<m;i++){
    for(j=0;j<p;j++) temp_vector.push_back(U[i][j]);
    U1.push_back(temp_vector);
    temp_vector.clear();
  }
  for(i=0;i<n;i++){
    for(j=0;j<p;j++) temp_vector.push_back(VT[j][i]);
    V1.push_back(temp_vector);
    temp_vector.clear();
  } 
  for(i=0;i<n;i++) 
    for(j=0;j<m;j++) {
      AI[i][j] = 0.0;
      for(k=0;k<p;k++) AI[i][j]=AI[i][j] + V1[i][k]/ A[k][k] * U1[j][k];
    }

  if (plane ==0 ) {
    for(i=0;i<m;i++) read[i]=linename.Cell[ bpm_index[i]]->X[0]; 
    for(i=0;i<n;i++) {
      dkick[i]=0;
      for(j=0;j<m;j++) dkick[i]= dkick[i] - AI[i][j]*read[j]; 
      linename.Cell[ kicker_index[i] ]->SetP( "HKICK", linename.Cell[ kicker_index[i] ]->GetP("HKICK")+ dkick[i] ); 
    }
  }
  else{
    for(i=0;i<m;i++) read[i]=linename.Cell[ bpm_index[i]]->X[2]; 
    for(i=0;i<n;i++) {
      dkick[i]=0;
      for(j=0;j<m;j++) dkick[i]= dkick[i] - AI[i][j]*read[j]; 
      linename.Cell[ kicker_index[i] ]->SetP( "VKICK", linename.Cell[ kicker_index[i] ]->GetP("VKICK")+ dkick[i]); 
    }
  }
}

void Local_Three_Bump(Line linename, int plane,  const char *corr1,   const char *corr2,   const char *corr3, double kick1)
{
  int loc1, loc2, loc3;
  double beta1, beta2, beta3, scale2, scale3, phi_21,phi_31, phi_32,  kick2, kick3;

  loc1=Get_Index(linename,corr1, 1);
  loc2=Get_Index(linename,corr2, 1);
  loc3=Get_Index(linename,corr3, 1);
  if(plane == 1) {
    beta1 =  linename.Cell[loc1]->Beta1;
    beta2 =  linename.Cell[loc2]->Beta1;
    beta3 =  linename.Cell[loc3]->Beta1;
    phi_21 = ( linename.Cell[loc2]->Mu1 - linename.Cell[loc1]->Mu1 )  * 2.0 * PI ; 
    phi_31 = ( linename.Cell[loc3]->Mu1 - linename.Cell[loc1]->Mu1 )  * 2.0 * PI ; 
    phi_32 = ( linename.Cell[loc3]->Mu1 - linename.Cell[loc2]->Mu1 )  * 2.0 * PI ; 
    scale2=-sqrt(beta1/beta2)*sin(phi_31)/sin(phi_32);
    scale3=-sqrt(beta1/beta3)*sin(phi_21)/sin(-phi_32); 
    kick2=scale2*kick1; 
    kick3=scale3*kick1;
    Set_KL(linename,corr1,"HKICK", kick1);
    Set_KL(linename,corr2,"HKICK", kick2);   
    Set_KL(linename,corr3,"HKICK", kick3);
 }
  else if ( plane == 0 ){
    beta1 =  linename.Cell[loc1]->Beta2;
    beta2 =  linename.Cell[loc2]->Beta2;
    beta3 =  linename.Cell[loc3]->Beta2;
    phi_21 = ( linename.Cell[loc2]->Mu2 - linename.Cell[loc1]->Mu2 )  * 2.0 * PI ; 
    phi_31 = ( linename.Cell[loc3]->Mu2 - linename.Cell[loc1]->Mu2 )  * 2.0 * PI ; 
    phi_32 = ( linename.Cell[loc3]->Mu2 - linename.Cell[loc2]->Mu2 )  * 2.0 * PI ; 
    kick1 = 0.2e-03;
    scale2=-sqrt(beta1/beta2)*sin(phi_31)/sin(phi_32);
    scale3=-sqrt(beta1/beta3)*sin(phi_21)/sin(-phi_32); 
    kick2=scale2*kick1; 
    kick3=scale3*kick1;
    Set_KL(linename,corr1,"VKICK", kick1);
    Set_KL(linename,corr2,"VKICK", kick2);   
    Set_KL(linename,corr3,"VKICK", kick3);
 }
  else{
  }
  if(false){
    Cal_Twiss(linename, 0);
    Print_Twiss(linename,"./twiss");
    exit(0);
  }
}

double RMS_Leakage_Orbit( Line linename, int plane, int i1, int i2 )
{
  int i, count;
  double mean, sum, rms;
  
  sum=0.;
  count=0;
  for(i=0;i<linename.Ncell;i++){
    if( (i-i1)*(i-i2) > 0 )  {
      count++;
      sum += linename.Cell[i]->X[plane*2];
    }
  }
  mean=sum/count;
  
  sum=0.;
  count=0;
  for(i=0;i<linename.Ncell;i++) {
    if( (i-i1)*(i-i2) > 0 ){
      count++;
      sum += (linename.Cell[i]->X[plane*2]-mean) *  (linename.Cell[i]->X[plane*2]-mean) ;
    }
  }
  rms=sqrt ( sum / count );
  
  return rms; 
}

void Correct_Orbit_SlidingBump1(Line & linename, int m, int n, vector<int> bpm_index, vector<int> kicker_index, int plane)
{
  int i,j;
  double  reading , exceed=0.5e-03;
  int     loc1, loc2, loc3, loc_bpm;
  double  s1, s2, s3, phi_31, phi_21, phi_32, dphi_bpm, beta1,  beta2,  beta3, kick1, kick2, kick3, scale2, scale3;
  int     dir;
  double distance;
  
  for(i=0;i<kicker_index.size();i++ ){

    if(i== kicker_index.size()-1 ){
      loc1=kicker_index[kicker_index.size()-1];
      loc2=kicker_index[0];
      loc3=kicker_index[1];
    }
    else if (i== kicker_index.size()-2 ){
      loc1=kicker_index[kicker_index.size()-2];
      loc2=kicker_index[kicker_index.size()-1];
      loc3=kicker_index[0];
    }
    else{
      loc1=kicker_index[i];
      loc2=kicker_index[i+1];
      loc3=kicker_index[i+2];
    }
    s1 = linename.Cell[loc1]->S;
    s2 = linename.Cell[loc2]->S;
    s3 = linename.Cell[loc3]->S;

    if( s2 > s1 and s2 > s3 and s1 > s3  ) {
      loc_bpm=  loc2 ;
      dir=0;
      if ( plane ==0 ) {
	dphi_bpm=( linename.Tune1 - linename.Cell[loc_bpm]->Mu1 +  linename.Cell[loc3]->Mu1 )   * 2 * PI;
      }
      else{
	dphi_bpm=( linename.Tune2 - linename.Cell[loc_bpm]->Mu2 +  linename.Cell[loc3]->Mu2 )   * 2 * PI;
      }
    }
    else if ( s1> s2 and s1 > s3 and s3 > s2) {
      loc_bpm= loc2 ;
      dir=1;     
      if ( plane ==0 ) {
	dphi_bpm=( linename.Cell[loc_bpm]->Mu1 + linename.Tune1 -  linename.Cell[loc1]->Mu1 )   * 2 * PI;
      }
      else{
	dphi_bpm=( linename.Cell[loc_bpm]->Mu2 + linename.Tune2 -  linename.Cell[loc1]->Mu2 )   * 2 * PI;
      }
    }
    else if ( s3 > s2 and s2 > s1 ) {
      distance=100;
      loc_bpm=0;
      for(j=0;j<bpm_index.size();j++){
	if( abs( linename.Cell[ bpm_index[j] ]->S - linename.Cell[loc2]->S )  < distance ){
	  distance=  abs( linename.Cell[ bpm_index[j] ]->S - linename.Cell[loc2]->S );
	  loc_bpm   =   bpm_index[j];
	}
      }
      if(  linename.Cell[loc_bpm ]->S > linename.Cell[loc2]->S ) {
	dir=0; 
	if (plane == 0 ) {dphi_bpm= linename.Cell[loc3]->Mu1*2*PI - linename.Cell[ loc_bpm ]->Mu1 * 2.0 * PI;}
	else{ dphi_bpm= linename.Cell[loc3]->Mu2*2*PI - linename.Cell[ loc_bpm ]->Mu2 * 2.0 * PI;}
      }
      else{  
	dir=1; 
	if(plane==0 ){ dphi_bpm=  linename.Cell[ loc_bpm ]->Mu1 * 2.0 * PI - linename.Cell[loc1]->Mu1*2*PI; }
        else{ dphi_bpm=  linename.Cell[ loc_bpm ]->Mu2 * 2.0 * PI - linename.Cell[loc1]->Mu2*2*PI; }
      }
    }
    else{ cout<<"something wrong in sliding bump."<<endl; exit(1); }

    if ( plane == 0 ) {
      reading= linename.Cell[ loc_bpm ]->X[0];
      if(abs(reading) > exceed ) { 
	beta1= linename.Cell[loc1]->Beta1 ; 
	beta2= linename.Cell[loc2]->Beta1 ;  
	beta3= linename.Cell[loc3]->Beta1 ; 
	phi_21 = ( linename.Cell[loc2]->Mu1 - linename.Cell[loc1]->Mu1 )  * 2.0 * PI ; 
	phi_31 = ( linename.Cell[loc3]->Mu1 - linename.Cell[loc1]->Mu1 )  * 2.0 * PI ; 
	phi_32 = ( linename.Cell[loc3]->Mu1 - linename.Cell[loc2]->Mu1 )  * 2.0 * PI ; 
        if ( phi_21 < 0 ) phi_21 = phi_21 + linename.Tune1*2*PI;
        if ( phi_31 < 0 ) phi_31 = phi_31 + linename.Tune1*2*PI;
        if ( phi_32 < 0 ) phi_32 = phi_32 + linename.Tune1*2*PI;
	if( dir ==1 ) {
	  kick1=-reading
	    / sqrt(beta1* linename.Cell[ loc_bpm ]->Beta1 ) / sin( dphi_bpm  );
	  scale2=-sqrt(beta1/beta2)*sin(phi_31)/sin(phi_32);
	  scale3=-sqrt(beta1/beta3)*sin(phi_21)/sin(-phi_32); 
	  kick2=scale2*kick1;  kick3=scale3 * kick1;
	}
	else{
	  kick3=-reading 
	    / sqrt(beta3* linename.Cell[ loc_bpm ]->Beta1 ) / sin( dphi_bpm ); 
	  scale2=-sqrt(beta1/beta2)*sin(phi_31)/sin(phi_32);
	  scale3=-sqrt(beta1/beta3)*sin(phi_21)/sin(-phi_32); 
	  kick1=kick3/scale3; kick2=scale2*kick1;
	}
	linename.Cell[loc1]->SetP("HKICK", linename.Cell[loc1]->GetP("HKICK")+ kick1);
	linename.Cell[loc2]->SetP("HKICK", linename.Cell[loc2]->GetP("HKICK")+ kick2);
	linename.Cell[loc3]->SetP("HKICK", linename.Cell[loc3]->GetP("HKICK")+ kick3);
      }
    }
    else{
      reading= linename.Cell[ loc_bpm ]->X[2];
      if(abs(reading) > exceed ) { 
	beta1= linename.Cell[loc1]->Beta2 ; 
	beta2= linename.Cell[loc2]->Beta2 ;  
	beta3= linename.Cell[loc3]->Beta2 ; 
	phi_21 = ( linename.Cell[loc2]->Mu2 - linename.Cell[loc1]->Mu2 )  * 2.0 * PI ; 
	phi_31 = ( linename.Cell[loc3]->Mu2 - linename.Cell[loc1]->Mu2 )  * 2.0 * PI ; 
	phi_32 = ( linename.Cell[loc3]->Mu2 - linename.Cell[loc2]->Mu2 )  * 2.0 * PI ; 
        if ( phi_21 < 0 ) phi_21 = phi_21 + linename.Tune2*2*PI;
        if ( phi_31 < 0 ) phi_31 = phi_31 + linename.Tune2*2*PI;
        if ( phi_32 < 0 ) phi_32 = phi_32 + linename.Tune2*2*PI;
	if( dir ==1 ) {
	  kick1=-reading
	    / sqrt(beta1* linename.Cell[ loc_bpm ]->Beta2 ) / sin( dphi_bpm  );
	  scale2=-sqrt(beta1/beta2)*sin(phi_31)/sin(phi_32);
	  scale3=-sqrt(beta1/beta3)*sin(phi_21)/sin(-phi_32); 
	  kick2=scale2*kick1;  kick3=scale3 * kick1;
	}
	else{
	  kick3=-reading 
	    / sqrt(beta3* linename.Cell[ loc_bpm ]->Beta2 ) / sin( dphi_bpm ); 
	  scale2=-sqrt(beta1/beta2)*sin(phi_31)/sin(phi_32);
	  scale3=-sqrt(beta1/beta3)*sin(phi_21)/sin(-phi_32); 
	  kick1=kick3/scale3; kick2=scale2*kick1;
	}
	linename.Cell[loc1]->SetP("VKICK", linename.Cell[loc1]->GetP("VKICK")+ kick1);
	linename.Cell[loc2]->SetP("VKICK", linename.Cell[loc2]->GetP("VKICK")+ kick2);
	linename.Cell[loc3]->SetP("VKICK", linename.Cell[loc3]->GetP("VKICK")+ kick3);
      }
    }
  }
}

void Correct_Orbit_SlidingBump2(Line & linename, int m, int n, vector<int> bpm_index, vector<int> kicker_index, int plane)
{
  int i;
  double  reading , exceed=0.5e-03;
  int     loc1, loc2, loc3, loc_bpm;
  double  phi_31, phi_21, phi_32, dphi_bpm, beta1,  beta2,  beta3, kick1, kick2, kick3, scale2, scale3;
  
  for(i=0;i<kicker_index.size();i++ ){

    if(i== kicker_index.size()-1 ){
      loc1=kicker_index[kicker_index.size()-1];
      loc2=kicker_index[0];
      loc3=kicker_index[1];
    }
    else if (i== kicker_index.size()-2 ){
      loc1=kicker_index[kicker_index.size()-2];
      loc2=kicker_index[kicker_index.size()-1];
      loc3=kicker_index[0];
    }
    else{
      loc1=kicker_index[i];
      loc2=kicker_index[i+1];
      loc3=kicker_index[i+2];
    }
    
    loc_bpm= loc2;

    if ( plane == 0 ) {
      reading= linename.Cell[ loc_bpm ]->X[0];
      if(abs(reading) > exceed ) { 
	beta1= linename.Cell[loc1]->Beta1 ; 
	beta2= linename.Cell[loc2]->Beta1 ;  
	beta3= linename.Cell[loc3]->Beta1 ; 
	phi_21 = ( linename.Cell[loc2]->Mu1 - linename.Cell[loc1]->Mu1 )  * 2.0 * PI ; 
	phi_31 = ( linename.Cell[loc3]->Mu1 - linename.Cell[loc1]->Mu1 )  * 2.0 * PI ; 
	phi_32 = ( linename.Cell[loc3]->Mu1 - linename.Cell[loc2]->Mu1 )  * 2.0 * PI ; 
	dphi_bpm=( linename.Cell[loc_bpm]->Mu1 -  linename.Cell[loc1]->Mu1 )   * 2 * PI;
        if ( phi_21 < 0 ) phi_21 = phi_21 + linename.Tune1*2*PI;
        if ( phi_31 < 0 ) phi_31 = phi_31 + linename.Tune1*2*PI;
        if ( phi_32 < 0 ) phi_32 = phi_32 + linename.Tune1*2*PI;
	if ( dphi_bpm < 0 ) dphi_bpm = dphi_bpm + linename.Tune1*2*PI;
	kick1=-reading
	  / sqrt(beta1* linename.Cell[ loc_bpm ]->Beta1 ) / sin( dphi_bpm  );
	scale2=-sqrt(beta1/beta2)*sin(phi_31)/sin(phi_32);
	scale3=-sqrt(beta1/beta3)*sin(phi_21)/sin(-phi_32); 
	kick2=scale2*kick1;  kick3=scale3 * kick1;
	linename.Cell[loc1]->SetP("HKICK", linename.Cell[loc1]->GetP("HKICK")+ kick1);
	linename.Cell[loc2]->SetP("HKICK", linename.Cell[loc2]->GetP("HKICK")+ kick2);
	linename.Cell[loc3]->SetP("HKICK", linename.Cell[loc3]->GetP("HKICK")+ kick3);
      }
    }
    else{
      reading= linename.Cell[ loc_bpm ]->X[2];
      if(abs(reading) > exceed ) { 
	beta1= linename.Cell[loc1]->Beta2 ; 
	beta2= linename.Cell[loc2]->Beta2 ;  
	beta3= linename.Cell[loc3]->Beta2 ; 
	phi_21 = ( linename.Cell[loc2]->Mu2 - linename.Cell[loc1]->Mu2 )  * 2.0 * PI ; 
	phi_31 = ( linename.Cell[loc3]->Mu2 - linename.Cell[loc1]->Mu2 )  * 2.0 * PI ; 
	phi_32 = ( linename.Cell[loc3]->Mu2 - linename.Cell[loc2]->Mu2 )  * 2.0 * PI ; 
	dphi_bpm=( linename.Cell[loc_bpm]->Mu2  -  linename.Cell[loc1]->Mu2 )   * 2 * PI;
        if ( phi_21 < 0 ) phi_21 = phi_21 + linename.Tune2*2*PI;
        if ( phi_31 < 0 ) phi_31 = phi_31 + linename.Tune2*2*PI;
        if ( phi_32 < 0 ) phi_32 = phi_32 + linename.Tune2*2*PI;
	if ( dphi_bpm < 0 ) dphi_bpm = dphi_bpm + linename.Tune2*2*PI;
	kick1=-reading
	  / sqrt(beta1* linename.Cell[ loc_bpm ]->Beta2 ) / sin( dphi_bpm  );
	scale2=-sqrt(beta1/beta2)*sin(phi_31)/sin(phi_32);
	scale3=-sqrt(beta1/beta3)*sin(phi_21)/sin(-phi_32); 
	kick2=scale2*kick1;  kick3=scale3 * kick1;
	linename.Cell[loc1]->SetP("VKICK", linename.Cell[loc1]->GetP("VKICK")+ kick1);
	linename.Cell[loc2]->SetP("VKICK", linename.Cell[loc2]->GetP("VKICK")+ kick2);
	linename.Cell[loc3]->SetP("VKICK", linename.Cell[loc3]->GetP("VKICK")+ kick3);
      }
    }
  }
}

void Orbit_Status( Line linename, vector<int> bpm_index, int plane, double & orbit_mean, double & orbit_max, double & orbit_rms )
{
  int i;
  double temp, sum;
  double max, mean, rms;
  
  max=0;
  for(i=0;i<bpm_index.size();i++){
    temp=abs( linename.Cell[ bpm_index[i] ]->X[plane*2] );
    if( temp  > max ){
      max = temp;
    }
  }

  sum=0.;
  for(i=0;i<bpm_index.size();i++) sum += linename.Cell[ bpm_index[i] ]->X[plane*2];
  mean=sum/bpm_index.size();

  sum=0.;
  for(i=0;i<bpm_index.size();i++) sum += (linename.Cell[ bpm_index[i] ]->X[plane*2]-mean) *  (linename.Cell[ bpm_index[i] ]->X[plane*2]-mean) ;
  rms=sqrt ( sum / bpm_index.size() );

  orbit_mean = mean;
  orbit_max  = max;
  orbit_rms  = rms;
  
}


void Cal_Spin_Orbit(Line & linename, double deltap)
{
  int i,j,iter=0;
  double x0[9],x01[9],x1[9];
  double dx[3];
  double d=1.0e-06;
  double mat[2][2];
  double chi;
  double codeps=1e-12;
  int Max_iter=20;
  
  Cal_Orbit_Num(linename,deltap);

  for (i=0;i<6;i++) x0[i]=linename.Cell[linename.Ncell-1]->X[i];
  x0[6]=0.; x0[7]=1.0; x0[8]=0.;

  do{
    iter++;
    
    for(j=0;j<9;j++) x01[j]=x0[j];
    for(j=0;j<linename.Ncell;j++)  linename.Cell[j]->sPass(x01);
    
    for(j=0;j<9;j++) x1[j]=x0[j];
    x1[6]=x0[6]+d;
    for(j=0;j<linename.Ncell;j++) linename.Cell[j]->sPass(x1);
    mat[0][0]=(x1[6]-x01[6])/d;
    mat[1][0]=(x1[8]-x01[8])/d;

    for(j=0;j<9;j++) x1[j]=x0[j];
    x1[8]=x0[8]+d;
    for(j=0;j<linename.Ncell;j++) linename.Cell[j]->sPass(x1);
    mat[0][1]=(x1[6]-x01[6])/d;
    mat[1][1]=(x1[8]-x01[8])/d;
    
    dx[0]= x0[6]-x01[6];
    dx[1]= x0[8]-x01[8];
    chi=dx[0]*dx[0] + dx[1]*dx[1];
    chi=sqrt(chi/2.);

    if( chi > codeps ){
      for(i=0;i<2;i++) mat[i][i]=mat[i][i]-1.0000001;
      mat_inv(&mat[0][0],2); 
      x0[6]=x0[6] + mat[0][0]*dx[0]  + mat[0][1]*dx[1];
      x0[8]=x0[8] + mat[1][0]*dx[0]  + mat[1][1]*dx[1];
      x0[7]= sqrt(1.- x0[6]*x0[6] - x0[8]*x0[8]);
    } 
    
  }while (chi>codeps && iter <Max_iter );
  
  if(iter == Max_iter-1 ) 
    {
      cout<<"Failed to find Spin COD."<<endl;
      exit(1);
    }
  else
    {
      for(j=0;j<linename.Ncell; j++) {
	linename.Cell[j]->sPass(x0);
        chi=sqrt(x0[6]*x0[6] + x0[7]*x0[7] + x0[8]*x0[8] );
	for(i=0;i<3;i++) linename.Cell[j]->n0[i]=x0[6+i]/chi;
      } 
    }
}

void Cal_Spin_Tune(Line & linename, double deltap)
{
  int i,j;
  double n[3], m[3], l[3];
  double x[9];
  int Nturn=500;
  double sample[Nturn];
  double tune;  

  Cal_Spin_Orbit(linename, deltap);

  for(i=0;i<3;i++) n[i]  = linename.Cell[linename.Ncell-1]->n0[i];

  m[2]=0.;
  m[0]=1.0/sqrt( 1. + n[0]*n[0] /  n[1]*n[1]) ;
  m[1]=-n[0]*m[0]/n[1];
  vector_cross_product(n, m, l);
  
  for(i=0;i<6;i++) x[i]=linename.Cell[linename.Ncell-1]->X[i];
  x[6]=m[0];  x[7]=m[1];  x[8]=m[2];    

  for(j=0;j<Nturn;j++){
    for(i=0;i<linename.Ncell;i++) linename.Cell[i]->sPass(x);
    sample[j] = x[6];
  }
  
  FineTuneFinder(Nturn, sample, 0.001, 0.499,tune);
  linename.SpinTune=tune;
}

void Set_Initial_Spin(Line & linename, int Nturn, double n0[3], double x[9])
{
  int j,k;
  double b=0, b1=0., b2=0., b3=0.;
  double x1[9], x2[9], x3[9];
  
  for(j=0;j<6;j++) x1[j]=  x[j];   
  for(j=0;j<6;j++) x2[j]=  x[j];   
  for(j=0;j<6;j++) x3[j]=  x[j];
  x1[6]=1.; x1[7]=0.; x1[8]=0.;
  x2[6]=0.; x2[7]=1.; x2[8]=0.;  
  x3[6]=0.; x3[7]=0.; x3[8]=1.; 
  
  for(j=0;j<Nturn;j++){
    for(k=0;k<linename.Ncell;k++) linename.Cell[k]->sPass(x1);
    b1=b1+x1[6]*n0[0] +x1[7]*n0[1]+x1[8]*n0[2];
  }
  
  for(j=0;j<Nturn;j++){
    for(k=0;k<linename.Ncell;k++) linename.Cell[k]->sPass(x2);
    b2=b2+x2[6]*n0[0] +x2[7]*n0[1]+x2[8]*n0[2];
  }

  for(j=0;j<Nturn;j++){
    for(k=0;k<linename.Ncell;k++) linename.Cell[k]->sPass(x3);
    b3=b3+x3[6]*n0[0] +x3[7]*n0[1]+x3[8]*n0[2];
  }

  b=sqrt(b1*b1+b2*b2+b3*b3);
  b1=b1/b;   b2=b2/b;   b3=b3/b;
  x[6]=b1;   x[7]=b2;   x[8]=b3;
}

void Cal_Spin_Resonance(Line & linename, double Ggamma1, double Ggamma2, double  tunex, double tuney, double tunes)
{
  int i,k;
  double real,  imag, amp,  theta;
  double Gr1=int(Ggamma1) , Gr2=int(Ggamma2), Gr;
  double k1l, bety, muy ;
  double emit = 10.e-6 ;
  
  cout<<"Gr1  is  : "<<Gr1<<endl;
  cout<<"Gr2  is  : "<<Gr2<<endl;
   
  for(k=Gr1;k<Gr2+1;k++){
    
    Gr = k  - tuney ;
    real=0; imag=0; theta = 0 ;
    for(i=0;i<linename.Ncell;i++){
      if(linename.Cell[i]->TYPE==string("SBEND") || linename.Cell[i]->TYPE==string("GMULT") ||  linename.Cell[i]->TYPE==string("SMULT") ){
	theta= theta + linename.Cell[i]->GetP("ANGLE");
      }
      if(linename.Cell[i]->TYPE==string("QUAD") || linename.Cell[i]->TYPE==string("MULT") ||
         linename.Cell[i]->TYPE==string("GMULT")|| linename.Cell[i]->TYPE==string("SMULT")  ){
	k1l=  linename.Cell[i]->GetP("K1L");
	bety= linename.Cell[i]->Beta2;
	muy=  linename.Cell[i]->Mu2 * 2 * PI;
	real=  real + k1l * cos( Gr * theta + muy  ) * sqrt(bety);
	imag = imag + k1l * sin( Gr * theta + muy  ) * sqrt(bety);
      }
    }
    real= real * sqrt(emit/(Gr/GP.G) ) *(1+Gr ) / 4 / PI;
    imag= imag * sqrt(emit/(Gr/GP.G) ) *(1+Gr ) / 4 / PI;
    amp = sqrt(real * real + imag*imag );
    cout<<Gr<<"   "<<amp<<endl; 

    Gr = k  + tuney ;
    real=0; imag=0; theta = 0 ;
    for(i=0;i<linename.Ncell;i++){
      if(linename.Cell[i]->TYPE==string("SBEND") || linename.Cell[i]->TYPE==string("GMULT") || linename.Cell[i]->TYPE==string("SMULT")  ){
	theta= theta + linename.Cell[i]->GetP("ANGLE");
      }
      if(linename.Cell[i]->TYPE==string("QUAD") || linename.Cell[i]->TYPE==string("GMULT") || linename.Cell[i]->TYPE==string("SMULT") ){
	k1l= linename.Cell[i]->GetP("K1L");
	bety= linename.Cell[i]->Beta2;
	muy= linename.Cell[i]->Mu2 * 2 * PI;
	real= real + k1l * cos( Gr * theta - muy  ) * sqrt(bety);
	imag = imag + k1l * sin( Gr * theta - muy  ) * sqrt(bety);
      }
    }
    real= real * sqrt(emit/(Gr/GP.G) ) *(1+Gr ) / 4/ PI;
    imag= imag * sqrt(emit/(Gr/GP.G) ) *(1+Gr ) / 4/ PI;
    amp = sqrt(real * real + imag*imag );
    cout<<Gr<<"   "<<amp<<endl; 
  }
}

void Cal_SRLoss_U0rad( Line & linename)
//---calculate electron radiation loss analytically, gamma >> 1
{
  int i,j;
  double l1, angle1, href1;
  double sum = 0;

  double x0[6], L, Angle, K0L, actual_l, actual_angle; // for GSBEND

  if( GP.radiate == false ) {
    GP.U0rad =0.0;
  }
  else {
    for(i=0;i<linename.Ncell;i++){
      if(linename.Cell[i]->TYPE==string("SBEND") ){
	l1     = linename.Cell[i]->L;
	angle1 = linename.Cell[i]->GetP("ANGLE");
	href1  = angle1 / l1 ;
	sum    += Cr  * pow( (GP.energy /1000.), 4.0)  * href1 * ( angle1 /2./PI)  ;
      }
      else if ( linename.Cell[i]->TYPE==string("GSBEND") ){
	L=  linename.Cell[i]->L;
	Angle = linename.Cell[i]->GetP("ANGLE");
	K0L = linename.Cell[i]->GetP("K0L");
	for(j=0;j<6;j++) x0[j] = linename.Cell[i-1]->X[j];
	Cal_gsbend_l_angle(x0, L, Angle, K0L, actual_l, actual_angle);
	href1  = actual_angle / actual_l;
	angle1 = actual_angle;
	sum    += Cr  * pow( (GP.energy /1000.), 4.0)  * href1 * ( angle1 /2./PI)  ;
      }
    }
    GP.U0rad = sum * 1000. ;  
  }
  
}

void Cal_SRLoss_U0rad_Track( Line & linename)
//---calculate electron radiation loss analytically, gamma >> 1
{
  int i,j;
  double l1, angle1, href1;
  double sum = 0;

  double x0[6], L, Angle, K0L, actual_l, actual_angle; // for GSBEND

  if( GP.radiate == false ) {
    GP.U0rad =0.0;
  }
  else {
    for(i=0;i<linename.Ncell;i++){
      if(linename.Cell[i]->TYPE==string("SBEND") ){
	
	L=  linename.Cell[i]->L;
	Angle = linename.Cell[i]->GetP("ANGLE");
	K0L = Angle;
	for(j=0;j<6;j++) x0[j] = linename.Cell[i-1]->X[j];
	Cal_gsbend_l_angle(x0, L, Angle, K0L, actual_l, actual_angle);
	href1  = actual_angle / actual_l;
	angle1 = actual_angle;
	sum    += Cr  * pow( (GP.energy /1000.), 4.0)  * href1 * ( angle1 /2./PI)  ;
      }
      else if ( linename.Cell[i]->TYPE==string("GSBEND") ){
	L=  linename.Cell[i]->L;
	Angle = linename.Cell[i]->GetP("ANGLE");
	K0L = linename.Cell[i]->GetP("K0L");
	for(j=0;j<6;j++) x0[j] = linename.Cell[i-1]->X[j];
	Cal_gsbend_l_angle(x0, L, Angle, K0L, actual_l, actual_angle);
	href1  = actual_angle / actual_l;
	angle1 = actual_angle;
	sum    += Cr  * pow( (GP.energy /1000.), 4.0)  * href1 * ( angle1 /2./PI)  ;
      }
    }
    GP.U0rad = sum * 1000. ;  
  }
  
}

void  Cal_SRLoss_Particle(Line & linename, double  x[], double & loss)
 //---calculate electron radiation loss through element-by-element tracking 
{
  int i,j;
  double pt0;
  
  pt0=x[pt_];
  for(i=0;i<linename.Ncell;i++) {
    if(linename.Cell[i]->TYPE !=string("RFCAV")){
      linename.Cell[i]->Pass(x);
    }
  }
  
  loss = (x[pt_]-pt0) * GP.beta * GP.energy;
}


//=============================
//
//      Particle Tracking 
//
//=============================

void Track(Line & linename, double x[], int nturn, int & stable, int & lost_turn, int & lost_post)
{
  int j,k;
  //-----quick check 
  if( abs(x[0]) > 0.1 ||  abs(x[2]) > 0.1 || stable ==0 ) {
    stable = 0;    lost_turn= 0;    lost_post = 0;  return; }
  
  //----now we do tracking
  for(GP.turn=0; GP.turn < nturn; GP.turn++) {
    for(j=0;j<linename.Ncell;j++) {
      linename.Cell[j]->Pass(x);
      if( abs(x[0]) > linename.Cell[j]->APx or abs(x[2]) > linename.Cell[j]->APy or  isnan(x[0]) or isnan(x[2]) ) {
      	stable=0;  lost_turn= GP.turn;  lost_post=j; return ;
      }
      //each element to output
      //cout<<linename.Cell[j]->S<<"   ";
      //for(k=0;k<6;k++) cout<<x[k]<<" ";
      //cout<<endl;
    } 
    // each turn to outputgnuplot
    // for(k=0;k<6;k++) cout<<x[k]<<" ";
    // cout<<endl;
  }
}

void Track_spin(Line & linename, double x[], int nturn, int & stable, int & lost_turn, int & lost_post)
{
  int j;
  //-----quick check 
  if( abs(x[0]) > 0.1 ||  abs(x[2]) > 0.1 || stable ==0 ) {
    stable = 0;    lost_turn= 0;    lost_post = 0;  return; }
  
  //----now we do tracking
  for(GP.turn=0; GP.turn < nturn; GP.turn++) {
    for(j=0;j<linename.Ncell;j++) {
      linename.Cell[j]->sPass(x);
      if( abs(x[0]) > linename.Cell[j]->APx or abs(x[2]) > linename.Cell[j]->APy or  isnan(x[0]) or isnan(x[2])  ) {
      	stable=0;  lost_turn= GP.turn;  lost_post=j; return ;
      }
      //each element to output
      //cout<<linename.Cell[j]->S<<"   ";
      //for(k=0;k<6;k++) cout<<x[k]<<" ";
      //cout<<endl;
    } 
    // each turn to output
    // for(k=0;k<6;k++) cout<<x[k]<<" ";
    // cout<<endl;
  }
}

void Track_tbt(Line & linename, double x[], int nturn, double  x_tbt[], int & stable, int & lost_turn, int & lost_post)
{
  int j;

  //-----quick check 
  if( abs(x[0]) > 0.1 ||  abs(x[2]) > 0.1 || stable ==0 ) {
    stable = 0;     lost_turn= 0;     lost_post = 0;     return;   }

  //----now we track
  for(GP.turn=0;GP.turn<nturn; GP.turn++) {
    for(j=0;j<6;j++) x_tbt[GP.turn*6+j]=x[j];
    
    //for(j=0;j<6;j++) cout<<x[j]<<"  ";
    //cout<<endl;
    
    for(j=0;j<linename.Ncell;j++) {
      linename.Cell[j]->Pass(x);
      if( abs(x[0]) > linename.Cell[j]->APx or abs(x[2]) > linename.Cell[j]->APy or  isnan(x[0]) or isnan(x[2])   ) { 
	stable=0;  lost_turn= GP.turn;  lost_post=j; return ;
      } } }
}

void Track_tbt_spin(Line & linename, double x[], int nturn, double  x_tbt[], int & stable, int & lost_turn, int & lost_post)
{
  int j;

  //-----quick check 
  if( abs(x[0]) > 0.1 ||  abs(x[2]) > 0.1 || stable ==0 ) {
    stable = 0;     lost_turn= 0;     lost_post = 0;     return;   }

  //----now we track
  for(GP.turn=0;GP.turn<nturn; GP.turn++) {
    for(j=0;j<9;j++) x_tbt[GP.turn*9+j]=x[j];
    
    //for(j=0;j<6;j++) cout<<x[j]<<"  ";
    //cout<<endl;
    
    for(j=0;j<linename.Ncell;j++) {
      linename.Cell[j]->sPass(x);
      if( abs(x[0]) > linename.Cell[j]->APx or abs(x[2]) > linename.Cell[j]->APy or  isnan(x[0]) or isnan(x[2])  ) { 
	stable=0;  lost_turn= GP.turn;  lost_post=j; return ;
      } } }
}

void Track_tbt(Line & linename, double x[], int nturn, double  x_tbt[], int &bpm_index, int & stable, int & lost_turn, int & lost_post)
{
  int j,k;

  //-----quick check 
  if( abs(x[0]) > 0.1 ||  abs(x[2]) > 0.1 || stable ==0 ) {
    stable = 0;     lost_turn= 0;     lost_post = 0;     return;   }

  //----now we track
  for(GP.turn=0;GP.turn<nturn; GP.turn++) {
    for(j=0;j<linename.Ncell;j++) {
      linename.Cell[j]->Pass(x);
      if(j == bpm_index ){
	for(k=0;k<6;k++) x_tbt[GP.turn*6+k]=x[k];
	//for(k=0;k<6;k++) cout<<x[k]<<"  ";
	//cout<<endl;
      }
      if( abs(x[0]) > linename.Cell[j]->APx or abs(x[2]) > linename.Cell[j]->APy or  isnan(x[0]) or isnan(x[2])  ) { 
	stable=0;  lost_turn= GP.turn;  lost_post=j; return ;
      } } }
}

void Track_tbt_spin(Line & linename, double x[], int nturn, double  x_tbt[], int &bpm_index, int & stable, int & lost_turn, int & lost_post)
{
  int j,k;

  //-----quick check 
  if( abs(x[0]) > 0.1 ||  abs(x[2]) > 0.1 || stable ==0 ) {
    stable = 0;     lost_turn= 0;     lost_post = 0;     return;   }

  //----now we track
  for(GP.turn=0;GP.turn<nturn; GP.turn++) {
    for(j=0;j<linename.Ncell;j++) {
      linename.Cell[j]->sPass(x);
      if(j == bpm_index ){
	for(k=0;k<9;k++) x_tbt[GP.turn*9+k]=x[k];
	//for(k=0;k<6;k++) cout<<x[k]<<"  ";
	//cout<<endl;
      }
      if( abs(x[0]) > linename.Cell[j]->APx or abs(x[2]) > linename.Cell[j]->APy or  isnan(x[0]) or isnan(x[2])  ) { 
	stable=0;  lost_turn= GP.turn;  lost_post=j; return ;
      } } }
}

void  Cal_dz0_OnMomentumPart(Line & linename, double x0[], double dz0[])
//----calculate dz0 for on-momentum particle whose z not always zero for some special  cases  
{
  int j;
  
  for(j=0;j<linename.Ncell;j++) {
    x0[4] =  0.;  x0[5] =0.;  
    linename.Cell[j]->Pass(x0);
    dz0[j] = x0[4];
  }
  
}

void Track_wrt_OnMomentumPart(Line & linename, double x[], double dz0[],  int nturn, int & stable, int & lost_turn, int & lost_post)
//----tracking  w.r.t. on-momentum particle whose z not always zero for some special  cases  
{
  int j,k;
  //-----quick check 
  if( abs(x[0]) > 0.1 ||  abs(x[2]) > 0.1 || stable ==0 ) {
    stable = 0;    lost_turn= 0;    lost_post = 0;  return; }
  
  //----now we do tracking
  for(GP.turn=0; GP.turn < nturn; GP.turn++) {
    for(j=0;j<linename.Ncell;j++) {
      linename.Cell[j]->Pass(x);
      x[4]=x[4] - dz0[j];
      if( abs(x[0]) > linename.Cell[j]->APx or abs(x[2]) > linename.Cell[j]->APy or  isnan(x[0]) or isnan(x[2]) ) {
      	stable=0;  lost_turn= GP.turn;  lost_post=j; return ;
      }
      //each element to output
      //cout<<linename.Cell[j]->S<<"   ";
      //for(k=0;k<6;k++) cout<<x[k]<<" ";
      //cout<<endl;
    } 
    // each turn to outputgnuplot
    // for(k=0;k<6;k++) cout<<x[k]<<" ";
    // cout<<endl;
  }
  
}

//------find tune, chroms from particle tracking
void Cal_Tunes_Track(Line & linename, double x[], int Nturn, double & tune1,  double & tune2)
{
  int k;
  int stable = 1, lost_turn=0, lost_post=0;
  double  x_tbt[Nturn*6], xtemp[Nturn];

  if(x[0] < 1.0e-9 ) x[0]=x[0]+1.0e-9;
  if(x[2] < 1.0e-9 ) x[2]=x[2]+1.0e-9; 

  stable = 1; lost_turn=0; lost_post=0;
  Track_tbt(linename, x, Nturn, x_tbt, stable, lost_turn, lost_post);

  if(stable==1) {
    for(k=0;k<Nturn;k++) xtemp[k]=x_tbt[k*6+0];
    FineTuneFinder(Nturn, xtemp, tune_low, tune_high, tune1);
    for(k=0;k<Nturn;k++) xtemp[k]=x_tbt[k*6+2];
    FineTuneFinder(Nturn, xtemp,  tune_low, tune_high,tune2);
  }
  
}

void Cal_Tunes_Track(Line & linename, double deltap0)
{
  int k;
  double x[6];

  int Nturn=1024;
  int stable = 1, lost_turn=0, lost_post=0;
  double  x_tbt[Nturn*6], xtemp[1024];
  double tunex1, tuney1;
  
  Cal_Orbit_Num(linename, deltap0);

  x[0]=linename.Cell[linename.Ncell-1]->X[0] ;
  x[1]=linename.Cell[linename.Ncell-1]->X[1] ;
  x[2]=linename.Cell[linename.Ncell-1]->X[2] ;
  x[3]=linename.Cell[linename.Ncell-1]->X[3] ;
  x[4]=0. ;
  x[5]= deltap0;
  
  if(x[0] < 1.0e-9 ) x[0]=x[0]+1.0e-9;
  if(x[2] < 1.0e-9 ) x[2]=x[2]+1.0e-9; 

  stable = 1; lost_turn=0; lost_post=0;
  Track_tbt(linename, x, Nturn, x_tbt, stable, lost_turn, lost_post);

  if(stable==1) {
    for(k=0;k<1024;k++) xtemp[k]=x_tbt[k*6+0];
    FineTuneFinder(1024, xtemp, tune_low, tune_high, tunex1);
    for(k=0;k<1024;k++) xtemp[k]=x_tbt[k*6+2];
    FineTuneFinder(1024, xtemp, tune_low, tune_high, tuney1);
  }
  linename.Tune1= tunex1;
  linename.Tune2= tuney1;
}

void Cal_Chrom_Track(Line & linename)
{
  double deltap;
  double qx0, qy0, qxp, qyp, qxm, qym;

  deltap=0.0003;
  Cal_Tunes_Track(linename, deltap);
  qxp=linename.Tune1;
  qyp=linename.Tune2;
  
  deltap=-0.0003;
  Cal_Tunes_Track(linename, deltap);
  qxm=linename.Tune1;
  qym=linename.Tune2;

  deltap=0.00;
  Cal_Tunes_Track(linename, deltap);
  qx0=linename.Tune1;
  qy0=linename.Tune2;
  
  linename.Chromx1=(qxp-qxm)/2./0.0003;
  linename.Chromy1=(qyp-qym)/2./0.0003;
  linename.Chromx2=(qxp+qxm-2*qx0)/2./0.0003/0.0003;
  linename.Chromy2=(qyp+qym-2*qy0)/2./0.0003/0.0003;
}

//----track and save tbt data, fft results
void Track_tbt_FFT( Line & linename, double x0[], int Nturn, int m)
{
  int i;
  double x_tbt[6*Nturn];
  int stable=1, lost_turn=0, lost_post=0;
  double  famp[Nturn], fphi[Nturn], xfft[Nturn], yfft[Nturn], zfft[Nturn];
  fstream  f1; 
 
  Track_tbt(linename, x0, Nturn,x_tbt, stable, lost_turn, lost_post);
  f1.open("./x_tbt.dat",ios::out);
  for(i=0;i<Nturn;i++)  f1<< i<<"  "<<x_tbt[i*6+0] <<" "<<x_tbt[i*6+1]<<" "<<x_tbt[i*6+2]
			     <<"  "<<x_tbt[i*6+3] <<" "<<x_tbt[i*6+4]<<" "<<x_tbt[i*6+5]<<endl;
  f1.close();
  
  if(m!=0 ) {  
    for (i=0;i<Nturn;i++){
      famp[i]=x_tbt[6*i+0];  fphi[i]=0;
    }
    fft(m,famp,fphi);
    for(i=0;i<Nturn;i++){
      xfft[i]=sqrt(famp[i]*famp[i] + fphi[i]*fphi[i]);
    }
    
    for (i=0;i<Nturn;i++){
      famp[i]=x_tbt[6*i+2];  fphi[i]=0;
    }
    fft(m,famp,fphi);
    for(i=0;i<Nturn;i++){
      yfft[i]= sqrt(famp[i]*famp[i] + fphi[i]*fphi[i]);
    }
    
    for (i=0;i<Nturn;i++){
      famp[i]=x_tbt[6*i+4];  fphi[i]=0;
    }
    fft(m,famp,fphi);
    for(i=0;i<Nturn;i++){
      zfft[i]= sqrt(famp[i]*famp[i] + fphi[i]*fphi[i]);
    }

    f1.open("./fft.dat",ios::out);
    for(i=0;i<Nturn;i++)  f1<< 1.0*i/Nturn<<"  "<<xfft[i] <<"  "<<yfft[i] <<"  "<<zfft[i]<<endl;
    f1.close();
  }
}


void Track_tbt_tune_footprint( Line & linename, double deltap0, double emitx, double emity) 
// I prefer RF off for this , emitx, emity are unnormalized emittances
{
  double betax0, alfax0, betay0, alfay0;

  int i,j,k;
  double nsigma, nsigmax, nsigmay;
  double xn, pxn, yn, pyn, x[6];

  int Nturn=2048;
  int stable = 1, lost_turn=0, lost_post=0;
  double  x_tbt[Nturn*6];
  double xtemp[1024];
  double tunex1, tunex2, tuney1, tuney2;

  fstream  f1;
  f1.open("./footprint-output.dat", ios::out);
  f1.close();

  Cal_Twiss(linename, deltap0);
  betax0= linename.Cell[linename.Ncell-1]->Beta1; 
  betay0= linename.Cell[linename.Ncell-1]->Beta2;
  alfax0= linename.Cell[linename.Ncell-1]->Alfa1; 
  alfay0= linename.Cell[linename.Ncell-1]->Alfa2;

  Cal_Orbit_Num(linename,deltap0);

  for( i=0;i<5;i++)
    for(j=0;j<6;j++){

      nsigma=  j ;
      if(nsigma == 0 ) nsigma=0.1;
      nsigmax= nsigma * cos((i+1)*15 *PI/180);
      nsigmay= nsigma * sin((i+1)*15 *PI/180);

      xn=  sqrt(emitx) * nsigmax ; 
      pxn= 0.;
      yn=  sqrt(emity) * nsigmay ; 
      pyn= 0;
      
      x[0]= sqrt(betax0) * xn ;                                          
      x[1]=-alfax0 * xn / sqrt( betax0 ) + 1.0* pxn /  sqrt( betax0 );   
      x[2]= sqrt(betay0) * yn ;                                          
      x[3]=-alfay0 * yn / sqrt( betay0 ) + 1.0* pyn /  sqrt( betay0 );   
      x[4]= 0. ;
      x[5]=  deltap0;
      x[0]= x[0] + linename.Cell[linename.Ncell-1]->X[0]  ;
      x[1]= x[1] + linename.Cell[linename.Ncell-1]->X[1]  ;
      x[2]= x[2] + linename.Cell[linename.Ncell-1]->X[2]  ;
      x[3]= x[3] + linename.Cell[linename.Ncell-1]->X[3]  ; 

      if(x[0] == 0. ) x[0]=x[0]+1.0e-09;
      if(x[2] == 0. ) x[2]=x[2]+1.0e-09; 
      
      stable = 1; lost_turn=0; lost_post=0;
      Track_tbt(linename, x, Nturn, x_tbt, stable, lost_turn, lost_post);
      if(stable==1) {
	for(k=0;k<1024;k++) xtemp[k]=x_tbt[k*6+0];
	FineTuneFinder(1024,xtemp,tune_low, tune_high,tunex1);
	for(k=0;k<1024;k++) xtemp[k]=x_tbt[k*6+2];
	FineTuneFinder(1024,xtemp,tune_low, tune_high,tuney1);
	for(k=0;k<1024;k++) xtemp[k]=x_tbt[(k+1023)*6+0];
	FineTuneFinder(1024,xtemp,tune_low, tune_high,tunex2);
	for(k=0;k<1024;k++) xtemp[k]=x_tbt[(k+1023)*6+2];
	FineTuneFinder(1024,xtemp,tune_low, tune_high,tuney2);
	f1.open("./footprint-output.dat",  ios::out | ios::app);
	f1<<(i+1)*15<<" "<<j<<" "
	  <<setw(25)<<setprecision(18)<<scientific<<tunex1<<" "
	  <<setw(25)<<setprecision(18)<<scientific<<tunex2<<" "
	  <<setw(25)<<setprecision(18)<<scientific<<tuney1<<" "
	  <<setw(25)<<setprecision(18)<<scientific<<tuney2<<"  "
	  <<setw(25)<<setprecision(18)<<scientific<<tunex1-tunex2<<" "
	  <<setw(25)<<setprecision(18)<<scientific<<tuney1-tuney2<<endl;
	f1.close();
      }
    }
}

//-------frequency map, tune footprint etc from tracking
void Track_tbt_FMA( Line & linename, double deltap0, double sigmax0, double sigmay0 ) // I prefer RF off for this 
{
  int i,j,k;
  int Nturn=2048;
  double sigma_step=0.1;
  int sigma0=0, sigma1=6,  nsigma= int((sigma1-sigma0 ) *1.0 / sigma_step  )+1  ;  

  double x[6];
  double  x_tbt[Nturn*6];
  int stable = 1, lost_turn=0, lost_post=0;
  
  double xtemp[1024];
  double tunex1, tunex2, tuney1, tuney2;

  fstream  f1;
  f1.open("./FMA-output.dat", ios::out);
  f1.close();
  
  Cal_Orbit_Num(linename, deltap0);

  for(i=0; i<nsigma;  i++) {
    for(j=0; j<nsigma; j++ ) {
      if(  sqrt( 1.0*i*i +1.0*j*j) <= (6./sigma_step) ) {
        
        x[0] = i*sigma_step*sigmax0 + 0. ; //  for dispersive orbit can be added here
        x[1] = 0.;
        x[2] = j*sigma_step*sigmay0;
        x[3] = 0.;
        x[4] = 0.;
        x[5] = deltap0;
	x[0]= x[0] + linename.Cell[linename.Ncell-1]->X[0]  ;
	x[1]= x[1] + linename.Cell[linename.Ncell-1]->X[1]  ;
	x[2]= x[2] + linename.Cell[linename.Ncell-1]->X[2]  ;
	x[3]= x[3] + linename.Cell[linename.Ncell-1]->X[3]  ; 
        
        if(x[0] < 1.0e-09 ) x[0]=x[0]+1.0e-09;
        if(x[2] < 1.0e-09 ) x[2]=x[2]+1.0e-09; 

        stable = 1; lost_turn=0; lost_post=0;
        Track_tbt(linename, x, Nturn, x_tbt, stable, lost_turn, lost_post);
        if(stable==1) {
          for(k=0;k<1024;k++) xtemp[k]=x_tbt[k*6+0];
          FineTuneFinder(1024,xtemp,tune_low, tune_high, tunex1);
          for(k=0;k<1024;k++) xtemp[k]=x_tbt[k*6+2];
	  FineTuneFinder(1024,xtemp,tune_low, tune_high, tuney1);
          for(k=0;k<1024;k++) xtemp[k]=x_tbt[(k+1023)*6+0];
	  FineTuneFinder(1024,xtemp,tune_low, tune_high, tunex2);
          for(k=0;k<1024;k++) xtemp[k]=x_tbt[(k+1023)*6+2];
	  FineTuneFinder(1024,xtemp,tune_low, tune_high, tuney2);
          f1.open("./FMA-output.dat",  ios::out | ios::app);
          f1<<i*sigma_step<<" "<<j*sigma_step<<" "
            <<setw(25)<<setprecision(18)<<scientific<<tunex1<<" "
            <<setw(25)<<setprecision(18)<<scientific<<tunex2<<" "
            <<setw(25)<<setprecision(18)<<scientific<<tuney1<<" "
            <<setw(25)<<setprecision(18)<<scientific<<tuney2<<"  "
            <<setw(25)<<setprecision(18)<<scientific<<tunex1-tunex2<<" "
            <<setw(25)<<setprecision(18)<<scientific<<tuney1-tuney2<<endl;
          f1.close();
        }
       } 
    } 
  }
}

void Track_tbt_Lyapunov( Line & linename, double deltap0, double sigmax0, double sigmay0 )
{
  int i,j,k;

  int    Nturn=2048;
  double sigma_step=0.05;
  int    sigma0=0,  sigma1=6,  nsigma= int(  (sigma1-sigma0 ) *1.0 / sigma_step  )  ; 

  double x1[6], x2[6];
  double x1_tbt[Nturn*6], x2_tbt[Nturn*6];
  int    stable1, stable2, lost_turn=0, lost_post=0;
  double distance0=1.0e-06, distance[Nturn], Lay[Nturn];
  double sum, LayAvg;

  fstream  f1;
  f1.open("./Lmap-output.dat", ios::out);
  f1.close();
  
  for(i=0; i<nsigma;  i++) {
    for(j=0; j<nsigma; j++ ) {
      if(  sqrt( 1.0*i*i +1.0*j*j) <= (6./sigma_step) ) {
	
        x1[0] = i*sigma_step*sigmax0;
        x1[1] = 0.;
        x1[2] = j*sigma_step*sigmay0;
        x1[3] = 0.;
        x1[4] = 0.;
        x1[5] = deltap0;
	
        if(x1[0] == 0. ) x1[0]=x1[0]+1.0e-09;
        if(x1[2] == 0. ) x1[2]=x1[2]+1.0e-09; 
        x2[0] = x1[0]+0.707e-06;
        x2[1] = 0.;
        x2[2] = x1[2]+0.707e-06;
        x2[3] = 0.;
        x2[4] = 0.;
        x2[5] = deltap0;
	
        stable1=1;  lost_turn=0, lost_post=0;
        Track_tbt(linename, x1, Nturn, x1_tbt, stable1, lost_turn, lost_post);
	stable2=1;  lost_turn=0, lost_post=0;
	Track_tbt(linename, x2, Nturn, x2_tbt, stable2, lost_turn, lost_post);
        if(stable1 * stable2 == 1 ) {

          for (k=0;k<Nturn;k++){
            distance[k]=sqrt(  ( x1_tbt[k*6+0]-x2_tbt[k*6+0])*(x1_tbt[k*6+0]-x2_tbt[k*6+0])
                               +(x1_tbt[k*6+1]-x2_tbt[k*6+1])*(x1_tbt[k*6+1]-x2_tbt[k*6+1])
                               +(x1_tbt[k*6+2]-x2_tbt[k*6+2])*(x1_tbt[k*6+2]-x2_tbt[k*6+2])
                               +(x1_tbt[k*6+3]-x2_tbt[k*6+3])*(x1_tbt[k*6+3]-x2_tbt[k*6+3]) );
            Lay[k]=log(distance[k]/distance0)/(k+1);
          } 
	  
          sum=0.0;
          for(k=0;k<25;k++) sum=sum+Lay[Nturn-1-k];
          LayAvg=sum/25;
	  
	  f1.open("./Lmap-output.dat", ios::out| ios::app);
	  f1<<i*sigma_step<<"  "<<j*sigma_step<<setw(25)<<setprecision(18)<<scientific<<LayAvg<<endl;
	  f1.close();
        }

      }
    }  
  } 
}

void Track_tbt_ActionDiff( Line & linename, double nsigma, double phase, double emitx0, double emity0, int nturn, int nstep, int npart)
{
  int i,j,k;

  Cal_Twiss(linename, 0);
  double betx=linename.Cell[ linename.Ncell-1]->Beta1;
  double alfx=linename.Cell[ linename.Ncell-1]->Alfa1;
  double bety=linename.Cell[ linename.Ncell-1]->Beta2;
  double alfy=linename.Cell[ linename.Ncell-1]->Alfa2;

  double nsigmax, nsigmay;
  double xn, pxn, yn, pyn, TwoJx0, TwoJy0, phix0, phiy0 ;

  int stable=1, lost_turn=0, lost_post=0, ncount=0;
  double x[npart][6];
  double xtrack[6];
  double x_tbt[nturn*6];
  
  double temp_x[nturn], temp_px[nturn], beta, alpha, gamma,emit;
  double TwoJxPart[npart], TwoJyPart[npart];
  double meanx,meany, rmsx, rmsy;
  
  fstream  f1;
  f1.open("./Diff-output.dat", ios::out);
  f1.close();
  
  nsigmax = nsigma * cos(phase);
  nsigmay = nsigma * sin(phase);
  TwoJx0  = ( nsigmax *  nsigmax ) *  emitx0 ;   
  TwoJy0  = ( nsigmay *  nsigmay ) *  emity0 ; 
  for(i=0; i<npart; i++ ) {
    phix0=rnd(seed)*360;
    phiy0=rnd(seed)*360;
    xn  = sqrt(TwoJx0)*cos(phix0*PI/180);
    pxn =-sqrt(TwoJx0)*sin(phix0*PI/180);
    x[i][0]= xn * sqrt( betx) ;
    x[i][1]= -alfx* xn / sqrt( betx) + pxn / sqrt(betx);
    yn  = sqrt(TwoJy0)*cos(phiy0*PI/180);
    pyn =-sqrt(TwoJy0)*sin(phiy0*PI/180);
    x[i][2]= yn * sqrt( bety) ;
    x[i][3]= -alfy* yn / sqrt( bety) + pyn / sqrt(bety);
    x[i][4]=0.;
    x[i][5]=0.;
  }

  for(j=0;j<nstep;j++) {
    
    ncount=0;
    for(i=0;i<npart;i++) {
      stable=1;
      for(k=0;k<6;k++) xtrack[k]= x[i][k]; 
      if(xtrack[0] == 0. ) xtrack[0]=xtrack[0]+1.0e-09;
      if(xtrack[2] == 0. ) xtrack[2]=xtrack[2]+1.0e-09; 
      Track_tbt(linename, xtrack, nturn, x_tbt, stable, lost_turn, lost_post);
      for(k=0;k<6;k++) x[i][k] = xtrack[k]; 
      if(stable == 1) {
        ncount++;

	for(k=0;k<nturn;k++){
	  temp_x[k] = x_tbt[k*6+0];
	  temp_px[k]= x_tbt[k*6+1];
	}
        Cal_Twiss_Emit(temp_x, temp_px, nturn, beta, alpha, gamma, emit);   TwoJxPart[i]=emit;
        for(k=0;k<nturn;k++){
	  temp_x[k] = x_tbt[k*6+2];
          temp_px[k]= x_tbt[k*6+3]; 
	}
        Cal_Twiss_Emit(temp_x, temp_px, nturn, beta, alpha, gamma, emit);   TwoJyPart[i]=emit;
      }
    }
    Cal_Mean_RMS( TwoJxPart, ncount, meanx, rmsx);
    Cal_Mean_RMS( TwoJyPart, ncount, meany, rmsy);
    
    f1.open("./Diff-output.dat",  ios::out | ios::app);
    f1<< (j+1)*nturn<<"   "<<ncount<<"    "<<meanx<<"  "<<rmsx<<"     "<<meany<<"   "<<rmsy<<endl; 
    f1.close();
  }
}

//------dynamic aperture tracking
void Track_DA( Line & linename, int nturn, double deltap0, double sigmax0, double sigmay0, const char* filename ) // Track with Cell[i]->Pass(), twin particle tracking, 10 phase angles
{
  int i;
  double angle, sigma_s, sigma_e, sigma, sigma1, sigma2, sigma_max=50.0;  
  int    stable, stable1=1, stable2=1,  lost_turn=0, lost_post=0; 
  FILE *f2;
  f2=fopen(filename,"w");
  fclose(f2);
  double x[6];

  for(i=0; i<10;i++){
    
    angle=( 90/11.)*(i+1)* PI/180;
    sigma_s=0.;
    sigma_e=sigma_max;
    
    do {
      sigma = (sigma_s + sigma_e )/2.0;
      sigma1= sigma;
      sigma2= sigma + 0.05;

      x[0]= sigma1 * cos(angle) * sigmax0;
      x[1]= 0. ;
      x[2]= sigma1 * sin(angle) * sigmay0;
      x[3]= 0. ;
      x[4]= 0.;
      x[5]= deltap0;
      stable1=1; lost_turn=0; lost_post=0;
      Track(linename, x, nturn, stable1, lost_turn, lost_post);
      
      x[0]= sigma2 * cos(angle) * sigmax0;
      x[1]= 0. ;
      x[2]= sigma2 * sin(angle) * sigmay0;
      x[3]= 0. ;
      x[4]= 0.;
      x[5]= deltap0;
      stable2=1; lost_turn=0; lost_post=0;
      Track(linename, x, nturn, stable2, lost_turn, lost_post);
      
      if ( stable1* stable2 ==1 ) {
	stable =1 ; }
      else{
	stable=0;
      }
      if( stable ==0 ) 	sigma_e= sigma;
      if( stable ==1 ) 	sigma_s= sigma;
    } while ( abs( sigma_e- sigma_s) > 0.1 );

    f2=fopen(filename,"a");
    fprintf(f2,"%f %f \n", ( 90/11.)*(i+1),  sigma);
    fclose(f2);
  }
} 

void Track_DA1( Line & linename, int nturn, double deltap0, double sigmax0, double sigmay0, const char * filename ) // Track with Cell[i]->Pass(), twin particle tracking,  5 phase angles
{
  int i;
  double angle, sigma_s, sigma_e, sigma, sigma1, sigma2, sigma_max=50.0;  
  int    stable, stable1=1, stable2=1,  lost_turn=0, lost_post=0; 
  FILE *f2;
  f2=fopen(filename,"w");
  fclose(f2);
  double x[6];
  
  for(i=0; i<5;i++){
    
    angle=15*(i+1)* PI/180;
    sigma_s=0.;
    sigma_e=sigma_max;
    
    do {
      sigma = (sigma_s + sigma_e )/2.0;
      sigma1= sigma;
      sigma2= sigma + 0.05;

      x[0]= sigma1 * cos(angle) * sigmax0;
      x[1]= 0. ;
      x[2]= sigma1 * sin(angle) * sigmay0;
      x[3]= 0. ;
      x[4]= 0.;
      x[5]= deltap0;
      stable1=1; lost_turn=0; lost_post=0;
      Track(linename, x, nturn, stable1, lost_turn, lost_post);
      
      x[0]= sigma2 * cos(angle) * sigmax0;
      x[1]= 0. ;
      x[2]= sigma2 * sin(angle) * sigmay0;
      x[3]= 0. ;
      x[4]= 0.;
      x[5]= deltap0;
      stable2=1; lost_turn=0; lost_post=0;
      Track(linename, x, nturn, stable2, lost_turn, lost_post);
      
      if ( stable1* stable2 ==1 ) {
	stable =1 ; }
      else{
	stable=0;
      }
      if( stable ==0 ) 	sigma_e= sigma;
      if( stable ==1 ) 	sigma_s= sigma;
    } while ( abs( sigma_e- sigma_s) > 0.1  );

    f2=fopen(filename,"a");
    fprintf(f2,"%f %f \n", 15.*(i+1),  sigma);
    fclose(f2);
  }
} 

void Track_DA2( Line & linename, int nturn, double deltap0, double sigmax0, double sigmay0 , const char * filename) // Track with Cell[i]->Pass(), single particle tracking
{
  int i;
  double angle, sigma_s, sigma_e, sigma, sigma_max=50.0;  
  int    stable=1, lost_turn=0, lost_post=0; 
  FILE *f2;
  f2=fopen(filename,"w");
  fclose(f2);
  double x[6];

  for(i=0; i<10;i++){
    
    angle=( 90/11.)*(i+1)* PI/180;
    sigma_s=0.;
    sigma_e=sigma_max;
    
    do {
      sigma = (sigma_s + sigma_e )/2.0;
      x[0]= sigma * cos(angle) * sigmax0;
      x[1]= 0. ;
      x[2]= sigma * sin(angle) * sigmay0;
      x[3]= 0. ;
      x[4]= 0.;
      x[5]= deltap0;
      stable=1; lost_turn=0; lost_post=0;
      Track(linename, x, nturn, stable, lost_turn, lost_post);

      if( stable ==0 ) 	sigma_e= sigma;
      if( stable ==1 ) 	sigma_s= sigma;
    } while ( abs( sigma_e- sigma_s) > 0.1  );

    f2=fopen(filename,"a");
    fprintf(f2,"%f %f \n", ( 90/11.)*(i+1),  sigma);
    fclose(f2);
  }
} 


void Track_DA3( Line & linename, int nturn, double deltap0, double sigmax0, double sigmay0, const char*  filename ) // Track with Cell[i]->Pass(), single particle tracking
{
  int i;
  double angle, sigma_s, sigma_e, sigma, sigma_max=30.0;  
  int    stable=1, lost_turn=0, lost_post=0; 
  FILE *f2;
  f2=fopen(filename,"w");
  fclose(f2);
  double x[6];

  Cal_Orbit_Num(linename, 0);

  for(i=0; i<5;i++){
    
    angle=15*(i+1)* PI/180;
    sigma_s=0.;
    sigma_e=sigma_max;
    
    do {
      sigma = (sigma_s + sigma_e )/2.0;
      x[0]= sigma * cos(angle) * sigmax0;
      x[1]= 0. ;
      x[2]= sigma * sin(angle) * sigmay0;
      x[3]= 0. ;
      x[4]= 0.;
      x[5]= deltap0;
      stable=1; lost_turn=0; lost_post=0;
      Track(linename, x, nturn, stable, lost_turn, lost_post);

      if( stable ==0 ) 	sigma_e= sigma;
      if( stable ==1 ) 	sigma_s= sigma;
    } while ( abs( sigma_e- sigma_s) > 0.1  );

    f2=fopen(filename ,"a");
    fprintf(f2,"%f %f \n", 15.*(i+1),  sigma);
    fclose(f2);
  }
} 

void Track_DA_Uniform( Line & linename, int nturn, double deltap0, double angle, double sigmax0, double sigmay0, const char* filename ) // Track with Cell[i]->Pass(), twin particle tracking
{
  int i;
  double sigma_s, sigma_e, sigma, sigma1, sigma2, sigma3, sigma_max=30.0;  
  int    stable, stable1=1, lost_turn=0, lost_post=0; 
  FILE *f2;
  f2=fopen(filename,"w");
  fclose(f2);
  double x[6];

  for(i=0; i<150;i++){
    
    sigma = 1. + i*0.2;
    x[0]= sigma * cos(angle) * sigmax0;
    x[1]= 0. ;
    x[2]= sigma * sin(angle) * sigmay0;
    x[3]= 0. ;
    x[4]= 0.;
    x[5]= deltap0;
    stable=1; lost_turn=0; lost_post=0;
    Track(linename, x, nturn, stable, lost_turn, lost_post);
    
    f2=fopen(filename,"a");
    fprintf(f2,"%f %f %f \n", angle, sigma, 1.0*stable);
    fclose(f2);

    if(stable ==0 )  exit(0);
    
  }
}

void Track_Radial_old(Line & linename, double x[], int nturn, int & stable, int & lost_turn, int & lost_post)
{
  int i, j, k;
  double x0[6];
  
  //-----quick check 
  if( abs(x[0]) > 0.1 ||  abs(x[2]) > 0.1 || stable ==0 ) {
    stable = 0;    lost_turn= 0;    lost_post = 0;  return;
  }
  
  //----now we do tracking
  for(GP.turn=0; GP.turn < nturn; GP.turn++) {
    
    for(k=0;k<6;k++) x0[k] = 0.;   // linename.Cell[ linename.Ncell-1 ]->X[k];
    for(j=0;j<linename.Ncell;j++) {
      x0[4] =  0.;  x0[5] =0.;
      linename.Cell[j]->Pass(x);  
      linename.Cell[j]->Pass(x0);
      x[4] = x[4] - x0[4];
      if( abs(x[0]) > linename.Cell[j]->APx or abs(x[2]) > linename.Cell[j]->APy  or  isnan(x[0]) or isnan(x[2])  ) {
	stable=0;  lost_turn= GP.turn;  lost_post=j; return ;
      }
      //each element to output
      //cout<<linename.Cell[j]->S<<"   ";
      //for(k=0;k<6;k++) cout<<x[k]<<" ";
      //cout<<endl;
    } 
    //each turn to output
    if(GP.turn  % 100==0) {
      for(k=0;k<6;k++) cout<<x[k]<<" ";
      cout<<endl;
     }
  }

}

void Track_Radial(Line & linename, double x[], int nturn, int & stable, int & lost_turn, int & lost_post)
{
  int i, j, k;
  double x0[6];
  double x04[ linename.Ncell ];
  
  //-----quick check 
  if( abs(x[0]) > 0.1 ||  abs(x[2]) > 0.1 || stable ==0 ) {
    stable = 0;    lost_turn= 0;    lost_post = 0;  return;
  }

  //-----get  on momentum particle's  x0[4]
  for(k=0;k<6;k++) x0[k] = 0.;   // linename.Cell[ linename.Ncell-1 ]->X[k];
  for(j=0;j<linename.Ncell;j++) {
    x0[4] =  0.;  x0[5] =0.;
    linename.Cell[j]->Pass(x0);
    x04[j] = x0[4];
  }
  
  //----now we do tracking
  for(GP.turn=0; GP.turn < nturn; GP.turn++) {
    for(j=0;j<linename.Ncell;j++) {
      linename.Cell[j]->Pass(x);  
      x[4] = x[4] - x04[j];
      if( abs(x[0]) > linename.Cell[j]->APx or abs(x[2]) > linename.Cell[j]->APy  or  isnan(x[0]) or isnan(x[2])  ) {
	stable=0;  lost_turn= GP.turn;  lost_post=j; return ;
      }
      //each element to output
      //cout<<linename.Cell[j]->S<<"   ";
      //for(k=0;k<6;k++) cout<<x[k]<<" ";
      //cout<<endl;
    } 
    //each turn to output
    //if(GP.turn  % 100==0) {
    //  for(k=0;k<6;k++) cout<<x[k]<<" ";
    //  cout<<endl;
    // }
  }

}

void Track_DA_Radial( Line & linename, int nturn, double deltap0, double sigmax0, double sigmay0, const char* filename ) // Track with Cell[i]->Pass(), twin particle tracking, 10 phase angles
{
  int i;
  double angle, sigma_s, sigma_e, sigma, sigma1, sigma2, sigma_max=15.0;  
  int    stable, stable1=1, stable2=1,  lost_turn=0, lost_post=0; 
  FILE *f2;
  f2=fopen(filename,"w");
  fclose(f2);
  double x[6];

  for(i=0; i<10;i++){
    
    angle=( 90/11.)*(i+1)* PI/180;
    sigma_s=0.;
    sigma_e=sigma_max;
    
    do {
      sigma = (sigma_s + sigma_e )/2.0;
      sigma1= sigma;
      sigma2= sigma + 0.05;

      x[0]= sigma1 * cos(angle) * sigmax0;
      x[1]= 0. ;
      x[2]= sigma1 * sin(angle) * sigmay0;
      x[3]= 0. ;
      x[4]= 0.;
      x[5]= deltap0;
      stable1=1; lost_turn=0; lost_post=0;
      Track_Radial(linename, x, nturn, stable1, lost_turn, lost_post);
      
      x[0]= sigma2 * cos(angle) * sigmax0;
      x[1]= 0. ;
      x[2]= sigma2 * sin(angle) * sigmay0;
      x[3]= 0. ;
      x[4]= 0.;
      x[5]= deltap0;
      stable2=1; lost_turn=0; lost_post=0;
      Track_Radial(linename, x, nturn, stable2, lost_turn, lost_post);
      
      if ( stable1* stable2 ==1 ) {
	stable =1 ; }
      else{
	stable=0;
      }
      if( stable ==0 ) 	sigma_e= sigma;
      if( stable ==1 ) 	sigma_s= sigma;
    } while ( abs( sigma_e- sigma_s) > 0.1  );

    f2=fopen(filename,"a");
    fprintf(f2,"%f %f \n", ( 90/11.)*(i+1),  sigma);
    fclose(f2);
  }
} 

void Track_DA_Radial_Uniform( Line & linename, int nturn, double deltap0, double angle, double sigmax0, double sigmay0, const char* filename ) // Track with Cell[i]->Pass(), twin particle tracking
{
  int i,k;
  double sigma_s, sigma_e, sigma, sigma1, sigma2, sigma3, sigma_max=100.0;  
  int    stable, stable1=1, lost_turn=0, lost_post=0; 
  FILE *f2;
  f2=fopen(filename,"w");
  fclose(f2);
  double x[6];

  for(i=0; i<150;i++){

    sigma = 1.0 + i*0.2;
    x[0]= sigma * cos(angle) * sigmax0;
    x[1]= 0.;
    x[2]= sigma * sin(angle) * sigmay0;
    x[3]= 0.;
    x[4]= 0.;
    x[5]= deltap0;
    stable=1; lost_turn=0; lost_post=0;
    Track_Radial(linename, x, nturn, stable, lost_turn, lost_post);
    
    f2=fopen(filename,"a");
    fprintf(f2,"%f %f %f \n", angle, sigma, 1.0*stable);
    fclose(f2);
    if(stable ==0 )  exit(0);
    
  }
} 

void Track_Spin(Line & linename, double x[9], int nturn, int & stable, int & lost_turn, int & lost_post)
{
  int j;
  //-----quick check 
  if( abs(x[0]) > 0.1 ||  abs(x[2]) > 0.1 || stable ==0 ) {
    stable = 0;    lost_turn= 0;    lost_post = 0;    return; }
  
  //----now we do tracking
  for(GP.turn=0; GP.turn < nturn; GP.turn++) {
    for(j=0;j<linename.Ncell;j++) {
      linename.Cell[j]->sPass(x);
      if( abs(x[0]) > linename.Cell[j]->APx or abs(x[2]) > linename.Cell[j]->APy or  isnan(x[0]) or isnan(x[2])  ) { 
      	stable=0;  lost_turn= GP.turn;  lost_post=j; return ;
      }
      //each element to output
      //cout<<linename.Cell[j]->S<<"   ";
      //for(k=0;k<8;k++) cout<<x[k]<<" ";
      //cout<<endl;
    } 
    // each turn to output
    // for(k=0;k<8;k++) cout<<x[k]<<" ";
    // cout<<endl;
  }
}

/*-----

int nelement;
vector <int>     type;
vector <string>  element_name;  
vector <double>  nelement_dx, nelement_dy, nelement_dt, nelement_APx, nelement_APy;
int Ndrift, Nbend, Nquad, Nsext, Nmult, Ngmult;
int Nmatrix, Nrf, Nhobb, Nelens, Nrotat, Nsnake, Nerhicbb, Ncrabrf;

vector <double> drift_p;
vector <double> bend_p;
vector <double> quad_p;
vector <double> sext_p;
vector <double> mult_p;
vector <double> gmult_p;
vector <double> matrix_p;
vector <double> rf_p;
vector <double> hobb_p;
vector <double> elens_p;
vector <double> rotat_p;
vector <double> snake_p;
vector <double> erhicbb_p;
vector <double> crabrf_p;

void Track_Fast_Prepare(Line & linename)
{
  int i,j,k;
  char name[125];

  //-----------do some statistics
  nelement=linename.Ncell;
  Ndrift=0;
  Nbend=0;
  Nquad=0;
  Nsext=0;
  Nmult=0;
  Ngmult=0;
  Nmatrix=0;
  Nrf=0;
  Nhobb=0;
  Nelens=0;
  Nrotat=0;
  Nsnake=0;
  Nerhicbb=0;
  Ncrabrf=0;

  for(i=0; i<linename.Ncell; i++){
    if(linename.Cell[i]->TYPE==string("DRIFT") ) {
      Ndrift++;
    }
    else if(linename.Cell[i]->TYPE==string("SBEND") ) {
      Nbend++;
    }
    else if(linename.Cell[i]->TYPE==string("QUAD") ){
      Nquad++;
    }
    else if(linename.Cell[i]->TYPE==string("SEXT") ){
      Nsext++;
    }
    else if(linename.Cell[i]->TYPE==string("MULT") ){
      Nmult++;
    }
    else if(linename.Cell[i]->TYPE==string("GMULT") ){
      Ngmult++;
    }
    else if(linename.Cell[i]->TYPE==string("MATRIX") ){
      Nmatrix++;
    }
    else if(linename.Cell[i]->TYPE==string("RFCAV") ){
      Nrf++;
    }
    else if(linename.Cell[i]->TYPE==string("BEAMBEAM") ){
      Nhobb++;
    }
    else if(linename.Cell[i]->TYPE==string("ELENS") ){
      Nelens++;
    }
    else if(linename.Cell[i]->TYPE==string("ROTAT") ){
      Nrotat++;
    }
    else if(linename.Cell[i]->TYPE==string("SNAKE") ){
      Nsnake++;
    }
    else if(linename.Cell[i]->TYPE==string("ERHICBB") ){
      Nerhicbb++;
    }
    else if(linename.Cell[i]->TYPE==string("CRABRF") ){
      Ncrabrf++;
    }
    else{ 
      cout<<"Element "<<linename.Cell[i]->NAME<<" , "<<linename.Cell[i]->TYPE<<" "<<"changed to DRIFT. "<<endl;
      Ndrift++;
    }
  }
  
  if( true ) {
    cout<<"---------------------------------"<<endl;
    cout<<"Statistics: "<<endl;
    cout<<"There are "<<Ndrift<<" DRIFT"<<endl;
    cout<<"There are "<<Nbend<< " SBEDN"<<endl;
    cout<<"There are "<<Nquad<< " QUAD"<<endl;
    cout<<"There are "<<Nsext<< " SEXT"<<endl;
    cout<<"There are "<<Nmult<< " MULT"<<endl;
    cout<<"There are "<<Ngmult<<" GMULT"<<endl;
    cout<<"There are "<<Nmatrix<<" MATRIX"<<endl;
    cout<<"There are "<<Nrf<< " RFCAV"<<endl;
    cout<<"There are "<<Nhobb<< " BEAMBEAM"<<endl;
    cout<<"There are "<<Nelens<<" ELENS"<<endl;
    cout<<"There are "<<Nrotat<<" ROTAT"<<endl;
    cout<<"There are "<<Nsnake<<" SNAKE"<<endl;
    cout<<"There are "<<Nerhicbb<<" ERHICBB"<<endl;
    cout<<"There are "<<Ncrabrf<<" CRABRF"<<endl;

    cout<<"Totally "<<linename.Ncell<<" elements"<<endl;
  }
 
  //--- store strengths
  for(i=0; i<nelement; i++){
    element_name.push_back(linename.Cell[i]->NAME );
    nelement_dx.push_back(linename.Cell[i]->DX );
    nelement_dy.push_back(linename.Cell[i]->DY );
    nelement_dt.push_back(linename.Cell[i]->DPSI );
    nelement_APx.push_back(linename.Cell[i]->APx );
    nelement_APy.push_back(linename.Cell[i]->APy ); 

    if(linename.Cell[i]->TYPE==string("DRIFT") ) {
      type.push_back( 0 );
      drift_p.push_back(linename.Cell[i]->L);
    }
    else if(linename.Cell[i]->TYPE==string("SBEND") ) {
      type.push_back( 1  );
      bend_p.push_back( linename.Cell[i]->L );
      bend_p.push_back( linename.Cell[i]->GetP("ANGLE") );
      bend_p.push_back( linename.Cell[i]->GetP("E1") );
      bend_p.push_back( linename.Cell[i]->GetP("E2") );
      bend_p.push_back( linename.Cell[i]->GetP("Nint") );
    }
    else if(linename.Cell[i]->TYPE==string("QUAD") ){
      type.push_back( 2 );
      quad_p.push_back( linename.Cell[i]->L );
      quad_p.push_back( linename.Cell[i]->GetP("K1L") );
      quad_p.push_back( linename.Cell[i]->GetP("K1SL") );
      quad_p.push_back( linename.Cell[i]->GetP("Nint") );
    }
    else if(linename.Cell[i]->TYPE==string("SEXT") ){
      type.push_back( 3 );
      sext_p.push_back(linename.Cell[i]->L);
      sext_p.push_back(linename.Cell[i]->GetP("K2L"));
      sext_p.push_back(linename.Cell[i]->GetP("K2SL"));
      sext_p.push_back(linename.Cell[i]->GetP("Nint"));
    }
    else if(linename.Cell[i]->TYPE==string("MULT") ){
      type.push_back( 4 );
      mult_p.push_back(linename.Cell[i]->L);
      for(j=0;j<11;j++){
	sprintf(name, "K%dL",j);
        mult_p.push_back(linename.Cell[i]->GetP(name));
      }
      for(j=0;j<11;j++){
        sprintf(name, "K%dSL",j);
        mult_p.push_back(linename.Cell[i]->GetP(name));
      }
      mult_p.push_back(int(linename.Cell[i]->GetP("Nint")));
      mult_p.push_back(int(linename.Cell[i]->GetP("Norder")));
    }
  else if(linename.Cell[i]->TYPE==string("GMULT") ){
      type.push_back( 5 );
      gmult_p.push_back(linename.Cell[i]->L);
      gmult_p.push_back( linename.Cell[i]->GetP("ANGLE") );
      gmult_p.push_back( linename.Cell[i]->GetP("E1") );
      gmult_p.push_back( linename.Cell[i]->GetP("E2") );
      for(j=0;j<11;j++){
	sprintf(name, "K%dL",j);
        gmult_p.push_back(linename.Cell[i]->GetP(name));
      }
      for(j=0;j<11;j++){
        sprintf(name, "K%dSL",j);
        gmult_p.push_back(linename.Cell[i]->GetP(name));
      }
      gmult_p.push_back(int(linename.Cell[i]->GetP("Nint")));
      gmult_p.push_back(int(linename.Cell[i]->GetP("Norder")));
    }
   else if(linename.Cell[i]->TYPE==string("MATRIX") ){
      type.push_back( 6 );
      matrix_p.push_back( linename.Cell[i]->L); 
      matrix_p.push_back(  linename.Cell[i]->GetP("XCO_IN_X") );
      matrix_p.push_back(  linename.Cell[i]->GetP("XCO_IN_PX") );
      matrix_p.push_back(  linename.Cell[i]->GetP("XCO_IN_Y") );
      matrix_p.push_back(  linename.Cell[i]->GetP("XCO_IN_PY") );
      matrix_p.push_back(  linename.Cell[i]->GetP("XCO_IN_Z") );
      matrix_p.push_back(  linename.Cell[i]->GetP("XCO_IN_DELTA") );
      matrix_p.push_back(  linename.Cell[i]->GetP("XCO_OUT_X") );
      matrix_p.push_back(  linename.Cell[i]->GetP("XCO_OUT_PX") );
      matrix_p.push_back(  linename.Cell[i]->GetP("XCO_OUT_Y") );
      matrix_p.push_back(  linename.Cell[i]->GetP("XCO_OUT_PY") );
      matrix_p.push_back(  linename.Cell[i]->GetP("XCO_OUT_Z") );
      matrix_p.push_back(  linename.Cell[i]->GetP("XCO_OUT_DELTA") );
      for(j=0;j<6;j++)
	for(k=0;k<6;k++)
	  {
	    sprintf(name, "M%d%d",j+1,k+1);
	    matrix_p.push_back( linename.Cell[i]->GetP(name) ); 
	  }
    }
    else if(linename.Cell[i]->TYPE==string("RFCAV") ){
      type.push_back( 7 );
      rf_p.push_back( linename.Cell[i]->L );
      rf_p.push_back( linename.Cell[i]->GetP("VRF") );
      rf_p.push_back( linename.Cell[i]->GetP("FRF") );
      rf_p.push_back( linename.Cell[i]->GetP("PHASE0") );
    }
    else if(linename.Cell[i]->TYPE==string("BEAMBEAM") ){
      type.push_back( 8 );
      hobb_p.push_back( linename.Cell[i]->GetP("TREATMENT") );
      hobb_p.push_back( linename.Cell[i]->GetP("NP") );
      hobb_p.push_back( linename.Cell[i]->GetP("BBSCALE") );
      hobb_p.push_back( linename.Cell[i]->GetP("SIGMAL") );
      hobb_p.push_back( linename.Cell[i]->GetP("NSLICE") );
      hobb_p.push_back( linename.Cell[i]->GetP("EMITX") );
      hobb_p.push_back( linename.Cell[i]->GetP("BETAX") );
      hobb_p.push_back( linename.Cell[i]->GetP("ALFAX") );
      hobb_p.push_back( linename.Cell[i]->GetP("EMITY") );
      hobb_p.push_back( linename.Cell[i]->GetP("BETAY") );
      hobb_p.push_back( linename.Cell[i]->GetP("ALFAY") );
    }
    else if(linename.Cell[i]->TYPE==string("ELENS") ){
      type.push_back( 9  );
      elens_p.push_back( linename.Cell[i]->L );   
      elens_p.push_back( linename.Cell[i]->GetP("NE") );
      elens_p.push_back( linename.Cell[i]->GetP("BBSCALE") );
      elens_p.push_back( linename.Cell[i]->GetP("BETAE") );
      elens_p.push_back( linename.Cell[i]->GetP("NSLICE") );
      elens_p.push_back( linename.Cell[i]->GetP("SIGMAX") );
      elens_p.push_back( linename.Cell[i]->GetP("SIGMAY") );
    }
    else if(linename.Cell[i]->TYPE==string("ROTAT") ){
      type.push_back( 10 );
      rotat_p.push_back( linename.Cell[i]->L );
      rotat_p.push_back( linename.Cell[i]->GetP("NX") );
      rotat_p.push_back( linename.Cell[i]->GetP("NY") );
      rotat_p.push_back( linename.Cell[i]->GetP("NS") );
      rotat_p.push_back( linename.Cell[i]->GetP("ANGLE") );
    }
    else if(linename.Cell[i]->TYPE==string("SNAKE") ){
      type.push_back(11 );
      rotat_p.push_back( linename.Cell[i]->L );
      rotat_p.push_back( linename.Cell[i]->GetP("NX") );
      rotat_p.push_back( linename.Cell[i]->GetP("NY") );
      rotat_p.push_back( linename.Cell[i]->GetP("NS") );
      rotat_p.push_back( linename.Cell[i]->GetP("ANGLE") );
    }
    else if(linename.Cell[i]->TYPE==string("ERHICBB") ){
      type.push_back(12 );
      rotat_p.push_back( linename.Cell[i]->GetP("NE") );
      rotat_p.push_back( linename.Cell[i]->GetP("BBSCALE") );
      rotat_p.push_back( linename.Cell[i]->GetP("SIGMAX") );
      rotat_p.push_back( linename.Cell[i]->GetP("SIGMAY") );
    }
    else if(linename.Cell[i]->TYPE==string("CRABRF") ){
      type.push_back(13 );
      crabrf_p.push_back( linename.Cell[i]->L );
      crabrf_p.push_back( linename.Cell[i]->GetP("VRF") );
      crabrf_p.push_back( linename.Cell[i]->GetP("FRF") );
      crabrf_p.push_back( linename.Cell[i]->GetP("PHASE0") );
    }
    else{ 
      type.push_back( 0 );
      drift_p.push_back( linename.Cell[i]->L  );
    }
  }

}

void Track_Fast(Line & linename, double x[6], int nturn, int & stable, int & lost_turn, int & lost_post)
{
  int i,j;
  double  knl[11], knsl[11];
  double  m66[36], xco_in[6], xco_out[6];

  //-----quick check before track
  if( abs(x[0]) > 0.1 ||  abs(x[2]) > 0.1 || stable ==0 ) {
    stable = 0;      lost_turn= 0;       lost_post = 0;      return ;   }
  
  //---now we do track
  for(GP.turn=0;GP.turn<nturn; GP.turn++) {
    
    Ndrift=0;
    Nbend=0;
    Nquad=0;
    Nsext=0;
    Nmult=0;
    Ngmult=0;
    Nmatrix=0;     
    Nrf=0;  
    Nhobb=0;
    Nelens=0;
    Nrotat=0;
    Nsnake=0;
    Nerhicbb=0;
    Ncrabrf=0;
    
    for(i=0; i<nelement; i++){
      
      //--alignmnet error: can be used for each element 
      //GtoL(x, nelement_dx[i],nelement_dy[i],nelement_dt[i]); 
      
      if(type[i]==0 ) {
	DRIFT_Pass(x, drift_p[Ndrift]);
	Ndrift++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2])  ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==1 ) {
	SBEND_Pass(x,bend_p[Nbend*5+0],int(bend_p[Nbend*5+4]), bend_p[Nbend*5+1], bend_p[Nbend*5+2],bend_p[Nbend*5+3]); 
	Nbend++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2])  ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==2 ){
	QUAD_Pass(x, quad_p[Nquad*4+0], int(quad_p[Nquad*4+3]), quad_p[Nquad*4+1], quad_p[Nquad*4+2]); 
	Nquad++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2])  ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==3 ){
	SEXT_Pass(x, sext_p[Nsext*4+0], int(sext_p[Nsext*4+3]), sext_p[Nsext*4+1],  sext_p[Nsext*4+2]);
	Nsext++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2]) ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==4 ){
	for(j=0;j<11;j++) {
	  knl[j]=mult_p[Nmult*25+1+j]; 
	  knsl[j]=mult_p[Nmult*25+1+11+j];
	}
	MULT_Pass(x, mult_p[Nmult*25+0], int(mult_p[Nmult*25+23]), int(mult_p[Nmult*25+24]), knl, knsl);
	Nmult++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2]) ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==5 ){
	for(j=0;j<11;j++) {
	  knl[j]=gmult_p[Ngmult*28+4+j]; 
	  knsl[j]=gmult_p[Ngmult*28+4+11+j];
	}
	GMULT_Pass(x, gmult_p[Ngmult*28+0], int(gmult_p[Ngmult*28+26]),int(gmult_p[Ngmult*28+27]), gmult_p[Ngmult*28+1] ,gmult_p[Ngmult*28+2], gmult_p[Ngmult*28+3], knl, knsl);
	Ngmult++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2])  ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==6 ){
	for(j=0;j<6;j++) xco_in[j] = matrix_p[Nmatrix*49+1+j];
	for(j=0;j<6;j++) xco_out[j]= matrix_p[Nmatrix*49+1+6+j];
	for(j=0;j<36;j++) {
	  m66[j]=matrix_p[Nmatrix*49+13+j]; 
	}
	MATRIX_Pass(x, matrix_p[Nmatrix*49+0],xco_in, xco_out, m66);
	Nmatrix++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2])  ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if (type[i]==7 ){
        RFCAV_Pass(x, rf_p[Nrf*4+0], rf_p[Nrf*4+1], rf_p[Nrf*4+2],rf_p[Nrf*4+3] );
	Nrf++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2])  ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==8 ){
	BEAMBEAM_Pass(x, int(hobb_p[Nhobb*11+0]), hobb_p[Nhobb*11+1], hobb_p[Nhobb*11+2], hobb_p[Nhobb*11+3], int(hobb_p[Nhobb*11+4]), 
		      hobb_p[Nhobb*11+5], hobb_p[Nhobb*11+6], hobb_p[Nhobb*11+7], hobb_p[Nhobb*11+8], hobb_p[Nhobb*11+9],hobb_p[Nhobb*11+10] );
	Nhobb++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2])  ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==9 ){
        ELENS_Pass(x, elens_p[Nelens*7+0], elens_p[Nelens*7+1], elens_p[Nelens*7+2],
		   elens_p[Nelens*7+3], int(elens_p[Nelens*7+4]), elens_p[Nelens*7+5], elens_p[Nelens*7+6]);
	Nelens++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2])  ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==10 ){
        DRIFT_Pass(x,rotat_p[Nrotat*5+0]);
	Nrotat++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2])  ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==11 ){
        DRIFT_Pass(x,snake_p[Nsnake*5+0]);
	Nsnake++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2]) ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==12 ){
	ERHICBB_Pass(x, GP.gamma, erhicbb_p[Nerhicbb*4 +0], erhicbb_p[Nerhicbb*4 +1]);
        Nerhicbb++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i]or  isnan(x[0]) or isnan(x[2])  ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==13 ){
	CRABRF_Pass(x,crabrf_p[Ncrabrf*4 +0],crabrf_p[Ncrabrf*4 +1],crabrf_p[Ncrabrf*4 +2],crabrf_p[Ncrabrf*4 +3]);
        Ncrabrf++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2]) ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else{ }
      
      //LtoG(x, nelement_dx[i],nelement_dy[i],nelement_dt[i]); 
      
      //print out each element 
      //cout<<i<<"   ";
      //for(j=0;j<6;j++) cout<<x[j]<<" ";
      //cout<<endl;
    }

    //print out turn-by-turn data
    //for(j=0;j<6;j++) cout<<x[j]<<" ";
    //cout<<endl;
  }
  
}

//----track fast and sums of x^2, y^2, z^2
void Track_Fast_Emit(Line & linename, double x[6], int nturn, int & stable, int & lost_turn, int & lost_post, double & sum_x2, double & sum_y2, double & sum_z2)
{
  int i,j;
  double  knl[11], knsl[11];
  double  m66[36], xco_in[6], xco_out[6];

  sum_x2=0.;  sum_y2=0.;  sum_z2=0.;

  //-----quick check before track
  if( abs(x[0]) > 0.1 ||  abs(x[2]) > 0.1 || stable ==0 ) {
    stable = 0;      lost_turn= 0;       lost_post = 0;      return ;   }
  
  //---now we do track
  for(GP.turn=0;GP.turn<nturn; GP.turn++) {
    
    Ndrift=0;
    Nbend=0;
    Nquad=0;
    Nsext=0;
    Nmult=0;
    Ngmult=0;
    Nmatrix=0;     
    Nrf=0;  
    Nhobb=0;
    Nelens=0;
    Nrotat=0;
    Nsnake=0;
    Nerhicbb=0;
    Ncrabrf=0;
    
    for(i=0; i<nelement; i++){
      
      //--alignmnet error: can be used for each element 
      //GtoL(x, nelement_dx[i],nelement_dy[i],nelement_dt[i]); 
      
      if(type[i]==0 ) {
	DRIFT_Pass(x, drift_p[Ndrift]);
	Ndrift++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2]) ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==1 ) {
	SBEND_Pass(x,bend_p[Nbend*5+0],int(bend_p[Nbend*5+4]), bend_p[Nbend*5+1], bend_p[Nbend*5+2],bend_p[Nbend*5+3]); 
	Nbend++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2]) ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==2 ){
	QUAD_Pass(x, quad_p[Nquad*4+0], int(quad_p[Nquad*4+3]), quad_p[Nquad*4+1], quad_p[Nquad*4+2]); 
	Nquad++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2]) ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==3 ){
	SEXT_Pass(x, sext_p[Nsext*4+0], int(sext_p[Nsext*4+3]), sext_p[Nsext*4+1],  sext_p[Nsext*4+2]);
	Nsext++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2]) ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==4 ){
	for(j=0;j<11;j++) {
	  knl[j]=mult_p[Nmult*25+1+j]; 
	  knsl[j]=mult_p[Nmult*25+1+11+j];
	}
	MULT_Pass(x, mult_p[Nmult*25+0], int(mult_p[Nmult*25+23]), int(mult_p[Nmult*25+24]), knl, knsl);
	Nmult++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2]) ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==5 ){
	for(j=0;j<11;j++) {
	  knl[j]=gmult_p[Ngmult*28+4+j]; 
	  knsl[j]=gmult_p[Ngmult*28+4+11+j];
	}
	GMULT_Pass(x, gmult_p[Ngmult*28+0], int(gmult_p[Ngmult*28+26]),int(gmult_p[Ngmult*28+27]), gmult_p[Ngmult*28+1] ,gmult_p[Ngmult*28+2], gmult_p[Ngmult*28+3], knl, knsl);
	Ngmult++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2]) ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==6 ){
	for(j=0;j<6;j++) xco_in[j] = matrix_p[Nmatrix*49+1+j];
	for(j=0;j<6;j++) xco_out[j]= matrix_p[Nmatrix*49+1+6+j];
	for(j=0;j<36;j++) {
	  m66[j]=matrix_p[Nmatrix*49+13+j]; 
	}
	MATRIX_Pass(x, matrix_p[Nmatrix*49+0],xco_in, xco_out, m66);
	Nmatrix++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2]) ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if (type[i]==7 ){
        RFCAV_Pass(x, rf_p[Nrf*4+0], rf_p[Nrf*4+1], rf_p[Nrf*4+2],rf_p[Nrf*4+3] );
	Nrf++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2]) ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==8 ){
	BEAMBEAM_Pass(x, int(hobb_p[Nhobb*11+0]), hobb_p[Nhobb*11+1], hobb_p[Nhobb*11+2], hobb_p[Nhobb*11+3], int(hobb_p[Nhobb*11+4]), 
		      hobb_p[Nhobb*11+5], hobb_p[Nhobb*11+6], hobb_p[Nhobb*11+7], hobb_p[Nhobb*11+8], hobb_p[Nhobb*11+9],hobb_p[Nhobb*11+10] );
	Nhobb++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2]) ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==9 ){
        ELENS_Pass(x, elens_p[Nelens*7+0], elens_p[Nelens*7+1], elens_p[Nelens*7+2],
		   elens_p[Nelens*7+3], int(elens_p[Nelens*7+4]), elens_p[Nelens*7+5], elens_p[Nelens*7+6]);
	Nelens++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2]) ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==10 ){
        DRIFT_Pass(x,rotat_p[Nrotat*5+0]);
	Nrotat++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2]) ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==11 ){
        DRIFT_Pass(x,snake_p[Nsnake*5+0]);
	Nsnake++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2]) ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==12 ){
	ERHICBB_Pass(x, GP.gamma, erhicbb_p[Nerhicbb*4 +0], erhicbb_p[Nerhicbb*4 +1]);
        Nerhicbb++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2]) ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==13 ){
	CRABRF_Pass(x,crabrf_p[Ncrabrf*4 +0],crabrf_p[Ncrabrf*4 +1],crabrf_p[Ncrabrf*4 +2],crabrf_p[Ncrabrf*4 +3]);
        Ncrabrf++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2]) ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else{ }
      
      //LtoG(x, nelement_dx[i],nelement_dy[i],nelement_dt[i]); 
      
      //print out each element 
      //cout<<i<<"   ";
      //for(j=0;j<6;j++) cout<<x[j]<<" ";
      //cout<<endl;
    }
    sum_x2 += x[0]*x[0];
    sum_y2 += x[2]*x[2]; 
    sum_z2 += x[4]*x[4]; 
    //print out turn-by-turn data
    //for(j=0;j<6;j++) cout<<x[j]<<" ";
    //cout<<endl;
  }

}

void Track_Fast_DA( Line & linename, int nturn, double deltap0, double sigmax0, double sigmay0, char * filename ) // Track without class, twin particles
{
  int i;
  double angle, sigma_s, sigma_e, sigma, sigma1, sigma2, sigma_max=100.0;  
  int    stable, stable1=1, stable2=1,  lost_turn=0, lost_post=0; 
  FILE *f2;

  f2=fopen(filename,"w");
  fclose(f2);
  double x[6];

  for(i=0; i<10;i++){

    angle=( 90/11.)*(i+1)* PI/180;
    sigma_s=0.;
    sigma_e=sigma_max;
    
    do {
      sigma = (sigma_s + sigma_e )/2.0;
      sigma1= sigma;
      sigma2= sigma + 0.05;

      x[0]= sigma1 * cos(angle) * sigmax0;
      x[1]= 0. ;
      x[2]= sigma1 * sin(angle) * sigmay0;
      x[3]= 0. ;
      x[4]= 0.;
      x[5]= deltap0;

      stable1=1; lost_turn=0; lost_post=0;
      Track_Fast(linename, x, nturn, stable1, lost_turn, lost_post);
      
      x[0]= sigma2 * cos(angle) * sigmax0;
      x[1]= 0. ;
      x[2]= sigma2 * sin(angle) * sigmay0;
      x[3]= 0. ;
      x[4]= 0.;
      x[5]= deltap0;

      stable2=1; lost_turn=0; lost_post=0;
      Track_Fast(linename, x, nturn, stable2, lost_turn, lost_post);
      
      if ( stable1* stable2 ==1 ) {
	stable =1 ; }
      else{
	stable=0;
      }
      if( stable ==0 ) 	sigma_e= sigma;
      if( stable ==1 ) 	sigma_s= sigma;
    } while ( abs( sigma_e- sigma_s) > 0.1  );

    f2=fopen(filename,"a");
    fprintf(f2,"%f %f \n", ( 90/11.)*(i+1),  sigma);
    fclose(f2);
  }
} 

void Track_Fast_DA1( Line & linename, int nturn, double deltap0, double sigmax0, double sigmay0, char* filename ) // Track without class, twin particles
{
  int i;
  double angle, sigma_s, sigma_e, sigma, sigma1, sigma2, sigma_max=100.0;  
  int    stable, stable1=1, stable2=1,  lost_turn=0, lost_post=0; 
  FILE *f2;

  f2=fopen(filename,"w");
  fclose(f2);
  double x[6];
  
  Cal_Orbit_Num(linename, 0);

  for(i=0; i<5;i++){

    angle=15*(i+1)* PI/180;
    sigma_s=0.;
    sigma_e=sigma_max;
    
    do {
      sigma = (sigma_s + sigma_e )/2.0;
      sigma1= sigma;
      sigma2= sigma + 0.05;

      x[0]= sigma1 * cos(angle) * sigmax0;
      x[1]= 0. ;
      x[2]= sigma1 * sin(angle) * sigmay0;
      x[3]= 0. ;
      x[4]= 0.;
      x[5]= deltap0;

      stable1=1; lost_turn=0; lost_post=0;
      Track_Fast(linename, x, nturn, stable1, lost_turn, lost_post);
      
      x[0]= sigma2 * cos(angle) * sigmax0;
      x[1]= 0. ;
      x[2]= sigma2 * sin(angle) * sigmay0;
      x[3]= 0. ;
      x[4]= 0.;
      x[5]= deltap0;

      stable2=1; lost_turn=0; lost_post=0;
      Track_Fast(linename, x, nturn, stable2, lost_turn, lost_post);
      
      if ( stable1* stable2 ==1 ) {
	stable =1 ; }
      else{
	stable=0;
      }
      if( stable ==0 ) 	sigma_e= sigma;
      if( stable ==1 ) 	sigma_s= sigma;
    } while ( abs( sigma_e- sigma_s) > 0.1  );

    f2=fopen(filename,"a");
    fprintf(f2,"%f %f \n", 15.*(i+1),  sigma);
    fclose(f2);
  }
} 

void Track_Fast_DA2( Line & linename, int nturn, double deltap0, double sigmax0, double sigmay0, char* filename ) // Track without class, single particle
{
  int i;
  double angle, sigma_s, sigma_e, sigma, sigma_max=100.0;  
  int    stable, lost_turn=0, lost_post=0; 
  FILE *f2;

  f2=fopen(filename,"w");
  fclose(f2);
  double x[6];

  Cal_Orbit_Num(linename, 0);

  for(i=0; i<5;i++){

    angle=15*(i+1)* PI/180;
    sigma_s=0.;
    sigma_e=sigma_max;
    
    do {
      sigma = (sigma_s + sigma_e )/2.0;
      x[0]= sigma * cos(angle) * sigmax0;
      x[1]= 0. ;
      x[2]= sigma * sin(angle) * sigmay0;
      x[3]= 0. ;
      x[4]= 0.;
      x[5]= deltap0;

      stable=1; lost_turn=0; lost_post=0;
      Track_Fast(linename, x, nturn, stable, lost_turn, lost_post);
      
      if( stable ==0 ) 	sigma_e= sigma;
      if( stable ==1 ) 	sigma_s= sigma;
    } while ( abs( sigma_e- sigma_s) > 0.1  );

    f2=fopen(filename,"a");
    fprintf(f2,"%f %f \n", 15.*(i+1),  sigma);
    fclose(f2);
  }
} 

void Track_Fast_Spin(Line & linename, double x[], int nturn, int & stable, int & lost_turn, int & lost_post)
{
  int i,j;
  double  knl[11], knsl[11];
  double  m66[36], xco_in[6], xco_out[6];
  double  axis[3];

  //-----quick check before track
  if( abs(x[0]) > 0.1 ||  abs(x[2]) > 0.1 || stable ==0 ) {
    stable = 0;      lost_turn= 0;       lost_post = 0;      return ;   }
  
  //---now we do track
  for(GP.turn=0;GP.turn<nturn; GP.turn++) {
    
    Ndrift=0;
    Nbend=0;
    Nquad=0;
    Nsext=0;
    Nmult=0;
    Ngmult=0;
    Nmatrix=0;     
    Nrf=0;  
    Nhobb=0;
    Nelens=0;
    Nrotat=0;
    Nsnake=0;
    Nerhicbb=0;
    Ncrabrf=0;
    
    for(i=0; i<nelement; i++){
      
      //--alignmnet error: can be used for each element 
      //GtoL(x, nelement_dx[i],nelement_dy[i],nelement_dt[i]); 
      
      if(type[i]==0 ) {
	DRIFT_Pass(x, drift_p[Ndrift]);
	Ndrift++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2]) ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==1 ) {
	SBEND_sPass(x,bend_p[Nbend*5+0],int(bend_p[Nbend*5+4]), bend_p[Nbend*5+1], bend_p[Nbend*5+2],bend_p[Nbend*5+3]); 
	Nbend++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2]) ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==2 ){
	QUAD_sPass(x, quad_p[Nquad*4+0], int(quad_p[Nquad*4+3]), quad_p[Nquad*4+1], quad_p[Nquad*4+2]); 
	Nquad++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2])  ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==3 ){
	SEXT_sPass(x, sext_p[Nsext*4+0], int(sext_p[Nsext*4+3]), sext_p[Nsext*4+1],  sext_p[Nsext*4+2]);
	Nsext++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2])  ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==4 ){
	for(j=0;j<11;j++) {
	  knl[j]=mult_p[Nmult*25+1+j]; 
	  knsl[j]=mult_p[Nmult*25+1+11+j];
	}
	MULT_sPass(x, mult_p[Nmult*25+0], int(mult_p[Nmult*25+23]), int(mult_p[Nmult*25+24]), knl, knsl);
	Nmult++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2]) ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==5 ){
	for(j=0;j<11;j++) {
	  knl[j]=gmult_p[Nmult*28+4+j]; 
	  knsl[j]=gmult_p[Nmult*28+4+11+j];
	}
	GMULT_sPass(x, gmult_p[Nmult*28+0], int(gmult_p[Nmult*28+26]),int(gmult_p[Nmult*28+27]), gmult_p[Nmult*28+1] ,gmult_p[Nmult*28+2], gmult_p[Nmult*28+3], knl, knsl);
	Ngmult++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2]) ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==6 ){
	for(j=0;j<6;j++) xco_in[j] = matrix_p[Nmatrix*49+1+j];
	for(j=0;j<6;j++) xco_out[j]= matrix_p[Nmatrix*49+1+6+j];
	for(j=0;j<36;j++) {
	  m66[j]=matrix_p[Nmatrix*49+13+j]; 
	}
	MATRIX_Pass(x, matrix_p[Nmatrix*49+0],xco_in, xco_out, m66);
	Nmatrix++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2]) ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if (type[i]==7 ){
        RFCAV_Pass(x, rf_p[Nrf*4+0], rf_p[Nrf*4+1], rf_p[Nrf*4+2],rf_p[Nrf*4+3] );
	Nrf++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2]) ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==8 ){
	BEAMBEAM_Pass(x, int(hobb_p[Nhobb*11+0]), hobb_p[Nhobb*11+1], hobb_p[Nhobb*11+2], int(hobb_p[Nhobb*11+3]), hobb_p[Nhobb*11+4], 
		      hobb_p[Nhobb*11+5], hobb_p[Nhobb*11+6], hobb_p[Nhobb*11+7], hobb_p[Nhobb*11+8], hobb_p[Nhobb*11+9],hobb_p[Nhobb*11+10] );
	Nhobb++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2]) ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==9 ){
	ELENS_Pass(x, elens_p[Nelens*7+0], elens_p[Nelens*7+1], elens_p[Nelens*7+2],
		   elens_p[Nelens*7+3], int(elens_p[Nelens*7+4]), elens_p[Nelens*7+5], elens_p[Nelens*7+6]);
	Nelens++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2])  ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==10 ){
        for(j=0;j<3;j++) axis[j]=rotat_p[Nrotat*5+1+j];
        ROTAT_sPass(x,rotat_p[Nrotat*5+0],axis,rotat_p[Nrotat*5+4]);
	Nrotat++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2]) ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==10 ){
        for(j=0;j<3;j++) axis[j]=snake_p[Nsnake*5+1+j];
        SNAKE_sPass(x,snake_p[Nsnake*5+0],axis,rotat_p[Nsnake*5+4]);
	Nsnake++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2])  ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==12 ){
	ERHICBB_Pass(x, GP.gamma, erhicbb_p[Nerhicbb*4 +0], erhicbb_p[Nerhicbb*4 +1]);
        Nerhicbb++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2]) ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==13 ){
	CRABRF_Pass(x,crabrf_p[Ncrabrf*4 +0],crabrf_p[Ncrabrf*4 +1],crabrf_p[Ncrabrf*4 +2],crabrf_p[Ncrabrf*4 +3]);
        Ncrabrf++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] or  isnan(x[0]) or isnan(x[2]) ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else{ }
      
      //LtoG(x, nelement_dx[i],nelement_dy[i],nelement_dt[i]); 
      
      //print out each element 
      // cout<<i<<"   ";
      // for(j=0;j<9;j++) cout<<x[j]<<" ";
      // cout<<endl;
    }
    
    // print out turn-by-turn data
    //for(j=0;j<9;j++) cout<<x[j]<<" ";
    //cout<<endl;
  }
  
}

--*/

void BeamBeam4D_SoftGaussian(double Np1, int Npart_1, double PartDist_1[], 
                             double Np2, int Npart_2, double PartDist_2[], 
                             double gamma1, double gamma2, double bbscale1, double bbscale2)
{
  int i, j;
  double mean_1[2],rms_1[2],mean_2[2],rms_2[2];
  double dx, dy, dpx, dpy;
  
  Cal_Mean_RMS(Npart_1, PartDist_1,  mean_1, rms_1);
  Cal_Mean_RMS(Npart_2, PartDist_2,  mean_2, rms_2);

  for(i=0;i<Npart_1;i++){
    dx= PartDist_1[6*i+0] - mean_2[0];
    dy= PartDist_1[6*i+2] - mean_2[1];
    BB4D(dx, dy, gamma1, Np2, rms_2[0], rms_2[1], dpx,  dpy);
    PartDist_1[6*i+1] +=  dpx * bbscale1; 
    PartDist_1[6*i+3] +=  dpy * bbscale1;
  } 

  for(i=0;i<Npart_2;i++){
    dx= PartDist_2[6*i+0] - mean_1[0];
    dy= PartDist_2[6*i+2] - mean_1[1];
    BB4D(dx, dy, gamma2, Np1, rms_1[0], rms_1[1], dpx,  dpy);
    PartDist_2[6*i+1] +=  dpx * bbscale2; 
    PartDist_2[6*i+3] +=  dpy * bbscale2;
  }
}

void BeamBeam6D_SoftGaussian(double Np1, int Mp1, double x1[], double Np2, int Mp2, double x2[], 
                             double gamma1, double gamma2, double bbscale1, double bbscale2,
			     int Nslice1, int Nslice2)
//  bunch slicing will be done in this function
{
  int    i,j,k;

  int    Nslice=Nslice1 * Nslice2;
  int    Index1[Nslice], Index2[Nslice];
  double Z1[Nslice1], Z2[Nslice2], S[Nslice],  Z10[Nslice];
  double temp[6], temp1;
  int    istart, iend, itemp;

  int    Mp1_coll=Mp1/Nslice1, Mp2_coll=Mp2/Nslice2;
  double xc1[6*Mp1_coll], xc2[6*Mp2_coll];
  double mean1[2], rms1[2], mean2[2], rms2[2], emit1[2], emit2[2], beta2[2], alfa2[2], gama2[2],  beta1[2], alfa1[2], gama1[2];
  double dx, dy, dpx, dpy;

  char   filename[128];
  fstream f1;

  //----sorting particles in each  bunch
  //cout<<"sorting particles, begin. "<<endl;
  for(i=0;i<Mp1;i++){
    for(j=i+1;j<Mp1;j++){
      if(x1[j*6+4]>x1[i*6+4]){
	for(k=0;k<6;k++) temp[k]   = x1[j*6+k];
        for(k=0;k<6;k++) x1[j*6+k] = x1[i*6+k];
	for(k=0;k<6;k++) x1[i*6+k] = temp[k];
      }
    }
  }

  for(i=0;i<Mp2;i++){
    for(j=i+1;j<Mp2;j++){
      if(x2[j*6+4]>x2[i*6+4]){
	for(k=0;k<6;k++) temp[k]   = x2[j*6+k];
        for(k=0;k<6;k++) x2[j*6+k] = x2[i*6+k];
	for(k=0;k<6;k++) x2[i*6+k] = temp[k];
      }
    }
  }
  //cout<<"sorting particles, end. "<<endl;
  
  /*for checking purpose*/
  if(false){ 
    
    f1.open("./bunch1.dat",ios::out);
    for(j=0;j<Mp1;j++) {
      for(k=0;k<6;k++) f1<<x1[j*6+k]<<" ";
      f1<<endl;
    }
    f1.close(); 
    for(i=0;i<Nslice1;i++){
      sprintf(filename, "./bunch1_slice%d.dat",i);
      f1.open(filename,ios::out);
      istart=Mp1_coll*i;
      iend=Mp1_coll*(i+1);
      for(j=istart;j<iend;j++) {
	for(k=0;k<6;k++) f1<<x1[j*6+k]<<" ";
	f1<<endl;
      }
      f1.close();
    }

    f1.open("./bunch2.dat",ios::out);
    for(j=0;j<Mp2;j++) {
      for(k=0;k<6;k++) f1<<x2[j*6+k]<<" ";
      f1<<endl;
    }
    f1.close(); 
    for(i=0;i<Nslice2;i++){
      sprintf(filename, "./bunch2_slice%d.dat",i);
      f1.open(filename,ios::out);
      istart=Mp2_coll*i;
      iend=Mp2_coll*(i+1);
      for(j=istart;j<iend;j++) {
	for(k=0;k<6;k++) f1<<x2[j*6+k]<<" ";
	f1<<endl;
      }
      f1.close();
    }

  }
  
  /*----check purpose----*/
  if(false){
    f1.open("./bunch1_stat.dat",ios::out);
    for(i=0;i<Nslice1;i++){
      for(j=0;j<Mp1_coll;j++){
	istart=Mp1_coll * i ;
	for(k=0;k<6;k++)  xc1[j*6+k] = x1[(istart+j)*6+k]; 
      }
      Cal_Mean_RMS(Mp1_coll, xc1, mean1, rms1);
      f1<<i<<"  "<<mean1[0]<<"  "<<mean1[1]<<"  "<<rms1[0]<<"  "<<rms1[1]<<endl;
    }
    f1.close();
    
    f1.open("./bunch2_stat.dat",ios::out);
    for(i=0;i<Nslice2;i++){
      for(j=0;j<Mp2_coll;j++){
	istart=Mp2_coll * i ;
	for(k=0;k<6;k++)  xc2[j*6+k] = x2[(istart+j)*6+k]; 
      }
      Cal_Mean_RMS(Mp2_coll, xc2, mean1, rms1);
      f1<<i<<"  "<<mean1[0]<<"  "<<mean1[1]<<"  "<<rms1[0]<<"  "<<rms1[1]<<endl;
    }
    f1.close();
  }
  
  //---sorting slices for collision

  for(i=0;i<Nslice1;i++){
    temp1=0;
    istart=i*Mp1_coll;  
    iend  =(i+1)*Mp1_coll;
    for(j=istart;j<iend;j++){
      temp1 =temp1 + x1[6*j+4]; 
    }
    Z1[i]=temp1 /Mp1_coll;
  }

  for(i=0;i<Nslice2;i++){
    temp1=0;
    istart=i*Mp2_coll;  
    iend  =(i+1)*Mp2_coll;
    for(j=istart;j<iend;j++){
      temp1 =temp1 + x2[6*j+4]; 
    }
    Z2[i]=temp1 /Mp2_coll;
  }

  /*---check purpose----*/
  if(true){ 
    f1.open("./bunch1_Z.dat",ios::out);
    for(j=0;j<Nslice1;j++) {
      f1<<j<<"  "<<Z1[j]<<endl;
    }
    f1.close(); 
    f1.open("./bunch2_Z.dat",ios::out);
    for(j=0;j<Nslice2;j++) {
      f1<<j<<"  "<<Z2[j]<<endl;
    }
    f1.close(); 
  }

  for(i=0;i<Nslice1;i++){
    for(j=0;j<Nslice2;j++){
      Index1[i*Nslice2 + j] = i;
      Index2[i*Nslice2 + j] = j;
      S[i*Nslice2 + j]= (Z1[i] -Z2[j])/2.;
      Z10[i*Nslice2 + j]= -(Z1[i] +Z2[j])/2.;
    }
  }

  for(i=0;i<Nslice;i++){
    for(j=i+1;j<Nslice;j++){
      if(Z10[j] < Z10[i] ){ 
	temp1 =  Z10[j];  Z10[j] = Z10[i];  Z10[i]=temp1;
	
	temp1 =  S[j]  ;  S[j]   = S[i]  ;  S[i]  =temp1;
        itemp =  Index1[j];  Index1[j] = Index1[i];  Index1[i] =itemp;
        itemp =  Index2[j];  Index2[j] = Index2[i];  Index2[i] =itemp; 
      }
    }
  }

  /*----check purpose----*/
  if(false){
    f1.open("./collide_order.dat",ios::out);
    for(i=0;i<Nslice;i++){
      f1<<i<<"   "<<Index1[i]<<"  "<<Index2[i]<<"   "<<Z1[Index1[i]]<<"   "<<Z2[Index2[i]]<<"  "<<S[i]<<"   "<<Z10[i]<<endl;
    }
    f1.close();
  }

  //-----collide slice to slice : be careful both beams should have same x^, y^ directions

  for(i=0;i<Nslice;i++){

    istart= Index1[i] * Mp1_coll;
    for(j=0;j<Mp1_coll;j++){
      for(k=0;k<6;k++)  xc1[j*6+k] = x1[(istart+j)*6+k]; 
    }
    istart= Index2[i] * Mp2_coll;
    for(j=0;j<Mp2_coll;j++){
      for(k=0;k<6;k++)  xc2[j*6+k] = x2[(istart+j)*6+k]; 
    }
    
    for(j=0;j<Mp1_coll;j++){
      xc1[j*6+0] = xc1[j*6+0]  +  xc1[j*6+1] * S[i];
      xc1[j*6+2] = xc1[j*6+2]  +  xc1[j*6+3] * S[i];
    }
    for(j=0;j<Mp2_coll;j++){
      xc2[j*6+0] = xc2[j*6+0]  -  xc2[j*6+1] * S[i];
      xc2[j*6+2] = xc2[j*6+2]  -  xc2[j*6+3] * S[i];
    }
    
    Cal_Mean_RMS(Mp1_coll, xc1, mean1, rms1);
    Cal_Mean_RMS(Mp2_coll, xc2, mean2, rms2);
    Cal_Emit(Mp1_coll, xc1, emit1, beta1, alfa1, gama1);
    Cal_Emit(Mp2_coll, xc2, emit2, beta2, alfa2, gama2);
      
    if(false){
      if( Index2[i] == 2 ) {
	//cout<<S[i]<<"  "<<mean1[0]<<"  "<< mean1[1]<<"  "<<rms1[0]<<"  "<<rms1[1]<<"  "<<emit1[0]<<"   "<<emit1[1]<<endl;
        cout<<S[i]<<"  "<<mean2[0]<<"  "<< mean2[1]<<"  "<<rms2[0]<<"  "<<rms2[1]<<"  "<<emit2[0]<<"   "<<emit2[1]<<endl;
      }
    }

    for(j=0;j<Mp1_coll;j++){
      dx= xc1[6*j+0] - mean2[0];
      dy= xc1[6*j+2] - mean2[1];
      BB4D(dx, dy, gamma1, Np2 / Nslice2, rms2[0], rms2[1], dpx,  dpy);
      xc1[6*j+1] = xc1[6*j+1] + dpx * bbscale1; 
      xc1[6*j+3] = xc1[6*j+3] + dpy * bbscale1; 
    }
    for(j=0;j<Mp2_coll;j++){
      dx= xc2[6*j+0] - mean1[0];
      dy= xc2[6*j+2] - mean1[1];
      BB4D(dx, dy, gamma2, Np1 / Nslice1, rms1[0], rms1[1], dpx,  dpy);
      xc2[6*j+1] = xc2[6*j+1] + dpx * bbscale2; 
      xc2[6*j+3] = xc2[6*j+3] + dpy * bbscale2; 
    }


    for(j=0;j<Mp1_coll;j++){
      xc1[j*6+0] = xc1[j*6+0]  -  xc1[j*6+1] * S[i];
      xc1[j*6+2] = xc1[j*6+2]  -  xc1[j*6+3] * S[i];
    }
    for(j=0;j<Mp2_coll;j++){
      xc2[j*6+0] = xc2[j*6+0]  +  xc2[j*6+1] * S[i];
      xc2[j*6+2] = xc2[j*6+2]  +  xc2[j*6+3] * S[i];
    }

    Cal_Mean_RMS(Mp1_coll, xc1, mean1, rms1);
    Cal_Mean_RMS(Mp2_coll, xc2, mean2, rms2);
    Cal_Emit(Mp1_coll, xc1, emit1, beta1, alfa1, gama1);
    Cal_Emit(Mp2_coll, xc2, emit2, beta2, alfa2, gama2);
    
    if(false){
      if( Index2[i] == 2 ) {
	//cout<<S[i]<<"  "<<mean1[0]<<"  "<< mean1[1]<<"  "<<rms1[0]<<"  "<<rms1[1]<<"  "<<emit1[0]<<"   "<<emit1[1]<<endl;
        cout<<S[i]<<"  "<<mean2[0]<<"  "<< mean2[1]<<"  "<<rms2[0]<<"  "<<rms2[1]<<"  "<<emit2[0]<<"   "<<emit2[1]<<endl;
      }
    }

    istart= Index1[i] * Mp1_coll;
    for(j=0;j<Mp1_coll;j++){
      for(k=0;k<6;k++) x1[(istart+j)*6+k] =  xc1[j*6+k]; 
    }
    istart= Index2[i] * Mp2_coll;
    for(j=0;j<Mp2_coll;j++){
      for(k=0;k<6;k++)  x2[(istart+j)*6+k] = xc2[j*6+k]; 
    }

  }
  
}

void BeamBeam6D_SoftGaussian(double Np1, int Mp1, double x1[], double Np2, int Mp2, double x2[], 
			     double gamma1, double gamma2, double bbscale1, double bbscale2,
			     int Nslice1, int Nslice2, double zb1[], double zb2[])
//----bondary of bunch slices are pre-defined with zb1[], zb2[]
{
  int    i, j, k, i1;

  double Z1_l[Nslice1],  Z2_l[Nslice2];
  double Z1_g[Nslice1],  Z2_g[Nslice2];
  int    Mp1_l[Nslice1], Mp2_l[Nslice2]; 
  int    Nslice=Nslice1 * Nslice2;
  int    itemp, Index1_g[Nslice], Index2_g[Nslice];
  double S_g[Nslice], Z10_g[Nslice];

  int    ncount;
  int    Mp1_coll, Mp2_coll;
  int    dimxc1= int(2*Mp1/Nslice1 ), dimxc2= int(2*Mp2/Nslice2);
  double *xc1 =new double [dimxc1*6], *xc2 = new double [dimxc2*6];
  double sum[2], mean1[2], rms1[2], mean2[2], rms2[2];
  
  double  temp;
  char    filename[128];
  fstream f1;

  //---- slicing bunches--------
  
  for(i=0;i<Nslice1;i++) {
    Mp1_l[i]=0 ;     Mp2_l[i]=0 ;  
    Z1_l[i] =0.;     Z2_l[i] =0.;  
    Z1_g[i] =0.;     Z2_g[i] =0.;
  }    
  
  for(i=0;i<Nslice1;i++){
    temp=0; 
    for(j=0;j<Mp1;j++){
      if( x1[j*6+4] < zb1[i]  and  x1[j*6+4] >= zb1[i+1] ) {
	temp =temp + x1[6*j+4];
	Mp1_l[i]++;
      }
    }
    Z1_l[i]=temp;
  }
  
  for(i=0;i<Nslice2;i++){
    temp=0; 
    for(j=0;j<Mp2;j++){
      if( x2[j*6+4] < zb2[i]  and  x2[j*6+4] >= zb2[i+1] ) {
	temp =temp + x2[6*j+4]; 
	Mp2_l[i]++;
      }
    }
    Z2_l[i]=temp;
  }
  
  for(i=0;i<Nslice1;i++) Z1_g[i]= Z1_l[i]/Mp1_l[i];
  for(i=0;i<Nslice2;i++) Z2_g[i]= Z2_l[i]/Mp2_l[i];
  
  //------set up collision-----
  
  for(i=0;i<Nslice1;i++){
    for(j=0;j<Nslice2;j++){
      Index1_g[i*Nslice2 + j] = i;
      Index2_g[i*Nslice2 + j] = j;
      S_g[i*Nslice2 + j]= (Z1_g[i] - Z2_g[j])/2.;
      Z10_g[i*Nslice2 + j]= -(Z1_g[i] + Z2_g[j])/2.;
    }
  }
  
  for(i=0;i<Nslice;i++){
    for(j=i+1;j<Nslice;j++){
	if(Z10_g[j] < Z10_g[i] ){ 
	  temp =  Z10_g[j];  Z10_g[j] = Z10_g[i];  Z10_g[i]=temp;
	  temp =  S_g[j]  ;  S_g[j]   = S_g[i]  ;  S_g[i]  =temp;
	  itemp =  Index1_g[j];  Index1_g[j] = Index1_g[i];  Index1_g[i] =itemp;
	  itemp =  Index2_g[j];  Index2_g[j] = Index2_g[i];  Index2_g[i] =itemp; 
	}
    }
  }

  //------ colliding slice to slice---
  
    for(i=0;i<Nslice;i++){

      /*----pick out----*/
      
      ncount=0;
      for(j=0;j<Mp1;j++){
        if( x1[j*6+4] < zb1[ Index1_g[i] ]  and  x1[j*6+4] >= zb1[ Index1_g[i]+1 ]  ) {
	  for(k=0;k<6;k++) xc1[ncount*6+k]  =   x1[j*6+k];
	  ncount++;
	}
      }
      Mp1_coll = ncount;
      
      ncount=0;
      for(j=0;j<Mp2;j++){
        if( x2[j*6+4] < zb2[ Index2_g[i] ]  and  x2[j*6+4] >= zb2[  Index2_g[i]+1 ]  ) {
	  for(k=0;k<6;k++) xc2[ncount*6+k]  =   x2[j*6+k];
	  ncount++;
	}
      }
      Mp2_coll = ncount;
      
      if(Mp1_coll > 2*Mp1/Nslice1 or Mp1_coll < Mp1/Nslice1/2 or  Mp2_coll > 2*Mp2/Nslice2 or Mp2_coll < Mp2/Nslice2/2 ) {
	cout<<"Unbalanced population for predefined slices. exit. "<<endl;
        cout<<endl;
      }

      /*-----beam-beam interaction-------*/
      
      for(j=0;j<Mp1_coll;j++){
	xc1[j*6+0] = xc1[j*6+0]  +  xc1[j*6+1] * S_g[i];
	xc1[j*6+2] = xc1[j*6+2]  +  xc1[j*6+3] * S_g[i];
      }
      for(j=0;j<Mp2_coll;j++){
	xc2[j*6+0] = xc2[j*6+0]  -  xc2[j*6+1] * S_g[i];
	xc2[j*6+2] = xc2[j*6+2]  -  xc2[j*6+3] * S_g[i];
      }
      BeamBeam4D_SoftGaussian(Np1*Mp1_coll/Mp1, Mp1_coll, xc1, Np2*Mp2_coll/Mp2, Mp2_coll, xc2, gamma1, gamma2, bbscale1, bbscale2);
      for(j=0;j<Mp1_coll;j++){
	xc1[j*6+0] = xc1[j*6+0]  -  xc1[j*6+1] * S_g[i];
	xc1[j*6+2] = xc1[j*6+2]  -  xc1[j*6+3] * S_g[i];
      }
      for(j=0;j<Mp2_coll;j++){
	xc2[j*6+0] = xc2[j*6+0]  +  xc2[j*6+1] * S_g[i];
	xc2[j*6+2] = xc2[j*6+2]  +  xc2[j*6+3] * S_g[i];
      }

      /*---stack back------*/
      
      ncount=0;
      for(j=0;j<Mp1;j++){
	if( x1[j*6+4] < zb1[ Index1_g[i] ]  and  x1[j*6+4] >= zb1[  Index1_g[i]+1 ]  ) {
	  for(k=0;k<6;k++) x1[j*6+k] = xc1[ncount*6+k];
	  ncount++;
	}
      }
      
      ncount=0;
      for(j=0;j<Mp2;j++){
	if( x2[j*6+4] < zb2[ Index2_g[i] ]  and  x2[j*6+4] >= zb2[  Index2_g[i]+1 ]  ) {
	  for(k=0;k<6;k++) x2[j*6+k] = xc2[ncount*6+k]; 
	  ncount++;
	}
      }
      
    }
    
    delete [] xc1;
    delete [] xc2;   
}

void BeamBeam4D_PIC(double Np1, int Mp1, double xc1[], double Np2, int Mp2, double xc2[], 
                    double gamma1, double gamma2, double bbscale1, double bbscale2)
//----solve Posisson Equation with Green Function + FFT
{
  int i, j, i1, j1;
  int     Nx=64, Ny=64;
  int     nsigma=15;
  double  x0, y0,  hx, hy, rmsx, rmsy, xoffset, yoffset;
  double  xp,yp, ex,ey;
  double  src[Nx*Ny], phi[Nx*Ny], Ex[Nx*Ny], Ey[Nx*Ny];
  double  grn_c_r[2*Nx*2*Ny], grn_c_i[2*Nx*2*Ny];

  //---set particles frame
  double mean1[2], rms1[2],  mean2[2], rms2[2];
  
  Cal_Mean_RMS(Mp1, xc1, mean1, rms1);
  Cal_Mean_RMS(Mp2, xc2, mean2, rms2);
  
  if(rms1[0]  > rms2[0]) {
    rmsx=rms2[0];
  }
  else{
    rmsx=rms1[0];
  }
  if(rms1[1]  > rms2[1]) {
    rmsy=rms2[1];
  }
  else{
    rmsy=rms1[1];
  }
  
  hx= nsigma * rmsx / Nx ;   hy= nsigma * rmsy / Ny ;

  //----slice 1  to slice 2:
  x0=mean1[0] - Nx * hx / 2 ;
  y0=mean1[1] - Ny * hy / 2 ;
  xoffset= mean2[0] - mean1[0];
  yoffset= mean2[1] - mean1[1];
  
  for(i1=0; i1<Nx*Ny;i1++){
    src[i1]=0;   phi[i1]=0;
  }
  
  for(i1=0;i1<Mp1;i1++){
    xp = xc1[6*i1+0];  yp = xc1[6*i1+2];
    if( abs(xp-mean1[0]) < (Nx/2 - 3)*hx  and   abs(yp-mean1[1]) < (Ny/2 -3 )*hy ) {
      Charge2Grid(xp,yp, x0, y0, Nx, Ny, hx,  hy, src);
    }
  }
  Cal_Effective_Green_Function(xoffset, yoffset, Nx, Ny, hx, hy, grn_c_r, grn_c_i);
  Cal_Potential_Convolution(Nx, Ny, hx, hy, src, grn_c_r, grn_c_i, phi);
  Cal_Electric_Field(Nx, Ny, hx, hy, phi, Ex, Ey);
  
  for(i1=0;i1<Mp2;i1++){
    xp=xc2[6*i1+0];  yp=xc2[6*i1+2];
    if( abs(xp-mean2[0]) < (Nx/2 - 3)*hx  and   abs(yp-mean2[1]) < (Ny/2 -3 )*hy ) {
      Interpolate_Electric_Field(Nx, Ny, hx, hy, x0, y0, Ex, Ey, xp-xoffset, yp-yoffset, ex, ey);
      xc2[6*i1+1] = xc2[6*i1+1]  + ex * 1.60217662e-19 * bbscale2 /( gamma2  * 938.27201e6 ) * (Np1/Mp1) ;
      xc2[6*i1+3] = xc2[6*i1+3]  + ey * 1.60217662e-19 * bbscale2 /( gamma2  * 938.27201e6 ) * (Np1/Mp1) ;
    }

  }
  
  //----slice 2  to slice 1:
  x0=mean2[0] - Nx * hx / 2 ;
  y0=mean2[1] - Ny * hy / 2 ;
  xoffset= mean1[0] - mean2[0];
  yoffset= mean1[1] - mean2[1];
  
  for(i1=0; i1<Nx*Ny;i1++){
    src[i1]=0;   phi[i1]=0;
  }
  
  for(i1=0;i1<Mp2;i1++){
    xp = xc2[6*i1+0];  yp = xc2[6*i1+2];
    if( abs(xp-mean2[0]) < (Nx/2 - 3)*hx  and   abs(yp-mean2[1]) < (Ny/2 -3 )*hy ) {
      Charge2Grid(xp,yp, x0, y0, Nx, Ny, hx,  hy, src);
    }
  }
  Cal_Effective_Green_Function(xoffset, yoffset, Nx, Ny, hx, hy, grn_c_r, grn_c_i);
  Cal_Potential_Convolution(Nx, Ny, hx, hy, src, grn_c_r, grn_c_i, phi);
  Cal_Electric_Field(Nx, Ny, hx, hy, phi, Ex, Ey);
  
  for(i1=0;i1<Mp1;i1++){
    xp=xc1[6*i1+0];  yp=xc1[6*i1+2];
    if( abs(xp-mean1[0]) < (Nx/2 - 3)*hx  and   abs(yp-mean1[1]) < (Ny/2 -3 )*hy ) {
      Interpolate_Electric_Field(Nx, Ny, hx, hy, x0, y0, Ex, Ey, xp-xoffset, yp-yoffset, ex, ey);
      xc1[6*i1+1] = xc1[6*i1+1]  + ex * 1.60217662e-19 * bbscale1 /(gamma1 * 938.27201e6)  * (Np2/Mp2) ;
      xc1[6*i1+3] = xc1[6*i1+3]  + ey * 1.60217662e-19 * bbscale1 /(gamma1 * 938.27201e6)  * (Np2/Mp2) ;
    }
  }
  
}

void BeamBeam6D_PIC(double Np1, int Mp1, double x1[], double Np2, int Mp2, double x2[], 
                    double gamma1, double gamma2, double bbscale1, double bbscale2,
		    int Nslice1, int Nslice2)
// solve Posisson Equation with Green Function + FFT
// sorting in slices to be done in this function 
{
  int    i,j,k;

  int    Nslice=Nslice1 * Nslice2;
  int    Index1[Nslice], Index2[Nslice];
  double Z1[Nslice1], Z2[Nslice2], S[Nslice],  Z10[Nslice];
  double temp[6], temp1;
  int    istart, iend, itemp;

  int    Mp1_coll=Mp1/Nslice1, Mp2_coll=Mp2/Nslice2;
  double xc1[6*Mp1_coll], xc2[6*Mp2_coll];
  double mean1[2], rms1[2], mean2[2], rms2[2], emit1[2], emit2[2], beta2[2], alfa2[2], gama2[2],  beta1[2], alfa1[2], gama1[2];
  double dx, dy, dpx, dpy;

  char   filename[128];
  fstream f1;

  //----sorting particles in each  bunch
  //cout<<"sorting particles, begin. "<<endl;
  for(i=0;i<Mp1;i++){
    for(j=i+1;j<Mp1;j++){
      if(x1[j*6+4]>x1[i*6+4]){
	for(k=0;k<6;k++) temp[k]   = x1[j*6+k];
        for(k=0;k<6;k++) x1[j*6+k] = x1[i*6+k];
	for(k=0;k<6;k++) x1[i*6+k] = temp[k];
      }
    }
  }

  for(i=0;i<Mp2;i++){
    for(j=i+1;j<Mp2;j++){
      if(x2[j*6+4]>x2[i*6+4]){
	for(k=0;k<6;k++) temp[k]   = x2[j*6+k];
        for(k=0;k<6;k++) x2[j*6+k] = x2[i*6+k];
	for(k=0;k<6;k++) x2[i*6+k] = temp[k];
      }
    }
  }
  //cout<<"sorting particles, end. "<<endl;
  
  /*for checking purpose*/
  if(true){ 
    
    f1.open("./bunch1.dat",ios::out);
    for(j=0;j<Mp1;j++) {
      for(k=0;k<6;k++) f1<<x1[j*6+k]<<" ";
      f1<<endl;
    }
    f1.close(); 
    for(i=0;i<Nslice1;i++){
      sprintf(filename, "./bunch1_slice%d.dat",i);
      f1.open(filename,ios::out);
      istart=Mp1_coll*i;
      iend=Mp1_coll*(i+1);
      for(j=istart;j<iend;j++) {
	for(k=0;k<6;k++) f1<<x1[j*6+k]<<" ";
	f1<<endl;
      }
      f1.close();
    }

    f1.open("./bunch2.dat",ios::out);
    for(j=0;j<Mp2;j++) {
      for(k=0;k<6;k++) f1<<x2[j*6+k]<<" ";
      f1<<endl;
    }
    f1.close(); 
    for(i=0;i<Nslice2;i++){
      sprintf(filename, "./bunch2_slice%d.dat",i);
      f1.open(filename,ios::out);
      istart=Mp2_coll*i;
      iend=Mp2_coll*(i+1);
      for(j=istart;j<iend;j++) {
	for(k=0;k<6;k++) f1<<x2[j*6+k]<<" ";
	f1<<endl;
      }
      f1.close();
    }
  }
  
  /*----check purpose----*/
  if(false){
    f1.open("./bunch1_stat.dat",ios::out);
    for(i=0;i<Nslice1;i++){
      for(j=0;j<Mp1_coll;j++){
	istart=Mp1_coll * i ;
	for(k=0;k<6;k++)  xc1[j*6+k] = x1[(istart+j)*6+k]; 
      }
      Cal_Mean_RMS(Mp1_coll, xc1, mean1, rms1);
      f1<<i<<"  "<<mean1[0]<<"  "<<mean1[1]<<"  "<<rms1[0]<<"  "<<rms1[1]<<endl;
    }
    f1.close();
    
    f1.open("./bunch2_stat.dat",ios::out);
    for(i=0;i<Nslice2;i++){
      for(j=0;j<Mp2_coll;j++){
	istart=Mp2_coll * i ;
	for(k=0;k<6;k++)  xc2[j*6+k] = x2[(istart+j)*6+k]; 
      }
      Cal_Mean_RMS(Mp2_coll, xc2, mean1, rms1);
      f1<<i<<"  "<<mean1[0]<<"  "<<mean1[1]<<"  "<<rms1[0]<<"  "<<rms1[1]<<endl;
    }
    f1.close();
  }
  
  //---sorting slices for collision
  for(i=0;i<Nslice1;i++){
    temp1=0;
    istart=i*Mp1_coll;  
    iend  =(i+1)*Mp1_coll;
    for(j=istart;j<iend;j++){
      temp1 =temp1 + x1[6*j+4]; 
    }
    Z1[i]=temp1 /Mp1_coll;
  }

  for(i=0;i<Nslice2;i++){
    temp1=0;
    istart=i*Mp2_coll;  
    iend  =(i+1)*Mp2_coll;
    for(j=istart;j<iend;j++){
      temp1 =temp1 + x2[6*j+4]; 
    }
    Z2[i]=temp1 /Mp2_coll;
  }

  /*---check purpose----*/
  if(false){ 
    f1.open("./bunch1_Z.dat",ios::out);
    for(j=0;j<Nslice1;j++) {
      f1<<j<<"  "<<Z1[j]<<endl;
    }
    f1.close(); 
    f1.open("./bunch2_Z.dat",ios::out);
    for(j=0;j<Nslice2;j++) {
      f1<<j<<"  "<<Z2[j]<<endl;
    }
    f1.close(); 
  }

  for(i=0;i<Nslice1;i++){
    for(j=0;j<Nslice2;j++){
      Index1[i*Nslice2 + j] = i;
      Index2[i*Nslice2 + j] = j;
      S[i*Nslice2 + j]= (Z1[i] -Z2[j])/2.;
      Z10[i*Nslice2 + j]= -(Z1[i] +Z2[j])/2.;
    }
  }

  for(i=0;i<Nslice;i++){
    for(j=i+1;j<Nslice;j++){
      if(Z10[j] < Z10[i] ){ 
	temp1 =  Z10[j];  Z10[j] = Z10[i];  Z10[i]=temp1;
	
	temp1 =  S[j]  ;  S[j]   = S[i]  ;  S[i]  =temp1;
        itemp =  Index1[j];  Index1[j] = Index1[i];  Index1[i] =itemp;
        itemp =  Index2[j];  Index2[j] = Index2[i];  Index2[i] =itemp; 
      }
    }
  }

  //-----collide slice to slice : be careful both beams should have same x^, y^ directions

  for(i=0;i<Nslice;i++){

    istart= Index1[i] * Mp1_coll;
    for(j=0;j<Mp1_coll;j++){
      for(k=0;k<6;k++)  xc1[j*6+k] = x1[(istart+j)*6+k]; 
    }
    istart= Index2[i] * Mp2_coll;
    for(j=0;j<Mp2_coll;j++){
      for(k=0;k<6;k++)  xc2[j*6+k] = x2[(istart+j)*6+k]; 
    }
    
    for(j=0;j<Mp1_coll;j++){
      xc1[j*6+0] = xc1[j*6+0]  +  xc1[j*6+1] * S[i];
      xc1[j*6+2] = xc1[j*6+2]  +  xc1[j*6+3] * S[i];
    }
    for(j=0;j<Mp2_coll;j++){
      xc2[j*6+0] = xc2[j*6+0]  -  xc2[j*6+1] * S[i];
      xc2[j*6+2] = xc2[j*6+2]  -  xc2[j*6+3] * S[i];
    }
    
    Cal_Mean_RMS(Mp1_coll, xc1, mean1, rms1);
    Cal_Mean_RMS(Mp2_coll, xc2, mean2, rms2);
    Cal_Emit(Mp1_coll, xc1, emit1, beta1, alfa1, gama1);
    Cal_Emit(Mp2_coll, xc2, emit2, beta2, alfa2, gama2);
      
    if(false){
      if( Index2[i] == 2 ) {
	//cout<<S[i]<<"  "<<mean1[0]<<"  "<< mean1[1]<<"  "<<rms1[0]<<"  "<<rms1[1]<<"  "<<emit1[0]<<"   "<<emit1[1]<<endl;
        cout<<S[i]<<"  "<<mean2[0]<<"  "<< mean2[1]<<"  "<<rms2[0]<<"  "<<rms2[1]<<"  "<<emit2[0]<<"   "<<emit2[1]<<endl;
      }
    }

    BeamBeam4D_PIC(Np1/Nslice1, Mp1_coll, xc1, Np2/Nslice2, Mp2_coll, xc2, gamma1, gamma2, bbscale1, bbscale2);

    for(j=0;j<Mp1_coll;j++){
      xc1[j*6+0] = xc1[j*6+0]  -  xc1[j*6+1] * S[i];
      xc1[j*6+2] = xc1[j*6+2]  -  xc1[j*6+3] * S[i];
    }
    for(j=0;j<Mp2_coll;j++){
      xc2[j*6+0] = xc2[j*6+0]  +  xc2[j*6+1] * S[i];
      xc2[j*6+2] = xc2[j*6+2]  +  xc2[j*6+3] * S[i];
    }

    if(false){

      Cal_Mean_RMS(Mp1_coll, xc1, mean1, rms1);
      Cal_Mean_RMS(Mp2_coll, xc2, mean2, rms2);
      Cal_Emit(Mp1_coll, xc1, emit1, beta1, alfa1, gama1);
      Cal_Emit(Mp2_coll, xc2, emit2, beta2, alfa2, gama2);
      
      if( Index2[i] == 2 ) {
	//cout<<S[i]<<"  "<<mean1[0]<<"  "<< mean1[1]<<"  "<<rms1[0]<<"  "<<rms1[1]<<"  "<<emit1[0]<<"   "<<emit1[1]<<endl;
        cout<<S[i]<<"  "<<mean2[0]<<"  "<< mean2[1]<<"  "<<rms2[0]<<"  "<<rms2[1]<<"  "<<emit2[0]<<"   "<<emit2[1]<<endl;
      }
    }

    istart= Index1[i] * Mp1_coll;
    for(j=0;j<Mp1_coll;j++){
      for(k=0;k<6;k++) x1[(istart+j)*6+k] =  xc1[j*6+k]; 
    }
    istart= Index2[i] * Mp2_coll;
    for(j=0;j<Mp2_coll;j++){
      for(k=0;k<6;k++)  x2[(istart+j)*6+k] = xc2[j*6+k]; 
    }

  }
  
}

double  Lumi_Cal_PIC_4D(double Np1, int Mp1, double Np2, int Mp2, double x1[], double x2[], int Ncoll, double freq)
{
  int i,j,k;
  int    Nx=64, Ny=64;
  double src1[Nx*Ny],   src2[Nx*Ny];
  double x0,  y0, hx, hy;
  double xmax, xmin, ymax, ymin;  
  fstream f1;
  double lumi;
  
  if(false){
    f1.open("./bunch1_coll.dat",ios::out);
    for(j=0;j<Mp1;j++){
      for(k=0;k<6;k++) f1<<x1[j*6+k]<<" ";
      f1<<endl;
    }
    f1.close();
    f1.open("./bunch2_coll.dat",ios::out);
    for(j=0;j<Mp2;j++){
      for(k=0;k<6;k++) f1<<x2[j*6+k]<<" ";
      f1<<endl;
    }
    f1.close();
  }

  xmin=1000.;   xmax=-1000.;
  ymin=1000.;   ymax=-1000.;
  for(i=0;i<Mp1;i++){
    if(x1[6*i+0] < xmin ) xmin=x1[6*i+0] ;
    if(x1[6*i+0] > xmax ) xmax=x1[6*i+0] ;
    if(x1[6*i+2] < ymin ) ymin=x1[6*i+2] ;
    if(x1[6*i+2] > ymax ) ymax=x1[6*i+2] ;
  }
  for(i=0;i<Mp2;i++){
    if(x2[6*i+0] < xmin ) xmin=x2[6*i+0] ;
    if(x2[6*i+0] > xmax ) xmax=x2[6*i+0] ;
    if(x2[6*i+2] < ymin ) ymin=x2[6*i+2] ;
    if(x2[6*i+2] > ymax ) ymax=x2[6*i+2] ;
  }

  x0=xmin;  y0=ymin;
  hx=(xmax-xmin)/Nx; hy=(ymax- ymin)/Ny;
  x0=xmin-4*hx;  y0=ymin-4*hy;
  hx=(xmax-xmin+8*hx)/Nx; hy=(ymax- ymin+8*hy)/Ny;

  for(i=0;i<Nx*Ny;i++) {
    src1[i]=0; 
    src2[i]=0; 
  }
  
  for(i=0;i<Mp1;i++) {
    Charge2Grid(x1[ 6*i +0], x1[ 6*i +2], x0, y0, Nx, Ny, hx, hy, src1);
  }

  for(i=0;i<Mp2;i++) {
    Charge2Grid(x2[ 6*i +0], x2[ 6*i +2], x0, y0, Nx, Ny, hx, hy, src2);
  }

  /*----for check purposes---*/
  if(false){
    f1.open("./distr_grid_1.dat",ios::out);
    for(i=0;i<Nx;i++) 
      for(j=0;j<Ny;j++) {
	f1<<i<<"  "<<j<<"   "<<src1[i*Ny+j]<<endl;
      }
    f1.close();
    
    f1.open("./distr_grid_2.dat",ios::out);
    for(i=0;i<Nx;i++) 
      for(j=0;j<Ny;j++) {
	f1<<i<<"  "<<j<<"   "<<src2[i*Ny+j]<<endl;
      }
    f1.close();
  }
  
  lumi=0;
  for(i=0;i<Nx*Ny;i++) {
    lumi = lumi + src1[i] * src2[i]; 
  }
  lumi= Ncoll * lumi * (Np1/Mp1) * (Np2/Mp2) * freq / hx / hy / 10000;
  
  //cout<<"I am here, local lumi = "<<lumi<<endl;

  return lumi;
} 

double  Lumi_Cal_PIC_6D(double Np1, int Mp1, double Np2, int Mp2, int Nslice1, int Nslice2, double x1[], double x2[], 
                        double  gamma1, double gamma2, double bbscale1, double bbscale2, int Ncoll, double freq)
// solve Posisson Equation with Green Function + FFT
// sorting in slices to be done in this function   
{
  int i,j,k;

  int    Nslice=Nslice1 * Nslice2;
  int    Index1[Nslice], Index2[Nslice];
  double Z1[Nslice1], Z2[Nslice2], S[Nslice],  Z10[Nslice];
  double temp[6], temp1;
  int    istart, iend, itemp;

  int    Mp1_coll=Mp1/Nslice1, Mp2_coll=Mp2/Nslice2;
  double xc1[6*Mp1_coll], xc2[6*Mp2_coll];
  double mean1[2], rms1[2], mean2[2], rms2[2], emit1[2], emit2[2];
  
  double dx, dy, dpx, dpy;

  char    filename[128];
  fstream f1;
  double  lumi=0;

  //----sorting particles in each  bunch
  cout<<"sorting particles, begin. "<<endl;
  for(i=0;i<Mp1;i++){
    for(j=i+1;j<Mp1;j++){
      if(x1[j*6+4]>x1[i*6+4]){
	for(k=0;k<6;k++) temp[k]   = x1[j*6+k];
        for(k=0;k<6;k++) x1[j*6+k] = x1[i*6+k];
	for(k=0;k<6;k++) x1[i*6+k] = temp[k];
      }
    }
  }

  for(i=0;i<Mp2;i++){
    for(j=i+1;j<Mp2;j++){
      if(x2[j*6+4]>x2[i*6+4]){
	for(k=0;k<6;k++) temp[k]   = x2[j*6+k];
        for(k=0;k<6;k++) x2[j*6+k] = x2[i*6+k];
	for(k=0;k<6;k++) x2[i*6+k] = temp[k];
      }
    }
  }
  cout<<"sorting particles, end. "<<endl;
  
  /*for checking purpose*/
  if(false){ 
    
    f1.open("./bunch1.dat",ios::out);
    for(j=0;j<Mp1;j++) {
      for(k=0;k<6;k++) f1<<x1[j*6+k]<<" ";
      f1<<endl;
    }
    f1.close(); 
    for(i=0;i<Nslice1;i++){
      sprintf(filename, "./bunch1_slice%d.dat",i);
      f1.open(filename,ios::out);
      istart=Mp1_coll*i;
      iend=Mp1_coll*(i+1);
      for(j=istart;j<iend;j++) {
	for(k=0;k<6;k++) f1<<x1[j*6+k]<<" ";
	f1<<endl;
      }
      f1.close();
    }

    f1.open("./bunch2.dat",ios::out);
    for(j=0;j<Mp2;j++) {
      for(k=0;k<6;k++) f1<<x2[j*6+k]<<" ";
      f1<<endl;
    }
    f1.close(); 
    for(i=0;i<Nslice2;i++){
      sprintf(filename, "./bunch2_slice%d.dat",i);
      f1.open(filename,ios::out);
      istart=Mp2_coll*i;
      iend=Mp2_coll*(i+1);
      for(j=istart;j<iend;j++) {
	for(k=0;k<6;k++) f1<<x2[j*6+k]<<" ";
	f1<<endl;
      }
      f1.close();
    }

  }

  /*for checking purpose*/
  if(false){
    f1.open("./bunch1_stat.dat",ios::out);
    for(i=0;i<Nslice1;i++){
      for(j=0;j<Mp1_coll;j++){
	istart=Mp1_coll * i ;
	for(k=0;k<6;k++)  xc1[j*6+k] = x1[(istart+j)*6+k]; 
      }
      Cal_Mean_RMS(Mp1_coll, xc1, mean1, rms1);
      f1<<i<<"  "<<mean1[0]<<"  "<<mean1[1]<<"  "<<rms1[0]<<"  "<<rms1[1]<<endl;
    }
    f1.close();
    
    f1.open("./bunch2_stat.dat",ios::out);
    for(i=0;i<Nslice2;i++){
      for(j=0;j<Mp2_coll;j++){
	istart=Mp2_coll * i ;
	for(k=0;k<6;k++)  xc2[j*6+k] = x2[(istart+j)*6+k]; 
      }
      Cal_Mean_RMS(Mp2_coll, xc2, mean2, rms2);
      f1<<i<<"  "<<mean2[0]<<"  "<<mean2[1]<<"  "<<rms2[0]<<"  "<<rms2[1]<<endl;
    }
    f1.close();
  }
  
  //---sorting slices for collision
  for(i=0;i<Nslice1;i++){
    temp1=0;
    istart=i*Mp1_coll;  
    iend  =(i+1)*Mp1_coll;
    for(j=istart;j<iend;j++){
      temp1 =temp1 + x1[6*j+4]; 
    }
    Z1[i]=temp1 /Mp1_coll;
  }

  for(i=0;i<Nslice2;i++){
    temp1=0;
    istart=i*Mp2_coll;  
    iend  =(i+1)*Mp2_coll;
    for(j=istart;j<iend;j++){
      temp1 =temp1 + x2[6*j+4]; 
    }
    Z2[i]=temp1 /Mp2_coll;
  }

  /*---check purpose----*/
  if(false){ 
    f1.open("./bunch1_Z.dat",ios::out);
    for(j=0;j<Nslice1;j++) {
      f1<<j<<"  "<<Z1[j]<<endl;
    }
    f1.close(); 
    f1.open("./bunch2_Z.dat",ios::out);
    for(j=0;j<Nslice2;j++) {
      f1<<j<<"  "<<Z2[j]<<endl;
    }
    f1.close(); 
  }

  for(i=0;i<Nslice1;i++){
    for(j=0;j<Nslice2;j++){
      Index1[i*Nslice2 + j] = i;
      Index2[i*Nslice2 + j] = j;
      S[i*Nslice2 + j]= (Z1[i] -Z2[j])/2.;
      Z10[i*Nslice2 + j]= -(Z1[i] +Z2[j])/2.;
    }
  }

  for(i=0;i<Nslice;i++){
    for(j=i+1;j<Nslice;j++){
      if(Z10[j] < Z10[i] ){ 
	temp1 =  Z10[j];  Z10[j] = Z10[i];  Z10[i]=temp1;
	temp1 =  S[j]  ;  S[j]   = S[i]  ;  S[i]  =temp1;
        itemp =  Index1[j];  Index1[j] = Index1[i];  Index1[i] =itemp;
        itemp =  Index2[j];  Index2[j] = Index2[i];  Index2[i] =itemp; 
      }
    }
  }

  /*----check purpose----*/
  if(true){
    f1.open("./collide_order.dat",ios::out);
    for(i=0;i<Nslice;i++){
      f1<<i<<"   "<<Index1[i]<<"  "<<Index2[i]<<"   "<<Z1[Index1[i]]<<"   "<<Z2[Index2[i]]<<"  "<<S[i]<<"   "<<Z10[i]<<endl;
    }
    f1.close();
  }

  //-----collide slice to slice : be careful both beams should have same x^, y^ directions

  for(i=0;i<Nslice;i++){

    istart= Index1[i] * Mp1_coll;
    for(j=0;j<Mp1_coll;j++){
      for(k=0;k<6;k++)  xc1[j*6+k] = x1[(istart+j)*6+k]; 
    }
    istart= Index2[i] * Mp2_coll;
    for(j=0;j<Mp2_coll;j++){
      for(k=0;k<6;k++)  xc2[j*6+k] = x2[(istart+j)*6+k]; 
    }
    
    for(j=0;j<Mp1_coll;j++){
      xc1[j*6+0] = xc1[j*6+0]  +  xc1[j*6+1] * S[i];
      xc1[j*6+2] = xc1[j*6+2]  +  xc1[j*6+3] * S[i];
    }
    for(j=0;j<Mp2_coll;j++){
      xc2[j*6+0] = xc2[j*6+0]  -  xc2[j*6+1] * S[i];
      xc2[j*6+2] = xc2[j*6+2]  -  xc2[j*6+3] * S[i];
    }
    
    lumi=lumi+ Lumi_Cal_PIC_4D(Np1/Nslice1,  Mp1_coll,  Np2/Nslice2, Mp2_coll, xc1, xc2, Ncoll, freq);
    
    Cal_Mean_RMS(Mp1_coll, xc1, mean1, rms1);
    Cal_Mean_RMS(Mp2_coll, xc2, mean2, rms2);   
    
   if(true){
      for(j=0;j<Mp1_coll;j++){
	dx= xc1[6*j+0] - mean2[0];
	dy= xc1[6*j+2] - mean2[1];
	BB4D(dx, dy, gamma1, Np2 / Nslice2, rms2[0], rms2[1], dpx,  dpy);
	xc1[6*j+1] = xc1[6*j+1] + dpx * bbscale1; 
	xc1[6*j+3] = xc1[6*j+3] + dpy * bbscale1; 
      }
      for(j=0;j<Mp2_coll;j++){
	dx= xc2[6*j+0] - mean1[0];
	dy= xc2[6*j+2] - mean1[1];
	BB4D(dx, dy, gamma2, Np1 / Nslice1, rms1[0], rms1[1], dpx,  dpy);
	xc2[6*j+1] = xc2[6*j+1] + dpx * bbscale2; 
	xc2[6*j+3] = xc2[6*j+3] + dpy * bbscale2; 
      }
    }
    
    for(j=0;j<Mp1_coll;j++){
      xc1[j*6+0] = xc1[j*6+0]  -  xc1[j*6+1] * S[i];
      xc1[j*6+2] = xc1[j*6+2]  -  xc1[j*6+3] * S[i];
    }
    for(j=0;j<Mp2_coll;j++){
      xc2[j*6+0] = xc2[j*6+0]  +  xc2[j*6+1] * S[i];
      xc2[j*6+2] = xc2[j*6+2]  +  xc2[j*6+3] * S[i];
    }
    
    istart= Index1[i] * Mp1_coll;
    for(j=0;j<Mp1_coll;j++){
      for(k=0;k<6;k++) x1[(istart+j)*6+k] =  xc1[j*6+k]; 
    }
    istart= Index2[i] * Mp2_coll;
    for(j=0;j<Mp2_coll;j++){
      for(k=0;k<6;k++)  x2[(istart+j)*6+k] = xc2[j*6+k]; 
    }

  }

  return lumi;
} 

double  Lumi_Cal_PIC_6D_Angle(double Np1, int Mp1, double Np2, int Mp2, double theta, int Nslice1, int Nslice2, double x1[], double x2[],
                              double  gamma1, double gamma2,  double bbscale1, double bbscale2, int Ncoll, double freq)
// theta  : half crossing angle
// solve Posisson Equation with Green Function + FFT
// sorting in slices to be done in this function   
{
  int i,j,k;

  double xtrack[6];

  int    Nslice=Nslice1 * Nslice2;
  int    Index1[Nslice], Index2[Nslice];
  double Z1[Nslice1], Z2[Nslice2], S[Nslice],  Z10[Nslice];
  double temp[6], temp1;
  int    istart, iend, itemp;

  int    Mp1_coll=Mp1/Nslice1, Mp2_coll=Mp2/Nslice2;
  double xc1[6*Mp1_coll], xc2[6*Mp2_coll];
  double mean1[2], rms1[2], mean2[2], rms2[2], emit1[2], emit2[2], beta2[2], alfa2[2], gama2[2], beta1[2], alfa1[2], gama1[2];

  double dx, dy, dpx, dpy;

  char    filename[128];
  fstream f1;
  double  part, lumi=0;


  //---Lorentz transfer to head-on collision
  for(i=0;i<Mp1;i++){
    for(j=0;j<6;j++) xtrack[j]=x1[i*6 + j];
    Lorentz_Transfer(xtrack, theta);
    for(j=0;j<6;j++) x1[i*6 + j]=xtrack[j];
  } 
  
  for(i=0;i<Mp2;i++){
    for(j=0;j<6;j++) xtrack[j]=x2[i*6 + j];
    Lorentz_Transfer(xtrack, theta);
    for(j=0;j<6;j++) x2[i*6 + j]=xtrack[j];
  } 
  
  //----sorting particles in each  bunch
  cout<<"sorting particles, begin. "<<endl;
  for(i=0;i<Mp1;i++){
    for(j=i+1;j<Mp1;j++){
      if(x1[j*6+4]>x1[i*6+4]){
	for(k=0;k<6;k++) temp[k]   = x1[j*6+k];
        for(k=0;k<6;k++) x1[j*6+k] = x1[i*6+k];
	for(k=0;k<6;k++) x1[i*6+k] = temp[k];
      }
    }
  }

  for(i=0;i<Mp2;i++){
    for(j=i+1;j<Mp2;j++){
      if(x2[j*6+4]>x2[i*6+4]){
	for(k=0;k<6;k++) temp[k]   = x2[j*6+k];
        for(k=0;k<6;k++) x2[j*6+k] = x2[i*6+k];
	for(k=0;k<6;k++) x2[i*6+k] = temp[k];
      }
    }
  }
  cout<<"sorting particles, end. "<<endl;
  
  /*for checking purpose*/
  if(false){ 
    
    f1.open("./bunch1.dat",ios::out);
    for(j=0;j<Mp1;j++) {
      for(k=0;k<6;k++) f1<<x1[j*6+k]<<" ";
      f1<<endl;
    }
    f1.close(); 
    for(i=0;i<Nslice1;i++){
      sprintf(filename, "./bunch1_slice%d.dat",i);
      f1.open(filename,ios::out);
      istart=Mp1_coll*i;
      iend=Mp1_coll*(i+1);
      for(j=istart;j<iend;j++) {
	for(k=0;k<6;k++) f1<<x1[j*6+k]<<" ";
	f1<<endl;
      }
      f1.close();
    }

    f1.open("./bunch2.dat",ios::out);
    for(j=0;j<Mp2;j++) {
      for(k=0;k<6;k++) f1<<x2[j*6+k]<<" ";
      f1<<endl;
    }
    f1.close(); 
    for(i=0;i<Nslice2;i++){
      sprintf(filename, "./bunch2_slice%d.dat",i);
      f1.open(filename,ios::out);
      istart=Mp2_coll*i;
      iend=Mp2_coll*(i+1);
      for(j=istart;j<iend;j++) {
	for(k=0;k<6;k++) f1<<x2[j*6+k]<<" ";
	f1<<endl;
      }
      f1.close();
    }

  }
 
  /*----check purpose----*/
  if(true){
    f1.open("./bunch1_stat.dat",ios::out);
    for(i=0;i<Nslice1;i++){
      for(j=0;j<Mp1_coll;j++){
	istart=Mp1_coll * i ;
	for(k=0;k<6;k++)  xc1[j*6+k] = x1[(istart+j)*6+k]; 
      }
      Cal_Mean_RMS(Mp1_coll, xc1, mean1, rms1);
      f1<<i<<"  "<<mean1[0]<<"  "<<mean1[1]<<"  "<<rms1[0]<<"  "<<rms1[1]<<endl;
    }
    f1.close();
    
    f1.open("./bunch2_stat.dat",ios::out);
    for(i=0;i<Nslice2;i++){
      for(j=0;j<Mp2_coll;j++){
	istart=Mp2_coll * i ;
	for(k=0;k<6;k++)  xc2[j*6+k] = x2[(istart+j)*6+k]; 
      }
      Cal_Mean_RMS(Mp2_coll, xc2, mean2, rms2);
      f1<<i<<"  "<<mean2[0]<<"  "<<mean2[1]<<"  "<<rms2[0]<<"  "<<rms2[1]<<endl;
    }
    f1.close();
  }
  
  //---sorting slices for collision

  for(i=0;i<Nslice1;i++){
    temp1=0;
    istart=i*Mp1_coll;  
    iend  =(i+1)*Mp1_coll;
    for(j=istart;j<iend;j++){
      temp1 =temp1 + x1[6*j+4]; 
    }
    Z1[i]=temp1 /Mp1_coll;
  }

  for(i=0;i<Nslice2;i++){
    temp1=0;
    istart=i*Mp2_coll;  
    iend  =(i+1)*Mp2_coll;
    for(j=istart;j<iend;j++){
      temp1 =temp1 + x2[6*j+4]; 
    }
    Z2[i]=temp1 /Mp2_coll;
  }

  /*---check purpose----*/
  if(true){ 
    f1.open("./bunch1_Z.dat",ios::out);
    for(j=0;j<Nslice1;j++) {
      f1<<j<<"  "<<Z1[j]<<endl;
    }
    f1.close(); 
    f1.open("./bunch2_Z.dat",ios::out);
    for(j=0;j<Nslice2;j++) {
      f1<<j<<"  "<<Z2[j]<<endl;
    }
    f1.close(); 
  }

  for(i=0;i<Nslice1;i++){
    for(j=0;j<Nslice2;j++){
      Index1[i*Nslice2 + j] = i;
      Index2[i*Nslice2 + j] = j;
      S[i*Nslice2 + j]= (Z1[i] -Z2[j])/2.;
      Z10[i*Nslice2 + j]= -(Z1[i] +Z2[j])/2.;
    }
  }

  for(i=0;i<Nslice;i++){
    for(j=i+1;j<Nslice;j++){
      if(Z10[j] < Z10[i] ){ 
	temp1 =  Z10[j];  Z10[j] = Z10[i];  Z10[i]=temp1;
	
	temp1 =  S[j]  ;  S[j]   = S[i]  ;  S[i]  =temp1;
        itemp =  Index1[j];  Index1[j] = Index1[i];  Index1[i] =itemp;
        itemp =  Index2[j];  Index2[j] = Index2[i];  Index2[i] =itemp; 
      }
    }
  }

  /*----check purpose----*/
  if(true){
    f1.open("./collide_order.dat",ios::out);
    for(i=0;i<Nslice;i++){
      f1<<i<<"   "<<Index1[i]<<"  "<<Index2[i]<<"   "<<Z1[Index1[i]]<<"   "<<Z2[Index2[i]]<<"  "<<S[i]<<"   "<<Z10[i]<<endl;
    }
    f1.close();
  }

  //-----collide slice to slice : be careful both beams should have same x^, y^ directions

  for(i=0;i<Nslice;i++){

    istart= Index1[i] * Mp1_coll;
    for(j=0;j<Mp1_coll;j++){
      for(k=0;k<6;k++)  xc1[j*6+k] = x1[(istart+j)*6+k]; 
    }
    istart= Index2[i] * Mp2_coll;
    for(j=0;j<Mp2_coll;j++){
      for(k=0;k<6;k++)  xc2[j*6+k] = x2[(istart+j)*6+k]; 
    }
    
    Cal_Mean_RMS(Mp1_coll, xc1, mean1, rms1);
    Cal_Mean_RMS(Mp2_coll, xc2, mean2, rms2);
    Cal_Emit(Mp1_coll, xc1, emit1, beta1, alfa1,gama1);
    Cal_Emit(Mp2_coll, xc2, emit2, beta2, alfa2,gama2);
    if(false){
      if( Index2[i] == 2 ) {
	//cout<<S[i]<<"  "<<mean1[0]<<"  "<< mean1[1]<<"  "<<rms1[0]<<"  "<<rms1[1]<<"  "<<emit1[0]<<"   "<<emit1[1]<<endl;
        cout<<S[i]<<"  "<<mean2[0]<<"  "<< mean2[1]<<"  "<<rms2[0]<<"  "<<rms2[1]<<"  "<<emit2[0]<<"   "<<emit2[1]<<endl;
      }
    }

    for(j=0;j<Mp1_coll;j++){
      xc1[j*6+0] = xc1[j*6+0]  +  xc1[j*6+1] * S[i];
      xc1[j*6+2] = xc1[j*6+2]  +  xc1[j*6+3] * S[i];

    }
    for(j=0;j<Mp2_coll;j++){
      xc2[j*6+0] = xc2[j*6+0]  -  xc2[j*6+1] * S[i];
      xc2[j*6+2] = xc2[j*6+2]  -  xc2[j*6+3] * S[i];
    }
    
    part= Lumi_Cal_PIC_4D(Np1/Nslice1,  Mp1_coll,  Np2/Nslice2, Mp2_coll, xc1, xc2, Ncoll, freq);
    lumi=lumi+part;

    if(false){
      cout<<"Collision "<<Index1[i]<<", "<<Index2[i]<<"  "<<S[i]<<"   :   "<<part<<endl;
    }

    Cal_Mean_RMS(Mp1_coll, xc1, mean1, rms1);
    Cal_Mean_RMS(Mp2_coll, xc2, mean2, rms2);  

    if(true){
      for(j=0;j<Mp1_coll;j++){
	dx= xc1[6*j+0] - mean2[0];
	dy= xc1[6*j+2] - mean2[1];
	BB4D(dx, dy, gamma1, Np2 / Nslice2, rms2[0], rms2[1], dpx,  dpy);
	xc1[6*j+1] = xc1[6*j+1] + dpx * bbscale1; 
	xc1[6*j+3] = xc1[6*j+3] + dpy * bbscale1; 
      }
      for(j=0;j<Mp2_coll;j++){
	dx= xc2[6*j+0] - mean1[0];
	dy= xc2[6*j+2] - mean1[1];
	BB4D(dx, dy, gamma2, Np1 / Nslice1, rms1[0], rms1[1], dpx,  dpy);
	xc2[6*j+1] = xc2[6*j+1] + dpx * bbscale2; 
	xc2[6*j+3] = xc2[6*j+3] + dpy * bbscale2; 
      }
    }
    
    for(j=0;j<Mp1_coll;j++){
      xc1[j*6+0] = xc1[j*6+0]  -  xc1[j*6+1] * S[i];
      xc1[j*6+2] = xc1[j*6+2]  -  xc1[j*6+3] * S[i];
    }
    for(j=0;j<Mp2_coll;j++){
      xc2[j*6+0] = xc2[j*6+0]  +  xc2[j*6+1] * S[i];
      xc2[j*6+2] = xc2[j*6+2]  +  xc2[j*6+3] * S[i];
    }
    
    istart= Index1[i] * Mp1_coll;
    for(j=0;j<Mp1_coll;j++){
      for(k=0;k<6;k++) x1[(istart+j)*6+k] =  xc1[j*6+k]; 
    }
    istart= Index2[i] * Mp2_coll;
    for(j=0;j<Mp2_coll;j++){
      for(k=0;k<6;k++)  x2[(istart+j)*6+k] = xc2[j*6+k]; 
    }
    
  }

  return  lumi;

} 

void BeamBeam4D_PIC_NonOffset(double Np1, int Mp1, double xc1[], double Np2, int Mp2, double xc2[],
			      double gamma1, double gamma2, double bbscale1, double bbscale2,
			      int Nx, int Ny, double hx, double hy, double grn_c_r[], double  grn_c_i[] )
// solve Posisson Equation with Green Function + FFT
// green fucntion pre-calculated  
{
  int     i, j, i1;
  double  x0, y0, xoffset, yoffset;
  double  xp, yp, ex,ey;
  double  src[Nx*Ny], phi[Nx*Ny];

  x0= -Nx * hx / 2 ;   y0= -Ny * hy / 2 ;
  xoffset= 0. ;   yoffset= 0. ;

  //----slice 1  to slice 2:
  for(i1=0; i1<Nx*Ny;i1++){
    src[i1]=0;   phi[i1]=0;
  }
  for(i1=0;i1<Mp1;i1++){
    xp = xc1[6*i1+0];  yp = xc1[6*i1+2];
    if( abs(xp) < (Nx/2 - 3)*hx  and   abs(yp) < (Ny/2 -3 )*hy ) {
      Charge2Grid(xp, yp, x0, y0, Nx, Ny, hx, hy, src);
    }
  }
  Cal_Potential_Convolution(Nx, Ny, hx, hy, src, grn_c_r, grn_c_i, phi);
  for(i1=0;i1<Mp2;i1++){
    xp=xc2[6*i1+0];  yp=xc2[6*i1+2];
    if( abs(xp) < (Nx/2 - 3)*hx  and   abs(yp) < (Ny/2 -3 )*hy ) {
      Interpolate_Electric_Field_1Potential(Nx, Ny, hx, hy, x0, y0,  phi, xp, yp, ex, ey);
      xc2[6*i1+1] = xc2[6*i1+1]  + ex * 1.60217662e-19 * bbscale2 /( gamma2  * 938.27201e6 ) * (Np1/Mp1) ;
      xc2[6*i1+3] = xc2[6*i1+3]  + ey * 1.60217662e-19 * bbscale2 /( gamma2  * 938.27201e6 ) * (Np1/Mp1) ;
    }
  }
  
  //----slice 2  to slice 1:
  for(i1=0; i1<Nx*Ny;i1++){
    src[i1]=0;   phi[i1]=0;
  }
  for(i1=0;i1<Mp2;i1++){
    xp = xc2[6*i1+0];  yp = xc2[6*i1+2];
    if( abs(xp) < (Nx/2 - 3)*hx  and   abs(yp) < (Ny/2 -3 )*hy ) {
      Charge2Grid(xp, yp, x0, y0, Nx, Ny, hx, hy, src);
    }
  }
  Cal_Potential_Convolution(Nx, Ny, hx, hy, src, grn_c_r, grn_c_i, phi);
  for(i1=0;i1<Mp1;i1++){
    xp=xc1[6*i1+0];  yp=xc1[6*i1+2];
    if( abs(xp) < (Nx/2 - 3)*hx  and  abs(yp) < (Ny/2 -3 )*hy ) {
      Interpolate_Electric_Field_1Potential(Nx, Ny, hx, hy, x0, y0,  phi, xp, yp, ex, ey);
      xc1[6*i1+1] = xc1[6*i1+1]  + ex * 1.60217662e-19 * bbscale1 /(gamma1 * 938.27201e6)  * (Np2/Mp2) ;
      xc1[6*i1+3] = xc1[6*i1+3]  + ey * 1.60217662e-19 * bbscale1 /(gamma1 * 938.27201e6)  * (Np2/Mp2) ;
    }
  }

}

void BeamBeam6D_PIC_NonOffset(double Np1, int Mp1, double x1[], double Np2, int Mp2, double x2[], 
			      double gamma1, double gamma2, double bbscale1, double bbscale2,
			      int Nslice1, int Nslice2, double zb1[], double zb2[], int Nx, int Ny, double hx, double hy, double grn_c_r[], double  grn_c_i[] )
// solve Posisson Equation with Green Function + FFT
// green fucntion and slice bounar are pre-calculated  
{
  int    i, j, k, i1;

  double Z1_l[Nslice1],  Z2_l[Nslice2];
  double Z1_g[Nslice1],  Z2_g[Nslice2];
  int    Mp1_l[Nslice1], Mp2_l[Nslice2]; 
  int    Nslice=Nslice1 * Nslice2;
  int    itemp, Index1_g[Nslice], Index2_g[Nslice];
  double S_g[Nslice], Z10_g[Nslice];

  int    ncount;
  int    Mp1_coll, Mp2_coll;
  int    dimxc1= int(2*Mp1/Nslice1 ), dimxc2= int(2*Mp2/Nslice2);
  double *xc1 =new double [dimxc1*6], *xc2 = new double [dimxc2*6];
  double sum[2], mean1[2], rms1[2], mean2[2], rms2[2];
  
  double  temp;
  char    filename[128];
  fstream f1;

  //---- slicing bunches--------
  
  for(i=0;i<Nslice1;i++) {
    Mp1_l[i]=0 ;     Mp2_l[i]=0 ;  
    Z1_l[i] =0.;     Z2_l[i] =0.;  
    Z1_g[i] =0.;     Z2_g[i] =0.;
  }    
  
  for(i=0;i<Nslice1;i++){
    temp=0; 
    for(j=0;j<Mp1;j++){
      if( x1[j*6+4] < zb1[i]  and  x1[j*6+4] >= zb1[i+1] ) {
	temp =temp + x1[6*j+4];
	Mp1_l[i]++;
      }
    }
    Z1_l[i]=temp / Mp1_l[i] ;
    Z1_g[i]= Z1_l[i];
  }
  
  for(i=0;i<Nslice2;i++){
    temp=0; 
    for(j=0;j<Mp2;j++){
      if( x2[j*6+4] < zb2[i]  and  x2[j*6+4] >= zb2[i+1] ) {
	temp =temp + x2[6*j+4]; 
	Mp2_l[i]++;
      }
    }
    Z2_l[i]=temp / Mp2_l[i];
    Z2_g[i]= Z2_l[i];
  }
  
  //------set up collision-----
  
  for(i=0;i<Nslice1;i++){
    for(j=0;j<Nslice2;j++){
      Index1_g[i*Nslice2 + j] = i;
      Index2_g[i*Nslice2 + j] = j;
      S_g[i*Nslice2 + j]= (Z1_g[i] - Z2_g[j])/2.;
      Z10_g[i*Nslice2 + j]= -(Z1_g[i] + Z2_g[j])/2.;
    }
  }
  
  for(i=0;i<Nslice;i++){
    for(j=i+1;j<Nslice;j++){
	if(Z10_g[j] < Z10_g[i] ){ 
	  temp =  Z10_g[j];  Z10_g[j] = Z10_g[i];  Z10_g[i]=temp;
	  temp =  S_g[j]  ;  S_g[j]   = S_g[i]  ;  S_g[i]  =temp;
	  itemp =  Index1_g[j];  Index1_g[j] = Index1_g[i];  Index1_g[i] =itemp;
	  itemp =  Index2_g[j];  Index2_g[j] = Index2_g[i];  Index2_g[i] =itemp; 
	}
    }
  }

  //------ colliding slice to slice---
  
    for(i=0;i<Nslice;i++){

      /*----pick out----*/
      
      ncount=0;
      for(j=0;j<Mp1;j++){
        if( x1[j*6+4] < zb1[ Index1_g[i] ]  and  x1[j*6+4] >= zb1[ Index1_g[i]+1 ]  ) {
	  for(k=0;k<6;k++) xc1[ncount*6+k]  =   x1[j*6+k];
	  ncount++;
	}
      }
      Mp1_coll = ncount;
      
      ncount=0;
      for(j=0;j<Mp2;j++){
        if( x2[j*6+4] < zb2[ Index2_g[i] ]  and  x2[j*6+4] >= zb2[  Index2_g[i]+1 ]  ) {
	  for(k=0;k<6;k++) xc2[ncount*6+k]  =   x2[j*6+k];
	  ncount++;
	}
      }
      Mp2_coll = ncount;
      
      if(Mp1_coll > 2*Mp1/Nslice1 or Mp1_coll < Mp1/Nslice1/2 or  Mp2_coll > 2*Mp2/Nslice2 or Mp2_coll < Mp2/Nslice2/2 ) {
	cout<<"Unbalanced population for predefined slices. exit. "<<endl;
        cout<<endl;
      }

      /*-----beam-beam interaction-------*/
      
      for(j=0;j<Mp1_coll;j++){
	xc1[j*6+0] = xc1[j*6+0]  +  xc1[j*6+1] * S_g[i];
	xc1[j*6+2] = xc1[j*6+2]  +  xc1[j*6+3] * S_g[i];
      }
      for(j=0;j<Mp2_coll;j++){
	xc2[j*6+0] = xc2[j*6+0]  -  xc2[j*6+1] * S_g[i];
	xc2[j*6+2] = xc2[j*6+2]  -  xc2[j*6+3] * S_g[i];
      }
      BeamBeam4D_PIC_NonOffset(Np1*Mp1_coll/Mp1, Mp1_coll, xc1, Np2*Mp2_coll/Mp2, Mp2_coll, xc2, gamma1, gamma2, bbscale1, bbscale2, Nx, Ny, hx, hy, grn_c_r,grn_c_i);
      for(j=0;j<Mp1_coll;j++){
	xc1[j*6+0] = xc1[j*6+0]  -  xc1[j*6+1] * S_g[i];
	xc1[j*6+2] = xc1[j*6+2]  -  xc1[j*6+3] * S_g[i];
      }
      for(j=0;j<Mp2_coll;j++){
	xc2[j*6+0] = xc2[j*6+0]  +  xc2[j*6+1] * S_g[i];
	xc2[j*6+2] = xc2[j*6+2]  +  xc2[j*6+3] * S_g[i];
      }

      /*---stack back------*/
      
      ncount=0;
      for(j=0;j<Mp1;j++){
	if( x1[j*6+4] < zb1[ Index1_g[i] ]  and  x1[j*6+4] >= zb1[  Index1_g[i]+1 ]  ) {
	  for(k=0;k<6;k++) x1[j*6+k] = xc1[ncount*6+k];
	  ncount++;
	}
      }
      
      ncount=0;
      for(j=0;j<Mp2;j++){
	if( x2[j*6+4] < zb2[ Index2_g[i] ]  and  x2[j*6+4] >= zb2[  Index2_g[i]+1 ]  ) {
	  for(k=0;k<6;k++) x2[j*6+k] = xc2[ncount*6+k]; 
	  ncount++;
	}
      }
      
    }
    
    delete [] xc1;
    delete [] xc2;   
}

void BeamBeam6D_PIC_Interpolation(double Np1, int Mp1, double x1[], double Np2, int Mp2, double x2[], 
				  double gamma1, double gamma2, double bbscale1, double bbscale2,
				  int Nslice1, int Nslice2, double zb1[], double zb2[], int Nx, int Ny, double hx, double hy, double grn_c_r[], double  grn_c_i[] )
// solve Posisson Equation with Green Function + FFT
// green fucntion and slice bounar are pre-calculated  
{
  int    i, j, k, i1;

  double Z1_l[Nslice1],  Z2_l[Nslice2];
  double Z1_g[Nslice1],  Z2_g[Nslice2];
  int    Mp1_l[Nslice1], Mp2_l[Nslice2]; 
  int    Nslice=Nslice1 * Nslice2;
  int    itemp, Index1_g[Nslice], Index2_g[Nslice];
  double S_g[Nslice],      Z10_g[Nslice];

  int    ncount;
  int    Mp1_coll, Mp2_coll;
  int    dimxc1= int( 2*Mp1/Nslice1 ), dimxc2= int(2*Mp2/Nslice2);
  double *xc1 =new double [dimxc1*6], *xc2 = new double [dimxc2*6];
  double sum[2], mean1[2], rms1[2], mean2[2], rms2[2];
  
  double  x0, y0, xoffset, yoffset;
  double  xp, yp, ex,ey;
  double  temp, zfront, zback, sfront, sback, sp;
  double  src[Nx*Ny], phi1[Nx*Ny], phi2[Nx*Ny];

  char  filename[128];
  fstream f1;

  //-----slicing two bunches-----
  
  for(i=0;i<Nslice1;i++) {
    Mp1_l[i]=0 ;     Mp2_l[i]=0 ;  
    Z1_l[i] =0.;     Z2_l[i] =0.;  
    Z1_g[i] =0.;     Z2_g[i] =0.;
  }    
  
  for(i=0;i<Nslice1;i++){
    temp=0; 
    for(j=0;j<Mp1;j++){
      if( x1[j*6+4] < zb1[i]  and  x1[j*6+4] >= zb1[i+1] ) {
	temp =temp + x1[6*j+4];
	Mp1_l[i]++;
      }
    }
    Z1_l[i]=temp;
  }
  
  for(i=0;i<Nslice2;i++){
    temp=0; 
    for(j=0;j<Mp2;j++){
      if( x2[j*6+4] < zb2[i]  and  x2[j*6+4] >= zb2[i+1] ) {
	temp =temp + x2[6*j+4]; 
	Mp2_l[i]++;
      }
    }
    Z2_l[i]=temp;
  }
  
  for(i=0;i<Nslice1;i++) Z1_g[i]= Z1_l[i]/Mp1_l[i];
  for(i=0;i<Nslice2;i++) Z2_g[i]= Z2_l[i]/Mp2_l[i];

  //----collision order--------

  for(i=0;i<Nslice1;i++){
    for(j=0;j<Nslice2;j++){
      Index1_g[i*Nslice2 + j] = i;
      Index2_g[i*Nslice2 + j] = j;
      S_g[i*Nslice2 + j]= (Z1_g[i] - Z2_g[j])/2.;
      Z10_g[i*Nslice2 + j]= -(Z1_g[i] + Z2_g[j])/2.;
    }
  }
  
  for(i=0;i<Nslice;i++){
    for(j=i+1;j<Nslice;j++){
	if(Z10_g[j] < Z10_g[i] ){ 
	  temp =  Z10_g[j];  Z10_g[j] = Z10_g[i];  Z10_g[i]=temp;
	  temp =  S_g[j]  ;  S_g[j]   = S_g[i]  ;  S_g[i]  =temp;
	  itemp =  Index1_g[j];  Index1_g[j] = Index1_g[i];  Index1_g[i] =itemp;
	  itemp =  Index2_g[j];  Index2_g[j] = Index2_g[i];  Index2_g[i] =itemp; 
	}
    }
  }

  //------ colliding slice to slice---
  
    for(i=0;i<Nslice;i++){
      
      /*----pick out ------*/
      
      ncount=0;
      for(j=0;j<Mp1;j++){
        if( x1[j*6+4] < zb1[ Index1_g[i] ]  and  x1[j*6+4] >= zb1[ Index1_g[i]+1 ]  ) {
	  for(k=0;k<6;k++) xc1[ncount*6+k]  =   x1[j*6+k];
	  ncount++;
	}
      }
      Mp1_coll = ncount;
      
      ncount=0;
      for(j=0;j<Mp2;j++){
        if( x2[j*6+4] < zb2[ Index2_g[i] ]  and  x2[j*6+4] >= zb2[  Index2_g[i]+1 ]  ) {
	  for(k=0;k<6;k++) xc2[ncount*6+k]  =   x2[j*6+k];
	  ncount++;
	}
      }
      Mp2_coll = ncount;
      
      if(Mp1_coll > 2*Mp1/Nslice1 or Mp1_coll < Mp1/Nslice1/2 or  Mp2_coll > 2*Mp2/Nslice2 or Mp2_coll < Mp2/Nslice2/2 ) {
	cout<<"Unbalanced population for predefined slices. exit. "<<endl;
        cout<<endl;
      }
      
      //------slice 1  to slice 2 

      zfront=-100;
      zback = 100;
      for(j=0;j<Mp2_coll;j++){
	if( xc2[j*6+4] > zfront) zfront = xc2[j*6+4] ;
	if( xc2[j*6+4] < zback)  zback  = xc2[j*6+4] ;
      }     
      sfront= ( zfront - Z1_g[ Index1_g[i] ] )/2;
      sback = ( zback  - Z1_g[ Index1_g[i] ] )/2;
      
      x0= -Nx * hx / 2 ;   y0= -Ny * hy / 2 ;
      xoffset= 0. ;   yoffset= 0. ;
      
      for(j=0;j<Mp1_coll;j++){
	xc1[j*6+0] = xc1[j*6+0]  -  xc1[j*6+1] * sfront;
	xc1[j*6+2] = xc1[j*6+2]  -  xc1[j*6+3] * sfront;
      }
      for(i1=0; i1<Nx*Ny;i1++){
	src[i1]=0; phi1[i1]=0;
      }
      for(i1=0;i1<Mp1_coll;i1++){
	xp = xc1[6*i1+0];  yp = xc1[6*i1+2];
	if( abs(xp) < (Nx/2 - 3)*hx  and abs(yp) < (Ny/2 -3 )*hy ) {
	  Charge2Grid(xp,yp, x0, y0, Nx, Ny, hx, hy, src);
	}
      }
      Cal_Potential_Convolution(Nx, Ny, hx, hy, src, grn_c_r, grn_c_i, phi1);
      for(j=0;j<Mp1_coll;j++){
	xc1[j*6+0] = xc1[j*6+0]  +  xc1[j*6+1] * sfront;
	xc1[j*6+2] = xc1[j*6+2]  +  xc1[j*6+3] * sfront;
      }
      
      for(j=0;j<Mp1_coll;j++){
	xc1[j*6+0] = xc1[j*6+0]  -  xc1[j*6+1] * sback;
	xc1[j*6+2] = xc1[j*6+2]  -  xc1[j*6+3] * sback;
      }
      for(i1=0; i1<Nx*Ny;i1++){
	src[i1]=0; phi2[i1]=0;
      }
      for(i1=0;i1<Mp1_coll;i1++){
	xp = xc1[6*i1+0];  yp = xc1[6*i1+2];
	if( abs(xp) < (Nx/2 - 3)*hx  and abs(yp) < (Ny/2 -3 )*hy ) {
	  Charge2Grid(xp,yp, x0, y0, Nx, Ny, hx, hy, src);
	}
      }
      Cal_Potential_Convolution(Nx, Ny, hx, hy, src, grn_c_r, grn_c_i, phi2);
      for(j=0;j<Mp1_coll;j++){
	xc1[j*6+0] = xc1[j*6+0]  +  xc1[j*6+1] * sback;
	xc1[j*6+2] = xc1[j*6+2]  +  xc1[j*6+3] * sback;
      }

      for(j=0;j<Mp2_coll;j++){
	sp=  (xc2[6*j+4] - Z1_g[ Index1_g[i] ] )/2;
	xc2[6*j+0] = xc2[6*j+0]  +  xc2[6*j+1] * sp;
	xc2[6*j+2] = xc2[6*j+2]  +  xc2[6*j+3] * sp;
	xp=xc2[6*j+0]; yp=xc2[6*j+2];
	if( abs(xp) < (Nx/2 - 3)*hx  and abs(yp) < (Ny/2 -3 )*hy ) {
	  Interpolate_Electric_Field_2Potential(Nx, Ny, hx, hy, x0, y0, sfront, sback, phi1, phi2, xp-xoffset, yp-yoffset, sp, ex, ey);
	  xc2[6*j+1] = xc2[6*j+1]  + ex * 1.60217662e-19 * bbscale2 /( gamma2  * 938.27201e6 ) * (Np1/Mp1) ;
	  xc2[6*j+3] = xc2[6*j+3]  + ey * 1.60217662e-19 * bbscale2 /( gamma2  * 938.27201e6 ) * (Np1/Mp1) ;
	}
	xc2[6*j+0] = xc2[6*j+0]  -  xc2[6*j+1] * sp;
	xc2[6*j+2] = xc2[6*j+2]  -  xc2[6*j+3] * sp;
      }

      //----slice 2  to slice 1
      
      zfront=-100;
      zback = 100;
      for(j=0;j<Mp1_coll;j++){
	if( xc1[j*6+4] > zfront) zfront = xc1[j*6+4] ;
	if( xc1[j*6+4] < zback)  zback  = xc1[j*6+4] ;
      }     
      sfront= ( zfront - Z2_g[ Index2_g[i] ] )/2;
      sback = ( zback  - Z2_g[ Index2_g[i] ] )/2;
      
      for(j=0;j<Mp2_coll;j++){
	xc2[j*6+0] = xc2[j*6+0]  -  xc2[j*6+1] * sfront;
	xc2[j*6+2] = xc2[j*6+2]  -  xc2[j*6+3] * sfront;
      }
      for(i1=0; i1<Nx*Ny;i1++){
	src[i1]=0; phi1[i1]=0;
      }
      for(i1=0;i1<Mp2_coll;i1++){
	xp = xc2[6*i1+0];  yp = xc2[6*i1+2];
	if( abs(xp) < (Nx/2 - 3)*hx  and abs(yp) < (Ny/2 -3 )*hy ) {
	  Charge2Grid(xp,yp, x0, y0, Nx, Ny, hx, hy, src);
	}
      }
      Cal_Potential_Convolution(Nx, Ny, hx, hy, src, grn_c_r, grn_c_i, phi1);
      for(j=0;j<Mp2_coll;j++){
	xc2[j*6+0] = xc2[j*6+0]  +  xc2[j*6+1] * sfront;
	xc2[j*6+2] = xc2[j*6+2]  +  xc2[j*6+3] * sfront;
      }
      
      for(j=0;j<Mp2_coll;j++){
	xc2[j*6+0] = xc2[j*6+0]  -  xc2[j*6+1] * sback;
	xc2[j*6+2] = xc2[j*6+2]  -  xc2[j*6+3] * sback;
      }
      for(i1=0; i1<Nx*Ny;i1++){
	  src[i1]=0; phi2[i1]=0;
      }
      for(i1=0;i1<Mp2_coll;i1++){
	xp = xc2[6*i1+0];  yp = xc2[6*i1+2];
	if( abs(xp) < (Nx/2 - 3)*hx  and abs(yp) < (Ny/2 -3 )*hy ) {
	  Charge2Grid(xp,yp, x0, y0, Nx, Ny, hx, hy, src);
	}
      }
      Cal_Potential_Convolution(Nx, Ny, hx, hy, src, grn_c_r, grn_c_i, phi2);
      for(j=0;j<Mp2_coll;j++){
	xc2[j*6+0] = xc2[j*6+0]  +  xc2[j*6+1] * sback;
	xc2[j*6+2] = xc2[j*6+2]  +  xc2[j*6+3] * sback;
      }

      for(j=0;j<Mp1_coll;j++){
	sp=  (  xc1[j*6+4]  - Z2_g[ Index2_g[i] ])/2;
	xc1[6*j+0] = xc1[j*6+0]  +  xc1[j*6+1] * sp;   
	xc1[6*j+2] = xc1[j*6+2]  +  xc1[j*6+3] * sp;
	xp=xc1[6*j+0];  yp=xc1[6*j+2];
	if(  abs(xp) < (Nx/2 - 3)*hx  and abs(yp) < (Ny/2 -3 )*hy ) {
	  Interpolate_Electric_Field_2Potential(Nx, Ny, hx, hy, x0, y0, sfront, sback,  phi1, phi2, xp-xoffset, yp-yoffset, sp, ex, ey);
	  xc1[6*j+1] = xc1[6*j+1]  + ex * 1.60217662e-19 * bbscale1 /(gamma1 * 938.27201e6)  * (Np2/Mp2) ;
	  xc1[6*j+3] = xc1[6*j+3]  + ey * 1.60217662e-19 * bbscale1 /(gamma1 * 938.27201e6)  * (Np2/Mp2) ;
	}
	xc1[6*j+0] = xc1[6*j+0]  -  xc1[6*j+1] * sp;
	xc1[6*j+2] = xc1[6*j+2]  -  xc1[6*j+3] * sp;
      }   
      
      //---stack back ----
      
      ncount=0;
      for(j=0;j<Mp1;j++){
	if( x1[j*6+4] < zb1[ Index1_g[i] ]  and  x1[j*6+4] >= zb1[  Index1_g[i]+1 ]  ) {
	  for(k=0;k<6;k++) x1[j*6+k] = xc1[ncount*6+k];
	  ncount++;
	}
      }
      
      ncount=0;
      for(j=0;j<Mp2;j++){
	if( x2[j*6+4] < zb2[ Index2_g[i] ]  and  x2[j*6+4] >= zb2[  Index2_g[i]+1 ]  ) {
	  for(k=0;k<6;k++) x2[j*6+k] = xc2[ncount*6+k]; 
	  ncount++;
	}
      }

    }

    delete [] xc1;
    delete [] xc2;

}


#include "simtrack.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace  std; 

gsl_rng            *gsl_r3;
const gsl_rng_type *gslT;

double  bn_r = 0, an_r = 0;       //  random  field  errors
double  bn_s = 10, an_s = 10;     //  static  field  errors
double  DAseed = 1234;            //  seed  for rndom errors
double  DAangle  =  15 ;          // DA search  direction, degree
double  DAdeltap = 20e-04;          // particle's  dp/p0
double  DAturn = 100000;          // particle's  dp/p0

//-----------------supporting  file---------------------------

void Track_Radial_Test(Line & linename, double x[], int nturn, int & stable, int & lost_turn, int & lost_post)
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
      if( abs(x[0]) > linename.Cell[j]->APx or abs(x[2]) > linename.Cell[j]->APy
	  or  isnan(x[0]) or isnan(x[2])  ) {
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

void Track_DA_Radial_Uniform_Test( Line & linename, int nturn, double deltap0, double angle, double sigmax0, double sigmay0, const char* filename ) 
{
  int i,k;
  double sigma_s, sigma_e, sigma, sigma1, sigma2, sigma3, sigma_max=100.0;  
  int    stable, stable1=1, lost_turn=0, lost_post=0; 
  FILE *f2;
  f2=fopen(filename,"w");
  fclose(f2);
  double x[6];

  for(i=0; i<100;i++){

    sigma = 3.0 + i*0.2;
    x[0]= sigma * cos(angle) * sigmax0;
    x[1]= 0.;
    x[2]= sigma * sin(angle) * sigmay0;
    x[3]= 0.;
    x[4]= 0.;
    x[5]= deltap0;
    stable=1; lost_turn=0; lost_post=0;
    Track_Radial_Test(linename, x, nturn, stable, lost_turn, lost_post);
    
    f2=fopen(filename,"a");
    fprintf(f2,"%f %f %f \n", angle, sigma, 1.0*stable);
    fclose(f2);
    if(stable ==0 )  exit(0);
    
  }
} 

void Assign_IR_FieldErrors(Line & linename)
{

  int     i,j,k;
  Element * temp_element;
  int     nfactor = 1;
  double  rad;

  fstream f1;
  string  line;
  string  name_temp;
  double  rad_temp;
  vector <string> namlist;
  vector <double> radlist;

  //------check  random number generators
  if(false){
    for(i=0;i<100;i++)
      // cout<<gsl_ran_flat(gsl_r3, -1.0, 1.0)<<endl;
      cout<<gsl_rng_uniform(gsl_r3)<<endl;
    exit(0);
  }
  
  /*---read  in the known magnet information---*/
  
  f1.open("./magnet_radius.dat",ios::in);
  if(!f1)
    {
      cout<<"error in opening file ./magnet_radius.dat."<<endl;
      exit(0);
    }
  
  while( getline(f1, line, '\n') ){
    //cout<<line<<endl;
    istringstream ss(line);
    ss>>name_temp>>rad_temp;
    namlist.push_back( name_temp );
    radlist.push_back( rad_temp );
  }
  
  f1.close();
  
  for(i=0;i<namlist.size();i++){
    cout<<namlist[i]<<"  "<<radlist[i]<<endl;
  }
  
  /*---addd field errors  to IR dipoles---*/
  
  if(true){
    
    cout<<".....adding field errors to IR Sbends"<<endl;
    
    string  temp_name;
    double  temp_l, temp_angle, temp_e1, temp_e2;
    double  knl[11], knsl[11];

    i=0;
    while(i < linename.Ncell-1 ) {
      
      if( ( linename.Cell[i]->S <100 or linename.Cell[i]->S > 3733 )  and  linename.Cell[i]->TYPE == string("SBEND")  ){
	
	temp_name =  linename.Cell[i]->NAME;
        temp_l =     linename.Cell[i]->L;
	temp_angle = linename.Cell[i]->GetP("ANGLE");
	temp_e1    = linename.Cell[i]->GetP("E1");
	temp_e2    = linename.Cell[i]->GetP("E2");

	rad= 50.;
	for(j=0;j<namlist.size();j++){
	  if(temp_name == namlist[j] ){
	    rad= radlist[j];
	    cout<<temp_name <<"  set to  "<<rad<<endl;
	  }
	}

	for(j=0;j<11;j++) {
	  knl[j]=0; knsl[j]=0;
	}
	
	nfactor =1 ;
        for(j=2;j<11;j++) {
	  nfactor= nfactor*j;
	  //knl[j]  =  temp_angle * 0.0001 * (bn_s + bn_r * gsl_ran_gaussian ( gsl_r3, 1.0) ) * nfactor / pow( 60 * 0.001,  j); 
	  //knsl[j] =  temp_angle * 0.0001 * (an_s + an_r * gsl_ran_gaussian ( gsl_r3, 1.0) ) * nfactor / pow( 60 * 0.001,  j);
	  knl[j]  =  temp_angle * 0.0001 * (bn_s + bn_r * gsl_ran_flat(gsl_r3, -1.0, 1.0) ) * nfactor / pow( rad * 0.001,  j); 
	  knsl[j] =  temp_angle * 0.0001 * (an_s + an_r * gsl_ran_flat(gsl_r3, -1.0, 1.0) ) * nfactor / pow( rad * 0.001,  j);
	}
	 
	linename.Delete(i);
        temp_element= new SBENDMULT(temp_name.c_str(),temp_l, temp_angle, temp_e1, temp_e2, knl, knsl);
	linename.Insert(i,temp_element);
	
	cout<<" ...field  errors  assigned  to "<<temp_name<<" with Rref "<<rad<<endl;
	
      }
      
      i++;
    }
    
  }
  

  if(true){
    
    cout<<".....adding field errors to offset  IR sbends"<<endl;
    
    string  temp_name;
    double  temp_l, temp_angle;
    double  knl[11], knsl[11];

    i=0;
    while(i < linename.Ncell-1 ) {
      
      if( ( linename.Cell[i]->S <100 or linename.Cell[i]->S > 3733 )  and
	  linename.Cell[i]->TYPE==string("MULT") and linename.Cell[i]->L != 0. ){
	
	temp_name =  linename.Cell[i]->NAME;
        temp_l =     linename.Cell[i]->L;
	temp_angle = linename.Cell[i]->GetP("K0L");

	rad= 50.;
	for(j=0;j<namlist.size();j++){
	  if(temp_name == namlist[j] ){
	    rad= radlist[j];
	    cout<<temp_name <<"  set to  "<<rad<<endl;
	  }
	}

	for(j=0;j<11;j++) {
	  knl[j]=0; knsl[j]=0;
	}

	knl[0]=temp_angle;

	nfactor =1 ;
        for(j=2;j<11;j++) {
	  nfactor= nfactor*j;
	  //knl[j]  =  temp_angle * 0.0001 * (bn_s + bn_r * gsl_ran_gaussian ( gsl_r3, 1.0) ) * nfactor / pow( 60 * 0.001,  j); 
	  //knsl[j] =  temp_angle * 0.0001 * (an_s + an_r * gsl_ran_gaussian ( gsl_r3, 1.0) ) * nfactor / pow( 60 * 0.001,  j);
	  knl[j]  =  temp_angle * 0.0001 * (bn_s + bn_r * gsl_ran_flat(gsl_r3, -1.0, 1.0) )* nfactor / pow( rad * 0.001,  j); 
	  knsl[j] =  temp_angle * 0.0001 * (an_s + an_r * gsl_ran_flat(gsl_r3, -1.0, 1.0) ) * nfactor / pow( rad * 0.001,  j);
	  
	}
	 
	linename.Delete(i);
        temp_element= new MULT(temp_name.c_str(),temp_l, knl, knsl);
	linename.Insert(i,temp_element);
	
	cout<<" ...field  errors  assigned  to "<<temp_name<<" with Rref "<<rad<<endl;
	
      }
      
      i++;
    }
    
  }
  
  
  /*---add field errors to quadrupoles---*/
  
  if(true){

    cout<<".....adding field errors to IR quadrupoles"<<endl;
    
    string  temp_name;
    double  temp_l, temp_k1l, temp_k1sl;
    int     temp_nint;
    double  knl[11], knsl[11];

    i=0;
    while(i < linename.Ncell-1 ) {
      
      if( ( linename.Cell[i]->S <100 or linename.Cell[i]->S > 3733 )  and  linename.Cell[i]->TYPE == string("QUAD")  ){

	temp_name  =  linename.Cell[i]->NAME;
        temp_l     =  linename.Cell[i]->L;
	temp_k1l   =  linename.Cell[i]->GetP("K1L");
	temp_k1sl  =  linename.Cell[i]->GetP("K1SL");

	rad= 40.;
	for(j=0;j<namlist.size();j++){
	  if(temp_name == namlist[j] ){
	    rad= radlist[j];
	    cout<<temp_name <<"  set to  "<<rad<<endl;
	  }
	}
	
	for(j=0;j<11;j++) {
	  knl[j]=0; knsl[j]=0;
	}
   
	knl[1]=  temp_k1l ;  knsl[1]=  temp_k1sl ;
	
	nfactor=1;
	for(j=2;j<11;j++) {
	  nfactor= nfactor*j;
	  //knl[j]=  temp_k1l * (40*0.001)  * 0.0001 * (bn_s + bn_r * gsl_ran_gaussian ( gsl_r3, 1.0) ) * nfactor / pow( 40 * 0.001,  j);
	  //knsl[j]= temp_k1l * (40*0.001)  * 0.0001 * (an_s + an_r * gsl_ran_gaussian ( gsl_r3, 1.0) ) * nfactor / pow( 40 * 0.001,  j);
	  knl[j]=  temp_k1l * (rad*0.001)  * 0.0001 * (bn_s + bn_r * gsl_ran_flat(gsl_r3, -1.0, 1.0) ) * nfactor / pow( rad * 0.001,  j);
	  knsl[j]= temp_k1l * (rad*0.001)  * 0.0001 * (an_s + an_r * gsl_ran_flat(gsl_r3, -1.0, 1.0) ) * nfactor / pow( rad * 0.001,  j);
	} 

        linename.Delete(i);
        temp_element= new MULT(temp_name.c_str(), temp_l, knl, knsl);
	linename.Insert(i, temp_element);
	
	cout<<" ...field  errors  assigned  to "<<temp_name<<" with Rref "<<rad<<endl;
	
      }
      
      i++;
    }
    
  }

}

void Fit_Tune_HSR(Line & linename, double q1, double q2 )
{
  int i, count=0;
  double tunex0, tuney0, tunex1=q1, tuney1=q2, dtunex, dtuney;
  double qf_k1l_0, qd_k1l_0;
  double dk1l_qf, dtunex_qf,  dtuney_qf, dk1l_qd, dtunex_qd,  dtuney_qd;
  double scale_qf, scale_qd;
  vector<int> qf_index, qd_index ;
  char * qf_name;
  char * qd_name;

  /*------sort QF  and QD  families-------*/
  qf_k1l_0 =  linename.Cell[ Get_Index(linename, "YI7_QF15", 1) ]->GetP("K1L");
  qd_k1l_0 =  linename.Cell[ Get_Index(linename, "YI7_QD16", 1) ]->GetP("K1L");

  for(i=0;i<linename.Ncell;i++){
    if( linename.Cell[i]->TYPE==string("QUAD") ){
      if ( abs( linename.Cell[i]->GetP("K1L") -  qf_k1l_0 ) < 1.0e-5 ){
	qf_index.push_back(i);
      }
      else if ( abs( linename.Cell[i]->GetP("K1L") -  qd_k1l_0 ) < 1.0e-5 ){
	qd_index.push_back(i);
      }
      else{
      }
    }
  }
  
  if( qf_index.size() == 0 or qd_index.size() == 0 ){
    cout<<" Problem encountered  during tune fitting. exit."<<endl;
    exit(0);
  }

  /*--------tune fitting----------*/
  
  Cal_Twiss(linename,0.0);
  tunex0=linename.Tune1;
  tuney0=linename.Tune2;

  while( (tunex1-tunex0)*(tunex1-tunex0) + (tuney1-tuney0)*(tuney1-tuney0)  > 1.0e-8 and count<30 ) {
    
    qf_k1l_0= linename.Cell[ qf_index[0] ]->GetP("K1L");
    qd_k1l_0= linename.Cell[ qd_index[0] ]->GetP("K1L");
    
    dk1l_qf=  qf_k1l_0 * 0.001;
    for(i=0;i<qf_index.size();i++) {
      linename.Cell[ qf_index[i] ]->SetP("K1L", qf_k1l_0 + dk1l_qf );
    }
    Cal_Twiss(linename,0.0);
    dtunex_qf=linename.Tune1 - tunex0;
    dtuney_qf=linename.Tune2 - tuney0;
    for(i=0;i<qf_index.size();i++) {
      linename.Cell[ qf_index[i] ]->SetP("K1L", qf_k1l_0 );
    }
    
    dk1l_qd=  qd_k1l_0 * 0.001;
    for(i=0;i<qd_index.size();i++) {
      linename.Cell[ qd_index[i] ]->SetP("K1L", qd_k1l_0 + dk1l_qd );
    }
    Cal_Twiss(linename,0.0);
    dtunex_qd=linename.Tune1 - tunex0;
    dtuney_qd=linename.Tune2 - tuney0;
    for(i=0;i<qd_index.size();i++) {
      linename.Cell[ qd_index[i] ]->SetP("K1L", qd_k1l_0 );
    }
    
    dtunex=tunex1- tunex0;
    dtuney=tuney1- tuney0;
    LinearEquations(dtunex_qf, dtunex_qd,dtunex, dtuney_qf, dtuney_qd, dtuney, scale_qf, scale_qd);
    for(i=0;i<qf_index.size();i++) {
      linename.Cell[ qf_index[i] ]->SetP("K1L", qf_k1l_0 +  dk1l_qf * scale_qf );
    }
    for(i=0;i<qd_index.size();i++) {
      linename.Cell[ qd_index[i] ]->SetP("K1L", qd_k1l_0 +  dk1l_qd * scale_qd );
    }
    Cal_Twiss(linename,0.0);
    tunex0=linename.Tune1;
    tuney0=linename.Tune2;

    count++;
    cout<<"iteration "<<count<<" : "<<tunex0<<"  "<<tuney0<<endl;
    
  }

}


//-------ORBIT correction
void Correct_Orbit_SVD1(Line & linename, int m, int n, vector<int> bpm_index, vector<int> kicker_index, int plane)
{
  if( n > m ) {
    cout<<" Numer of BPM m should be larger than the number of correctors n."<<endl;
    exit(0);
  }
  
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
  
  Cal_Twiss(linename,0.);
}

void Add_IR_BPMs_Correctors(Line & linename )
{
  int  i, loc_temp;
  Element * temp_element;
  
   /*----add h correctors------*/
  
  loc_temp = Get_Index(linename,"Q3PR",1);
  temp_element= new KICKER("HKICK1",0,0,0);
  linename.Insert(loc_temp, temp_element);
  
  loc_temp = Get_Index(linename,"Q2PR",1);
  temp_element= new KICKER("HKICK2",0,0,0);
  linename.Insert(loc_temp, temp_element);
  
  loc_temp = Get_Index(linename,"B0APF",1);
  temp_element= new KICKER("HKICK3",0,0,0);
  linename.Insert(loc_temp, temp_element);
  
  loc_temp = Get_Index(linename,"B1PF",1);
  temp_element= new KICKER("HKICK4",0,0,0);
  linename.Insert(loc_temp, temp_element);
  
  /*----add v correctors------*/
  
  loc_temp = Get_Index(linename,"Q3PR",1);
  temp_element= new KICKER("VKICK1",0,0,0);
  linename.Insert(loc_temp, temp_element);
  
  loc_temp = Get_Index(linename,"Q1APR",1);
  temp_element= new KICKER("VKICK2",0,0,0);
  linename.Insert(loc_temp, temp_element);

  loc_temp = Get_Index(linename,"Q1APF",1);
  temp_element= new KICKER("VKICK3",0,0,0);
  linename.Insert(loc_temp, temp_element);
  
  loc_temp = Get_Index(linename,"Q3PF",1);
  temp_element= new KICKER("VKICK4",0,0,0);
  linename.Insert(loc_temp, temp_element);

  /*-----add bpm----------*/

  loc_temp = Get_Index(linename,"IP6",1);
  temp_element= new BPM("BPM1",0);
  linename.Insert(loc_temp, temp_element);

  loc_temp = Get_Index(linename,"O_CRAB_IP6R",1);
  temp_element= new BPM("BPM2",0);
  linename.Insert(loc_temp, temp_element);
  
  loc_temp = Get_Index(linename,"O_CRAB_IP6F",1);
  temp_element= new BPM("BPM3",0);
  linename.Insert(loc_temp, temp_element);

  loc_temp = Get_Index(linename,"Q1APR",1);
  temp_element= new BPM("BPM4",0);
  linename.Insert(loc_temp, temp_element);

}

void IR_Orbit_Correction_HSR( Line & linename )
{
  int i;
  vector<int>  bpm_index;
  vector<int>  cor_index;
 
  /*------H correction----*/

  bpm_index.push_back( Get_Index(linename,"BPM1",1) );
  bpm_index.push_back( Get_Index(linename,"BPM2",1) );  
  bpm_index.push_back( Get_Index(linename,"BPM3",1) );
  bpm_index.push_back( Get_Index(linename,"BPM4",1) );  

  cor_index.push_back(  Get_Index(linename,"HKICK1",1) );
  cor_index.push_back(  Get_Index(linename,"HKICK2",1) );
  cor_index.push_back(  Get_Index(linename,"HKICK3",1) );
  cor_index.push_back(  Get_Index(linename,"HKICK4",1) );

  Cal_Twiss(linename,0);
  cout<<"...before horizontal orbit correction:"<<endl;
  cout<<linename.Cell[ bpm_index[0]  ]->X[0]<<endl;
  cout<<linename.Cell[ bpm_index[1]  ]->X[0]<<endl;
  cout<<linename.Cell[ bpm_index[2]  ]->X[0]<<endl;
  cout<<linename.Cell[ bpm_index[3]  ]->X[0]<<endl; 

  Correct_Orbit_SVD1(linename, bpm_index.size(), cor_index.size(), bpm_index, cor_index, 0);

  Cal_Twiss(linename,0);
  cout<<"...after horizontal  orbit correction:"<<endl;
  cout<<linename.Cell[ bpm_index[0]  ]->X[0]<<endl;
  cout<<linename.Cell[ bpm_index[1]  ]->X[0]<<endl;
  cout<<linename.Cell[ bpm_index[2]  ]->X[0]<<endl;
  cout<<linename.Cell[ bpm_index[3]  ]->X[0]<<endl;

  /*------V correction----*/

  bpm_index.push_back( Get_Index(linename,"BPM1",1) );
  bpm_index.push_back( Get_Index(linename,"BPM2",1) );  
  bpm_index.push_back( Get_Index(linename,"BPM3",1) );
  bpm_index.push_back( Get_Index(linename,"BPM4",1) );  

  cor_index.push_back(  Get_Index(linename,"VKICK1",1) );
  cor_index.push_back(  Get_Index(linename,"VKICK2",1) );
  cor_index.push_back(  Get_Index(linename,"VKICK3",1) );
  cor_index.push_back(  Get_Index(linename,"VKICK4",1) );

  Cal_Twiss(linename,0);
  cout<<"...before vertical orbit correction:"<<endl;
  cout<<linename.Cell[ bpm_index[0]  ]->X[2]<<endl;
  cout<<linename.Cell[ bpm_index[1]  ]->X[2]<<endl;
  cout<<linename.Cell[ bpm_index[2]  ]->X[2]<<endl;
  cout<<linename.Cell[ bpm_index[3]  ]->X[2]<<endl; 

  Correct_Orbit_SVD1(linename, bpm_index.size(), cor_index.size(), bpm_index, cor_index, 1);

  Cal_Twiss(linename,0);
  cout<<"...after vertical orbit correction:"<<endl;
  cout<<linename.Cell[ bpm_index[0]  ]->X[2]<<endl;
  cout<<linename.Cell[ bpm_index[1]  ]->X[2]<<endl;
  cout<<linename.Cell[ bpm_index[2]  ]->X[2]<<endl;
  cout<<linename.Cell[ bpm_index[3]  ]->X[2]<<endl;

}

//------------------------------
//
//   EIC  HSR DA
//
//------------------------------

int main(int argc, char* argv[])
{

  //----general parameters
  int   i,j,k,l;
  Line  hsr;
  int   loc_ip6, loc_rf, loc_temp;
  Element * temp_element, * temp_element1, * temp_element2;
  double speed_light=2.9979e8;
  double  xoffset_cc1=0, yoffset_cc1=0, xoffset_cc2=0, yoffset_cc2=0, xoffset_ip6=0, yoffset_ip6=0 ;
  fstream  f1;
    
  //----set global parameters  for proton
  Set_RefPartEnergy(275000./938.27201, 938.27201, 1, 1);
  GP.H_expand = false;

  //-----parameters  from command line
  if (argc!=9) {
    cout<<"input format :  "<<endl;
    cout <<"DA_cal "<<" turn   deltap "<<" angle "<<" bn_static "<<" an_static "<<" bn_random "<<" an_random "<<"  seed "<< endl; 
    exit (EXIT_FAILURE);
  }

  DAturn=atof(argv[1]);
  DAdeltap=atof(argv[2]);
  DAangle=atof(argv[3]);
  bn_s=atof(argv[4]);
  an_s=atof(argv[5]);
  bn_r=atof(argv[6]);
  an_r=atof(argv[7]);
  DAseed=atof(argv[8]);

  cout<<"....read-in parameters: "<<endl;
  cout<<"  turns "<<"   deltap "<<" angle "<<" bn_static "<<" an_static "<<" bn_random "<<" an_random "<<"  seed  :"<< endl;
  cout<<"   "<<DAturn<<"    "<<DAdeltap<<"    "<<DAangle<<"    "<<bn_s<<"           "<<an_s<<"          "<<bn_r<<"          "<<an_r<<"    "<<DAseed<<endl;

  //----initiate  random  seeds
  gsl_rng_env_setup();
  gslT = gsl_rng_default;
  gsl_r3 = gsl_rng_alloc(gslT);
  gsl_rng_set(gsl_r3, DAseed);

  //---read in lattices
  Read_BMAD_Lattice(hsr, "./SimTrack_hsr_base_forIRerr.tab");
 
  //--twiss without any change
  if(true){
    cout<<"...Twiss as  it is read in: "<<endl;
    Cal_Optics(hsr);
    Print_Orbit(hsr,"./orbit0.dat");
    Print_Twiss(hsr,"./twiss0.dat");   
    Print_Optics_Summary(hsr);
  }

  //---- include  RF cavities
  if(true) {
    for(i=0; i<hsr.Ncell;i++){
      if( hsr.Cell[i]->TYPE=="RFCAV" ) {
	loc_rf=i;  break;
      } }
    cout<<" RF  at  : "<<loc_rf<<endl;
    temp_element= new RFCAV("RF197Mhz",0, 6., 2520 * speed_light/hsr.Length, 0.);
    hsr.Insert(loc_rf, temp_element);
    temp_element= new RFCAV("RF591Mhz",0, 20., 7560 * speed_light/hsr.Length, 0.);
    hsr.Insert(loc_rf, temp_element);  
  }

  //----asssign  IR field  errors
  Assign_IR_FieldErrors(hsr);

  if(true){
    cout<<"...After  installing IR field errors: "<<endl;
    Cal_Optics(hsr);
    Print_Optics_Summary(hsr);
  }

  //----do  IR orbit  correction
  Add_IR_BPMs_Correctors(hsr );
  IR_Orbit_Correction_HSR( hsr );
  
  //-----do  tune  matching
  cout<<"....tune  matching :"<<endl;
  Fit_Tune_HSR(hsr, 28.228, 27.210);

  if(true){
    cout<<"...After orbit/tune  re-matching: "<<endl;
    Cal_Optics(hsr);
    Print_Orbit(hsr,"orbit1.dat");
    Print_Twiss(hsr,"twiss1.dat");
    Print_Optics_Summary(hsr);
  }
  
  for(i=0; i<hsr.Ncell;i++){
    if( hsr.Cell[i]->NAME==string("IP6") ) {
      xoffset_ip6=hsr.Cell[i]->X[0];
      yoffset_ip6=hsr.Cell[i]->X[2];
    }
    if( hsr.Cell[i]->NAME==string("O_CRAB_IP6R") ) {
      xoffset_cc1=(hsr.Cell[i]->X[0] + hsr.Cell[i-1]->X[0])/2 ;
      yoffset_cc1=(hsr.Cell[i]->X[2] + hsr.Cell[i-1]->X[2])/2 ;
    }   
    if( hsr.Cell[i]->NAME==string("O_CRAB_IP6F") ) {
      xoffset_cc2=(hsr.Cell[i]->X[0] + hsr.Cell[i-1]->X[0])/2 ;
      yoffset_cc2=(hsr.Cell[i]->X[2] + hsr.Cell[i-1]->X[2])/2 ;
    }
  }
  cout<<"offsets at IP6: "<<xoffset_ip6<<"  "<<yoffset_ip6<<endl;
  cout<<"        at CC1: "<<xoffset_cc1<<"  "<<yoffset_cc1<<endl;
  cout<<"        at CC2: "<<xoffset_cc2<<"  "<<yoffset_cc2<<endl;
  
  //-----include  crab cavity and  its  multiples
  if(true) {
    
    double  mvolt=8.5;
    double  b1= 2.82e-3 , b2= 2.24e-3, b3=5.91e-1, b4= 0;
    double  a1= -9.77e-7, a2=-3.86e-5, a3=2.34e-3, a4= 0;
    double  scale= 25.816695 / mvolt;         //  since  measurement done at 8.5MV
    string line;
    
    f1.open("./ccmult.dat",ios::in);
    if(!f1)
      {
	cout<<"error in opening the file."<<endl;
	exit(0);
      }
    f1>>mvolt;
    f1>>b1>>a1;
    f1>>b2>>a2;    
    f1>>b3>>a3;
    f1>>b4>>a4;
    
    cout<<"....read in ccmult parameters."<<endl;
    cout<<mvolt<<" "<<b1<<" "<<a1<<" "<<b2<<" "<<a2<<" "<<b3<<" "<<a3<<" "<<b4<<" "<<a4<<endl;

    scale= 25.816695 / mvolt;
    
    for(i=0; i<hsr.Ncell;i++){
      if( hsr.Cell[i]->NAME=="O_CRAB_IP6R" ) {
	loc_rf=i; break;
      }
    }
    cout<<" First CC at :"<<loc_rf<<endl;
    
    hsr.Delete(loc_rf);
    temp_element2= new DRIFT("DRIFTCC",15.200000/2.0);
    hsr.Insert(loc_rf, temp_element2);
    //temp_element1= new CRABRF("CC1", 0, 25.816695,  2520 * speed_light/hsr.Length, 0.);
    temp_element1 = new CCMULT("CCMULT1", 0,  25.816695,  2520 * speed_light/hsr.Length, 0.,
    			       b1*scale, a1*scale, b2*scale, a2*scale,b3*scale, a3*scale, b4*scale, a4*scale );
    hsr.Insert(loc_rf, temp_element1);
    hsr.Cell[ loc_rf ]->DX = xoffset_cc1; 
    hsr.Cell[ loc_rf ]->DY = yoffset_cc1;
    temp_element2= new DRIFT("DRIFTCC",15.200000/2.0);
    hsr.Insert(loc_rf, temp_element2);
    
    
    for(i=0; i<hsr.Ncell;i++){
      if( hsr.Cell[i]->NAME=="O_CRAB_IP6F" ) {
	loc_rf=i;  break;
      }
    }
    cout<<" Second CC at :"<<loc_rf<<endl;
    
    hsr.Delete(loc_rf);
    temp_element2= new DRIFT("DRIFTCC", 15.060000/2.0);
    hsr.Insert(loc_rf, temp_element2);
    //temp_element1= new CRABRF("CC2", 0, 25.816695,  2520 * speed_light/hsr.Length, 0.);
    temp_element1 = new CCMULT("CCMULT2", 0,  25.816695,  2520 * speed_light/hsr.Length, 0.,
    			       b1*scale, a1*scale, b2*scale, a2*scale, b3*scale, a3*scale,b4*scale, a4*scale);
    hsr.Insert(loc_rf, temp_element1);
    hsr.Cell[ loc_rf ]->DX = xoffset_cc2; 
    hsr.Cell[ loc_rf ]->DY = yoffset_cc2;
    temp_element2= new DRIFT("DRIFTCC", 15.060000/2.0);
    hsr.Insert(loc_rf, temp_element2);
    
  }

  //------include crossing angle  collision
  if(true) {
    
    loc_ip6=0;
    for(i=0; i<hsr.Ncell;i++){
      if( hsr.Cell[i]->NAME==string("IP6") ) {
	loc_ip6=i;
	break;
      }
    }
    
    double  xoffset = hsr.Cell[ loc_ip6 ]->X[0];
    double  yoffset = hsr.Cell[ loc_ip6 ]->X[2];
    cout<<"BB IP  at : "<<loc_ip6<<endl;
    
    temp_element1= new LBT("LBT", 0.0125);
    temp_element2= new ILBT("ILBT",0.0125);    
    //temp_element = new BEAMBEAM("IP6", 6, 1.72e11, -1.0, 0.02, 5, 20e-9, 0.45, 0.,  1.3e-9, 0.056, 0.);
    temp_element = new BEAMBEAM("IP6", 6, 1.72e11, -1.0, 0.02, 5, 20e-9, 0.45, 0.,  1.3e-9, 0.056, 0.);
    
    hsr.Delete(loc_ip6);
    hsr.Insert(loc_ip6,temp_element2);
    
    hsr.Insert(loc_ip6,temp_element);
    //hsr.Cell[ loc_ip6 ]->DX = xoffset; 
    //hsr.Cell[ loc_ip6 ]->DY = yoffset;
    cout<<"...orbit offset at IP : "<<xoffset<<" "<<yoffset<<endl;
    
    hsr.Insert(loc_ip6,temp_element1);
    
  }

  //---last  check  optics
  if(true){
    cout<<"...Twiss  with  BB: "<<endl;
    //Cal_Optics(hsr);
    //GP.twiss_6d=  true;
    //Cal_Twiss_6D(hsr);
    Cal_Twiss(hsr,0);
    Print_Orbit(hsr,"orbit2.dat");
    Print_Twiss(hsr,"twiss2.dat");
    Print_Optics_Summary(hsr);
  }

  if(true) {
    for(i=0; i<hsr.Ncell;i++){
      if( hsr.Cell[i]->NAME==string("IP6") ) {
	loc_ip6=i;
	cout<<"...Offset at IP: "<<hsr.Cell[ loc_ip6 ]->X[0]<<"   "<<hsr.Cell[ loc_ip6 ]->X[2]<<endl;
      }
      if( hsr.Cell[i]->NAME==string("CCMULT1") ) {
	loc_ip6=i;
	cout<<"...Offset at CCMULT1: "<<hsr.Cell[ loc_ip6 ]->X[0]<<"   "<<hsr.Cell[ loc_ip6 ]->X[2]<<endl;
      }
      if( hsr.Cell[i]->NAME==string("CCMULT2") ) {
	loc_ip6=i;
	cout<<"...Offset at CCMULT2: "<<hsr.Cell[ loc_ip6 ]->X[0]<<"   "<<hsr.Cell[ loc_ip6 ]->X[2]<<endl;
      }      
    }
  }

  //----concatenate
  Clean_Up(hsr);
  
  //----calculate DA
  double emitx=11.3e-9,  emity=1.0e-9;
  double betx= 0.8,      alfx=  0;
  double bety= 0.072,    alfy=  0;
  double sigmax=sqrt(emitx*betx), sigmay=sqrt(emity*bety);
  cout<<"sigmax, sigmay: "<<sigmax<<"  "<<sigmay<<endl;

  //-----forradial shift lattice, uniform DA searching 

  Track_DA_Radial_Uniform_Test( hsr, DAturn, DAdeltap, DAangle*PI/180, sigmax, sigmay,"./DA-output.dat");
    
  //--the end.
  
  return(0);
  
}


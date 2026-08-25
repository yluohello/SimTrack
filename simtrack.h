#ifndef SIMTRACK_H
#define SIMTRACK_H
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_cdf.h>

using namespace std;

//================================  
//
//  Global  constants 
//
//================================

//-----global beam and tracking parameters
class  GlobalVariables
{
 public:
  GlobalVariables()
    {
      A= 1;                     // reference particle:   only used for longitudinal optics calculation
      Q= 1;                     // reference particle:   only used for RF acceleration and longitudinal optics calculation  
      m0     =1.6726216e-27;    // reference particle:   total rest mass,   unit in kg
      m      =1.6726216e-27;    // reference particle:   total mass,  unit in kg      
      q      =1.6021765e-19;    // reference particle:   total charge, unit in C
      gamma = 106.;             // reference particle:   gamma,  proton rest energy 938.27201MeV 
      beta  = 0.999955;         // reference particle:   volecity dived by light speed
      energy0=938.27201;        // reference particle:   total rest energy in MeV
      energy= 99456.8;          // reference  particle:  total energy in MeV
      p     = 1.;               // reference  particle:  total  momentum in SI unit
      cp    = 1.;               // reference  particle:  total momentum * c in SI unit
      brho =  331.738;          // reference particle's magnetic rigid= GP.p/GP.q, with sign ! 
                                // currentlt only used in Survey(), will be replaced after setting energy.
      G=1.7928474;              // spin G factor                         
      Gr=191.079706193;         // spin G.gamma
      step_deltap=0.0001;       // only used in numeric tune  calculations vs  deltap
      step_deltaz=0.02 ;        // only used in numeric tunecalculations  vs  z0     
      turn = 0;                 // for tracking turn control for modulating elements
      harm = 360;               // only used for longitudinal parameter calculation
      circumference = 3833.8451;// only used to FRF calculation on acceleration
      dgamma =0. ;              // GP.gamma's change in each turn, for the reference particle, only used for acceleration
      U0rad  =0. ;              // U0rad for reference particle 's SR energy loss, in unit of MeV, analytically calculated
      gammat =23.50631023;      // GP.gammat only for RF cavity on acceleration;
      H_expand = false;         // if true, ultrarelativistic approximation, beta=1. if not, exact Hamitonian for sbend, drift
      radiate  = false;         // if true, synchrotron radiation  will be turned on for sbends.
      //      twiss_6d = false;         // flag for 5D or 6D Twiss calculation, 0 -> 5d Twiss ( delta=const.), 1-> 6d Twiss.
      quad_fringe = false;      // if true, to incldue quadrupole fringe fields on both ends
    }
  int    A, Q, harm, turn;
  double m0, m, q, energy0, energy, p, cp, gamma, beta, brho;
  double G, Gr, step_deltap, step_deltaz;
  double dgamma, gammat;
  double circumference, U0rad;
  //bool   H_expand, radiate, twiss_6d, quad_fringe;
  bool   H_expand, radiate, quad_fringe;
};

extern GlobalVariables GP;
extern double PI;
extern double speed_light;

void   Print_GlobalVariables();
void   Set_RefPartEnergy(double gamma, double E0, double Q, double A);
void   Set_Spin_G(double g);

void   Cal_ParticleDelta(double pt, double &  delta1, double & gamma1, double & beta1);
double DeltaToPt(double delta);
double PtToDelta(double pt);

void   Cal_Momenta_Velocity(double x[6], double p[3], double v[3], double & gamma1, double & delta1, double & beta1);
void   Cal_Momenta_Velocity1(double x[6], double p[3], double v[3], double & gamma1, double & delta1, double & beta1);

double gamma_to_velocity(double gamma);
double velocity_to_gamma(double v[3]);

//===========================================
//
//          supporting functions
//
//===========================================

void    SplitString(const std::string& input, int& itemCount, std::vector<std::string>& items);
int     StringToInteger(const std::string & input );
double  StringToDouble(const std::string & input );
bool    StringInStringVector(const std::vector<std::string>& vec, const std::string& target);

int    fac( int i );
double rnd(double & r);
double gaussian(double u,double g, double & r);
double guassian_with_cut(double cut);
void   Init_GSL_Seeds();

double splint(double XA[], double YA[], int N, double X);
double Interpolate(double x[], double y[], int n, double t);
double Interpolate1(double x[], double y[], int n, double t);

void Generate_SectionMap_4D(double matrix[], double betax1, double alfax1, double betax2, double alfax2, 
 			                     double betay1, double alfay1, double betay2, double alfay2,
			                     double dphix,   double dphiy);
void Generate_SectionMap_6D(double matrix[], double betax1, double alfax1, double betax2, double alfax2, 
 			                     double betay1, double alfay1, double betay2, double alfay2,
			                     double dphix,   double dphiy);
void Generate_OneTurnMap_4D(double matrix[], double mux, double betax, double alfax, 
			                     double muy, double betay, double alfay );
void Generate_OneTurnMap_6D(double m66[],    double mux, double betax, double alfax, 
                                             double muy, double betay, double alfay, 
			                     double muz, double betaz, double alfaz );

void Matrix_Transfer(double matrix[36], double x[6]);
void Matrix_Transfer(double matrix[36], double x[], int  Np);

void Transfer_Ring(int Npart, double x[], double mux, double chromx, double betax, double alfax, double muy, double chromy, double betay, double alfay, double muz, double betaz, double alfaz);
//void Transfer_Ring2(int Npart, double x[], double mux, double chromx, double betax, double alfax, double muy, double chromy, double betay, double alfay, double muz, double betaz, double alfaz);
//void Transfer_Ring2c(int Npart, double x[], double mux, double chromx, double betax, double alfax, double muy, double chromy, double betay, double alfay, double muz, double betaz, double alfaz);

void   Cal_Mean(int Npart, double PartDist[], double mean[6] );
void   Cal_Max_Min_Mean_RMS( double x[], int n, double & max, double & min,  double & mean,  double & rms);
void   Cal_Mean_RMS(double x[], int n, double & mean, double & rms );
void   Cal_Mean_RMS(int Npart, double PartDist[], double mean[2], double rms[2]);
void   Cal_Mean_RMS_6D(int Npart,  double PartDist[],  int id[],  double mean[], double rms[]);

void   Cal_Twiss_Emit( double x[], double px[], double n, double & beta, double & alpha, double & gamma, double & emit);
void   Cal_Emit(int Npart, double PartDist[], double emit[2], double beta[2], double alfa[2], double gama[2]);
void   Cal_Emit_3D(int Npart, double PartDist[], double emit[], double beta[], double alfa[], double gama[]);

void   Cal_Num_Bin( double x[], int m, double y[], int n, double xmax, double xmin);
void   Cal_Num_Bin_Weight( double x[], int m, double weight[], double y[], int n, double xmax, double xmin);

void  Cal_Max_6D(int Mp2, double x2[],  double maxi[]);

void   llsq( int n, double x[], double y[], double & a, double  & b );
void   llsq_err( int n, double x[], double y[], double & a, double  & b ,  double  & ar,  double  & br, double & chi2, double  & r2);
void   pfit(double x[],double y[], int n, double a[], int m);
void   Linear_Fit(double x[],double y[], int n, double & term0, double & term1, double & term2, double & term3);

void   Generate_Eye_Matrix( double eye[], int m );
void   Print_Matrix(double  mat[], int m, int n);
void   mat_vect_mult(double mat[], int m,  double xin[], double xout[]);
void   mat_vect_mult(double mat[], int m,  int n, double xin[], double xout[]);
void   mat_vect_mult(double m66[],  double xin[6], double xout[6]);
void   mat_mult(double a[],  double b[], double c[], int m, int n, int k);
double mat_det(double a[],int n);
void   mat_transpose(double a[], double at[], int m, int n);
int    mat22_inv(double a[]);
int    mat_inv(double a[],int n);
void   mat_change_hessenberg(double a[], int n);
int    mat_root_hessenberg(double a[],int n,double u[],double v[],double eps,int jt);

template <class T> void vector_cross_product(T a[], T b[], T c[]);

void   LinearEquations(double a1, double b1, double c1, double a2, double b2, double c2, double & x, double & y  );
void   rotation(double & x1, double & y1, double x0, double y0, double dtheta);

void   EigenSolver(double Matrix[4][4] , double wr[4], double  wi[4], double vr[4][4], double vi[4][4]); 
void   EigenSolver_6D(double Matrix[6][6] , double wr[6], double  wi[6], double vr[6][6], double vi[6][6]); 

void   brmul(double a[], double b[] , int m, int n, int k, double c[] );
static void ppp(double *a,double *e,double *s,double *v,int m,int n);
static void sss(double fg[2],double cs[2]);
int    bmuav(double a[],int m,int n,double *u,double *v,double eps,int ka);
void   mat_inv_svd(int m, int  n, double A[]);

void   fft(int m, double*x, double*y);
void   ifft(int m, double*x, double*y);
void   FineTuneFinder(int n, double x[], double px[], double tune_lower, double tune_upper, double & peaktune);
void   FineTuneFinder(int n, double x[], double tune_lower, double tune_upper, double & peaktune);

void Track_RungeKutta_BField_1step(double m0, double q, double gamma, double r0[3], double v0[3], void (*GetB)(double r[3], double B[]), double dt);
void Track_RungeKutta_EMField_1step(double m0, double q, double gamma0, double r0[3], double p0[3], void (*GetB)(double r[3], double B[]), double dt);
void Track_RungeKutta_Spin_1step(double m0, double q, double gamma, double r0[3], double p0[3], double s0[3], void (*GetB)(double r[3], double B[]), double dt);
  
int Powerof2(int n,int *m,int *twopm);
int FFT1D(double x[],double y[], int m, int dir);
int FFT2D(double x[], double y[], int nx, int ny,int dir);
void  Cal_Grid_Weight(double x1, double y1, double x0, double y0, int Nx, int Ny, double hx,  double hy, int xmesh[], int ymesh[], double xweight[], double yweight[]);
void  Charge2Grid(double x, double y, double x0, double y0, int Nx, int Ny, double hx,  double hy, double src[]);
double int_green(double x, double y);
double calc_egreen(double x1, double x2, double y1, double y2, double hx, double hy);
void Cal_Effective_Green_Function(double xoffset, double yoffset, int Nx, int Ny, double hx, double hy, double grn_c_r[], double grn_c_i[] );
void Cal_Potential_Direct(double xoffset, double yoffset, int Nx, int Ny, double hx, double hy, double src[],  double phi[]);
void Cal_Potential_Convolution(int Nx, int Ny, double hx, double hy, double src[], double grn_c_r[],double grn_c_i[], double phi[]);
void Cal_Potential_Convolution2(int Nx, int Ny, double hx, double hy, double xoffset, double yoffset, double src[], double phi[]);
void Cal_Electric_Field(int Nx, int Ny, double hx, double hy, double phi[], double Ex[], double  Ey[]);
void Interpolate_Electric_Field(int Nx, int Ny, double hx, double hy, double x0, double y0,  double Ex[], double  Ey[], double x, double y, double & ex, double & ey);
void Cal_Electric_PIC(int Nx, int Ny, double x0, double y0, double xoffset, double yoffset, double hx, double hy, double src[], double Ex[], double  Ey[]);
void Interpolate_Electric_Field_1Potential(int Nx, int Ny, double hx, double hy, double x0, double y0,   double phi[] , double x, double y, double & ex, double & ey);
void Interpolate_Electric_Field_2Potential(int Nx, int Ny, double hx, double hy, double x0, double y0, double sfront, double sback,  double phi1[], double  phi2[],  double x, double y, double sp, double & ex, double & ey);
void Zboundary_slices_Gaussian(double zb[], int Nslice, double sigmal);

class tps
{
public:
  tps ();
  tps (double d);
  tps (double d[7]);
  double &operator[] (int i);
  friend tps operator-(tps x);
  friend ostream & operator<< (ostream & stream, tps x);
  friend tps operator+ (tps x, tps y);
  friend tps operator- (tps x, tps y);
  friend tps operator* (tps x, tps y);
  friend tps DAinv (tps x);
  friend tps operator/ (tps x, tps y);
  friend tps sqr (tps x);
  friend tps sqrt (tps x);
  friend tps sin (tps x);
  friend tps asin (tps x);
  friend tps cos (tps x);
  friend tps tan (tps x);
  friend tps atan (tps x);
  friend tps sinh (tps x);
  friend tps cosh (tps x);
  friend tps exp (tps x);
  friend tps ln (tps x);
private:
  double sample[7];
};

class linmap
{
public:
  linmap ();
  void identity ();
  linmap (double x[6]);
  tps & operator[] (int i);
  friend linmap operator+ (linmap x, linmap y);
  void print ();
 private:
  tps map0[6];
};

void Getmat (linmap map0, double x[36]);
void Getpos (linmap map0, double x[6]);

//=========================================
//
//         particle transfer functions
//
//=========================================

template<class T>  void  GtoL1(T x[], double DX, double DY, double DT);
template<class T>  void LtoG1(T x[], double DX, double DY, double DT);
void Cal_W_Matrix(double theta, double phi, double psi, double W[]);
void Cal_VE_WE_Matrx(double length, double angle, double psi, double VE[], double WE[]);
template<class T>  void  GtoL2(T x[], double  length, double angle, double dx, double dy, double ds, double theta, double phi, double psi);
template<class T>  void LtoG2(T x[], double length, double angle, double dx, double dy, double ds, double theta, double phi, double psi);
template<class T>  void  GtoL(T x[], double length, double angle, double dx, double dy, double ds, double theta, double phi, double psi ) ;
template<class T>  void LtoG(T x[], double length, double angle, double dx, double dy, double ds, double theta, double phi, double psi);

template <class T> void DRIFT_Pass(T x[], double L);
template <class T> void drift_polar_pass(T x[], double L, double href);
template <class T> void bend_kick_pass(T x[], double L, double href  );
template <class T> void general_bend_kick_pass(T x[], double L, double href, double hreal);
template <class T> void bend_kick_pass_exact(T x[], double L, double href);
//template <class T> void sbend_exact_pass_v0(T x[], double L, double Angle);
void sbend_exact_pass_v0(double x[], double L, double Angle);
template <class T> void sbend_exact_pass(T x[], double L, double Angle, double cosAngle, double sinAngle);
template <class T> void general_sbend_exact_pass(T x[], double L, double Angle, double cosAngle, double sinAngle, double b0);
void Cal_gsbend_l_angle(double x0[6], double L, double Angle, double K0L, double & actual_l, double & actual_angle);
template <class T> void quad_kick_pass(T x[], double k1l, double k1sl);
template <class T> void sext_kick_pass(T x[], double k2l, double k2sl);
template <class T> void oct_kick_pass(T x[], double k3l, double k3sl);
template <class T> void mult_kick_pass(T x[], int Norder, double KNL[11], double KNSL[11]);
template <class T> void bend_mult_kick_pass(T x[], double L, double href, int Norder, double KNL[11], double KNSL[11]);
template <class T> void general_bend_mult_kick_pass(T x[], double L, double href, double hreal, int Norder, double KNL[11], double KNSL[11]);

template <class T>  void cal_Bfield(T x[], int Nint, int Norder, double Angle, double KNL[11], double KNSL[11]);
template <class T>  void radiate(T x[], double L, double href, T BLbrho[]);
template <class T>  void radiate1(T x[], double L, double href, T BLbrho[]);
void radiate2(double x[], double L, double href, double BLbrho[]);

template <class T> void get_Axy_wiggler(T AxoBrho[], T AyoBrho[], double z, T x[], double b0, double kx, double kz, double phiz0);
template <class T> void Wiggler_Pass_Forest(T x[], double l, int nint, double b0, double kx, double kz, double phiz0);
template <class T> void get_Ax_Wu(T & AxoBrho, T & IntpAxpydx, double z, T x[], double b0, double kx, double kz, double phiz0);
template <class T> void get_Ay_Wu(T & AyoBrho, T & IntpAypxdy, double z, T x[], double b0, double kx, double kz, double phiz0);
template <class T> void get_wiggler_B(T x[], T B[], double b0, double kx, double kz, double phiz0);
template <class T> void Wiggler_Pass_Wu(T x[], double l, int nint, double b0, double kx, double kz, double phiz0);

void spin(double x[], double angle, double href, double BLbrho[]);

template <class T> void SBEND_Pass(T x[], double L, int Nint, double Angle, double E1, double E2);
template <class T> void GSBEND_Pass(T x[], double L, int Nint, double Angle, double K0L, double E1, double E2);
template <class T> void QUAD_Pass(T x[], double L, int Nint, double k1l, double k1sl);
template <class T> void QUAD_Pass_Radiate(T x[], double L, int Nint, double k1l, double k1sl);
template <class T> void SEXT_Pass(T x[], double L, int Nint, double k2l, double k2sl);
template <class T> void SEXT_Pass_Radiate(T x[], double L, int Nint, double k2l, double k2sl);
template <class T> void OCT_Pass(T x[], double L, int Nint, double k3l, double k3sl);
template <class T> void OCT_Pass_Radiate(T x[], double L, int Nint, double k3l, double k3sl);
template <class T> void MULT_Pass(T x[], double L, int Nint, int Norder, double KNL[11], double KNSL[11]);
template <class T> void MULT_Pass_Radiate(T x[], double L, int Nint, int Norder, double KNL[11], double KNSL[11]);
template <class T> void GMULT_Pass(T x[], double L, int Nint,  int Norder, double Angle,  double E1, double E2, double KNL[11], double KNSL[11]);
template <class T> void SBENDMULT_Pass(T x[], double L, int Nint,  int Norder, double Angle,  double E1, double E2, double KNL[11], double KNSL[11]);
template <class T> void SMULT_Pass(T x[], double L, int Nint,  int Norder, double Angle,  double E1, double E2, double KNL[11], double KNSL[11]);
template <class T> void GSBENDMULT_Pass(T x[], double L, int Nint,  int Norder,double Angle, double K0L, double E1, double E2, double KNL[11], double KNSL[11] ); 
template <class T> void SOLEN_Pass(T x[], double L, int Nint, double KS);
template <class T> void WIGGLER_Pass(T x[], double L, int Nint, double B0, double KX, double KZ, double PHIZ0);
template <class T> void MATRIX_Pass( T x[], double L, double XCO_IN[6], double XCO_OUT[6], double M66[36]);
template <class T> void KICK_Pass(T x[], double L, double HKICK, double VKICK);
template <class T> void ACMULT_Pass(T x[], double L, int Norder, double KL, int TTURNS, double PHI0);
template <class T> void ACDIP_Pass(T x[], double L, double HKICKMAX, double VKICKMAX, double NUD, double TURNS, double TURNE, double PHID);
void  RFCAV_Pass(double x[], double L, double VRF, double FRF, double PHASE0);
void  CRABRF_Pass(double x[], double L, double VRF, double FRF, double PHASE0);
void  CCMULT_Pass(double x[], double L, double VRF, double FRF, double PHASE0, double B1, double A1, double B2, double A2, double B3, double A3);
void  LBT_Pass(double x[], double theta);
void  ILBT_Pass(double x[], double theta);
template <class T> void DIFFUSE_Pass(T x[], double DIFF_X, double  DIFF_Y, double DIFF_DELTA);
template <class T> void COOLING_Pass(T x[], double ALPHA);
template <class T> void TRANS_Pass( T x[], double DX, double DY, double DS);
template <class T> void SROTAT_Pass( T x[], double PSI);
template <class T> void YROTAT_Pass( T x[], double THETA);
template <class T> void YROTAT_MADX_Pass( T x[], double THETA);
template <class T> void XROTAT_Pass( T x[], double PHI);
template <class T> void PATCH_Pass( T x[], double DX, double DY, double DS, double THETA);

void cerrf(double xx, double yy, double & wx, double & wy );
void BB4D(double x, double y, double gamma, double N, double sigmax, double sigmay, double & Dpx, double & Dpy);
void BB4D(double xi[6], double gamma, double N, double bbscale, double sigmax, double sigmay);
void BB4D(double x, double y, double Q1, double E1, double N2, double Q2, double sigmax, double sigmay, double & Dpx, double & Dpy);
void BB4D(double xi[6], double Q1, double E1, double N2, double Q2, double sigmax, double sigmay);
void BB6D(double x[], double gamma, double Np, double bbscale, double sigma_l, int N_slice, double emitx_rms,  double betax_star, double alfx_star,  double emity_rms,  double betay_star, double alfy_star);
void BB6D_Lumi(double x[], double gamma, double Np, double bbscale, double sigma_l, int N_slice, double emitx_rms,  double betax_star, double alfx_star, double emity_rms,  double betay_star, double alfy_star, double & lumi_part);
void  Lorentz_Transfer(double x[6], double theta);
void  Lorentz_Transfer_Inverse(double x[6], double theta);
void  Lorentz_Transfer(int Np, double x[], double theta);
void  Lorentz_Transfer_Inverse(int Np, double x[], double theta);
void  BB6D_Angle(double x[], double gamma, double Np, double bbscale, double sigma_l, double theta, int N_slice, double emitx_rms, double betax_star, double alfx_star, double emity_rms,  double betay_star, double alfy_star);
void BB6D_Angle_Lumi(double x[6], double gamma, double Np, double bbscale, double sigma_l, double theta, int N_slice, double emitx_rms,  double betax_star, double alfx_star, double emity_rms,  double betay_star, double alfy_star, double & lumi_part);

void elens_pass_round_Gaussian_topoff(double x[], double gamma, double Ne, double Le, double beta_e, int N_slice,  double sigmax, double sigmay) ;
void elens_pass_round_Gaussian_truncated(double x[], double gamma,  double Ne, double Le, double beta_e, int N_slice,  double sigmax, double sigmay) ;
void elens_pass_round_uniform(double x[], double gamma,  double Ne, double Le, double beta_e, int N_slice, double sigmax, double sigmay) ;

void BEAMBEAM_Pass(double x[], int TREATMENT, double NP, double BBSCALE, double SIGMAL, int NSLICE, double EMITX,  double BETAX, double ALFAX, double EMITY, double BETAY, double ALFAY);
void LRBB_Pass(double x[], double NP, double BBSCALE, double SEPX, double SEPY, double SIGMAX, double SIGMAY); 
void CRBB_Pass(double x[], double NP, double THETA, double BBSCALE, double SIGMAL, int NSLICE, double EMITX,  double BETAX, double ALFAX, double EMITY, double BETAY, double ALFAY);
void ELENS_Pass(double x[],  double Le, double Ne, double bbscale, double beta_e, int N_slice, double sigmax, double sigmay);
void HELENS_Pass(double x[],  double Le, double Ne, double bbscale, double beta_e, int N_slice, double rout, double rin);
void ERHICBB_Pass(double x[], double gamma, double Ne, double bbscale);

void SBEND_sPass(double x[], double L, int Nint, double Angle, double E1, double E2);
void QUAD_sPass(double x[], double L, int Nint, double k1l, double k1sl);
void SEXT_sPass(double x[], double L, int Nint, double k2l, double k2sl);
void OCT_sPass(double x[], double L, int Nint, double k3l, double k3sl);
void MULT_sPass(double x[], double L, int Nint, int Norder, double KNL[11], double KNSL[11]);
void GMULT_sPass(double x[], double L, int Nint, int Norder, double Angle, double E1, double E2, double KNL[11], double KNSL[11]);
void SBENDMULT_sPass(double x[], double L, int Nint, int Norder, double Angle, double E1, double E2, double KNL[11], double KNSL[11]);
void SMULT_sPass(double x[], double L, int Nint, int Norder, double Angle, double E1, double E2, double KNL[11], double KNSL[11]);
void SOLEN_sPass(double x[], double L, int Nint, double KS);
void KICK_sPass(double x[], double L, double HKICK, double VKICK);
void ACMULT_sPass(double x[], double L, int Norder, double KLMAX, double KSLMAX, int TTURNS, double PHI0);
void ACDIP_sPass(double x[], double L, double HKICKMAX, double VKICKMAX, double NUD, double TURNS, double TURNE, double PHID);
void ROTAT_sPass(double x[9], double L, double n[3], double angle);
void SNAKE_sPass(double x[9], double L, double n[3], double angle);

//===========================================
//
//       Element definition
//
//============================================

class Element
{
 public:
  Element(string name);
  virtual ~Element() = default;
  
  virtual void    SetP(const char *name, double value)=0;
  virtual double  GetP(const char *name)=0;
  virtual void    Pass(double x[6])=0;
  virtual void    DAPass(tps x[6])=0; 
  virtual void    sPass(double x[9])=0;
  
  string  NAME, TYPE, GROUP;
  double  S, L;
  double  DX, DY, DS;                   //  this line and next line: misalignment of elements
  double  DTHETA, DPHI, DPSI;          //   dpsi(s as axis), dphi(x as axis), dtheta(y as axis)
  double  X[6], T[36], M[36], A[36];
  double  Beta1, Alfa1, Beta2, Alfa2,  Beta3, Alfa3, Mu1, Mu2, Mu3;
  double  r, c11, c12, c21, c22;
  double  Etax, Etay, Etaxp, Etayp;  //  momentum dispersion
  double  Ksix, Ksiy, Ksixp,  Ksiyp;  // crab dispersion  
  double  APx, APy;
  double  n0[3];    
};

//---------------DRIFT-----------------------------------
class DRIFT: public Element
{
 public:
  DRIFT(string name, double l);
  void   SetP(const char *name, double value); 
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
};

//--------------------------SBEND------------------------------------
class SBEND: public Element
{
 public:
  SBEND(string name, double l, double angle, double e1, double e2);
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  double  ANGLE, E1,E2;   // positve ANGLE bends particle to negative Vx dirction
  int Nint;
};

//--------------------------GSBEND------------------------------------
class GSBEND: public Element
{
 public:
  GSBEND(string name, double l, double angle, double k0l, double e1, double e2);
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  double  ANGLE, K0L, E1,E2;   // positve ANGLE bends particle to negative x dirction
  int Nint;
};

//---------------------QUAD-----------------------------------------
class QUAD: public Element
{
 public:
  QUAD(string name, double l, double k1l, double k1sl );
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  double  K1L, K1SL;
  int Nint, Norder;
};

//--------------------------------SKEWQ------------------------------
class SKEWQ: public Element
{
 public:
  SKEWQ(string name, double l, double k1sl);
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  double  K1SL;
  int Nint, Norder;
};

//----------------------------SEXT-------------------------------------------
class SEXT: public Element
{
 public:
  SEXT(string name, double l, double k2l, double k2sl);
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  double  K2L, K2SL;
  int Nint, Norder;
};

//-------------------------OCT---------------------------------------------
class OCT: public Element
{
 public:
  OCT(string name, double l, double k3l, double k3sl);
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  double  K3L,K3SL;
  int Nint,Norder;
};

//-------------------------------MULT--------------------------------------
class MULT: public Element
{
 public:
  MULT(string name, double l, double knl[11], double knsl[11]);
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  double  KNL[11], KNSL[11];
  int Nint, Norder;
};

//-------------------------------GMULT--------------------------------------
class GMULT: public Element
{
 public:
  GMULT(string name, double l, double angle, double e1, double e2, double knl[11], double knsl[11]);
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  double  ANGLE, E1,E2, KNL[11], KNSL[11];
  int Nint, Norder;
};

//-------------------------------SBENDMULT,exactly same as  GMULT--------------------------------------
class SBENDMULT: public Element
{
 public:
  SBENDMULT(string name, double l, double angle, double e1, double e2, double knl[11], double knsl[11]);
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void  Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  double  ANGLE, E1,E2, KNL[11], KNSL[11];
  int Nint, Norder;
};

//-------------------------------SMULT--------------------------------------
class SMULT: public Element
{
 public:
  SMULT(string name, double l, double angle, double e1, double e2, double knl[11], double knsl[11]);
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void  Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  double  ANGLE, E1,E2, KNL[11], KNSL[11];
  int Nint, Norder;
};

//--------------------------GSBENDMULT------------------------------------
class GSBENDMULT: public Element
{
 public:
  GSBENDMULT(string name, double l, double angle, double k0l, double e1, double e2, double knl[11], double knsl[11]);
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  double  ANGLE, K0L, E1, E2, KNL[11], KNSL[11];   // positve ANGLE bends particle to negative x dirction
  int Nint, Norder;
};

//----------------------------SOLEN-----------------------------------
class SOLEN: public Element
{
 public:
  SOLEN(string name, double l, double ks);
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  double  KS, Nint;    // KS = Bs/(brho ),  the same as MAD definition
};

//----------------------------Wiggler-----------------------------------
class WIGGLER: public Element
{
 public:
  WIGGLER(string name, double l, int nint, double b0, double kx, double kz, double phiz0);
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  int     Nint;
  double  B0, KX, KZ, PHIZ0;    
};


//-----------------------MATRIX-----------------------------------
class MATRIX: public Element
{
 public:
  MATRIX(string name, double l, double xco_in[6], double xco_out[6], double m66[36]);
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
    double  M66[36];
    double  XCO_IN[6], XCO_OUT[6];
};

//---------------------------KICKER-----------------------------------------
class KICKER: public Element
{
 public:
  KICKER(string name, double l, double hkick, double vkick);
  void SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  double  HKICK, VKICK;
};

//----------------------------HKICKER--------------------------------------
class HKICKER: public Element
{
 public:
  HKICKER(string name, double l, double hkick);
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  double  HKICK;   // kick in x'
};

//------------------------------VKICKER------------------------------------------
class VKICKER: public Element
{
 public:
  VKICKER(string name, double l, double vkick);
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  double  VKICK;  // kick in y'
};

//---------------------------HACMULT ( horizontal single order ac multipole )-----------------
class HACMULT: public Element
{
 public:
  HACMULT(string name, double l, int norder, double kl, int tturns, double phi0 );
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  int    Norder, TTURNS;              //  T period in unit of turns
  double KLMAX, PHI0;         //    KLMAX is kick amplitude
};

//----------------------------VACMULT (vertical single order ac multipole )----------------
class VACMULT: public Element
{
 public:
  VACMULT(string name, double l, int norder, double kl, int tturns, double phi0 );
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  int    Norder, TTURNS;            // TTURNs period in numbers of turns, 
  double KSLMAX, PHI0;                 // KSLMAX is kick amplitude
};

//----------------------------HACDIP-------------------------------------------
class HACDIP: public Element
{
 public:
  HACDIP(string name, double l, double hkickmax, double nud, double phid, int turns, int turne );
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  int     TURNS, TURNE;          // kick strength ramping between  TURNS and TURNE turns
  double  HKICKMAX, NUD, PHID;   // HKICKMAX, maximum kick angle
};

//------------------------------VACDIP------------------------------------------
class VACDIP: public Element
{
 public:
  VACDIP(string name, double l, double vkickmax, double nud, double phid, int turns, int turne );
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  int     TURNS, TURNE;
  double  VKICKMAX, NUD, PHID;  // VKICKMAX, maximum kick angle
};

//---------------------------------BPM-------------------------------
class BPM: public Element
{
 public:
  BPM(string name, double l);
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]); 
};

//----------------------------HBPM-----------------------------------
class HBPM: public Element
{
 public:
  HBPM(string name, double l);
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
};

//---------------------------------VBPM------------------------------
class VBPM: public Element
{
 public:
  VBPM(string name, double l);
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
};

//-----------------------------MARKER--------------------------
class MARKER: public Element
{
 public:
  MARKER(string name, double l);
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
};

//----------------------------RFCAV--------------------------------
class RFCAV: public Element
{
 public:
  RFCAV(string name, double l, double vrf, double frf, double phase0);
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  double VRF, FRF, PHASE0;
};

//---------------------CRAB CAVITY--------------------------------
class CRABRF: public Element
{
 public:
  CRABRF(string name, double l, double vrf, double frf, double phase0);
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  double VRF, FRF, PHASE0;
};

//---------------------CC MULTIPOLES--------------------------------
class CCMULT: public Element
{
 public:
  CCMULT(string name, double l, double vrf, double frf, double phase0, double b1, double a1, double b2, double a2, double b3, double a3, double b4, double a4);
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  double VRF, FRF, PHASE0;
  double B1, A1, B2, A2, B3, A3, B4, A4;
};

//---------------------Lorentz Boost Transfer-------------------------------
class LBT: public Element
{
 public:
  LBT(string name, double theta);
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  double THETA;
};

//---------------------Inverse Lorentz Boost Transfer-------------------------------
class ILBT: public Element
{
 public:
  ILBT(string name, double theta);
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  double THETA;
};

//------------------head-on beam-beam------------------------------
class BEAMBEAM: public Element
{
 public:
  BEAMBEAM(string name, int treatment, double np, double bbscale, double sigmal, int nslice, double emitx,  
	   double betax, double alfax, double emity,  double betay, double alfay);
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  int NSLICE, TREATMENT;
  double  NP, BBSCALE, EMITX, EMITY, SIGMAL, BETAX, ALFAX, BETAY, ALFAY;
};

//----LRBB ( 4-D long-range BEAMBEAM )------------
class LRBB: public Element
{
 public:
  LRBB(string name, double np, double bbscale, double sepx, double sepy, double sigmax, double sigmay);
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  double  NP, BBSCALE, SEPX, SEPY, SIGMAX, SIGMAY;
};

//----CRBB ( 6-d BB with crossing angle )
class CRBB: public Element
{
 public:
  CRBB(string name, double np, double theta, double bbscale, double sigmal, int nslice, 
       double emitx,  double betax, double alfax, double emity, double betay, double alfay);
      // emitx, emity:  un-normalized rms emittance,  sigma=SQRT[ emitx * betax ]
      // theta: half crossing angle  
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  int     NSLICE;
  double  NP, THETA, BBSCALE,SIGMAL, EMITX, EMITY, BETAX, ALFAX, BETAY, ALFAY;
};

//-----------------------------ERHICBB----------------------------------------------------
class ERHICBB: public Element
{
 public:
  ERHICBB(string name, double ne, double bbscale, double sigmax, double sigmay);
 //  sigmax, sigmay for electron slice at IP:  only for tune-shift calculation purpose.
 //  electron beam only has one slice, its sigma changes around IP  
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  double NE, BBSCALE, SIGMAX, SIGMAY;
};

//----------------------E-LENS------------------------------
class ELENS: public Element
{
 public:
  ELENS(string name, double l, double ne, double bbscale, int nslice, double betae, double sigmax, double sigmay);
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  double NE, BBSCALE,BETAE, SIGMAX, SIGMAY;
  int    NSLICE;
};

//----------------------HOLLOW E-LENS------------------------------
class HELENS: public Element
{
 public:
  HELENS(string name, double l, double ne, double bbscale, int nslice, double betae, double rout, double rin);
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  double NE, BETAE, BBSCALE, ROUT, RIN;
  int    NSLICE;
};

//---------------ROTAT  (spine rotator)-----------------------------------
class ROTAT: public Element
{
 public:
  ROTAT(string name, double l, double n[3], double angle);
  void   SetP(const char *name, double value); 
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
private:
  double N[3], ANGLE;  // ANGLE > 0, clockwise w.r.t. axis
};

//---------------SNAKE (spin snake )-----------------------------------
class SNAKE: public Element
{
 public:
  SNAKE(string name, double l, double n[3], double angle);
  void   SetP(const char *name, double value); 
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
private:
  double N[3], ANGLE; // ANGLE > 0, clockwise w.r.t. axis
};

//---------------------------DIFFUSE--------------------------------
class DIFFUSE: public Element
{
 public:
  DIFFUSE(string name, double diff_x, double diff_y, double diff_delta );
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  double DIFF_X, DIFF_Y, DIFF_DELTA;  
};

//---------------------------COOLING--------------------------------
class COOLING: public Element
{
 public:
  COOLING(string name, double alpha );
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  double ALPHA;  
};

//--------------------TRANS (coordinate system change)-------------------------------
class TRANS: public Element
{
 public:
  TRANS(string name, double dx, double dy, double ds);
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  double DX, DY, DS;
};

//--------------------SROTATION (coordinate system change)-------------------------------
class SROTAT: public Element
{
 public:
  SROTAT(string name, double psi);
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  double PSI;
};

//--------------------YROTATION (coordinate system change)-------------------------------
class YROTAT: public Element
{
 public:
  YROTAT(string name, double theta);
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  double THETA;
};

//--------------------XROTATION (coordinate system change)-------------------------------
class XROTAT: public Element
{
 public:
  XROTAT(string name, double phi);
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps x[6]);
  void   sPass(double x[9]);
 private:
  double PHI;
};

//--------------------PATCH (coordinate system change)-------------------------------
class PATCH: public Element
{
 public:
  PATCH(string name, double dx, double dy, double ds, double theta);
  void   SetP(const char *name, double value);
  double GetP(const char *name);
  void   Pass(double x[6]);
  void   DAPass(tps  x[6]);
  void   sPass(double x[9]);
 private:
  double DX, DY, DS, THETA;
};

//=========================================
//         
//           Line definition
//
//=========================================

class Line
{
 public:
  Line();
  ~Line();
  
  void Update();
  void Append(Element * x);
  void Delete(int i);
  void Insert( int i, Element * temp);
  void Replace( int i, Element * temp);
  void Rewind( int i );
  void Invert();
  void Empty();
  
  vector  <Element *> Cell;
  double Length;            //---length of reference orbit
  long   Ncell;             //---number of elements 
  double Tune1, Tune2, Tune3;  //----tunes
  double SpinTune;
  double Chromx1, Chromy1, Chromx2, Chromy2, Chromx3, Chromy3;  //   chromaticities
  double frev0;             //---decided by Length
  double Vrf_tot;           //---total rf voltage
  double Orbit_Length ;     //---length of closed orbit length.

  double Alfa0, Alfa1, Alfa2, Gammat, Slip;  //  transition parameters
  double Qs;                //  longitudinal tune
  double Bucket_length;     //  in unit of ns 
  double Bucket_height;     //  (dp/p0)_max
  double Bucket_area;       //  in ( phi_rf, dE/(hw_rev) ) phase space per nucleon
  double Bunch_length;      //  in unit of m, 6*sigma_l, in unit of ns
  double Bunch_area;        //  in ( phi_rf, dE/(hw_rev) ) phase space per nucleon
  double Bunch_height;      //  (dp/p0)_max for the given bunch area
};

/*----
void Line_Append_Element(Line & linename, Element * new_element);
void Line_Append_Elements(Line & linename1, Line & linename2);

void Line_Delete_Element(Line & linename, int k);
void Line_Delete_Elements(Line & linename, int i1, int i2);

void Line_Insert_Element(Line & linename, int k, Element * new_element1);
void Line_Insert_Elements(Line & linename1, Line & linename2, int i1);

void Line_Replace_Element(Line & linename, int k, Element * new_element1);

void Line_Rewind(Line & linename,  int k );
void Line_Invert(Line & linename);
void Line_Repeat(Line & linename1, Line & linename2, int n);
void Line_Connect(Line & linename1, Line & linename2, Line & linename3);
void Line_PickOut_Segment(Line & linename1, Line & linename2, int i1, int i2);
--*/

int  Count_One_Name_Element(Line & linename, const char * name);
int  Count_One_Type_Elements(Line & linename, const char * name);

int  Get_Index(Line & linename, const char * name, int k);

void Print_One_Type_Elements(Line & linename, const char * name);
void Print_All_Element_Types(Line & linename);

void Print_Sbend_Parameters(Line & linename, const char* name, int  ith);
void Print_Quad_Parameters(Line & linename, const char* name, int  ith);
void Print_Sext_Parameters(Line & linename, const char* name, int  ith);
void Print_Mult_Parameters(Line & linename, const char* name, int  ith);

void Print_All_Sbend_Parameters(Line & linename);
void Print_All_Quad_Parameters(Line & linename);
void Print_All_Sext_Parameters(Line & linename);
void Print_All_Mult_Parameters(Line & linename);

void Get_Mult_Parameters( Line & linename, int i, double knl[11], double knsl[11] );
void Get_Element_Parameter( Line & linename, const char* name, int ith, const char* para );

void Print_Line_Elements(Line & linename, const char * name);

void Set_Integration_Steps(double bl, double ql, double gl );
void Change_Integration_Steps(double bl, double ql, double gl );
void Change_Nint_Bend(Line & linename, double ds);
void Change_Nint_Mult(Line & linename, double ds);

void Split_Drift(Line & linename, int i, int m );
void Split_SBend(Line & linename, int i, int m );
void Split_Quad(Line & linename, int i, int m );
void Split_Sext(Line & linename, int i, int m );
void Split_Mult(Line & linename, int i, int m );
void Split_GMULT(Line & linename, int i, int m );
void Split_SBENDMULT(Line & linename, int i, int m );
void Split_SMult(Line & linename, int i, int m );

void Split_Drift_All(Line & linename,int m );
void Split_SBend_All(Line & linename, int m );
void Split_Quad_All(Line & linename, int m );
void Split_Sext_All(Line & linename, int m );
void Split_Mult_All(Line & linename, int m );
void Split_GMULT_All(Line & linename, int m );
void Split_SBENDMULT_All(Line & linename, int m );
void Split_SMult_All(Line & linename, int m );

void Concat_Drift(Line & linename);
void Clean_Up(Line & linename);
void Make_Thin(Line & linename);

double Get_KL(Line & linename, const char * name, const char *  kl);

void   Set_KL(Line & linename, const char * name, const char *  kl, double strength);
void   Set_KL(Line & linename, int  index,  const char *  kl, double strength);
void   Set_KL(Line & linename, vector<int> index_elements,  const char *  kl, double strength);

void   Set_dKL(Line & linename, const char * name, const char * kl, double dstrength);
void   Set_dKL(Line & linename, int index, const char * kl, double dstrength);
void   Set_dKL(Line & linename, vector<int> index_elements,  const char *  kl, double strength);
  
void     Cal_CC_Voltage(double theta, double freq_cc, double beta_star, double beta_cc, double energy, double & vrf);
double   Cal_CC_Voltage_To_Kick(double Vcc, double freq_cc, double energy);
double   Cal_CC_Kick_To_Voltage(double kick, double freq_cc, double energy);

void  Read_MADX_Lattice(Line & linename, const char * filename);
void  Print_MADX_Lattice(Line & linename, const char * filename);
void  Print_MADX_Lattice_v2(Line & linename, const char * filename);

void  Read_BMAD_Lattice(Line & linename, const char * filename);

void  Tokenize(const string& str, vector<string>& tokens);
void  Read_LatticeSequence(Line & linename, const char * filename);
void  Print_LatticeSequence(Line & linename, const char * filename);

//================================
//
//     beam dynamics
//   
//================================

void Survey(Line & linename, double X0, double Y0, double Theta0, int dir, const char* filename);

void Cal_Orbit_Num(Line & linename, double deltap);
void Cal_Orbit_Num(Line & linename, double z, double deltap);
void Cal_OneTurnMap(Line & linename, double deltap);
void Cal_OneElementMap(Element * temp_element, double x[], double t66[]);
void Cal_ElementMap(Line & linename, double deltap);
void Get_ElementMap(Line & linename, int i1, double m66[]);
void Cal_SectionMap(Line & linename, int i1, int i2, double deltap, double t66[36] );
void Cal_SectionMap(Line & linename, int i1, int i2, double x0[], double t66[36] );

void Cal_A(Line & linename, double deltap);
void Trace_A(Line & linename, double deltap);
void Cal_Twiss(Line & linename, double deltap);

void Cal_OneTurnMap_Z(Line & linename, double z0);
void Cal_ElementMap_Z(Line & linename, double z0);
void Cal_Twiss_Z(Line & linename, double z0);
void Cal_Tune_vs_Z(Line & linename, const char *filename);
void Cal_Beta_Star_vs_Z(Line & linename, const char *filename);  

void Trace_Orbit(Line & linename, int istart, double x[]);
void Trace_Twiss(Line & linename, double deltap, double x[4], double Beta1, double Beta2, double Alfa1, double Alfa2, double c11, double c12, double c21, double c22);

void Cal_Orbit_Num_6D(Line & linename);
void Cal_OneTurnMap_6D(Line & linename);
void Cal_OneTurnMap_Element_6D(Line & linename, int k);
void Cal_ElementMap_6D(Line & linename);
void Cal_A_6D(Line & linename);
void Cal_A_6D_Rev(Line & linename); 
void Trace_A_6D(Line & linename);
void Cal_Twiss_6D(Line & linename);

void Cal_A_4D(double m[], double a[]);
void Cal_A_6D(double m[], double a[]);

void Cal_Chrom(Line & linename);
void Cal_Chrom_Dispersion_Num( Line & linename);

void Cal_Dispersion(Line & linename);
void Cal_Momentum_Dispersion(Line & linename);
void Trace_Momentum_Dispersion(Line & linename, int i1, int i2, double Eta1[4]);
void Cal_Timeflight_Dispersion(Line & linename);
void Trace_Timeflight_Dispersion(Line & linename, int i1,  int  i2, double Ksi1[4]);

void Cal_Optics(Line & linename); 

void Fit_Tune(Line & linename, double q1, double q2, const char * qf_name, const char * qd_name);
void Fit_Tune_RHICelens(Line & linename, double q1, double q2);

double  Get_RHIC_QF_K1L(Line & linename);
double  Get_RHIC_QD_K1L(Line & linename);
void    Set_RHIC_QF_dK1L(Line & linename, double dk1l);
void    Set_RHIC_QD_dK1L(Line & linename, double dk1l);
void    Fit_Tune_RHIC(Line & linename, double q1, double q2 );

void Fit_Chrom(Line & linename, double chrom1x_want, double chrom1y_want, const char * sf_name, const char * sd_name );
void Fit_Chrom_RHIC8fam(Line & linename, double chrom1x_want, double chrom1y_want );
void Fit_Chrom_RHIC(Line & linename, double chrom1x_want, double chrom1y_want );

void chrom_fit(double qx[],double qy[],double & chromx1,double & chromy1,double & chromx2,double & chromy2,double & chromx3,double & chromy3 );
void Cal_Chrom_Num( Line & linename);
void Correct_Chrom_Manual( Line & linename);

void Cal_Tune_vs_Deltap(Line & linename, const char *filename);
void Plot_Tune_vs_Deltap(Line & linename, const char* filename);
void Cal_Beta_Star_vs_Deltap(Line & linename, const char *filename);
void Plot_Beta_Star_vs_Deltap(Line & linename, const char* filename);

void Cal_Beta_vs_Deltap(Line & linename, const char *filename);
void Cal_Dispersion_vs_Deltap(Line & linename, const char *filename);
void Cal_Chromatic_Functions(Line & linename, const char* filename );
void Cal_W_Functions(Line & linename, const char *filename);

void Cal_Half_Integer_RDT(Line & linename,  const char* filename);
void Cal_Half_Integer_RDT_SextFamily(Line & linename,  const char* sextname);
void Cal_Q2_Source(Line & linename, const char* filename );
void Cal_Q2_Source_Section(Line & linename, const char* filename );

void Cal_Coupling_Coefficient( Line & linename );
void Cal_Coupling_Coefficient_Source( Line & linename, const char * filename);
void Cal_BeamSize_SigmaMatrx(Line & linename, const char* filename );
void Twiss_Propagation_Drift(double  L,  double beta0[], double alfa0[], double beta1[], double alfa1[], double dphase[] );
void Twiss_Propagation_Matrix(double  T[],  double beta0[], double alfa0[], double beta1[], double alfa1[], double dphase[] );
void Coupling_From_Solenoid(double L, double KS, int nstep, double beta0[], double alfa0[], double phasex0, double phasey0, double & creal, double & cimag);
void Cal_Coupling_Coefficient_Updated( Line & linename, double & creal, double & cimag );
void Cal_Coupling_One_Location(Line linename, int  i, double & creal, double & cimag );
void Cal_Coupling_Along_Ring(Line linename, const char* filename );
void Cal_Coupling_Along_Ring_Old(Line linename, const char* filename );
void Harmonic_Analysis(double xtemp[], int nturn, double tune, double & peak, double &phase );
void Meas_Coupling_Along_Ring(Line linename, int nbpm, int nturn, int m, const  char* filename) ;
void Meas_Coupling_Along_Ring_Old(Line linename , int nbpm, int nturn) ;

void SaveData_To_File(double xtemp[], int nturn, const char*  filename );
void FFT_Save_File(double xtemp[], int nturn, int m,  const char*  filename );
template <class T> void SOLEN_Pass_Old1(T x[], double L, int Nint, double KS);

void Extract_Orbits(Line linename, double holder[] );
void Make_Orbit_Difference(Line linename, double orbit0[], double orbit1[], double diff[] );
void Save_Orbit_Difference(Line linename, double orbit0[], double orbit1[], const char* filename);

void Cal_Sext_RDTs( Line & linename );
void Cal_3Qx_RDTs_Ring( Line & linename,const char * filename);
void Cal_3Qx_RDTs_Source( Line & linename );
void Cal_SkewSext_RDTs( Line & linename ) ;
void Cal_3Qy_RDTs_Ring( Line & linename,const char * filename);
void Cal_3Qy_RDTs_Source( Line & linename );
void Cal_Detuning_Sext( Line & linename );

void Cal_Detuning_Oct( Line & linename );
void Cal_Q2_Oct( Line & linename )  ;

double Cal_Pathlength( Line & linename, double deltap);
void   Cal_Gammat(Line & linename);
void   Cal_Orbit_Length(Line & linename, double deltap);
void   Cal_Qs(Line & linename);
void   Cal_Bucket_Area(Line & linename);
double RF_F_function(double phi_s, double phi_right, double phi_left);
void   Cal_Bunch_Area(Line & linename, double full_length);
void   Cal_Bunch_Height(Line & linename, double bunch_area);

void Print_Optics(Line & linename,  int  i);
void Print_Orbit(Line & linename, const char* filename);
void Print_Twiss(Line & linename, const char* filename);
void Print_Twiss_6D(Line & linename, const char* filename);
void Print_Twiss_Coupling(Line & linename, const char* filename);
void Print_Dispersion(Line & linename, const char* filename);
void Print_A_Matrix(Line & linename, const char* filename);

void Print_Optics_Summary(Line & linename);
void Print_Longitudinal_Summary( Line & linename);
void Plot_Twiss(Line & linename);
void Plot_Orbit(Line & linename);

void Print_PhaseAdvances(Line & linename, int index1, int index2);
void Add_Phaser(Line & linename, int loc, const char * name,  double mux, double muy);

void   Correct_Orbit_SVD(Line & linename, int m, int n, vector<int> bpm_index, vector<int> kicker_index, int plane);
void   Local_Three_Bump(Line linename, int plane,  const char *corr1,   const char *corr2,   const char *corr3, double kick1);
double RMS_Leakage_Orbit( Line linename, int plane, int i1, int i2 );
void   Correct_Orbit_SlidingBump1(Line & linename, int m, int n, vector<int> bpm_index, vector<int> kicker_index, int plane);
void   Correct_Orbit_SlidingBump2(Line & linename, int m, int n, vector<int> bpm_index, vector<int> kicker_index, int plane);
void   Orbit_Status( Line linename, vector<int> bpm_index, int plane, double &orbit_mean, double &orbit_max, double & orbit_rms );

void  Cal_Spin_Orbit(Line & linename, double deltap);
void  Cal_Spin_Tune(Line & linename, double deltap);
void  Cal_Spin_Resonance(Line & linename, double gamma1, double gamma2, double  tunex, double tuney, double tunes);
void  Set_Initial_Spin(Line & linename, int Nturn, double n0[3], double x0[9]);

void  Cal_SRLoss_U0rad( Line & linename);
void  Cal_SRLoss_U0rad_Track( Line & linename);
void  Cal_SRLoss_Particle(Line & linename, double x[], double &loss);

//=============================
//
//      Particle Tracking 
//
//=============================

void Track(Line & linename, double x[], int nturn, int & stable, int & lost_turn, int & lost_post);
void Track_spin(Line & linename, double x[], int nturn, int & stable, int & lost_turn, int & lost_post);
void Track_tbt(Line & linename, double x[], int nturn, double  x_tbt[], int & stable, int & lost_turn, int & lost_post);
void Track_tbt_spin(Line & linename, double x[], int nturn, double  x_tbt[], int & stable, int & lost_turn, int & lost_post);
void Track_tbt(Line & linename, double x[], int nturn, double  x_tbt[], int &bpm_index, int & stable, int & lost_turn, int & lost_post);
void Track_tbt_spin(Line & linename, double x[], int nturn, double  x_tbt[], int &bpm_index, int & stable, int & lost_turn, int & lost_post);

void Cal_dz0_OnMomentumPart(Line & linename, double x0[], double dz0[]);
void Track_wrt_OnMomentumPart(Line & linename, double x[], double dz0[],  int nturn, int & stable, int & lost_turn, int & lost_post);

void Cal_Tunes_Track(Line & linename, double x[], int Nturn, double & tune1,  double & tune2);
void Cal_Tunes_Track(Line & linename, double deltap0);
void Cal_Chrom_Track(Line & linename);

void Track_tbt_FFT( Line & linename, double x0[], int Nturn, int m);
void Track_tbt_tune_footprint( Line & linename, double deltap0, double emitx, double emity) ;
void Track_tbt_FMA( Line & linename, double deltap0, double sigmax0, double sigmay0 );
void Track_tbt_Lyapunov( Line & linename, double deltap0, double sigmax0, double sigmay0 );
void Track_tbt_ActionDiff( Line & linename, double nsigma, double phase, double emitx0, double emity0, int nturn, int nstep, int npart);

void Track_DA( Line & linename, int nturn, double deltap0, double sigmax0, double sigmay0,  const char* filename);
void Track_DA1( Line & linename, int nturn, double deltap0, double sigmax0, double sigmay0, const char* filename);
void Track_DA2( Line & linename, int nturn, double deltap0, double sigmax0, double sigmay0, const char* filename);
void Track_DA3( Line & linename, int nturn, double deltap0, double sigmax0, double sigmay0, const char* filename);
void Track_DA_Uniform( Line & linename, int nturn, double deltap0, double angle, double sigmax0, double sigmay0, const char* filename );

void Track_Radial(Line & linename, double x[], int nturn, int & stable, int & lost_turn, int & lost_post);
void Track_DA_Radial( Line & linename, int nturn, double deltap0, double sigmax0, double sigmay0,  const char* filename);
void Track_DA_Radial_Uniform( Line & linename, int nturn, double deltap0, double angle, double sigmax0, double sigmay0, const char* filename );

void Track_Spin(Line & linename, double x[], int nturn, int & stable, int & lost_turn, int & lost_post);

/*-----
void Track_Fast_Prepare(Line & linename);
void Track_Fast(Line & linename, double x[6], int nturn, int & stable, int & lost_turn, int & lost_post);
void Track_Fast_Emit(Line & linename, double x[6], int nturn, int & stable, int & lost_turn, int & lost_post, double & sum_x2, double & sum_y2, double & sum_z2);

void Track_Fast_DA( Line & linename, int nturn, double deltap0, double sigmax0, double sigmay0, char* filename );
void Track_Fast_DA1( Line & linename, int nturn, double deltap0, double sigmax0, double sigmay0, char* filename );
void Track_Fast_DA2( Line & linename, int nturn, double deltap0, double sigmax0, double sigmay0, char* filename );
void Track_Fast_Spin(Line & linename, double x[], int nturn, int & stable, int & lost_turn, int & lost_post);
----*/

void BeamBeam4D_SoftGaussian(double Np1, int Npart_1, double PartDist_1[], 
                             double Np2, int Npart_2, double PartDist_2[], 
                             double gamma1, double gamma2, double bbscale1, double bbscale2);
void BeamBeam6D_SoftGaussian(double Np1, int Mp1, double x1[], double Np2, int Mp2, double x2[], 
                             double gamma1, double gamma2, double bbscale1, double bbscale2,
			     int Nslice1, int Nslice2);
void BeamBeam6D_SoftGaussian(double Np1, int Mp1, double x1[], double Np2, int Mp2, double x2[], 
			     double gamma1, double gamma2, double bbscale1, double bbscale2,
			     int Nslice1, int Nslice2, double zb1[], double zb2[]);


void BeamBeam4D_PIC(double Np1, int Mp1, double xc1[], double Np2, int Mp2, double xc2[], 
                    double gamma1, double gamma2, double bbscale1, double bbscale2);
void BeamBeam6D_PIC(double Np1, int Mp1, double x1[], double Np2, int Mp2, double x2[], 
                    double gamma1, double gamma2, double bbscale1, double bbscale2,
		    int Nslice1, int Nslice2);

double  Lumi_Cal_PIC_4D(double Np1, int Mp1, double Np2, int Mp2, double x1[], double x2[], int Ncoll, double freq);
double  Lumi_Cal_PIC_6D(double Np1, int Mp1, double Np2, int Mp2, int Nslice1, int Nslice2, double x1[], double x2[], 
                        double  gamma1, double gamma2,  double bbscale1, double bbscale2, int Ncoll, double freq);
double  Lumi_Cal_PIC_6D_Angle(double Np1, int Mp1, double Np2, int Mp2, double theta, int Nslice1, int Nslice2, double x1[], double x2[],
                              double  gamma1, double gamma2,  double bbscale1, double bbscale2, int Ncoll, double freq);

void BeamBeam4D_PIC_NonOffset(double Np1, int Mp1, double xc1[], double Np2, int Mp2, double xc2[],
		    double gamma1, double gamma2, double bbscale1, double bbscale2,
		    int Nx, int Ny, double hx, double hy, double grn_c_r[], double  grn_c_i[] );
void BeamBeam6D_PIC_NonOffset(double Np1, int Mp1, double x1[], double Np2, int Mp2, double x2[], 
                              double gamma1, double gamma2, double bbscale1, double bbscale2,
			      int Nslice1, int Nslice2, double zb1[], double zb2[], int Nx, int Ny, double hx, double hy, double grn_c_r[], double  grn_c_i[] );
void BeamBeam6D_PIC_Interpolation(double Np1, int Mp1, double x1[], double Np2, int Mp2, double x2[], 
                                  double gamma1, double gamma2, double bbscale1, double bbscale2,
			          int Nslice1, int Nslice2, double zb1[], double zb2[], int Nx, int Ny, double hx, double hy, double grn_c_r[], double  grn_c_i[] );

#endif

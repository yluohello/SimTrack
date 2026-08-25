#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include "simtrack.h"

using namespace  std;

int main()
{
  int i,j,k;
  int loc_ip6, loc_ip8, loc_ip10, loc_elens,loc_rf;
  Line  rhic;
  Element * temp_element;

  //----set global parameters  for Au

  Set_RefPartEnergy( 100, 938.27201*197, 79, 197); 

  //---read in lattices
  Read_MADX_Lattice(rhic, "./parameters_input");
   
  //---check optics
  if(true){
    Cal_Twiss(rhic,0);
    Print_Orbit(rhic,"./orbit0.dat");
    Print_Twiss(rhic,"./twiss0.dat");
    Cal_Chrom(rhic);
    Cal_Dispersion(rhic);
    Print_Optics_Summary(rhic);
    
    Cal_Beta_Star_vs_Deltap(rhic, "./beta_vs_deltap.dat");
    Cal_Tune_vs_Deltap(rhic, "./tunes_vs_deltap.dat");
  }

  //---install RF
  if(true) {
    for(i=0; i<rhic.Ncell;i++){
      if( rhic.Cell[i]->TYPE=="RFCAV" ) {
	loc_rf=i;  break;
      } }
    rhic.Cell[loc_rf]->SetP("VRF",4.0);
    rhic.Cell[loc_rf]->SetP("FRF",78250.42279*GP.harm);
  }

  //************************************
  //  
  //      DA tracking 
  //
  //************************************
  
  int    nturn=1000000;
  double deltap =0.0015;
  double sigmax0=sqrt(2.5e-06 * 0.65  / GP.gamma);
  double sigmay0=sqrt(2.5e-06 * 0.65  / GP.gamma);  
  Track_DA( rhic, nturn, deltap, sigmax0, sigmay0, "DA-output.dat");

  //------the end

  return 0;

}

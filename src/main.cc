#include <TFile.h>

#include "fitter.h"

int main(){
  auto file = std::make_unique<TFile>( "/home/mikhail/bmn_run8/qa.recent_vf_2024_06_03.root", "READ" );
  TH2* obj{nullptr};
  file->GetObject( "h2_pq_mass2_tof400", obj );
  auto plot400 = M2PqPlot{};
  plot400.SetH2M2Pq( dynamic_cast<TH2D*>(obj->Clone("tof400")) );
  
  plot400.AddFitter( Particle().PdgCode(2212).M2Lo(0.6).M2Hi(1.4).PqLo(0.5).PqHi(8) );
  plot400.AddFitter( Particle().PdgCode(211).M2Lo(-0.3).M2Hi(0.4).PqLo(1.0).PqHi(2.5) );
  plot400.AddFitter( Particle().PdgCode(321).M2Lo(0.1).M2Hi(0.4).PqLo(1.0).PqHi(2.5) );
  plot400.AddFitter( Particle().PdgCode(1000010020).M2Lo(3.0).M2Hi(5.0).PqLo(1.0).PqHi(8.) );
  plot400.AddFitter( Particle().PdgCode(1000010030).M2Lo(7.0).M2Hi(9.0).PqLo(2.0).PqHi(6.) );
  plot400.AddFitter( Particle().PdgCode(1000020030).M2Lo(1.5).M2Hi(2.5).PqLo(1.0).PqHi(3.0) );
  
  plot400.Fit( 60, 1.0, 7 );

  plot400.Print( "file_out_400.pdf" );
  plot400.SaveApproximations( "identification_400.root" );
  plot400.Dump( "file_out_400.root" );

  file->GetObject( "h2_pq_mass2_tof700", obj );
  auto plot700 = M2PqPlot{};
  plot700.SetH2M2Pq( dynamic_cast<TH2D*>(obj->Clone("tof700")) );
  
  plot700.AddFitter( Particle().PdgCode(2212).M2Lo(0.6).M2Hi(1.4).PqLo(0.5).PqHi(8) );
  plot700.AddFitter( Particle().PdgCode(211).M2Lo(-0.3).M2Hi(0.4).PqLo(1.0).PqHi(2.5) );
  // plot700.AddFitter( Particle().PdgCode(321).M2Lo(0.1).M2Hi(0.4).PqLo(1.0).PqHi(2.5) );
  plot700.AddFitter( Particle().PdgCode(1000010020).M2Lo(3.0).M2Hi(5.0).PqLo(1.0).PqHi(8.) );
  plot700.AddFitter( Particle().PdgCode(1000010030).M2Lo(7.0).M2Hi(9.0).PqLo(2.0).PqHi(6.) );
  plot700.AddFitter( Particle().PdgCode(1000020030).M2Lo(1.5).M2Hi(2.5).PqLo(1.0).PqHi(3.0) );
  
  plot700.Fit( 60, 1.0, 7 );

  plot700.Print( "file_out_700.pdf" );
  plot700.SaveApproximations( "identification_700.root" );
  plot700.Dump( "file_out_700.root" );

  
  return 0;
}
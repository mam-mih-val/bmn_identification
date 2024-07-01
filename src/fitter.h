#ifndef FITTER_H
#define FITTER_H

#include <algorithm>
#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <TH1D.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TPDF.h>
#include <TFile.h>
#include <TCanvas.h>

struct Particle{
  Particle() = default;
  Particle& PdgCode( int pdg ){ pdg_code = pdg; return *this; }
  Particle& PqLo( double val ){ pq_lo = val; return *this; }
  Particle& PqHi( double val ){ pq_hi = val; return *this; }
  Particle& M2Lo( double val ){ m2_lo = val; return *this; }
  Particle& M2Hi( double val ){ m2_hi = val; return *this; }

  int pdg_code;
  double pq_lo;
  double pq_hi;
  double m2_lo;
  double m2_hi;
};

class Approximation{
public:
  Approximation( int pdg ) : pdg_code_(pdg), gaus_( new TF1( "func", "gaus" ) ) {}
  ~Approximation() = default;
  Approximation(Approximation&&) noexcept = default;
  Approximation& operator=(Approximation&&) noexcept = default;
  
  Approximation& SetMean( std::unique_ptr<TF1> func ){ f1_mean_ = std::move(func); return *this; }
  Approximation& SetSigma( std::unique_ptr<TF1> func ){ f1_sigma_ = std::move(func); return *this; }
  Approximation& SetAmplitude( std::unique_ptr<TF1> func ){ f1_amplitude_ = std::move(func); return *this; }
  double Eval( double pq, double m2 ) const {
    auto m = f1_mean_->Eval(pq);
    auto s = f1_sigma_->Eval(pq);
    auto a = f1_amplitude_->Eval(pq);
    gaus_->SetParameter(0, a);
    gaus_->SetParameter(1, m);
    gaus_->SetParameter(2, s);
    return gaus_->Eval(m2);
  }
  TF1* GetApproximation(double pq) const { 
    auto m = f1_mean_->Eval(pq);
    auto s = f1_sigma_->Eval(pq);
    auto a = f1_amplitude_->Eval(pq);

    auto func = new TF1( std::data( std::to_string( pdg_code_ ).append("_").append( std::to_string(pq) ) ), "gaus", -1, 10 );
    func->SetParameter( 0, a );
    func->SetParameter( 1, m );
    func->SetParameter( 2, s );

    return func;
  }
  int GetPdg() const { return pdg_code_; }
  void Save(){
    f1_mean_->Write();
    f1_sigma_->Write();
    f1_amplitude_->Write();
  }
private:
  int pdg_code_;
  std::unique_ptr<TF1> gaus_{};
  std::unique_ptr<TF1> f1_mean_;
  std::unique_ptr<TF1> f1_sigma_;
  std::unique_ptr<TF1> f1_amplitude_;
};

class ApproximationPlot{
public:
  ApproximationPlot() = default;
  TH2D* CalculatePurityPlot( int pdg_code ){
    auto result = new TH2D( std::to_string(pdg_code).c_str(), ";pq;m2", 500, 0, 10, 500, -1, 10 );
    for( size_t pq_bin=1; pq_bin <= result->GetXaxis()->GetNbins(); ++pq_bin ){
      auto pq = result->GetXaxis()->GetBinCenter(pq_bin);
      for( size_t m2_bin=1; m2_bin <= result->GetYaxis()->GetNbins(); ++m2_bin ){
        auto m2 = result->GetYaxis()->GetBinCenter(m2_bin);
        const auto& approx_pdg = std::find_if( approximations_.begin(), approximations_.end(), [pdg_code](Approximation& a){ return a.GetPdg() == pdg_code; } );
        auto val_pdg = approx_pdg->Eval(pq, m2);
        auto val_rest = 0.0;
        for( auto& a : approximations_) {
          val_rest+=a.Eval(pq, m2);
        }
        result->SetBinContent(pq_bin, m2_bin, val_pdg/val_rest);
      }
    }
    return result;
  }
  void AddApproximation( int pdg, Approximation approx ){ approximations_.emplace_back( std::move(approx) ); }
private:
  std::vector<Approximation> approximations_{};
};

class Fitter{
public:
  Fitter(const Particle& p) : 
    pdg_code_(p.pdg_code), 
    m2_lo_(p.m2_lo), 
    m2_hi_(p.m2_hi),
    pq_lo_(p.pq_lo),
    pq_hi_(p.pq_hi) {  }
  Fitter(int pdg_code, double min, double max) : pdg_code_(pdg_code), m2_lo_(min), m2_hi_(max) {  }
  void Fit( TH1D* h1_m2, double pq ){
    if( pq > pq_hi_ )
      return;
    if( pq < pq_lo_ )
      return;
    pq_.push_back(pq);
    auto fit_name  = std::string( std::to_string(pdg_code_) ).append( "_" ).append( std::to_string(pq) );
    fit_function_.emplace_back( new TF1( fit_name.c_str(), "gaus" ) );
    h1_m2->Fit( fit_function_.back().get(), "", "", m2_lo_, m2_hi_ );
    amplitude_.push_back( fit_function_.back()->GetParameter( 0 ) );
    amplitude_err_.push_back( fit_function_.back()->GetParError( 0 ) );
    
    mean_.push_back( fit_function_.back()->GetParameter( 1 ) );
    mean_err_.push_back( fit_function_.back()->GetParError( 1 ) );
    
    sigma_.push_back( fit_function_.back()->GetParameter( 2 ) );
    sigma_err_.push_back( fit_function_.back()->GetParError( 2 ) );
  }
  void Dump( ){
    UpdateGraphs();
    g1_mean_->Write( std::data( std::string("g1_mean_").append( std::to_string(pdg_code_) ) ) );
    g1_sigma_->Write( std::data( std::string("g1_sigma_").append( std::to_string(pdg_code_) ) ) );
    g1_amplitude_->Write( std::data( std::string("g1_amplitude_").append( std::to_string(pdg_code_) ) ) );
  }

  Approximation Approximate() {
    UpdateGraphs();
    auto f1_mean = std::make_unique<TF1>( std::string("mean_").append(std::to_string(pdg_code_)).c_str(), "pol1" );
    auto f1_sigma = std::make_unique<TF1>( std::string("sigma_").append( std::to_string(pdg_code_)).c_str(), "pol2" );
    auto f1_amplitude = std::make_unique<TF1>( std::string("amplitude_").append(std::to_string(pdg_code_)).c_str(), "[0] * pow(x, [1]) * exp( [2]*x )" );
    f1_amplitude->SetParLimits(1, 1, 10);
    f1_amplitude->SetParLimits(2, -999., 0);
    g1_mean_->Fit( f1_mean.get() );
    g1_sigma_->Fit( f1_sigma.get() );
    g1_amplitude_->Fit( f1_amplitude.get() );
    auto approx = Approximation(pdg_code_);
    approx.SetMean( std::move(f1_mean) );
    approx.SetSigma( std::move(f1_sigma) );
    approx.SetAmplitude( std::move(f1_amplitude) );
    return approx;
  }

  int GetPdgCode(){ return pdg_code_; }

private:
  void UpdateGraphs(){
    if( !g1_mean_ )
      g1_mean_.reset( new TGraphErrors( pq_.size(), pq_.data(), mean_.data(), nullptr, mean_err_.data() ) );
    if( !g1_sigma_ )
    g1_sigma_.reset( new TGraphErrors( pq_.size(), pq_.data(), sigma_.data(), nullptr, sigma_err_.data() ) );
    if( !g1_amplitude_ )
    g1_amplitude_.reset( new TGraphErrors( pq_.size(), pq_.data(), amplitude_.data(), nullptr, amplitude_err_.data() ) );
  }

  std::vector< std::unique_ptr<TF1> > fit_function_{};
  
  int pdg_code_{};
  double m2_lo_{0.};
  double m2_hi_{0.};
  double pq_lo_{0.};
  double pq_hi_{0.};

  std::vector<double> pq_{};

  std::vector<double> mean_{};
  std::vector<double> mean_err_{};
  
  std::vector<double> sigma_{};
  std::vector<double> sigma_err_{};
  
  std::vector<double> amplitude_{};
  std::vector<double> amplitude_err_{};
  
  std::unique_ptr<TGraphErrors> g1_mean_{};
  std::unique_ptr<TGraphErrors> g1_sigma_{};
  std::unique_ptr<TGraphErrors> g1_amplitude_{};
};

class M2PqPlot{
public:
  M2PqPlot() = default;
  ~M2PqPlot() = default;

  void SetH2M2Pq( TH2D* histo ){ h2_m2_pq_.reset(histo); }
  void AddFitter( Particle p ){ fitters_.emplace_back( p ); }
  void Fit( size_t n_points, double min, double max ){
    auto first_bin = h2_m2_pq_->GetXaxis()->FindBin( min );
    auto last_bin = h2_m2_pq_->GetXaxis()->FindBin( max );
    auto step = static_cast<size_t>( ( last_bin - first_bin ) / n_points );
    for( size_t bin = first_bin; bin < last_bin; bin+=step ){
      auto pq_lo = h2_m2_pq_->GetXaxis()->GetBinCenter( bin );
      auto pq_hi = h2_m2_pq_->GetXaxis()->GetBinCenter( bin+step );
      auto pq_center = ( pq_hi + pq_lo ) / 2.0;
      auto proj_name = std::string("h1_m2_pq_").append( std::to_string( pq_center ) );
      h1_m2_.emplace_back( h2_m2_pq_->ProjectionY( proj_name.c_str(), bin, bin+step ) );
      pq_slices_.emplace_back(pq_center);
      for( auto& fitter : fitters_ ){ fitter.Fit( h1_m2_.back().get(), pq_center ); }
    }
    for( auto& fitter : fitters_ ){ approximation_plot_.AddApproximation( fitter.GetPdgCode(), std::move(fitter.Approximate()) ); }
  }
  void Dump( std::string file_name ){
    auto file = new TFile( file_name.c_str(), "RECREATE" );
    file->cd();
    h2_m2_pq_->Write();
    for( auto& fitter : fitters_ ){ fitter.Dump(); }
    for( auto& h1 : h1_m2_ ){ h1->Write(); }

    for( auto& fitter : fitters_){
      auto hist = approximation_plot_.CalculatePurityPlot(fitter.GetPdgCode());
      hist->Write();
    }
    file->Close();
  }
  void SaveApproximations(std::string file_name){
    auto file = new TFile( file_name.c_str(), "RECREATE" );
    file->cd(); std::vector<Approximation> approximations{};
    for( auto& fitter : fitters_ ){ fitter.Approximate().Save(); }
    file->Close();
  }
  void Print( std::string file_name ){
    std::cout << "Here" << std::endl;
    std::vector<Approximation> approximations{};
    for( auto& fitter : fitters_ ){ approximations.push_back( fitter.Approximate() ); }
    auto pdf = new TPDF(file_name.data());
    int i=0;
    for( const auto& h1 : h1_m2_ ){
      h1->Draw();
      for( const auto& a : approximations ){
        auto f = a.GetApproximation(pq_slices_.at( i ));
        f->Draw("same");
      }
      gPad->SetLogy();
      gPad->Update();
      ++i;
    }
    pdf->Close();
  }

private:
  std::vector<Fitter> fitters_{};
  ApproximationPlot approximation_plot_{};
  std::unique_ptr<TH2D> h2_m2_pq_{};
  std::vector< double > pq_slices_{};
  std::vector< std::unique_ptr<TH1D> > h1_m2_{};
};


#endif // FITTER_H
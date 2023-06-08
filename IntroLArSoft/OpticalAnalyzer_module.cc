////////////////////////////////////////////////////////////////////////
// Class:       OpticalAnalyzer
// Plugin Type: analyzer (Unknown Unknown)
// File:        OpticalAnalyzer_module.cc
//
// Generated at Thu Aug  4 08:47:42 2022 by Rodrigo Alvarez Garrote using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

// include art
// #include "art/Framework/Core/EDAnalyzer.h"
// #include "art/Framework/Core/ModuleMacros.h"
// #include "art/Framework/Principal/Event.h"
// #include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
// #include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Persistency/Common/Ptr.h"
// #include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"

// include larsoft
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/CoreUtils/ServiceUtil.h"

// include ROOT
#include "TH1.h"
#include "TH1F.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TTimeStamp.h"
#include "TString.h"

namespace ANA {
  class OpticalAnalyzer;
}

class ANA::OpticalAnalyzer : public art::EDAnalyzer {
public:
  explicit OpticalAnalyzer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  OpticalAnalyzer(OpticalAnalyzer const&) = delete;
  OpticalAnalyzer(OpticalAnalyzer&&) = delete;
  OpticalAnalyzer& operator=(OpticalAnalyzer const&) = delete;
  OpticalAnalyzer& operator=(OpticalAnalyzer&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  void reconfigure(fhicl::ParameterSet const & p);

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  std::string fOpDetWaveformLabel;
  
  // Declare member data here.
  TTree *fMyAnaTree;
  TTree *fOpDetWaveforms;
  
  // MCTruth
  int f_pdg;
  // float f_E;
  double f_nueQ;
  double f_nueE;
  // float f_M; // Mass; from PDG unless overridden Should be in GeV

  double fx_i;
  double fy_i;
  double fz_i;

  double fx_f;
  double fy_f;
  double fz_f;
  int f_N;

  // std::map< int, TH1D* > averageWaveforms;
  std::vector<raw::ADC_Count_t> wvf;
  float fSampleFreq;
  // std::map< int, int   > waveformCount;
  // TH1D* averageWaveformAll;
  int eventCount;
};


ANA::OpticalAnalyzer::OpticalAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer(p){this->reconfigure(p);}  // ,
  // More initializers here.
  //......................................................
void ANA::OpticalAnalyzer::reconfigure(fhicl::ParameterSet const & p){
  fOpDetWaveformLabel = p.get<std::string> ("OpDetWaveformLabel");
} // Reconfigure
  
void ANA::OpticalAnalyzer::beginJob()
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  art::ServiceHandle<art::TFileService> tfs;
  //--------------------------------Gen MCTruth-----------------------------
  fMyAnaTree = tfs->make<TTree>("MCTruth", "MCTruth");
  fOpDetWaveforms = tfs->make<TTree>("OpDetWvfs", "OpDetWvfs");
  fMyAnaTree->Branch("Part#", &f_N);
  fMyAnaTree->Branch("pdg", &f_pdg);
  // fMyAnaTree->Branch("E", &f_E);
  fMyAnaTree->Branch("nueQ", &f_nueQ);
  fMyAnaTree->Branch("nueE", &f_nueE);
  // fMyAnaTree->Branch("M", &f_M);
  fMyAnaTree->Branch("x_i", &fx_i);
  fMyAnaTree->Branch("y_i", &fy_i);
  fMyAnaTree->Branch("z_i", &fz_i);
  fMyAnaTree->Branch("x_f", &fx_f);
  fMyAnaTree->Branch("y_f", &fy_f);
  fMyAnaTree->Branch("z_f", &fz_f);
  fOpDetWaveforms -> Branch("Sampling",  &fSampleFreq, "Sampling/F");
  fOpDetWaveforms -> Branch("OpDetWvfs",  &wvf);

} //beginJob

void ANA::OpticalAnalyzer::analyze(art::Event const& evt)
{
  // Implementation of required member function here.

  art::Handle<std::vector<simb::MCTruth>> mctruths;

  //--------------------------------MCTruth--------------------------------
  evt.getByLabel("marley", mctruths);
  int EvCounter = 0; 
  for (auto const &truth : *mctruths)
  { 
    std::cout << "\nProcessing event" << std::endl;
    int N = truth.NParticles();
    const simb::MCNeutrino &nue = truth.GetNeutrino();
    // std::ostream &operator << (std::ostream &output, const simb::MCNeutrino &nue);
    f_nueQ = sqrt(nue.QSqr());
    f_nueE = nue.Nu().E();
    std::cout << "MC Neutrino Energy: " << f_nueE << std::endl;
    const TLorentzVector &v4_i = nue.Nu().Position();
    const TLorentzVector &v4_f = nue.Nu().EndPosition();

    auto x_i = v4_i.X();
    auto y_i = v4_i.Y();
    auto z_i = v4_i.Z();

    auto x_f = v4_f.X();
    auto y_f = v4_f.Y();
    auto z_f = v4_f.Z();

    f_N  =   N;
    fx_i = x_i;
    fy_i = y_i;
    fz_i = z_i;
    fx_f = x_f;
    fy_f = y_f;
    fz_f = z_f;
      
    // for (int i = 0; i < N; ++i)
    // {
      // const simb::MCParticle &nu = truth.GetParticle(i);
      // float E = nu.E();
      // float M = nu.Mass(); // Mass; from PDG unless overridden Should be in GeV
      // const int pdg = nu.PdgCode();

      // if(E < 35){f_E = E;}
      // else{std::cout << "Excluded particle with energy " << E << " GeV" << std::endl;}
      // f_M = M;
      // f_pdg = pdg;
      // if(pdg == 11){std::cout << "Generated electron" << std::endl;}
      // if(pdg == 12){std::cout << "Generated neutrino_e" << std::endl;}
      // if(pdg == 22){std::cout << "Generated photon" << std::endl;}
      

      // my_vector[4][1][counter].push_back(pdg);
      // my_vector[4][2][counter].push_back(E);
      // my_vector[4][3][counter].push_back(x_i);
      // my_vector[4][4][counter].push_back(y_i);
      // my_vector[4][5][counter].push_back(z_i);
      // my_vector[4][6][counter].push_back(x_f);
      // my_vector[4][7][counter].push_back(y_f);
      // my_vector[4][8][counter].push_back(z_f);
    // }
    fMyAnaTree->Fill();
    EvCounter = EvCounter + 1; 
  }

  // Obtain parameters from TimeService
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
  fSampleFreq = clockData.OpticalClock().Frequency();
  // float fBaselineSubtract = 0.;
  // Assume starting at 0
  // fTimeBegin  = 0;
  // Get OpDetWaveforms from the event
  auto waveformHandle = evt.getHandle< std::vector< raw::OpDetWaveform >>(fOpDetWaveformLabel);

  // Access ART's TFileService, which will handle creating and writing
  // histograms for us
  art::ServiceHandle< art::TFileService > tfs;

  for (size_t i = 0; i < waveformHandle->size(); i++)
  {
    // This was probably required to overcome the "const" problem 
    // with OpDetPulse::Waveform()
    art::Ptr< raw::OpDetWaveform > waveformPtr(waveformHandle, i);
    raw::OpDetWaveform pulse = *waveformPtr;
    wvf = pulse.Waveform();
    fOpDetWaveforms->Fill();
    // int channel = pulse.ChannelNumber();

    // int fNticks = pulse.size();

    // std::cout << "Test" << std::endl;
    // Create the TH1 if it doesn't exist
    // auto waveform = averageWaveforms.find( channel );
    // if ( waveform == averageWaveforms.end() ) {
    // TString histName = TString::Format("avgwaveform_sample_%3E", fSampleFreq);
    // averageWaveforms[i] =  tfs->make< TH1D >(histName, ";t (us);", fNticks, 0, float(fNticks+1) / fSampleFreq);
    // TH1D *hist = tfs->make< TH1D >(histName, ";t (us);", fNticks, 0, float(fNticks+1) / fSampleFreq);
    // }
    // if (!averageWaveformAll) {
    //     averageWaveformAll =  tfs->make< TH1D >("avgwaveform_channel_all", ";t (us);", fNticks, 0, float(fNticks+1) / fSampleFreq);
    // }

    // Add this waveform to this histogram
    // for (size_t tick = 0; tick < pulse.size(); tick++) {
        // averageWaveforms[i]->Fill(double(tick)/fSampleFreq, pulse[tick] - fBaselineSubtract);
      // hist->Fill(double(tick)/fSampleFreq, pulse[tick] - fBaselineSubtract);
        // averageWaveformAll->Fill(double(tick)/fSampleFreq, pulse[tick] - fBaselineSubtract);
  }

    // Count number of waveforms on each channel
    // waveformCount[channel]++;
  
  // eventCount++;
}

void ANA::OpticalAnalyzer::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(ANA::OpticalAnalyzer)
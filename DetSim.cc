////////////////////////////////////////////////////////////////////////
// Class:       MyAnalysis
// Plugin Type: analyzer (art v3_05_01)
// File:        MyAnalysis_module.cc
//
// Generated at Tue Mar 21 10:26:24 2023 by Laura Pérez-Molina
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/PtrVector.h"

// Utility libraries
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h" // simb::MCParticle
#include "nusimdata/SimulationBase/MCTrajectory.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h" //top level information form Pandora inLArSoft in recob::PFParticle
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h" // DetSim Analysis
#include "larcorealg/Geometry/GeometryCore.h" // DetSim Analysis
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"// DetSim Analysis
#include "lardata/ArtDataHelper/TrackUtils.h"
#include "lardata/ArtDataHelper/HitCreator.h" // RawDigit
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larsim/MCCheater/BackTracker.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "larsim/MCCheater/BackTrackerService.h"

// ROOT includes.
#include "TH1.h"
#include "TH1F.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TTimeStamp.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TVector3.h"

// C++ includes
#include <TTree.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <vector>
#include <cmath>
#include <map>
#include <iterator> // std::begin(), std::end()
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <algorithm>

// #include "duneana/DAQSimAna/TriggerPrimitiveFinderTool.h"

#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dunereco/AnaUtils/DUNEAnaTrackUtils.h"
#include "dunereco/AnaUtils/DUNEAnaShowerUtils.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"


namespace ana { class MyAnalysis; }


class ana::MyAnalysis : public art::EDAnalyzer 
{
  public:
    explicit MyAnalysis(fhicl::ParameterSet const& p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    MyAnalysis(MyAnalysis const &) = delete;
    MyAnalysis(MyAnalysis &&) = delete;
    MyAnalysis& operator = (MyAnalysis const&) = delete;
    MyAnalysis& operator = (MyAnalysis&&) = delete;

    // Required functions.
    void analyze(art::Event const& e) override; //bulk of action
    void reconfigure(fhicl::ParameterSet const & p);

    // Selected optional functions.
    void beginJob() override;
    void endJob()   override;

    void reset(bool deepClean=false);

  private:

    // Declare member data here.
    std::vector<std::vector<std::string>> fProductsToDump;
    const geo::Geometry* fGeom;
    bool fRollUpUnsavedIDs;

    bool debug = true;

    //==================//
    // Tree information //
    //==================//
    TTree *fTruth;
    TTree *fReco;
    TTree *fDetSim;

    // Event Information //
    unsigned int fEvent;
    unsigned int fRun;
    unsigned int fSubRun;

    static const int kNMaxMCParticles = 1000000;
    static const int kNViews          = 3;
    static const int MaxChannels = 30720;
    static const int MaxSamples  = 6000;

    // MC Particles //
    unsigned int fNMCParticles;
    bool   fMCIsPrimary                [kNMaxMCParticles];
    int    fMCParticlePdgCode          [kNMaxMCParticles];
    double fMCParticleTrueEnergy       [kNMaxMCParticles]; 
    int    fMCParticleTrackID          [kNMaxMCParticles]; 
    int    fMCParticleParentTrackID    [kNMaxMCParticles]; 
    int    fMCParticleNTrajectoryPoint [kNMaxMCParticles]; 
    //Position
    double fMCParticleStartPositionX   [kNMaxMCParticles];
    double fMCParticleStartPositionY   [kNMaxMCParticles];
    double fMCParticleStartPositionZ   [kNMaxMCParticles];
    double fMCParticleStartPositionT   [kNMaxMCParticles];
    double fMCParticleEndPositionX     [kNMaxMCParticles];
    double fMCParticleEndPositionY     [kNMaxMCParticles];
    double fMCParticleEndPositionZ     [kNMaxMCParticles];
    double fMCParticleEndPositionT     [kNMaxMCParticles];
    //Momentum
    double fMCParticleStartMomentumX   [kNMaxMCParticles];
    double fMCParticleStartMomentumY   [kNMaxMCParticles];
    double fMCParticleStartMomentumZ   [kNMaxMCParticles];
    double fMCParticleStartMomentumE   [kNMaxMCParticles];
    double fMCParticleEndMomentumX     [kNMaxMCParticles];
    double fMCParticleEndMomentumY     [kNMaxMCParticles];
    double fMCParticleEndMomentumZ     [kNMaxMCParticles];
    double fMCParticleEndMomentumE     [kNMaxMCParticles];

    double fMCParticleVertexTime       [kNMaxMCParticles];
    double fMCParticleEndTime          [kNMaxMCParticles];
    int    fMCParticleNHits            [kNMaxMCParticles];
    int    fMCParticleNHitsView        [kNMaxMCParticles][kNViews];

    //Ancestor information
    std::map<int, int> TrackIDMap;
    std::map<int, int> ParentMap;
    std::map<int, int> AncestorMap;

    int TrackIDMap_keys   [kNMaxMCParticles];
    int TrackIDMap_values [kNMaxMCParticles];
    int ParentMap_keys    [kNMaxMCParticles];
    int ParentMap_values  [kNMaxMCParticles];
    int AncestorMap_keys  [kNMaxMCParticles];    
    int AncestorMap_values[kNMaxMCParticles];

    //SimEnergyDeposited //
    int fSimEdepTrackID[kNMaxMCParticles];
    int fSimEdepPDGCode[kNMaxMCParticles];
    double fSimEdepE   [kNMaxMCParticles];
    double fSimEdepX   [kNMaxMCParticles];
    double fSimEdepY   [kNMaxMCParticles];
    double fSimEdepZ   [kNMaxMCParticles];

    //Detector Hits Info
    int channel;
    int tdc;
    int adc;
    int view;
    int track_id;
    double energy;
    std::vector<int> track_ids;
    std::vector<double> energies;
}; //End class definition


ana::MyAnalysis::MyAnalysis(fhicl::ParameterSet const & p) : EDAnalyzer{p} //, 
   // More initializers here.
{
  art::ServiceHandle<art::TFileService> tfs;

  fTruth  = tfs->make<TTree>("Truth", "Truth");
  fReco   = tfs->make<TTree>("Reco",  "Reco");
  fDetSim = tfs->make<TTree>("DetSim","DetSim");

  debug = true;
  std::vector<std::vector<std::string>> ProductsToDump = p.get<std::vector<std::vector<std::string>>>("ProductsToDump");
  fProductsToDump = ProductsToDump;

  for (auto i : fProductsToDump)
  {
    auto i_label = i[0]; auto i_instance = i[1]; auto i_type = i[2]; auto o_label = i[3];
    TString out_label = o_label;

    // Set branches: Ifs over the different types of data
    if (i_type == "simb::MCParticle")
    {
      //Event branches//
      fTruth->Branch("Event",      &fEvent,     "Event/i");
      fTruth->Branch("Run",        &fRun,       "Run/i");
      fTruth->Branch("subRun",     &fSubRun,    "subRun/i");
      //MC truth branches//
      fTruth->Branch("nMCParticles",                &fNMCParticles,               "nMCParticles/i");
      fTruth->Branch("mcIsMCPrimary",               &fMCIsPrimary,                "MCIsPrimary[nMCParticles]/O");
      fTruth->Branch("mcParticlePdgCode",           &fMCParticlePdgCode,          "MCParticlePdgCode[nMCParticles]/I");
      fTruth->Branch("mcParticleTrueEnergy",        &fMCParticleTrueEnergy,       "MCParticleTrueEnergy[nMCParticles]/D");
      fTruth->Branch("mcParticleTrackID",           &fMCParticleTrackID,          "MCParticleTrackID[nMCParticles]/I");
      fTruth->Branch("mcParticleParentTrackID",     &fMCParticleParentTrackID,    "MCParticleParentTrackID[nMCParticles]/I");

      fTruth->Branch("mcParticleStartPositionX",    &fMCParticleStartPositionX,    "MCParticleStartPositionX[nMCParticles]/D");
      fTruth->Branch("mcParticleStartPositionY",    &fMCParticleStartPositionY,    "MCParticleStartPositionY[nMCParticles]/D");
      fTruth->Branch("mcParticleStartPositionZ",    &fMCParticleStartPositionZ,    "MCParticleStartPositionZ[nMCParticles]/D");
      fTruth->Branch("mcParticleStartPositionT",    &fMCParticleStartPositionT,    "MCParticleStartPositionT[nMCParticles]/D");
      fTruth->Branch("mcParticleStartMomentumX",    &fMCParticleStartMomentumX,    "MCParticleStartMomentumX[nMCParticles]/D");
      fTruth->Branch("mcParticleStartMomentumY",    &fMCParticleStartMomentumY,    "MCParticleStartMomentumY[nMCParticles]/D");
      fTruth->Branch("mcParticleStartMomentumZ",    &fMCParticleStartMomentumZ,    "MCParticleStartMomentumZ[nMCParticles]/D");
      fTruth->Branch("mcParticleStartMomentumE",    &fMCParticleStartMomentumE,    "MCParticleStartMomentumE[nMCParticles]/D");
      fTruth->Branch("mcParticleEndPositionX",      &fMCParticleEndPositionX,      "MCParticleEndPositionX[nMCParticles]/D");
      fTruth->Branch("mcParticleEndPositionY",      &fMCParticleEndPositionY,      "MCParticleEndPositionY[nMCParticles]/D");
      fTruth->Branch("mcParticleEndPositionZ",      &fMCParticleEndPositionZ,      "MCParticleEndPositionZ[nMCParticles]/D");
      fTruth->Branch("mcParticleEndPositionT",      &fMCParticleEndPositionT,      "MCParticleEndPositionT[nMCParticles]/D");
      fTruth->Branch("mcParticleEndMomentumX",      &fMCParticleEndMomentumX,      "MCParticleEndMomentumX[nMCParticles]/D");
      fTruth->Branch("mcParticleEndMomentumY",      &fMCParticleEndMomentumY,      "MCParticleEndMomentumY[nMCParticles]/D");
      fTruth->Branch("mcParticleEndMomentumZ",      &fMCParticleEndMomentumZ,      "MCParticleEndMomentumZ[nMCParticles]/D");
      fTruth->Branch("mcParticleEndMomentumE",      &fMCParticleEndMomentumE,      "MCParticleEndMomentumE[nMCParticles]/D");
      
      fTruth->Branch("mcParticleVertexTime",        &fMCParticleVertexTime,        "MCParticleVertexTime[nMCParticles]/D");
      fTruth->Branch("mcParticleEndTime",           &fMCParticleEndTime,           "MCParticleEndTime[nMCParticles]/D");
      
      // ANCESTOR INFORMATION //
      fTruth->Branch("TrackIDMap_values",  &TrackIDMap_values,  "TrackIDMap_values[nMCParticles]/I");
      fTruth->Branch("TrackIDMap_keys",    &TrackIDMap_keys,    "TrackIDMap_keys[nMCParticles]/I");
      fTruth->Branch("ParentMap_values",   &ParentMap_values ,  "ParentMap_values[nMCParticles]/I");
      fTruth->Branch("ParentMap_keys",     &ParentMap_keys,     "ParentMap_keys[nMCParticles]/I");   
      fTruth->Branch("AncestorMap_values", &AncestorMap_values, "AncestorMap_values[nMCParticles]/I");
      fTruth->Branch("AncestorMap_keys",   &AncestorMap_keys,   "AncestorMap_keys[nMCParticles]/I");
    }
    else if (i_type == "sim::SimEnergyDeposit")
    {
      // ENERGY DEPOSITED //
      fTruth->Branch("SimEdepTrackID",         &fSimEdepTrackID    ,"SimEdepTrackID[nMCParticles]/I");
      fTruth->Branch("SimEdepPDGCode",         &fSimEdepPDGCode    ,"SimEdepPDGCode[nMCParticles]/I");
      fTruth->Branch("SimEdepEnergy",          &fSimEdepE          ,"fSimEdepE[nMCParticles]/D");
      fTruth->Branch("SimEdepMiddlePositionX", &fSimEdepX ,"SimEdepX[nMCParticles]/D");
      fTruth->Branch("SimEdepMiddlePositionY", &fSimEdepY ,"SimEdepY[nMCParticles]/D");
      fTruth->Branch("SimEdepMiddlePositionZ", &fSimEdepZ ,"SimEdepZ[nMCParticles]/D");
    }
  } //for loop over products to dump

    fDetSim->Branch("Channels", &channel);
    fDetSim->Branch("TDCs", &tdc);
    fDetSim->Branch("ADCs", &adc);
    fDetSim->Branch("View", &view);
    fDetSim->Branch("TrackID", &track_id);
    fDetSim->Branch("TrackIDs", &track_ids);
    fDetSim->Branch("Energy", &energy);
    fDetSim->Branch("Energies", &energies);

  fRollUpUnsavedIDs    = p.get<bool>("RollUpUnsavedIDs"); 
  fGeom                = &*art::ServiceHandle<geo::Geometry>();
  // this->reconfigure(p);
} 


void ana::MyAnalysis::beginJob()
{
  // Implementation of optional member function here.
  // reset(true); //deep clean the variables

  // Call appropriate consumes<>() for any products to be retrieved by this module.
}//beginJob



void ana::MyAnalysis::analyze(const art::Event & evt)
{
  // reset(); //Don't deep clean

  struct detector_output
        {
            int channel   = {  };
            int tdc       = {  };
            int adc       = {  };
            int view      = {  };
            int track_id  = {-1};
            double energy = { 0};

            std::vector<int> track_ids   = {};
            std::vector<double> energies = {};

            detector_output(int _channel, int _tdc, int _adc, int _view, std::vector<sim::IDE> ides)
            { 
              channel = _channel;
              tdc     = _tdc;
              adc     = _adc;
              view    = _view;
              for (auto ide : ides)
              {
                track_ids.emplace_back(ide.trackID);
                energies.emplace_back(ide.energy);
                if (ide.energy > energy) 
                {
                  energy   = ide.energy;
                  track_id = ide.trackID;
                }
              }
            }
        };

  fEvent  = evt.id().event();
  fRun    = evt.id().run();
  fSubRun = evt.id().subRun();
  std::map<int, int> fTrackIDMap; // store out of the loop for the truth_tree a map of trackID to MCParticle index
  std::vector<detector_output> output_array; 


  std::cout << "=============== EVENT ID " << fEvent << " == RUN ID " << fRun << " == SUBRUN ID " << fSubRun << " ================" << std::endl;

  for (auto i : fProductsToDump)
  {
    auto i_label = i[0]; auto i_instance = i[1]; auto i_type = i[2]; auto o_label = i[3];
    if (debug){ std::cout << "---> Module_Label: " << i_label << ";\t Instance_Name: " << i_instance << ";\t Product_Type: " << i_type << ";\t Output_Name: " << o_label << std::endl; }
    
    if (std::find(trees.begin(), trees.end(), std::string("Truth")) != trees.end()) 
    {
      if (i_type == "simb::MCParticle")
      {
        if(!evt.isRealData())
        {
          art::ValidHandle<std::vector<simb::MCParticle>> mcParticles = evt.getValidHandle<std::vector<simb::MCParticle>>(i_label);

          if(mcParticles.isValid())
          {
            fNMCParticles = mcParticles->size();
            bool isMCPrimary(false);  
            for(unsigned int iMc = 0; iMc < fNMCParticles; iMc++)
            {
              const simb::MCParticle trueParticle = mcParticles->at(iMc);
              trueParticle.Process() == "primary"? isMCPrimary = true:isMCPrimary = false;
              fMCParticleTrueEnergy      [iMc] = trueParticle.E();
              fMCParticlePdgCode         [iMc] = trueParticle.PdgCode();
              fMCParticleTrackID         [iMc] = trueParticle.TrackId();
              fMCParticleVertexTime      [iMc] = trueParticle.T();
              fMCParticleEndTime         [iMc] = trueParticle.EndT();
              fMCParticleParentTrackID   [iMc] = trueParticle.Mother();
              fMCParticleNTrajectoryPoint[iMc] = trueParticle.NumberTrajectoryPoints();
              fMCIsPrimary               [iMc]    = isMCPrimary;
              // Position //
              fMCParticleStartPositionX  [iMc] = trueParticle.Position().X();
              fMCParticleStartPositionY  [iMc] = trueParticle.Position().Y();
              fMCParticleStartPositionZ  [iMc] = trueParticle.Position().Z();
              fMCParticleStartPositionT  [iMc] = trueParticle.Position().T();
              fMCParticleEndPositionX    [iMc] = trueParticle.EndPosition().X();
              fMCParticleEndPositionY    [iMc] = trueParticle.EndPosition().Y();
              fMCParticleEndPositionZ    [iMc] = trueParticle.EndPosition().Z();
              fMCParticleEndPositionT    [iMc] = trueParticle.EndPosition().T();
              // Momentum //
              fMCParticleStartMomentumX  [iMc] = trueParticle.Momentum().X();
              fMCParticleStartMomentumY  [iMc] = trueParticle.Momentum().Y();
              fMCParticleStartMomentumZ  [iMc] = trueParticle.Momentum().Z();
              fMCParticleStartMomentumE  [iMc] = trueParticle.Momentum().E();
              fMCParticleEndMomentumX    [iMc] = trueParticle.EndMomentum().X();
              fMCParticleEndMomentumY    [iMc] = trueParticle.EndMomentum().Y();
              fMCParticleEndMomentumZ    [iMc] = trueParticle.EndMomentum().Z();
              fMCParticleEndMomentumE    [iMc] = trueParticle.EndMomentum().E();
              
              ParentMap[trueParticle.TrackId()] = trueParticle.Mother();

              int mother = trueParticle.Mother(); int ancestry_level = 0; 
              if (mother != 0) 
              {
                while(ancestry_level < 20) 
                {
                  ancestry_level++ ;
                  int track_id = mother;
                  mother = ParentMap[track_id];
                  if (mother == 0) break;
                }
              }
              
              AncestorMap[trueParticle.TrackId()] = ancestry_level;
              TrackIDMap [trueParticle.TrackId()] = iMc;
              // std::cout<< "LEVEL " << AncestorMap[trueParticle.TrackId()]  << std::endl;

              TrackIDMap_keys   [iMc]  = trueParticle.TrackId();
              TrackIDMap_values [iMc]  = iMc;
              AncestorMap_keys  [iMc]  = trueParticle.TrackId();
              AncestorMap_values[iMc]  = ancestry_level;
              ParentMap_keys    [iMc]  = trueParticle.TrackId();
              ParentMap_values  [iMc]  = trueParticle.Mother();

            } //for loop filling MCInfo  

            fTrackIDMap = TrackIDMap; // to be used in other branches/trees

          } //if mcParticles.isValid()
        } //if !evt.isRealData()
      } //if i_type == "simb::MCTruth"
    } //if Truth is in the vector

    if (std::find(trees.begin(), trees.end(), std::string("DetSim")) != trees.end()) 
    {
        //================================================== Tree DetSim ==================================================//
        std::map<std::tuple<int, int, int, int>, sim::SimChannel> channels; //map with keys as a tuple of three integers and values of type sim::SimChannel
        art::Handle<std::vector<sim::SimChannel>> SimChaHandle;
        evt.getByLabel("tpcrawdecoder", "simpleSC", SimChaHandle);
        // evt.getByLabel(fDetSimLabel, fSimChaInstanceLabel, SimChaHandle);
        if(!SimChaHandle.isValid())
        {
          // std::cout<<"Unable to find std::vector<sim::SimChannel> with module label: " << fDetSimLabel << "; Instance: " << fSimChaInstanceLabel << std::endl;
          std::cout<<"Unable to find std::vector<sim::SimChannel> with module label: " << "tpcrawdecoder" << "; Instance: " << "simpleSC" << std::endl;
          return;
        } // !DetSimHandle
        
        art::Handle<std::vector<raw::RawDigit>> RawDigHandle;
        evt.getByLabel("tpcrawdecoder", "daq", RawDigHandle);
        // evt.getByLabel(fDetSimLabel, fRawDigInstanceLabel, RawDigHandle);
        if(!RawDigHandle.isValid())
        {
          // std::cout<<"Unable to find std::vector<raw::RawDigit> with module label: " << fDetSimLabel << "; Instance: " << fRawDigInstanceLabel << std::endl;
          std::cout<<"Unable to find std::vector<raw::RawDigit> with module label: " << "tpcrawdecoder" << "; Instance: " << "simpleSC" << std::endl;
          return;
        } // !RawDigHandle

        if(SimChaHandle.isValid() && RawDigHandle.isValid())
        {
          // tdc.clear(); adc.clear(); view.clear();
          channels.clear();

          std::vector<sim::SimChannel> sim_channels; 
          std::vector<raw::RawDigit> digits;

          for(auto sim_channel : *SimChaHandle) { sim_channels.push_back(sim_channel); }
          for(auto digit : *RawDigHandle) { digits.push_back(digit); }  // loop over all raw digits (i.e ADC counts for each channel for each time tick)

          int DigsSize = digits.size();                // number of digits [DUNE-FD: 30720]
          std::cout << "Digit has " << DigsSize << " size" << std::endl;
          std::cout << "Digit has " << digits[0].Samples() << " samples" << std::endl;
          std::cout << "Digit has " << (int)digits[0].GetPedestal() << " pedestal" << std::endl;

          for (int dig=0; dig<DigsSize; dig++)
          {
            // Processing each channel and tick.
            int channel  = digits[dig].Channel(); // channel number
            int NSamples = digits[dig].Samples();          // number of ADC samples (TDC ticks) [DUNE-FD: 6000]
            int pedestal = (int)digits[dig].GetPedestal(); // pedestal value [DUNE-FD: 2350]

            // uncompress the digits and remove the pedestal
            std::vector<short> uncompressed(NSamples); // uncompressed ADC values
            raw::Uncompress( digits[dig].ADCs(), uncompressed, pedestal, digits[dig].Compression()); 

            for (int tick=0; tick<NSamples; tick++)
            {
              geo::GeometryCore const* geometry_core = lar::providerFrom<geo::Geometry>(); // geometry provider --> which wire plane
              if(std::abs(uncompressed[tick]) > 20.0) // ADC threshold > 20
              {
                auto sim_channel = sim_channels[channel]; // sim channel
                
                std::vector<sim::IDE> ides = sim_channel.TrackIDsAndEnergies(tick, tick);
                for (auto const& ide : ides) 
                {
                  int particle_id = fTrackIDMap[ide.trackID];
                  view = geometry_core->View(channel);
                  std::tuple<int, int, int, int> key(channel, tick, particle_id, view);
                  channels[key] = sim_channel;
                  // Print the values being added to the channels map
                  // std::cout << "Added to channels map: Channel=" << channel
                  // << " Tick=" << tick << " ParticleID=" << particle_id << std::endl;
                }
              }
            }
          } // for (int dig=0; dig<DigsSize; dig++)

          if (channels.empty())  { std::cout << "WARNING channels empty"   << std::endl; }
          if (!channels.empty()) { std::cout << "GREAT channels not empty" << std::endl; }

          for (const auto& entry : channels) //iterates over each entry in the channels map
          {
            int channel, tick, particle_id, view;
            std::tie(channel, tick, particle_id, view) = entry.first;  //unpacks the tuple into the three variables
            const sim::SimChannel& sim_channel         = entry.second; //unpacks the sim channel into the variable

            output_array.emplace_back(channel, tick, particle_id, view, sim_channel.TrackIDsAndEnergies(tick, tick)); //creates a new entry in the output array
          }
        } // if(SimChaHandle.isValid() && RawDigHandle.isValid())
    } // if (tree == "DetSim")

  } // fProductsToDump loop


  fTruth->Fill();

  for (const auto& output : output_array) 
  {
    // Access the components of each detector_output struct
    channel   = output.channel;
    tdc       = output.tdc;
    adc       = output.adc;
    view      = output.view;
    track_id  = output.track_id;
    energy    = output.energy;
    track_ids = output.track_ids;
    energies  = output.energies;
    
    fDetSim->Fill();
    
    // Print or use the components as needed
    // std::cout << "Channel: " << channel << std::endl;
    // std::cout << "TDC: " << tdc << std::endl;
    // std::cout << "ADC: " << adc << std::endl;
    // std::cout << "Track ID: " << track_id << std::endl;
    // std::cout << "Energy: " << energy << std::endl;
    // // Print the track IDs
    // std::cout << "Track IDs: ";
    // for (int track : track_ids) { std::cout << track << " "; }
    // std::cout << std::endl;
    // // Print the energies
    // std::cout << "Energies: ";
    // for (double e : energies) { std::cout << e << " "; }
    // std::cout << std::endl;
  } // for (const auto& output : output_array)

}// analyze()



void ana::MyAnalysis::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(ana::MyAnalysis)
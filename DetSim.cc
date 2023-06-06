///////////////////////////////////////////////////////////////////////
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
    std::string fTruthLabel;
    std::string fDetSimLabel;
    std::string fRawDigInstanceLabel;
    std::string fSimChaInstanceLabel;
    bool fRollUpUnsavedIDs;

    //==================//
    // Tree information //
    //==================//
    TTree *fTree;
    TTree *fDetSim;

    //-------------------//
    // Event Information //
    //-------------------//
    unsigned int fEvent;
    unsigned int fRun;
    unsigned int fSubRun;

    static const int kNMaxMCParticles = 1000000;
    static const int kNMaxPFParticles = 2000;
    static const int kNMaxPFPClusters = 100;
    static const int kNViews          = 3;

    //--------------//
    // MC Particles //
    //--------------//
    unsigned int fNMCParticles;
    bool   fMCIsPrimary                [kNMaxMCParticles];
    int    fMCParticlePdgCode          [kNMaxMCParticles];
    double fMCParticleTrueEnergy       [kNMaxMCParticles]; 
    int    fMCParticleTrackID          [kNMaxMCParticles]; 
    int    fMCParticleParentTrackID    [kNMaxMCParticles]; 
    int    fMCParticleNTrajectoryPoint [kNMaxMCParticles]; 
    //////////
    //Position
    double fMCParticleStartPositionX   [kNMaxMCParticles];
    double fMCParticleStartPositionY   [kNMaxMCParticles];
    double fMCParticleStartPositionZ   [kNMaxMCParticles];
    double fMCParticleStartPositionT   [kNMaxMCParticles];
    double fMCParticleEndPositionX     [kNMaxMCParticles];
    double fMCParticleEndPositionY     [kNMaxMCParticles];
    double fMCParticleEndPositionZ     [kNMaxMCParticles];
    double fMCParticleEndPositionT     [kNMaxMCParticles];
    //////////
    //Momentum
    double fMCParticleStartMomentumX   [kNMaxMCParticles];
    double fMCParticleStartMomentumY   [kNMaxMCParticles];
    double fMCParticleStartMomentumZ   [kNMaxMCParticles];
    double fMCParticleStartMomentumE   [kNMaxMCParticles];
    double fMCParticleEndMomentumX     [kNMaxMCParticles];
    double fMCParticleEndMomentumY     [kNMaxMCParticles];
    double fMCParticleEndMomentumZ     [kNMaxMCParticles];
    double fMCParticleEndMomentumE     [kNMaxMCParticles];

    /////////////////////
    //Ancestor information
    std::map<int, int> TrackIDMap;
    std::map<int, int> ParentMap;
    std::map<int, int> AncestorMap;

    int TrackIDMap_keys [kNMaxMCParticles];
    int TrackIDMap_values [kNMaxMCParticles];
    int ParentMap_keys  [kNMaxMCParticles];
    int ParentMap_values  [kNMaxMCParticles];
    int AncestorMap_keys[kNMaxMCParticles];    
    int AncestorMap_values[kNMaxMCParticles];

    std::vector<std::vector<int>> channels;
    std::vector<std::vector<int>> tdc;
    std::vector<std::vector<int>> adc;
    std::vector<std::vector<int>> view;

}; //End class definition


ana::MyAnalysis::MyAnalysis(fhicl::ParameterSet const & p)
  : EDAnalyzer{p} //, 
   // More initializers here.
{
  this->reconfigure(p);
} 


void ana::MyAnalysis::reconfigure(fhicl::ParameterSet const& p)
{
  // assign fTruthLabel in the constructor using the parameter defined in MyAnalysis.fcl
  fTruthLabel          = p.get<std::string>("TruthLabel");
  fEdepInstanceLabel   = p.get<std::string>("InstanceName");
  fDetSimLabel         = p.get<std::string>("DetSimLabel");
  fRawDigInstanceLabel = p.get<std::string>("RawDigInstance");
  fSimChaInstanceLabel = p.get<std::string>("SimChaInstance");
  fRollUpUnsavedIDs    = p.get<bool>("RollUpUnsavedIDs"); 
  fGeom                = &*art::ServiceHandle<geo::Geometry>();
} // Reconfigure


void ana::MyAnalysis::beginJob()
{
  // Implementation of optional member function here.
  // reset(true); //deep clean the variables

  // Call appropriate consumes<>() for any products to be retrieved by this module.
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("Pandora_Output","Pandora_Output");
  fDetSim = tfs->make<TTree>("DetSim","DetSim");
  //////////////////
  //Event branches//
  //////////////////
  fTree->Branch("Event",      &fEvent,     "Event/i");
  fTree->Branch("Run",        &fRun,       "Run/i");
  fTree->Branch("subRun",     &fSubRun,    "subRun/i");
  /////////////////////
  //MC truth branches//
  /////////////////////
  fTree->Branch("nMCParticles",                &fNMCParticles,               "nMCParticles/i");
  fTree->Branch("mcIsMCPrimary",               &fMCIsPrimary,                "MCIsPrimary[nMCParticles]/O");
  fTree->Branch("mcParticlePdgCode",           &fMCParticlePdgCode,          "MCParticlePdgCode[nMCParticles]/I");
  fTree->Branch("mcParticleTrueEnergy",        &fMCParticleTrueEnergy,       "MCParticleTrueEnergy[nMCParticles]/D");
  fTree->Branch("mcParticleTrackID",           &fMCParticleTrackID,          "MCParticleTrackID[nMCParticles]/I");
  fTree->Branch("mcParticleParentTrackID",     &fMCParticleParentTrackID,    "MCParticleParentTrackID[nMCParticles]/I");
  // fTree->Branch("mcParticleNTrajectoryPoints", &fMCParticleNTrajectoryPoint, "MCParticleNTrajectoryPoint[nMCParticles]/I");

  fTree->Branch("mcParticleStartPositionX",    &fMCParticleStartPositionX,    "MCParticleStartPositionX[nMCParticles]/D");
  fTree->Branch("mcParticleStartPositionY",    &fMCParticleStartPositionY,    "MCParticleStartPositionY[nMCParticles]/D");
  fTree->Branch("mcParticleStartPositionZ",    &fMCParticleStartPositionZ,    "MCParticleStartPositionZ[nMCParticles]/D");
  fTree->Branch("mcParticleStartPositionT",    &fMCParticleStartPositionT,    "MCParticleStartPositionT[nMCParticles]/D");
  fTree->Branch("mcParticleStartMomentumX",    &fMCParticleStartMomentumX,    "MCParticleStartMomentumX[nMCParticles]/D");
  fTree->Branch("mcParticleStartMomentumY",    &fMCParticleStartMomentumY,    "MCParticleStartMomentumY[nMCParticles]/D");
  fTree->Branch("mcParticleStartMomentumZ",    &fMCParticleStartMomentumZ,    "MCParticleStartMomentumZ[nMCParticles]/D");
  fTree->Branch("mcParticleStartMomentumE",    &fMCParticleStartMomentumE,    "MCParticleStartMomentumE[nMCParticles]/D");
  fTree->Branch("mcParticleEndPositionX",      &fMCParticleEndPositionX,      "MCParticleEndPositionX[nMCParticles]/D");
  fTree->Branch("mcParticleEndPositionY",      &fMCParticleEndPositionY,      "MCParticleEndPositionY[nMCParticles]/D");
  fTree->Branch("mcParticleEndPositionZ",      &fMCParticleEndPositionZ,      "MCParticleEndPositionZ[nMCParticles]/D");
  fTree->Branch("mcParticleEndPositionT",      &fMCParticleEndPositionT,      "MCParticleEndPositionT[nMCParticles]/D");
  fTree->Branch("mcParticleEndMomentumX",      &fMCParticleEndMomentumX,      "MCParticleEndMomentumX[nMCParticles]/D");
  fTree->Branch("mcParticleEndMomentumY",      &fMCParticleEndMomentumY,      "MCParticleEndMomentumY[nMCParticles]/D");
  fTree->Branch("mcParticleEndMomentumZ",      &fMCParticleEndMomentumZ,      "MCParticleEndMomentumZ[nMCParticles]/D");
  fTree->Branch("mcParticleEndMomentumE",      &fMCParticleEndMomentumE,      "MCParticleEndMomentumE[nMCParticles]/D");
  
  fTree->Branch("TrackIDMap_values",  &TrackIDMap_values,  "TrackIDMap_values[nMCParticles]/I");
  fTree->Branch("TrackIDMap_keys",    &TrackIDMap_keys,    "TrackIDMap_keys[nMCParticles]/I");
  fTree->Branch("ParentMap_values",   &ParentMap_values ,  "ParentMap_values[nMCParticles]/I");
  fTree->Branch("ParentMap_keys",     &ParentMap_keys,     "ParentMap_keys[nMCParticles]/I");   
  fTree->Branch("AncestorMap_values", &AncestorMap_values, "AncestorMap_values[nMCParticles]/I");
  fTree->Branch("AncestorMap_keys",   &AncestorMap_keys,   "AncestorMap_keys[nMCParticles]/I");

  fDetSim->Branch("Channels", &channels);
  fDetSim->Branch("TDC",      &tdc     );    
  fDetSim->Branch("ADC",      &adc     );    
  fDetSim->Branch("View",     &view    );
}//beginJob



void ana::MyAnalysis::analyze(const art::Event & evt)
{

  // reset(); //Don't deep clean

  // const art::ServiceHandle<cheat::BackTrackerService> btServ;
  fEvent  = evt.id().event();
  fRun    = evt.id().run();
  fSubRun = evt.id().subRun();
  // auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt); //Contains all timing reference information for the detector. 

  std::cout << "=============== EVENT ID " << fEvent << " == RUN ID " << fRun << " == SUBRUN ID " << fSubRun << " ================" << std::endl;

  //Access the truth information//
  if(!evt.isRealData())
  {
    art::ValidHandle<std::vector<simb::MCParticle>> mcParticles = evt.getValidHandle<std::vector<simb::MCParticle>>(fTruthLabel);

    if(mcParticles.isValid())
    {
      fNMCParticles = mcParticles->size();
      bool isMCPrimary(false);  
      for(unsigned int iMc = 0; iMc < fNMCParticles; iMc++)
      {
        const simb::MCParticle trueParticle = mcParticles->at(iMc);
        fMCParticleTrueEnergy      [iMc] = trueParticle.E();
        fMCParticlePdgCode         [iMc] = trueParticle.PdgCode();
        fMCParticleTrackID         [iMc] = trueParticle.TrackId();
        fMCParticleVertexTime      [iMc] = trueParticle.T();
        fMCParticleEndTime         [iMc] = trueParticle.EndT();
        fMCParticleParentTrackID   [iMc] = trueParticle.Mother();
        fMCParticleNTrajectoryPoint[iMc] = trueParticle.NumberTrajectoryPoints();
        
        trueParticle.Process() == "primary"? isMCPrimary = true:isMCPrimary = false;

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

        TrackIDMap_keys[iMc]    = trueParticle.TrackId();
        TrackIDMap_values[iMc]  = iMc;
        AncestorMap_keys[iMc]   = trueParticle.TrackId();
        AncestorMap_values[iMc] = ancestry_level;
        ParentMap_keys[iMc]     = trueParticle.TrackId();
        ParentMap_values[iMc]   = trueParticle.Mother();

      } //for loop filling MCInfo  
    } //if mcParticles.isValid()
  } //if !evt.isRealData()


if (fDetSimLabel != "")
  {

    channels.clear(); tdc.clear(); adc.clear(); view.clear();
    art::Handle<std::vector<sim::SimChannel>> SimChaHandle;
    evt.getByLabel(fDetSimLabel, fSimChaInstanceLabel, SimChaHandle);
    if(!SimChaHandle.isValid())
    {
      std::cout<<"Unable to find std::vector<sim::SimChannel> with module label: " << fDetSimLabel << "; Instance: " << fSimChaInstanceLabel << std::endl;
      return;
    } // !DetSimHandle
    
    art::Handle<std::vector<raw::RawDigit>> RawDigHandle;
    evt.getByLabel(fDetSimLabel, fRawDigInstanceLabel, RawDigHandle);
    if(!RawDigHandle.isValid())
    {
      std::cout<<"Unable to find std::vector<raw::RawDigit> with module label: " << fDetSimLabel << "; Instance: " << fRawDigInstanceLabel << std::endl;
      return;
    } // !RawDigHandle

    if(SimChaHandle.isValid() && RawDigHandle.isValid())
    {
      // channels.clear(); tdc.clear(); adc.clear(); view.clear();
      std::vector<sim::SimChannel> sim_channels; 
      for(auto sim_channel : *SimChaHandle) { sim_channels.push_back(sim_channel); }

      for(auto digit : *RawDigHandle) // loop over all raw digits (i.e ADC counts for each channel for each time tick)
      // There are 15360 channels in the ProtoDUNE-SP detector.
      {
        int num_samples = digit.Samples();          // number of ADC samples (TDC ticks)
        int pedestal    = (int)digit.GetPedestal(); // pedestal value
        
        // uncompress the digits and remove the pedestal
        std::vector<short> uncompressed(num_samples); // uncompressed ADC values
        raw::Uncompress( digit.ADCs(), uncompressed, pedestal, digit.Compression()); 
        for (int ii = 0; ii < num_samples; ii++) { uncompressed[ii] -= pedestal; } // remove pedestal

        std::vector<int> channels_idparticle;
        for(int tdc_tick = 0; tdc_tick < num_samples; tdc_tick++)
        {
          auto channel     = digit.Channel();        // channel number
          auto sim_channel = sim_channels[channel];  // sim channel

          channels_idparticle.push_back(sim_channel.Channel()); // channel number
        }
        if (channels_idparticle.empty()) { std::cout << "GREAT" << std::endl; channels.push_back(std::vector<int>()); } 
        else { channels.push_back(channels_idparticle); }
      }
    }
  }

  if (!channels.empty()) 
  {
    std::cout << channels.size() << std::endl;
  }

  fDetSim->Fill();
  fTree->Fill();

}// analyze()



void ana::MyAnalysis::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(ana::MyAnalysis)
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
    static const int kNMaxPFParticles = 2000;
    static const int kNMaxPFPClusters = 100;
    static const int kNViews          = 3;
    static const int kevents = 100;

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


    // "Particle Flow" Particle (PFP) information 
    unsigned int fNPFParticles;
    int    fPFPID                             [kNMaxPFParticles];
    bool   fPFPIsPrimary                      [kNMaxPFParticles];
    int    fPFPTrueParticleMatchedID          [kNMaxPFParticles];
    int    fPFPTrueParticleMatchedPosition    [kNMaxPFParticles];
    int    fPFPParentID                       [kNMaxPFParticles];
    int    fPFPPdgCode                        [kNMaxPFParticles];
    int    fPFPNChildren                      [kNMaxPFParticles];
    int    fPFPNHits                          [kNMaxPFParticles];
    int    fPFPNHitsView                      [kNMaxPFParticles][kNViews];
    int    fPFPNSharedTrueParticleHits        [kNMaxPFParticles];
    int    fPFPNSharedTrueParticleHitsView    [kNMaxPFParticles][kNViews];
    int    fPFPTrueParticleMatchedIDView      [kNMaxPFParticles][kNViews];
    int    fPFPTrueParticleMatchedPositionView[kNMaxPFParticles][kNViews];
    double fPFPCompleteness                   [kNMaxMCParticles];
    double fPFPCompletenessView               [kNMaxMCParticles][kNViews];
    double fPFPPurity                         [kNMaxMCParticles];
    double fPFPPurityView                     [kNMaxMCParticles][kNViews];
    //////////
    //Clusters
    int    fPFPNClusters                      [kNMaxPFParticles];
    int    fPFPCluPlane                       [kNMaxPFParticles][kNMaxPFPClusters];
    int    fPFPCluView                        [kNMaxPFParticles][kNMaxPFPClusters];
    int    fPFPCluNHits                       [kNMaxPFParticles][kNMaxPFPClusters];
    double fPFPCluIntegral                    [kNMaxPFParticles][kNMaxPFPClusters];
    //////////
    //Tracks
    bool   fPFPIsTrack                        [kNMaxPFParticles];
    int    fPFPTrackID                        [kNMaxPFParticles];
    double fPFPTrackLength                    [kNMaxPFParticles];
    double fPFPTrackStartX                    [kNMaxPFParticles];
    double fPFPTrackStartY                    [kNMaxPFParticles];
    double fPFPTrackStartZ                    [kNMaxPFParticles];
    double fPFPTrackVertexX                   [kNMaxPFParticles];
    double fPFPTrackVertexY                   [kNMaxPFParticles];
    double fPFPTrackVertexZ                   [kNMaxPFParticles];
    double fPFPTrackEndX                      [kNMaxPFParticles];
    double fPFPTrackEndY                      [kNMaxPFParticles];
    double fPFPTrackEndZ                      [kNMaxPFParticles];
    double fPFPTrackTheta                     [kNMaxPFParticles];
    double fPFPTrackPhi                       [kNMaxPFParticles];
    double fPFPTrackZenithAngle               [kNMaxPFParticles];
    double fPFPTrackAzimuthAngle              [kNMaxPFParticles];
    double fPFPTrackStartDirectionX           [kNMaxPFParticles];
    double fPFPTrackStartDirectionY           [kNMaxPFParticles];
    double fPFPTrackStartDirectionZ           [kNMaxPFParticles];
    double fPFPTrackVertexDirectionX          [kNMaxPFParticles];
    double fPFPTrackVertexDirectionY          [kNMaxPFParticles];
    double fPFPTrackVertexDirectionZ          [kNMaxPFParticles];
    double fPFPTrackEndDirectionX             [kNMaxPFParticles];
    double fPFPTrackEndDirectionY             [kNMaxPFParticles];
    double fPFPTrackEndDirectionZ             [kNMaxPFParticles];
    float  fPFPTrackChi2                      [kNMaxPFParticles];
    int    fPFPTrackNdof                      [kNMaxPFParticles];
    //////////
    //Showers
    bool   fPFPIsShower                       [kNMaxPFParticles];
    int    fPFPShowerID                       [kNMaxPFParticles];
    int    fPFPShowerBestPlane                [kNMaxPFParticles];
    double fPFPShowerDirectionX               [kNMaxPFParticles];
    double fPFPShowerDirectionY               [kNMaxPFParticles];
    double fPFPShowerDirectionZ               [kNMaxPFParticles];
    double fPFPShowerDirectionErrX            [kNMaxPFParticles];
    double fPFPShowerDirectionErrY            [kNMaxPFParticles];
    double fPFPShowerDirectionErrZ            [kNMaxPFParticles];
    double fPFPShowerStartX                   [kNMaxPFParticles];
    double fPFPShowerStartY                   [kNMaxPFParticles];
    double fPFPShowerStartZ                   [kNMaxPFParticles];
    double fPFPShowerStartErrX                [kNMaxPFParticles];
    double fPFPShowerStartErrY                [kNMaxPFParticles];
    double fPFPShowerStartErrZ                [kNMaxPFParticles];
    double fPFPShowerLength                   [kNMaxPFParticles];
    double fPFPShowerOpenAngle                [kNMaxPFParticles];



    //Detector Hits Info
    // std::vector<Int_t> channels;
    // std::vector<Int_t> tdc;
    // std::vector<Int_t> adc;
    // std::vector<Int_t> view;
    
    std::vector<std::vector<int>> channels;
    std::vector<std::vector<int>> tdc;
    std::vector<std::vector<int>> adc;
    std::vector<std::vector<int>> view;

    // int channel[kNMaxMCParticles];
    // int tdc    [kNMaxMCParticles];
    // int adc    [kNMaxMCParticles];
    // int view   [kNMaxMCParticles];

}; //End class definition


ana::MyAnalysis::MyAnalysis(fhicl::ParameterSet const & p) : EDAnalyzer{p} //, 
   // More initializers here.
{
  art::ServiceHandle<art::TFileService> tfs;
  fTruth  = tfs->make<TTree>("Truth", "Truth");
  fReco   = tfs->make<TTree>("Reco",  "Reco");
  // fDetSim = tfs->make<TTree>("DetSim","DetSim");


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
      // fTruth->Branch("mcParticleNTrajectoryPoints", &fMCParticleNTrajectoryPoint, "MCParticleNTrajectoryPoint[nMCParticles]/I");

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
      // fTruth->Branch("mcParticleNHits",             &fMCParticleNHits,             "MCParticleNHits[nMCParticles]/I");
      // fTruth->Branch("mcParticleNHitsView",         &fMCParticleNHitsView,         "MCParticleNHitsView[nMCParticles][3]/I");
      
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
    // else if (i_type == "raw::RawDigit")
    // {
    //   // DETECTOR HITS //
    //   fDetSim->Branch("Channels", &channels);
    //   fDetSim->Branch("TDC",      &tdc     );    
    //   fDetSim->Branch("ADC",      &adc     );    
    //   fDetSim->Branch("View",     &view    );
    // // fTree->Branch("Channels",&channels,"Channels[nMCParticles]/I");
    // // fTree->Branch("TDC",     &tdc,     "TDC[nMCParticles]/I");
    // // fTree->Branch("ADC",     &adc,     "ADC[nMCParticles]/I");
    // // fTree->Branch("View",    &view,    "View[nMCParticles]/I");
    // } 

  } //for loop over products to dump

    //Event branches//
    fReco->Branch("Event",      &fEvent,     "Event/i");
    fReco->Branch("Run",        &fRun,       "Run/i");
    fReco->Branch("subRun",     &fSubRun,    "subRun/i");
    //PFP branches
    fReco->Branch("nPFParticles",                       &fNPFParticles,                      "nPFParticles/i");
    fReco->Branch("pfpTrueParticleMatchedID",           &fPFPTrueParticleMatchedID,           "PFPTrueParticleMatchedID[nPFParticles]/I");
    fReco->Branch("pfpTrueParticleMatchedPosition",     &fPFPTrueParticleMatchedPosition,     "PFPTrueParticleMatchedPosition[nPFParticles]/I");
    fReco->Branch("pfpIsPrimary",                       &fPFPIsPrimary,                       "PFPIsPrimary[nPFParticles]/O");
    fReco->Branch("pfpID",                              &fPFPID,                              "PFPID[nPFParticles]/I");
    fReco->Branch("pfpParentID",                        &fPFPParentID,                        "PFPParentID[nPFParticles]/I");
    fReco->Branch("pfpPdgCode",                         &fPFPPdgCode,                         "PFPPdgCode[nPFParticles]/I");
    fReco->Branch("pfpNChildren",                       &fPFPNChildren,                       "PFPNChildren[nPFParticles]/I");
    
    // fReco->Branch("pfpNHits",                           &fPFPNHits,                           "PFPNHits[nPFParticles]/I");
    // fReco->Branch("pfpNHitsView",                       &fPFPNHitsView,                       "PFPNHitsView[nPFParticles][3]/I");
    // fReco->Branch("pfpNSharedTrueParticleHits",         &fPFPNSharedTrueParticleHits,         "PFPNSharedTrueParticleHits[nPFParticles]/I");
    // fReco->Branch("pfpNSharedTrueParticleHitsView",     &fPFPNSharedTrueParticleHitsView,     "PFPNSharedTrueParticleHitsView[nPFParticles][3]/I");
    // fReco->Branch("pfpTrueParticleMatchedIDView",       &fPFPTrueParticleMatchedIDView,       "PFPTrueParticleMatchedIDView[nPFParticles][3]/I");
    // fReco->Branch("pfpTrueParticleMatchedPositionView", &fPFPTrueParticleMatchedPositionView, "PFPTrueParticleMatchedPositionView[nPFParticles][3]/I");
    
    fReco->Branch("pfpIsTrack",               &fPFPIsTrack,              "PFPIsTrack[nPFParticles]/O");
    fReco->Branch("pfpTrackID",               &fPFPTrackID,              "PFPNClusters[nPFParticles]/I");
    fReco->Branch("pfpTrackLength",           &fPFPTrackLength,          "PFPTrackLength[nPFParticles]/D");
    fReco->Branch("pfpTrackStartX",           &fPFPTrackStartX,          "PFPTrackStartX[nPFParticles]/D");
    fReco->Branch("pfpTrackStartY",           &fPFPTrackStartY,          "PFPTrackStartY[nPFParticles]/D");
    fReco->Branch("pfpTrackStartZ",           &fPFPTrackStartZ,          "PFPTrackStartZ[nPFParticles]/D");
    fReco->Branch("pfpTrackVertexX",          &fPFPTrackVertexX,         "PFPTrackVertexX[nPFParticles]/D");
    fReco->Branch("pfpTrackVertexY",          &fPFPTrackVertexY,         "PFPTrackVertexY[nPFParticles]/D");
    fReco->Branch("pfpTrackVertexZ",          &fPFPTrackVertexZ,         "PFPTrackVertexZ[nPFParticles]/D");
    fReco->Branch("pfpTrackEndX",             &fPFPTrackEndX,            "PFPTrackEndX[nPFParticles]/D");
    fReco->Branch("pfpTrackEndY",             &fPFPTrackEndY,            "PFPTrackEndY[nPFParticles]/D");
    fReco->Branch("pfpTrackEndZ",             &fPFPTrackEndZ,            "PFPTrackEndZ[nPFParticles]/D");
    fReco->Branch("pfpTrackTheta",            &fPFPTrackTheta,           "PFPTrackTheta[nPFParticles]/D");
    fReco->Branch("pfpTrackPhi",              &fPFPTrackPhi,             "PFPTrackPhi[nPFParticles]/D");
    fReco->Branch("pfpTrackZenithAngle",      &fPFPTrackZenithAngle,     "PFPTrackZenithAngle[nPFParticles]/D");
    fReco->Branch("pfpTrackAzimuthAngle",     &fPFPTrackAzimuthAngle,    "PFPTrackAzimuthAngle[nPFParticles]/D");
    fReco->Branch("pfpTrackStartDirectionX",  &fPFPTrackStartDirectionX, "PFPTrackStartDirectionX[nPFParticles]/D");
    fReco->Branch("pfpTrackStartDirectionY",  &fPFPTrackStartDirectionY, "PFPTrackStartDirectionY[nPFParticles]/D");
    fReco->Branch("pfpTrackStartDirectionZ",  &fPFPTrackStartDirectionZ, "PFPTrackStartDirectionZ[nPFParticles]/D");
    fReco->Branch("pfpTrackVertexDirectionX", &fPFPTrackVertexDirectionX ,"PFPTrackVertexDirectionX[nPFParticles]/D");
    fReco->Branch("pfpTrackVertexDirectionY", &fPFPTrackVertexDirectionY ,"PFPTrackVertexDirectionY[nPFParticles]/D");
    fReco->Branch("pfpTrackVertexDirectionZ", &fPFPTrackVertexDirectionZ ,"PFPTrackVertexDirectionZ[nPFParticles]/D");
    fReco->Branch("pfpTrackEndDirectionX",    &fPFPTrackEndDirectionX,    "PFPTrackEndDirectionX[nPFParticles]/D");
    fReco->Branch("pfpTrackEndDirectionY",    &fPFPTrackEndDirectionY,    "PFPTrackEndDirectionY[nPFParticles]/D");
    fReco->Branch("pfpTrackEndDirectionZ",    &fPFPTrackEndDirectionZ,    "PFPTrackEndDirectionZ[nPFParticles]/D");
    fReco->Branch("pfpTrackChi2",             &fPFPTrackChi2,             "PFPTrackChi2[nPFParticles]/F");
    fReco->Branch("pfpTrackStartNdof",        &fPFPTrackNdof,             "PFPTrackNdof[nPFParticles]/I");

    fReco->Branch("pfpNClusters",                       &fPFPNClusters,                       "PFPNClusters[nPFParticles]/I");
    fReco->Branch("pfpCluNHits",    fPFPCluNHits,    "PFPCluNHits[nPFParticles][100]/I");
    // fReco->Branch("pfpCluPlane",    fPFPCluPlane,    "PFPCluPlane[nPFParticles][100]/I");
    // fReco->Branch("pfpCluView",     fPFPCluView,     "PFPCluView[nPFParticles][100]/I");
    // fReco->Branch("pfpCluIntegral", fPFPCluIntegral, "PFPCluIntegral[nPFParticles][100]/D");

    fReco->Branch("pfpIsShower",            &fPFPIsShower,             "PFPIsShower[nPFParticles]/O");
    fReco->Branch("pfpShowerID",            &fPFPShowerID,            "PFPShowerID[nPFParticles]/I");
    fReco->Branch("pfpShowerBestPlane",     &fPFPShowerBestPlane,     "PFPShowerBestPlane[nPFParticles]/I");
    fReco->Branch("pfpShowerDirectionX",    &fPFPShowerDirectionX,    "PFPShowerDirectionX[nPFParticles]/D");
    fReco->Branch("pfpShowerDirectionY",    &fPFPShowerDirectionY,    "PFPShowerDirectionY[nPFParticles]/D");
    fReco->Branch("pfpShowerDirectionZ",    &fPFPShowerDirectionZ,    "PFPShowerDirectionZ[nPFParticles]/D");
    // fReco->Branch("pfpShowerDirectionErrX", &fPFPShowerDirectionErrX, "PFPShowerDirectionErrX[nPFParticles]/D");
    // fReco->Branch("pfpShowerDirectionErrY", &fPFPShowerDirectionErrY, "PFPShowerDirectionErrY[nPFParticles]/D");
    // fReco->Branch("pfpShowerDirectionErrZ", &fPFPShowerDirectionErrZ, "PFPShowerDirectionErrZ[nPFParticles]/D");
    fReco->Branch("pfpShowerStartX",        &fPFPShowerStartX,        "PFPShowerStartX[nPFParticles]/D");
    fReco->Branch("pfpShowerStartY",        &fPFPShowerStartY,        "PFPShowerStartY[nPFParticles]/D");
    fReco->Branch("pfpShowerStartZ",        &fPFPShowerStartZ,        "PFPShowerStartZ[nPFParticles]/D");
    // fReco->Branch("pfpShowerStartErrX",     &fPFPShowerStartErrX,     "PFPShowerStartErrX[nPFParticles]/D");
    // fReco->Branch("pfpShowerStartErrY",     &fPFPShowerStartErrY,     "PFPShowerStartErrY[nPFParticles]/D");
    // fReco->Branch("pfpShowerStartErrZ",     &fPFPShowerStartErrZ,     "PFPShowerStartErrZ[nPFParticles]/D");
    fReco->Branch("pfpShowerLength",        &fPFPShowerLength,        "PFPShowerLength[nPFParticles]/D");
    fReco->Branch("pfpShowerOpeningAngle",  &fPFPShowerOpenAngle,     "PFPShowerOpenAngle[nPFParticles]/D");

    fReco->Branch("pfpCompleteness",      &fPFPCompleteness,     "PFPCompleteness[nPFParticles]/D");
    fReco->Branch("pfpCompletenessView",  &fPFPCompletenessView, "PFPCompletenessView[nPFParticles][3]/D");
    fReco->Branch("pfpPurity",            &fPFPPurity,           "PFPPurity[nPFParticles]/D");
    fReco->Branch("pfpPurityView",        &fPFPPurityView,       "PFPPurityView[nPFParticles][3]/D");

  fRollUpUnsavedIDs    = p.get<bool>("RollUpUnsavedIDs"); 
  fGeom                = &*art::ServiceHandle<geo::Geometry>();
  // this->reconfigure(p);
} 


void ana::MyAnalysis::reconfigure(fhicl::ParameterSet const& p)
{
  // assign fTruthLabel in the constructor using the parameter defined in MyAnalysis.fcl
  // fTruthLabel          = p.get<std::string>("TruthLabel");
  // fEdepLabel           = p.get<std::string>("EdepLabel");
  // fEdepInstanceLabel   = p.get<std::string>("InstanceName");
  // fDetSimLabel         = p.get<std::string>("DetSimLabel");
  // fRawDigInstanceLabel = p.get<std::string>("RawDigInstance");
  // fSimChaInstanceLabel = p.get<std::string>("SimChaInstance");
  // fRollUpUnsavedIDs    = p.get<bool>("RollUpUnsavedIDs"); 
  // fGeom                = &*art::ServiceHandle<geo::Geometry>();
} // Reconfigure


void ana::MyAnalysis::beginJob()
{
  // Implementation of optional member function here.
  // reset(true); //deep clean the variables

  // Call appropriate consumes<>() for any products to be retrieved by this module.
}//beginJob



void ana::MyAnalysis::analyze(const art::Event & evt)
{
  // reset(); //Don't deep clean

  fEvent  = evt.id().event();
  fRun    = evt.id().run();
  fSubRun = evt.id().subRun();
  std::map<int, int> fTrackIDMap; // store out of the loop for the truth_tree a map of trackID to MCParticle index

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt); //Contains all timing reference information for the detector. 
  // const art::ServiceHandle<cheat::BackTrackerService> btServ;

  std::cout << "=============== EVENT ID " << fEvent << " == RUN ID " << fRun << " == SUBRUN ID " << fSubRun << " ================" << std::endl;

  for (auto i : fProductsToDump)
  {
    auto i_label = i[0]; auto i_instance = i[1]; auto i_type = i[2]; auto o_label = i[3];
    if (debug){ std::cout << "---> Module_Label: " << i_label << ";\t Instance_Name: " << i_instance << ";\t Product_Type: " << i_type << ";\t Output_Name: " << o_label << std::endl; }

    //================================================== Tree Truth ==================================================//
    //Access the truth information//
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


    // Access Deposited Energy 
    if (i_type == "sim::SimEnergyDeposit")
    {
      // fEdepLabel = i_label; 
      // auto EdepHandle = evt.getValidHandle<std::vector<sim::SimEnergyDeposit>>(fEdepLabel);
      // evt.getByLabel(fEdepLabel, fEdepInstanceLabel, EdepHandle);
      // art::ValidHandle<std::vector<sim::SimEnergyDeposit>> EdepHandle = evt.getValidHandle<std::vector<sim::SimEnergyDeposit>>(i_label);
      
      art::Handle<std::vector<sim::SimEnergyDeposit>> EdepHandle;
      evt.getByLabel(i_label, i_instance, EdepHandle);
      if(!EdepHandle.isValid()) { std::cout<<"Unable to find std::vector<sim::SimEnergyDeposit> with module label: " << i_label << std::endl; } // !EdepHandle
      if(EdepHandle.isValid())
      {
        for(auto edep : *EdepHandle) 
        {
          int particle_id = fTrackIDMap[edep.TrackID()];
          fSimEdepTrackID[particle_id] = edep.TrackID();
          fSimEdepPDGCode[particle_id] = edep.PdgCode();
          fSimEdepX      [particle_id] = edep.MidPointX();
          fSimEdepY      [particle_id] = edep.MidPointY();
          fSimEdepZ      [particle_id] = edep.MidPointZ();
          fSimEdepE      [particle_id] = edep.E();
        } // for(auto edep : *EdepHandle)
      } // EdepHandle.isValid()  


    } //if i_type == "sim::SimEnergyDeposit"

    //================================================== Tree DetSim ==================================================//

    // Access Detector Hits // ACTUALIZAR // MU DUNE-FD
    // if (i_type == "DetSim")
    // {
        // channels.clear(); tdc.clear(); adc.clear(); view.clear();
    //   art::Handle<std::vector<sim::SimChannel>> SimChaHandle;
    //   evt.getByLabel(fDetSimLabel, fSimChaInstanceLabel, SimChaHandle);
    //   if(!SimChaHandle.isValid())
    //   {
    //     std::cout<<"Unable to find std::vector<sim::SimChannel> with module label: " << fDetSimLabel << "; Instance: " << fSimChaInstanceLabel << std::endl;
    //     return;
    //   } // !DetSimHandle
      
    //   art::Handle<std::vector<raw::RawDigit>> RawDigHandle;
    //   evt.getByLabel(fDetSimLabel, fRawDigInstanceLabel, RawDigHandle);
    //   if(!RawDigHandle.isValid())
    //   {
    //     std::cout<<"Unable to find std::vector<raw::RawDigit> with module label: " << fDetSimLabel << "; Instance: " << fRawDigInstanceLabel << std::endl;
    //     return;
    //   } // !RawDigHandle

    //   if(SimChaHandle.isValid() && RawDigHandle.isValid())
    //   {
    //     // channels.clear(); tdc.clear(); adc.clear(); view.clear();

    //     std::vector<sim::SimChannel> sim_channels; 
    //     for(auto sim_channel : *SimChaHandle) { sim_channels.push_back(sim_channel); }

    //     // std::cout << "HERE: " << std::endl;
    //     // std::cout << sim_channels << std::endl;

    //     // channels.clear();

    //     for(auto digit : *RawDigHandle) // loop over all raw digits (i.e ADC counts for each channel for each time tick)
    //     // There are 15360 channels in the ProtoDUNE-SP detector.
    //     // RawDigit is the raw data from the detector. SimChannel is the simulated data from the detector.
    //     // RawDigit is the output of the electronics chain. SimChannel is the output of the simulation chain.
    //     // RawDigit is the ADC counts for each channel for each time tick. SimChannel is the number of electrons
    //     // that were deposited in each channel for each time tick.
    //     {
    //       int num_samples = digit.Samples();          // number of ADC samples (TDC ticks)
    //       int pedestal    = (int)digit.GetPedestal(); // pedestal value
          
    //       // uncompress the digits and remove the pedestal
    //       std::vector<short> uncompressed(num_samples); // uncompressed ADC values
    //       raw::Uncompress( digit.ADCs(), uncompressed, pedestal, digit.Compression()); 
    //       for (int ii = 0; ii < num_samples; ii++) { uncompressed[ii] -= pedestal; } // remove pedestal

    //       std::vector<int> channels_idparticle;
    //       for(int tdc_tick = 0; tdc_tick < num_samples; tdc_tick++)
    //       {
    //         auto channel     = digit.Channel();        // channel number
    //         auto sim_channel = sim_channels[channel];  // sim channel

    //         channels_idparticle.push_back(sim_channel.Channel()); // channel number
    //         // std::cout << channels_idparticle.size() << std::endl;
    //       // geo::GeometryCore const* geometry_core = lar::providerFrom<geo::Geometry>(); // geometry provider --> which wire plane
          
    //         // auto const& trackIDsAndEnergy = sim_channel.TrackIDsAndEnergies(l, l);
    //         // if(std::abs(uncompressed[tdc_tick]) > 20.0) // ADC threshold > 20
    //         // {
    //           // std::vector<sim::IDE> ides = sim_channel.TrackIDsAndEnergies(tdc_tick, tdc_tick);
    //           // for (auto const& ide : ides)
    //           // {
    //           //   int particle_id = TrackIDMap[ide.trackID];
    //           //   std::cout << "particle_id: " << particle_id << std::endl;
    //           //   // std::cout << "channel: " << sim_channel.Channel() << std::endl;

    //           //   // channels[particle_id].push_back(sim_channel.Channel()); // channel number
    //           //   // channels_idparticle[particle_id] = sim_channel.Channel(); 
    //           // }
    //       }
    //       if (channels_idparticle.empty()) { std::cout << "GREAT" << std::endl; channels.push_back(std::vector<int>()); } 
    //       else { channels.push_back(channels_idparticle); }
    //       // channels.push_back(channels_idparticle); // channel number



    //       // for(int tdc_tick = 0; tdc_tick < num_samples; tdc_tick++)
    //       // {
    //       //   // auto const& trackIDsAndEnergy = sim_channel.TrackIDsAndEnergies(l, l);
    //       //   if(std::abs(uncompressed[tdc_tick]) > 20.0) // ADC threshold > 20
    //       //   {
    //       //     std::vector<sim::IDE> ides = sim_channel.TrackIDsAndEnergies(tdc_tick, tdc_tick);
    //       //     // std::cout << "TrackID (sim): " << std::endl;
    //       //     for (auto const& ide : ides)
    //       //     {
    //       //       // std::cout << "TrackID (sim): " << ide.trackID << std::endl;
    //       //       int particle_id = TrackIDMap[ide.trackID];
    //       //       std::cout << "particle_id: " << particle_id << std::endl;
    //       //       channels[particle_id] = channels_element; 
    //       //       // std::cout<<"Channels: " << channels << std::endl;
    //       //       tdc[particle_id] = tdc_tick;
    //       //       adc[particle_id] =uncompressed[tdc_tick];
    //       //       // view[particle_id] = geometry_core->View(channels_element);
    //       //     }
    //       //   }
    //       // }
          




    //       // for (size_t i = 0; i < sim_channel.TDC_t().size(); ++i)
    //       // {
    //       //   int track_id = sim_channel.TDCIDEs_t()[i].trackID;
    //       //   int tdc      = sim_channel.TDC_t()[i];
    //       //   int adc      = sim_channel.Charge(tdc)[i];
    //       //   int channel  = sim_channel.Channel();
    //       //   int view     = geo->View(channel);
    //       //   channels.push_back(channel);
    //       //   tdc.push_back(tdc);
    //       //   adc.push_back(adc);
    //       //   view.push_back(view);
    //       // }
    //     }
    //   }
    // }
     
    // if (!channels.empty()) 
    // {
    //   std::cout << "GREAT channels not empty" << std::endl;
    //   std::cout << channels.size() << std::endl;
    // }


  } // fProductsToDump loop

    fTruth->Fill();
    // fDetSim->Fill();

    //================================================== Tree Reco ==================================================//

    // if (i_type == "recob::PFParticle")
    // {

      //Access the PFParticles from Pandora
      const std::vector<art::Ptr<recob::PFParticle>> pfparticleVect = dune_ana::DUNEAnaEventUtils::GetPFParticles(evt,"pandora");
      fNPFParticles = pfparticleVect.size();
      if(!fNPFParticles) { std::cout << "No PFParticles found!" << std::endl; }
      std::cout << "Number of PFParticles: " << fNPFParticles << std::endl;

      //Access the Clusters
      std::vector<art::Ptr<recob::Cluster>> clusterVect;
      auto clusterHandle = evt.getHandle<std::vector<recob::Cluster> >("pandora");
      if (clusterHandle){ art::fill_ptr_vector(clusterVect,clusterHandle); } 
      art::FindManyP<recob::Cluster> clusterParticleAssoc(pfparticleVect, evt, "pandora");

      //Access the Tracks
      auto trackHandle = evt.getHandle<std::vector<recob::Track> >("pandoraTrack");
      if (!trackHandle) { std::cout<<"Unable to find std::vector<recob::Track> with module label: " << "pandoraTrack" << std::endl; }
      std::vector<art::Ptr<recob::Track> > trackList;
      art::fill_ptr_vector(trackList, trackHandle);

      //std::map<int,int> pfpToMcMap;
      int iPfp(0);
      for(const art::Ptr<recob::PFParticle> &pfp: pfparticleVect) // loop over all PFParticles
      {
        fPFPID       [iPfp] = iPfp;
        fPFPIsPrimary[iPfp] = pfp->IsPrimary();
        fPFPPdgCode  [iPfp] = pfp->PdgCode();
        fPFPNChildren[iPfp] = pfp->NumDaughters();
        (pfp->IsPrimary())? fPFPParentID[iPfp] = -1 : fPFPParentID[iPfp] = pfp->Parent(); // if primary, parent ID = -1

        std::vector<art::Ptr<recob::Cluster>> pfpClusters = clusterParticleAssoc.at(pfp.key()); // get clusters associated to the PFParticle
        fPFPNClusters[iPfp] = pfpClusters.size(); // number of clusters associated to the PFParticle
        if(!pfpClusters.empty())
        {
          int iClu(0);
          for(const art::Ptr<recob::Cluster> &clu:pfpClusters)
          {
            fPFPCluPlane   [iPfp][iClu] = clu->Plane().asPlaneID().Plane;
            fPFPCluView    [iPfp][iClu] = clu->View();
            fPFPCluNHits   [iPfp][iClu] = clu->NHits();
            fPFPCluIntegral[iPfp][iClu] = clu->Integral();
            iClu++;
            if (iClu == kNMaxPFPClusters){ break; }
          }
        }


        std::vector<art::Ptr<recob::Hit>> pfpHits; 
        if(dune_ana::DUNEAnaPFParticleUtils::IsTrack(pfp, evt,"pandora","pandoraTrack")) // if PFParticle is a track
        {
          art::Ptr<recob::Track> track = dune_ana::DUNEAnaPFParticleUtils::GetTrack(pfp, evt, "pandora", "pandoraTrack"); 
          fPFPIsTrack              [iPfp] = true;
          fPFPTrackID              [iPfp] = track->ID();
          fPFPTrackLength          [iPfp] = track->Length();
          fPFPTrackStartX          [iPfp] = track->Start().X();
          fPFPTrackStartY          [iPfp] = track->Start().Y();
          fPFPTrackStartZ          [iPfp] = track->Start().Z();
          fPFPTrackVertexX         [iPfp] = track->Vertex().X();
          fPFPTrackVertexY         [iPfp] = track->Vertex().Y();
          fPFPTrackVertexZ         [iPfp] = track->Vertex().Z();
          fPFPTrackEndX            [iPfp] = track->End().X();
          fPFPTrackEndY            [iPfp] = track->End().Y();
          fPFPTrackEndZ            [iPfp] = track->End().Z();
          fPFPTrackTheta           [iPfp] = track->Theta();
          fPFPTrackPhi             [iPfp] = track->Phi();
          fPFPTrackZenithAngle     [iPfp] = track->ZenithAngle();
          fPFPTrackAzimuthAngle    [iPfp] = track->AzimuthAngle();
          fPFPTrackStartDirectionX [iPfp] = track->StartDirection().X();
          fPFPTrackStartDirectionY [iPfp] = track->StartDirection().Y();
          fPFPTrackStartDirectionZ [iPfp] = track->StartDirection().Z();
          fPFPTrackVertexDirectionX[iPfp] = track->VertexDirection().X();
          fPFPTrackVertexDirectionY[iPfp] = track->VertexDirection().Y();
          fPFPTrackVertexDirectionZ[iPfp] = track->VertexDirection().Z();
          fPFPTrackEndDirectionX   [iPfp] = track->EndDirection().X();
          fPFPTrackEndDirectionY   [iPfp] = track->EndDirection().Y();
          fPFPTrackEndDirectionZ   [iPfp] = track->EndDirection().Z();
          fPFPTrackChi2            [iPfp] = track->Chi2();
          fPFPTrackNdof            [iPfp] = track->Ndof();

          pfpHits = dune_ana::DUNEAnaTrackUtils::GetHits(track,evt,"pandoraTrack"); // get hits associated to the track
          
          std::vector<art::Ptr<recob::Hit> > pfpHitsView0 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pfp, evt, "pandora", 0);
          std::vector<art::Ptr<recob::Hit> > pfpHitsView1 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pfp, evt, "pandora", 1);
          std::vector<art::Ptr<recob::Hit> > pfpHitsView2 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pfp, evt, "pandora", 2);
          fPFPNHits    [iPfp]    = pfpHits.size();
          fPFPNHitsView[iPfp][0] = pfpHitsView0.size();
          fPFPNHitsView[iPfp][1] = pfpHitsView1.size();
          fPFPNHitsView[iPfp][2] = pfpHitsView2.size();

          if(!evt.isRealData()) // if MC
          {
            TruthMatchUtils::G4ID g4ID(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, pfpHits, fRollUpUnsavedIDs));
            if (TruthMatchUtils::Valid(g4ID))
            {
              fPFPTrueParticleMatchedID[iPfp] = g4ID; int pos(999999); 
              for(int unsigned ipos=0; ipos<fNMCParticles; ipos++) 
              {
                if(fMCParticleTrackID[ipos]==g4ID) {pos = ipos;} // get the position of the MCParticle in the MCParticle tree
              }
              fPFPTrueParticleMatchedPosition[iPfp] = pos;
            }

            TruthMatchUtils::G4ID g4IDView0(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,pfpHitsView0,fRollUpUnsavedIDs));
            if (TruthMatchUtils::Valid(g4ID))
            {
              fPFPTrueParticleMatchedIDView[iPfp][0] = g4IDView0;
              for(int unsigned ipos=0; ipos<fNMCParticles; ipos++) 
              {
                if (fMCParticleTrackID[ipos] != g4IDView0){ continue; }
                fPFPTrueParticleMatchedPositionView[iPfp][0] = ipos;
                break;
              }
            }

            TruthMatchUtils::G4ID g4IDView1(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, pfpHitsView1, fRollUpUnsavedIDs));
            if (TruthMatchUtils::Valid(g4ID))
            {
              fPFPTrueParticleMatchedIDView[iPfp][1] = g4IDView1;
              for(int unsigned ipos=0; ipos<fNMCParticles; ipos++) 
              {
                if (fMCParticleTrackID[ipos] != g4IDView1){ continue; }
                fPFPTrueParticleMatchedPositionView[iPfp][1] = ipos;
                break;
              }
            }

            TruthMatchUtils::G4ID g4IDView2(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,pfpHitsView2,fRollUpUnsavedIDs));
            if (TruthMatchUtils::Valid(g4ID))
            {
              fPFPTrueParticleMatchedIDView[iPfp][2] = g4IDView2;
              for(int unsigned ipos=0; ipos<fNMCParticles; ipos++) 
              {
                if (fMCParticleTrackID[ipos] != g4IDView2){ continue; }
                fPFPTrueParticleMatchedPositionView[iPfp][2] = ipos;
                break;
              }
            }
          } // if !evt.isRealData()
        } // if IsTrack

        if(dune_ana::DUNEAnaPFParticleUtils::IsShower(pfp, evt, "pandora", "pandoraShower")) // if PFParticle is a shower
        {
          art::Ptr<recob::Shower> shower = dune_ana::DUNEAnaPFParticleUtils::GetShower(pfp, evt, "pandora", "pandoraShower"); 
          fPFPIsShower           [iPfp] = true;
          fPFPShowerID           [iPfp] = shower->ID();
          fPFPShowerBestPlane    [iPfp] = shower->best_plane();
          fPFPShowerDirectionX   [iPfp] = shower->Direction().X();
          fPFPShowerDirectionY   [iPfp] = shower->Direction().Y();
          fPFPShowerDirectionZ   [iPfp] = shower->Direction().Z();
          fPFPShowerDirectionErrX[iPfp] = shower->DirectionErr().X();
          fPFPShowerDirectionErrY[iPfp] = shower->DirectionErr().Y();
          fPFPShowerDirectionErrZ[iPfp] = shower->DirectionErr().Z();
          fPFPShowerStartX       [iPfp] = shower->ShowerStart().X();
          fPFPShowerStartY       [iPfp] = shower->ShowerStart().Y();
          fPFPShowerStartZ       [iPfp] = shower->ShowerStart().Z();
          fPFPShowerStartErrX    [iPfp] = shower->ShowerStartErr().X();
          fPFPShowerStartErrY    [iPfp] = shower->ShowerStartErr().Y();
          fPFPShowerStartErrZ    [iPfp] = shower->ShowerStartErr().Z();
          fPFPShowerLength       [iPfp] = shower->Length();
          fPFPShowerOpenAngle    [iPfp] = shower->OpenAngle();

          pfpHits = dune_ana::DUNEAnaShowerUtils::GetHits(shower, evt, "pandoraShower");
          std::vector<art::Ptr<recob::Hit> > pfpHitsView0 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pfp, evt, "pandora", 0);
          std::vector<art::Ptr<recob::Hit> > pfpHitsView1 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pfp, evt, "pandora", 1);
          std::vector<art::Ptr<recob::Hit> > pfpHitsView2 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pfp, evt, "pandora", 2);
          fPFPNHits    [iPfp]    = pfpHits.size();
          fPFPNHitsView[iPfp][0] = pfpHitsView0.size();
          fPFPNHitsView[iPfp][1] = pfpHitsView1.size();
          fPFPNHitsView[iPfp][2] = pfpHitsView2.size();

          if(!evt.isRealData()) // if MC
          {
            TruthMatchUtils::G4ID g4ID(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, pfpHits, fRollUpUnsavedIDs));
            if (TruthMatchUtils::Valid(g4ID))
            {
              fPFPTrueParticleMatchedID[iPfp] = g4ID; int pos(999999); 
              for(int unsigned ipos=0; ipos<fNMCParticles; ipos++) 
              { 
                if(fMCParticleTrackID[ipos]==g4ID){ pos=ipos; }
              }
              fPFPTrueParticleMatchedPosition[iPfp] = pos;
            }

            TruthMatchUtils::G4ID g4IDView0(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, pfpHitsView0, fRollUpUnsavedIDs));
            if (TruthMatchUtils::Valid(g4ID))
            {
              fPFPTrueParticleMatchedIDView[iPfp][0] = g4IDView0;
              for(int unsigned ipos=0; ipos<fNMCParticles; ipos++) 
              {
                if (fMCParticleTrackID[ipos] != g4IDView0) { continue; }
                fPFPTrueParticleMatchedPositionView[iPfp][0] = ipos;
                break;
              }
            }

            TruthMatchUtils::G4ID g4IDView1(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, pfpHitsView1, fRollUpUnsavedIDs));
            if (TruthMatchUtils::Valid(g4ID))
            {
              fPFPTrueParticleMatchedIDView[iPfp][1] = g4IDView1;
              for(int unsigned ipos=0; ipos<fNMCParticles; ipos++) 
              {
                if (fMCParticleTrackID[ipos] != g4IDView1){ continue; }
                fPFPTrueParticleMatchedPositionView[iPfp][1] = ipos;
                break;
              }
            }

            TruthMatchUtils::G4ID g4IDView2(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, pfpHitsView2, fRollUpUnsavedIDs));
            if (TruthMatchUtils::Valid(g4ID))
            {
              fPFPTrueParticleMatchedIDView[iPfp][2] = g4IDView2;
              for(int unsigned ipos=0; ipos<fNMCParticles; ipos++) 
              {
                if (fMCParticleTrackID[ipos] != g4IDView2){ continue; }
                fPFPTrueParticleMatchedPositionView[iPfp][2] = ipos;
                break;
              }
            }
          } // if(!evt.isRealData())
        } // if(dune_ana::DUNEAnaPFParticleUtils::IsShower(pfp, evt, "pandora", "pandoraShower"))

        if(!evt.isRealData()) // if MC
        {
          //Fill shared MC particle to hits map for this specific PFP
          for (const auto& hit: pfpHits)
          {
            TruthMatchUtils::G4ID g4ID(TruthMatchUtils::TrueParticleID(clockData, hit, fRollUpUnsavedIDs));
            if (TruthMatchUtils::Valid(g4ID))
            {
              if(g4ID==fPFPTrueParticleMatchedID[iPfp])                  { ++fPFPNSharedTrueParticleHits    [iPfp]   ; }
              if(g4ID==fPFPTrueParticleMatchedID[iPfp] && hit->View()==0){ ++fPFPNSharedTrueParticleHitsView[iPfp][0]; }
              if(g4ID==fPFPTrueParticleMatchedID[iPfp] && hit->View()==1){ ++fPFPNSharedTrueParticleHitsView[iPfp][1]; }
              if(g4ID==fPFPTrueParticleMatchedID[iPfp] && hit->View()==2){ ++fPFPNSharedTrueParticleHitsView[iPfp][2]; }
            }
          }

          if(fPFPNHits    [iPfp]    > 0 && fPFPNHits    [iPfp]    < 999999) fPFPPurity    [iPfp]    = (float)fPFPNSharedTrueParticleHits    [iPfp]    / fPFPNHits    [iPfp];
          if(fPFPNHitsView[iPfp][0] > 0 && fPFPNHitsView[iPfp][0] < 999999) fPFPPurityView[iPfp][0] = (float)fPFPNSharedTrueParticleHitsView[iPfp][0] / fPFPNHitsView[iPfp][0];
          if(fPFPNHitsView[iPfp][1] > 0 && fPFPNHitsView[iPfp][1] < 999999) fPFPPurityView[iPfp][1] = (float)fPFPNSharedTrueParticleHitsView[iPfp][1] / fPFPNHitsView[iPfp][1];
          if(fPFPNHitsView[iPfp][2] > 0 && fPFPNHitsView[iPfp][2] < 999999) fPFPPurityView[iPfp][2] = (float)fPFPNSharedTrueParticleHitsView[iPfp][2] / fPFPNHitsView[iPfp][2];

          if(fPFPTrueParticleMatchedPosition[iPfp]<999999 && fMCParticleNHits    [fPFPTrueParticleMatchedPosition[iPfp]]    > 0 && fMCParticleNHits    [fPFPTrueParticleMatchedPosition[iPfp]]    < 999999) fPFPCompleteness    [iPfp]    = (float)fPFPNSharedTrueParticleHits    [iPfp]    / fMCParticleNHits    [fPFPTrueParticleMatchedPosition[iPfp]];
          if(fPFPTrueParticleMatchedPosition[iPfp]<999999 && fMCParticleNHitsView[fPFPTrueParticleMatchedPosition[iPfp]][0] > 0 && fMCParticleNHitsView[fPFPTrueParticleMatchedPosition[iPfp]][0] < 999999) fPFPCompletenessView[iPfp][0] = (float)fPFPNSharedTrueParticleHitsView[iPfp][0] / fMCParticleNHitsView[fPFPTrueParticleMatchedPosition[iPfp]][0];
          if(fPFPTrueParticleMatchedPosition[iPfp]<999999 && fMCParticleNHitsView[fPFPTrueParticleMatchedPosition[iPfp]][1] > 0 && fMCParticleNHitsView[fPFPTrueParticleMatchedPosition[iPfp]][1] < 999999) fPFPCompletenessView[iPfp][1] = (float)fPFPNSharedTrueParticleHitsView[iPfp][1] / fMCParticleNHitsView[fPFPTrueParticleMatchedPosition[iPfp]][1];
          if(fPFPTrueParticleMatchedPosition[iPfp]<999999 && fMCParticleNHitsView[fPFPTrueParticleMatchedPosition[iPfp]][2] > 0 && fMCParticleNHitsView[fPFPTrueParticleMatchedPosition[iPfp]][2] < 999999) fPFPCompletenessView[iPfp][2] = (float)fPFPNSharedTrueParticleHitsView[iPfp][2] / fMCParticleNHitsView[fPFPTrueParticleMatchedPosition[iPfp]][2];
        } // if(!evt.isRealData())

        iPfp++;

    } // for(const art::Ptr<recob::PFParticle> &pfp: pfparticleVect)

    fReco->Fill();
}// analyze()



void ana::MyAnalysis::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(ana::MyAnalysis)
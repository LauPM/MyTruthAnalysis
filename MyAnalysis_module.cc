////////////////////////////////////////////////////////////////////////
// Class:       MyAnalysis
// Plugin Type: analyzer (art v3_05_01)
// File:        MyAnalysis_module.cc
//
// Generated at Tue Mar 21 10:26:24 2023 by Laura Pérez-Molina
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
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
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/ArtDataHelper/TrackUtils.h"
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

#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dunereco/AnaUtils/DUNEAnaTrackUtils.h"
#include "dunereco/AnaUtils/DUNEAnaShowerUtils.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"


namespace ana 
{
  class MyAnalysis;
}


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
    std::string fEdepLabel;
    bool fRollUpUnsavedIDs;

    //==================//
    // Tree information //
    //==================//
    TTree *fTree;

    //-------------------//
    // Event Information //
    //-------------------//
    unsigned int fEventID;
    unsigned int fRunID;
    unsigned int fSubRunID;

    static const int kNMaxMCParticles = 100000;
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

    double fMCParticleVertexTime       [kNMaxMCParticles];
    double fMCParticleEndTime          [kNMaxMCParticles];
    int    fMCParticleNHits            [kNMaxMCParticles];
    int    fMCParticleNHitsView        [kNMaxMCParticles][kNViews];
    
    /////////////////////
    //SimEnergyDeposited
    int fSimEdepTrackID[kNMaxMCParticles];
    int fSimEdepPDGCode[kNMaxMCParticles];
    double fSimEdepE   [kNMaxMCParticles];
    double fSimEdepX   [kNMaxMCParticles];
    double fSimEdepY   [kNMaxMCParticles];
    double fSimEdepZ   [kNMaxMCParticles];

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
  fTruthLabel       = p.get<std::string>("TruthLabel");
  fEdepLabel      = p.get<std::string>("EdepLabel");
  fRollUpUnsavedIDs = p.get<bool>("RollUpUnsavedIDs"); 
  fGeom             = &*art::ServiceHandle<geo::Geometry>();
} // Reconfigure


void ana::MyAnalysis::beginJob()
{
  // Implementation of optional member function here.
  // reset(true); //deep clean the variables

  // Call appropriate consumes<>() for any products to be retrieved by this module.
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("Pandora_Output","Pandora_Output");
  //////////////////
  //Event branches//
  //////////////////
  fTree->Branch("eventID",      &fEventID,     "eventID/i");
  fTree->Branch("runID",        &fRunID,       "runID/i");
  fTree->Branch("subrunID",     &fSubRunID,    "subrunID/i");
  /////////////////////
  //MC truth branches//
  /////////////////////
  fTree->Branch("nMCParticles",                &fNMCParticles,               "nMCParticles/i");
  fTree->Branch("mcIsMCPrimary",               &fMCIsPrimary,                "MCIsPrimary[nMCParticles]/O");
  fTree->Branch("mcParticlePdgCode",           &fMCParticlePdgCode,          "MCParticlePdgCode[nMCParticles]/I");
  fTree->Branch("mcParticleTrueEnergy",        &fMCParticleTrueEnergy,       "MCParticleTrueEnergy[nMCParticles]/D");
  fTree->Branch("mcParticleTrackID",           &fMCParticleTrackID,          "MCParticleTrackID[nMCParticles]/I");
  fTree->Branch("mcParticleParentTrackID",     &fMCParticleParentTrackID,    "MCParticleParentTrackID[nMCParticles]/I");
  fTree->Branch("mcParticleNTrajectoryPoints", &fMCParticleNTrajectoryPoint, "MCParticleNTrajectoryPoint[nMCParticles]/I");

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
  
  fTree->Branch("mcParticleVertexTime",        &fMCParticleVertexTime,        "MCParticleVertexTime[nMCParticles]/D");
  fTree->Branch("mcParticleEndTime",           &fMCParticleEndTime,           "MCParticleEndTime[nMCParticles]/D");
  fTree->Branch("mcParticleNHits",             &fMCParticleNHits,             "MCParticleNHits[nMCParticles]/I");
  fTree->Branch("mcParticleNHitsView",         &fMCParticleNHitsView,         "MCParticleNHitsView[nMCParticles][3]/I");
  
  fTree->Branch("SimEdepTrackID",         &fSimEdepTrackID         ,"SimEdepTrackID[nMCParticles]/I");
  fTree->Branch("SimEdepPDGCode",         &fSimEdepPDGCode         ,"SimEdepPDGCode[nMCParticles]/I");
  fTree->Branch("SimEdepEnergy",          &fSimEdepE          ,"fSimEdepE[nMCParticles]/D");
  fTree->Branch("SimEdepMiddlePositionX", &fSimEdepX ,"SimEdepX[nMCParticles]/D");
  fTree->Branch("SimEdepMiddlePositionY", &fSimEdepY ,"SimEdepY[nMCParticles]/D");
  fTree->Branch("SimEdepMiddlePositionZ", &fSimEdepZ ,"SimEdepZ[nMCParticles]/D");

}//beginJob


void ana::MyAnalysis::analyze(const art::Event & evt)
{
  std::map<int, int> TrackIDMap;
  // reset(); //Don't deep clean
  const art::ServiceHandle<cheat::BackTrackerService> btServ;
  fEventID  = evt.id().event();
  fRunID    = evt.id().run();
  fSubRunID = evt.id().subRun();
  // auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt); //Contains all timing reference information for the detector. 

  std::cout << "=============== EVENT ID " << fEventID << " == RUN ID " << fRunID << " == SUBRUN ID " << fSubRunID << " ================" << std::endl;

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
        

        TrackIDMap[trueParticle.TrackId()] = iMc;

      } //for loop filling MCInfo  
    } //if mcParticles.isValid()

  } //if !evt.isRealData()

  // Access Deposited Energy
  if (fEdepLabel != "")
  {
    auto EdepHandle = evt.getHandle<std::vector<sim::SimEnergyDeposit>>(fEdepLabel);
    if(!EdepHandle.isValid())
    {
      std::cout<<"Unable to find std::vector<sim::SimEnergyDeposit> with module label: " << fEdepLabel << std::endl;
      return;
    } // !EdepHandle
    if(EdepHandle.isValid())
    {
      for(auto edep : *EdepHandle) 
      {
        int particle_id = TrackIDMap[edep.TrackID()];
        fSimEdepTrackID [particle_id] = edep.TrackID();
        fSimEdepPDGCode [particle_id] = edep.PdgCode();
        fSimEdepX[particle_id] = edep.MidPointX();
        fSimEdepY[particle_id] = edep.MidPointY();
        fSimEdepZ[particle_id] = edep.MidPointZ();
        fSimEdepE[particle_id] = edep.E();
      }
    }                     
  }

  fTree->Fill();

}// analyze()



void ana::MyAnalysis::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(ana::MyAnalysis)
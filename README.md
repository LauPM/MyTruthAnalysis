# MyTruthAnalysis

Analyzer that dumps:

* G4 stage particles information (Truth)
* Energy depositions
* DetSim [channels, ADCs, energy, tracks view]
* Reco Info

You can configure to run all the analyzer or create the trees you are interested in in MyAnalysis.fcl:
``` bash
TreesToWrite: [["Truth"], ["DetSim"], ["Reco"]]
ProductsToDump: [
                  ["largeant",      "",           "simb::MCParticle",         "G4_MCtruth"],  # Gen:    Montecarlo truth info from Generator Stage (particle PDG, original positions... etc)
                  ["IonAndScint",   "",           "sim::SimEnergyDeposit",    "EDep"],        # Edep:   Energy deposited in the detector
                  ["tpcrawdecoder", "daq",        "raw::RawDigit",            "RawDetSim"],   # DetSim: 
                  ["tpcrawdecoder", "simpleSC",   "sim::SimChannel",          "SimDetSim"],    # DetSim: 
                  ["pandora",       "",           "recob::PFParticle",        "Reco"],         #Reco
                  ["pandoraTrack",  "",           "recob::Track",             "Reco"],         #Reco
                  ["pandoraShower", "",           "recob::Shower",            "Reco"]         #Reco
                ]
 ```


The .fcl file is configured to dump the information in a ROOT file.

## LArSoft preparation

```bash
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
source localProducts_larsoft*/setup

setup larsoft v09_63_00  -q e20:prof  
mrbslp
setup

mrbsetenv
```

## How to use

```bash
lar -c run_MyAnalysis.fcl -s <input_file> -n <number_of_events>
```

## Output example
ROOT file with the analyzed Trees:
![Alt text](output.png)

## TODO


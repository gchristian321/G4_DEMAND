//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file DRAGON/src/DRAGONRunAction.cc
/// \brief Implementation of the DRAGON::DRAGONRunAction class

#include "G4Run.hh"                            //Geant4
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SDManager.hh"
#include "G4Ions.hh"
#include "G4DecayTable.hh"
#include "G4VDecayChannel.hh"
#include "G4HCtable.hh"
#include "G4RunManager.hh"
#include "G4IonTable.hh"
#include "G4RadioactiveDecay.hh"
#include "G4DigiManager.hh"
#include "G4EventManager.hh"


#include "DRAGONRunAction.hh"                  //Local
#include "DRAGONRun.hh"
#include "DRAGONPrimaryGeneratorAction.hh"
#include "DRAGONDetectorConstruction.hh"
#include "DRAGONRunActionMessenger.hh"
#include "DRAGONHistoManager.hh"
#include "DRAGONPhysicsList.hh"
#include "DRAGONSteppingAction.hh"
#include "DRAGONDigitizer.hh"
#include "DRAGONEventAction.hh"


namespace DRAGON
{

GlobalVariables* GlobalVariables::Instance = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DRAGONRunAction::DRAGONRunAction(DRAGONDetectorConstruction* det, DRAGONPrimaryGeneratorAction* prim, DRAGONPhysicsList* phys, DRAGONSteppingAction* step)
:fDetector(det),fPrimary(prim),fphys(phys),fstep(step),fHistoManager(nullptr),fRun(nullptr),idebug(0)
{
 uvinit();

 fHistoManager = new DRAGONHistoManager(fDetector,fPrimary,fphys);
 fRunMessenger = new DRAGONRunActionMessenger(this);
 
 if(fDetector){fDetector->SetRunAction(*this);}
 if(fPrimary){fPrimary->SetRunAction(*this);}
 if(fphys){fphys->SetRunAction(*this);}

 if(fHistoManager){fHistoManager->SetRunAction(*this);}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DRAGONRunAction::~DRAGONRunAction()
{
  delete fHistoManager;
  delete fRunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* DRAGONRunAction::GenerateRun()
{ 
  fRun = new DRAGONRun(); 
  return fRun;
}

void DRAGONRunAction::BeginOfRunAction(const G4Run* run)
{
    uginit();
    TIMINT = std::chrono::high_resolution_clock::now();
	
	auto eventAction = static_cast<DRAGONEventAction*>(G4EventManager::GetEventManager()->GetUserEventAction());
	
	auto DM = G4DigiManager::GetDMpointer();
    if (!DM->FindDigitizerModule("DSSD")) 
	   {
        auto digitizerDSSD = new DRAGONDigitizer("DSSD");
		digitizerDSSD->SetEventAction(eventAction);
		digitizerDSSD->SetDetConst(fDetector);
        DM->AddNewModule(digitizerDSSD);
        G4cout << "Registered digitizer: DSSD" << G4endl;
        }

    if (!DM->FindDigitizerModule("PMT")) 
	   {
        auto digitizerPMT = new DRAGONDigitizer("PMT");
		digitizerPMT->SetEventAction(eventAction);
		digitizerPMT->SetDetConst(fDetector);
        DM->AddNewModule(digitizerPMT);
        G4cout << "Registered digitizer: PMT" << G4endl;
        }

    fHistoManager->uhinit();

if (IsMaster()) {
    fHistoManager->ScanHitos();
	fHistoManager->PrintIDMap();
 return;
}


    if (!fDetector || !fphys) {
        G4cerr << "ERROR: Null pointer in RunAction - "
               << "fDetector: " << (fDetector ? "valid" : "null")
               << ", fphys: "   << (fphys ? "valid"   : "null") << G4endl;
        return;
    }

    if (!fPrimary) {
        G4cerr << "ERROR: fPrimary is null in RunAction (worker thread)" << G4endl;
        return;
    }

    if (!fstep) {
        G4cerr << "ERROR: fstep is null in RunAction" << G4endl;
        return;
    }

    fstep->Setirecoil(fphys->Getirecoil());
    fstep->Setires(fphys->Getires());
    fstep->Setalpha(fphys->Getalpha());
    fstep->Setn_detmate(fDetector->GetNDetMate());
    fstep->Setlkine(fphys->Getlkine());
    fstep->Setitckov(fphys->Getitckov());
    fstep->Settargtype(fDetector->GetTargType());
    fstep->Setresmass(fphys->Getresmass());
    fstep->Setrecoilmom(fphys->Getrecoilmom());
    fstep->Setresenerg(fphys->Getresenerg());
    fstep->Setprodm(fphys->Getprodm());
    fstep->Setbeammass(fphys->Getbeammass());
    fstep->Sete0recoil(fphys->Gete0recoil());
	fstep->Setlen_max(fDetector->Getlen_max());
	
    G4cout << "Run " << run->GetRunID() << ": SteppingAction parameters set successfully" << G4endl;
    std::cout << "ipart      = " << fstep->Getipart()      << std::endl;
    std::cout << "ires       = " << fstep->Getires()      << std::endl;
    std::cout << "irecoil    = " << fphys->Getirecoil()    << std::endl;
    std::cout << "alpha      = " << fphys->Getalpha()      << std::endl;
    std::cout << "n_detmate  = " << fDetector->GetNDetMate() << std::endl;
	std::cout << "Rrms  = " << fDetector->GetRrms() << std::endl;
	std::cout << "lkine      = " << fphys->Getlkine()      << std::endl;
    std::cout << "itckov     = " << fphys->Getitckov()     << std::endl;
    std::cout << "targtype   = " << fDetector->GetTargType() << std::endl;
    std::cout << "resmass    = " << fphys->Getresmass()    << std::endl;
    std::cout << "recoilmom  = " << fphys->Getrecoilmom()  << std::endl;
    std::cout << "resenerg   = " << fphys->Getresenerg()   << std::endl;
    std::cout << "prodm      = " << fphys->Getprodm()      << std::endl;
    std::cout << "beammass   = " << fphys->Getbeammass()   << std::endl;
    std::cout << "e0recoil   = " << fphys->Gete0recoil()   << std::endl;
	std::cout << "len_max   = " << fDetector->Getlen_max()  << std::endl;
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DRAGONRunAction::EndOfRunAction(const G4Run* run)
{
  uglast(run);   
}

void DRAGONRunAction::uvinit()
{
 //From uvinit.f 
 //C.
 //C.    *************** RUN DEFAULTS ***************
 //C.
 irevs = 0;
 idrun = 0; 
 idevt = 0;
 //C.
 //C.    *************** USER DEFAULTS ***************
 //C.
 //C.
 //C.
 //C Init of counters 
 nreact = 0;
 ntargexit = 0;
 nend = 0;
 nfcm2 = 0;
 ngascell = 0;
 nbeamout = 0;
 //C MT adds counters for the number of recoils that make it
 //C   to Q3 (Sext 1) and to Q8 (Quad 6).
 Num_Recoils_Q3 = 0;
 Num_Recoils_Q8 = 0;
 //C MT adds counter for the number of beam particles that reach
 //C   the end detector.
 Num_BeamPart_ENDV = 0;
 n_detector = 0;
 }

void DRAGONRunAction::uginit()
{
 if(idebug != 0)
   {
    PrintCustomParticles();
    PrintMaterials();
	fDetector->PrintVolumesInfo();
    PrintLogicalVolumes();
    PrintSensitiveDetectors();
    }
}


void DRAGONRunAction::PrintCustomParticles()
{
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4cout << "=== User-Defined Particles (Nuclei) ===" << G4endl;

    for (G4int i = 0; i < particleTable->size(); i++) {
        G4ParticleDefinition* particle = particleTable->GetParticle(i);
        if (particle && particle->GetParticleType() == "nucleus") {
            G4Ions* ion = dynamic_cast<G4Ions*>(particle);
            if (!ion) {
                G4cout << "Error: Particle " << particle->GetParticleName() << " is not a valid nucleus." << G4endl;
                continue;
            }

            G4cout << "=== Nucleus: " << particle->GetParticleName() << " ===" << G4endl;
            G4cout << "  Name: " << particle->GetParticleName() << G4endl;
            G4cout << "  PDG Code: " << particle->GetPDGEncoding() << G4endl;
            G4cout << "  Mass: " << particle->GetPDGMass()/GeV << " GeV" << G4endl;
            G4cout << "  Charge: " << particle->GetPDGCharge()/eplus << " e" << G4endl;
            G4cout << "  Spin: " << particle->GetPDGSpin() << G4endl;
            
            if (particle->GetPDGLifeTime() > 0) 
               {G4cout << "  Lifetime: " << particle->GetPDGLifeTime()/ns << " ns" << G4endl;} 
            else 
                {G4cout << "  Lifetime: Stable" << G4endl;}
            
            G4cout << "  Type: " << particle->GetParticleType() << G4endl;
            G4cout << "  Stable: " << (particle->GetPDGStable() ? "Yes" : "No") << G4endl;
            G4cout << "  Short-lived: " << (particle->IsShortLived() ? "Yes" : "No") << G4endl;
            G4cout << "  Baryon Number: " << particle->GetBaryonNumber() << G4endl;
            G4cout << "  Lepton Number: " << particle->GetLeptonNumber() << G4endl;
            G4cout << "  Decay Width: " << particle->GetPDGWidth()/MeV << " MeV" << G4endl;

            G4cout << "  Atomic Number (Z): " << ion->GetAtomicNumber() << G4endl;
            G4cout << "  Mass Number (A): " << ion->GetAtomicMass() << G4endl;
            G4cout << "  Excitation Energy: " << ion->GetExcitationEnergy()/MeV << " MeV" << G4endl;
            G4cout << "  Isomer Level: " << ion->GetIsomerLevel() << G4endl;

            G4DecayTable* decayTable = particle->GetDecayTable();
            if (decayTable) {
                G4cout << "  Decay Scheme:" << G4endl;
                for (G4int j = 0; j < decayTable->entries(); j++) {
                    G4VDecayChannel* channel = decayTable->GetDecayChannel(j);
                    if (channel) {
                        G4cout << "    Channel " << j + 1 << ": ";
                        G4cout << "Branching Ratio = " << channel->GetBR() * 100 << " %" << G4endl;
                        G4cout << "      Daughter Particles: ";
                        for (G4int k = 0; k < channel->GetNumberOfDaughters(); k++) {
                            G4cout << channel->GetDaughter(k)->GetParticleName();
                            if (k < channel->GetNumberOfDaughters() - 1) G4cout << ", ";
                        }
                        G4cout << G4endl;
                    }
                }
            } else {
                G4cout << "  Decay Scheme: Not defined (stable particle or no decay table)" << G4endl;
            }

            G4cout << "----------------------------------------" << G4endl;
        }
    }
}


void DRAGONRunAction::PrintIonProperties(G4ParticleDefinition* ion)
{
   if (!ion) {
        G4cout << "[ERROR] Null ion pointer." << G4endl;
        return;
    }

    if (!ion->GetPDGStable() && !ion->GetDecayTable()) {
        if (ion->GetParticleType() == "nucleus") {
            G4Ions* ionPtr = dynamic_cast<G4Ions*>(ion); 
            if (ionPtr) {
                // G4RadioactiveDecay* radioactiveDecay = new G4RadioactiveDecay();
                // G4DecayTable* decayTable = radioactiveDecay->LoadDecayTable(ionPtr);
                // ion->SetDecayTable(decayTable);
                // delete radioactiveDecay; 
            } else {
                G4cout << "[ERROR] Failed to cast ion to G4Ions for decay table loading." << G4endl;
            }
        } else {
            G4cout << "[WARNING] Ion is not a nucleus, cannot load decay table." << G4endl;
        }
    }

    G4cout << "===== Ion Properties =====" << G4endl;
    G4cout << "Name: " << ion->GetParticleName() << G4endl;
    G4cout << "Type: " << ion->GetParticleType() << G4endl;
    G4cout << "SubType: " << ion->GetParticleSubType() << G4endl;
    G4cout << "PDG Encoding: " << ion->GetPDGEncoding() << G4endl;
    G4cout << "Mass: " << G4BestUnit(ion->GetPDGMass(), "Energy") << G4endl;
    G4cout << "Charge: " << ion->GetPDGCharge() / eplus << " e" << G4endl;
    G4cout << "Spin: " << ion->GetPDGSpin() << G4endl;
    G4cout << "Magnetic Moment: " << ion->GetPDGMagneticMoment() << G4endl;
    G4cout << "Parity: " << ion->GetPDGiParity() << G4endl;
    G4cout << "Isospin: " << ion->GetPDGIsospin() << G4endl;
    G4cout << "Isospin3: " << ion->GetPDGIsospin3() << G4endl;
    G4cout << "Stable: " << (ion->GetPDGStable() ? "Yes" : "No") << G4endl;
    G4cout << "Lifetime: ";
    if (ion->GetPDGStable())
        G4cout << "Stable" << G4endl;
    else
        G4cout << G4BestUnit(ion->GetPDGLifeTime(), "Time") << G4endl;

    // Solo para núcleos
    if (ion->GetParticleType() == "nucleus") {
        G4cout << "Atomic Number (Z): " << ion->GetAtomicNumber() << G4endl;
        G4cout << "Atomic Mass (A): " << ion->GetAtomicMass() << G4endl;

        if (auto ionPtr = dynamic_cast<G4Ions*>(ion)) {
            G4cout << "Excitation Energy: " << G4BestUnit(ionPtr->GetExcitationEnergy(), "Energy") << G4endl;
        } else {
            G4cout << "Excitation Energy: [N/A - not a G4Ions]" << G4endl;
        }
    }

    // Mostrar tabla de decaimiento si existe
    if (!ion->GetPDGStable()) {
        G4DecayTable* decayTable = ion->GetDecayTable();
        if (decayTable) {
            G4cout << "--- Decay Table ---" << G4endl;
            G4int nDecays = decayTable->entries();
            for (G4int i = 0; i < nDecays; ++i) {
                G4VDecayChannel* channel = decayTable->GetDecayChannel(i);
                G4cout << "Mode " << i + 1 << ": "
                       << channel->GetVerboseLevel() << " "
                       << channel->GetBR() * 100.0 << "% -> ";
                for (G4int j = 0; j < channel->GetNumberOfDaughters(); ++j) {
                    G4cout << channel->GetDaughterName(j);
                    if (j < channel->GetNumberOfDaughters() - 1) G4cout << " + ";
                }
                G4cout << G4endl;
            }
        } else {
            G4cout << "[No decay table available]" << G4endl;
        }
    }

    G4cout << "==========================" << G4endl;
}



void DRAGONRunAction::PrintMaterials()
{
    G4MaterialTable* materialTable = G4Material::GetMaterialTable();
    G4cout << "=== Defined materials ===" << G4endl;

    for (const auto* material : *materialTable) {
        G4cout << *(material) << G4endl; 
    }
}


void DRAGONRunAction::PrintLogicalVolumes()
{
    G4LogicalVolumeStore* volumeStore = G4LogicalVolumeStore::GetInstance();
    G4PhysicalVolumeStore* physVolumeStore = G4PhysicalVolumeStore::GetInstance();
    G4cout << "=== Defined Logical Volumes ===" << G4endl;

    for (const auto* volume : *volumeStore) {
        G4cout << "Logical Volume: " << volume->GetName() << G4endl
               << "Material: " << volume->GetMaterial()->GetName() << G4endl
               << "Solid: " << volume->GetSolid()->GetName() << G4endl
               << "Volume: " << volume->GetSolid()->GetCubicVolume()/cm3 << " cm^3" 
               << G4endl;

        G4bool found = false;
        for (const auto* physVolume : *physVolumeStore) {
            if (physVolume->GetLogicalVolume() == volume) {
                found = true;
                G4cout << "  Physical Instance: " << physVolume->GetName() << G4endl;

                G4LogicalVolume* motherLog = physVolume->GetMotherLogical();
                G4cout << "    Mother Volume: " << (motherLog ? motherLog->GetName() : "None (World)") << G4endl;
                
                G4ThreeVector translation = physVolume->GetTranslation();
                G4cout << "    Relative Position: ("
                       << translation.x()/mm << ", "
                       << translation.y()/mm << ", "
                       << translation.z()/mm << ") mm" << G4endl;

                const G4RotationMatrix* rotation = physVolume->GetRotation();
                if (rotation) {
                   G4cout << "    Relative Orientation: Rotation present" << G4endl;
                   G4cout << "      Rotation Matrix:" << G4endl;
                   G4cout << "        [ " << rotation->xx() << ", " << rotation->xy() << ", " << rotation->xz() << " ]" << G4endl;
                   G4cout << "        [ " << rotation->yx() << ", " << rotation->yy() << ", " << rotation->yz() << " ]" << G4endl;
                   G4cout << "        [ " << rotation->zx() << ", " << rotation->zy() << ", " << rotation->zz() << " ]" << G4endl;
                   } else {
                G4cout << "    Relative Orientation: No rotation" << G4endl;
                }
  
            }
        }
        if (!found) {
            G4cout << "  No physical instances found for this logical volume." << G4endl;
        }
        G4cout << "----------------------------------------" << G4endl;
    }
}

void DRAGONRunAction::PrintSensitiveDetectors()
{
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();
    G4HCtable* hcTable = sdManager->GetHCtable();
    if (!hcTable) {
    G4cout << "  No sensitive detectors defined (no hit collection table)." << G4endl;
    return;
    }
    G4cout << "=== Defined Sensitive Detectors ===" << G4endl;

    G4int nCollections = hcTable->entries();
    for (G4int i = 0; i < nCollections; i++) {
        G4String colName = hcTable->GetHCname(i);
        G4cout << "Sensitive Detector Collection: " << colName << G4endl;
    }
    if (nCollections == 0) {
        G4cout << "  No sensitive detectors defined." << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DRAGONRunAction::SetIswit(const G4String& input) {
    std::istringstream iss(input);
    for (G4int i = 0; i < 10; ++i) {
        if (!(iss >> iswit[i])) {
            std::cout << "Invalid input string in method SetIswit. Expected 10 numbers.";
        }
	}
}

}

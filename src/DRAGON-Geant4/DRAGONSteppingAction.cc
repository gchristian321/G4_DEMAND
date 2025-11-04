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
/// \file DRAGON/src/DRAGONSteppingAction.cc
/// \brief Implementation of the DRAGON::DRAGONSteppingAction class

#include "G4Step.hh"                         //Geant4
#include "G4SteppingManager.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4DecayTable.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"
#include "G4RadioactiveDecay.hh"
#include "G4Decay.hh"
#include "G4SystemOfUnits.hh"  
 
#include "DRAGONSteppingAction.hh"           //Local
#include "DRAGONDetectorConstruction.hh"
#include "DRAGONRunAction.hh"
#include "DRAGONEventAction.hh"
#include "DRAGONSteppingActionMessenger.hh"
#include "DRAGONPhysicsList.hh"
#include "DRAGONPrimaryGeneratorAction.hh"
#include "DRAGONHistoManager.hh"

#include "geant3functions.hh"

#include "geom_dipole.hh"
#include "geom_edipol.hh"
#include "geom_mpole.hh"
#include "geom_sole.hh"


namespace DRAGON
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DRAGONSteppingAction::DRAGONSteppingAction(DRAGONEventAction* eventAction,DRAGONRunAction* runAction)
: fEventAction(eventAction),frunAction(runAction),max_step(10000),len_max(20000.),ipart(80)
{
 fSteppingActionMessenger = new DRAGONSteppingActionMessenger(this);
 
 uvinit();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DRAGONSteppingAction::uvinit()
{
 ntargexit = 0;
 nbeamout = 0;
}


void DRAGONSteppingAction::UserSteppingAction(const G4Step* step)
{


 //frunAction->GetHistoManager()->ScanHitos();
  
  G4Track* track = step->GetTrack();
  const G4ParticleDefinition* particle = track->GetParticleDefinition();
 
  //if(particle->GetParticleName() == "C12" && particle->GetPDGCharge() == 6 &&  track->GetDynamicParticle()->GetCharge() == 6)
  //   std::cout << "SE FORMO" << std::endl;
 
  std::cout << particle->GetParticleName() << std::endl;
  std::cout << particle->GetPDGCharge() << std::endl;
  std::cout << track->GetDynamicParticle()->GetCharge() << std::endl;
   
std::cout << "pppppppppppppppppppppppppppppppppppppp" << std::endl;
G4cout << "irecoil: " << irecoil << G4endl;
G4cout << "alpha: " << alpha << G4endl;
G4cout << "n_detmate: " << n_detmate << G4endl;
G4cout << "lkine: " << lkine << G4endl;
G4cout << "itckov: " << itckov << G4endl;
G4cout << "targtype: " << targtype << G4endl;
G4cout << "eres: " << frunAction->GetPrimGener()->Geteres() << G4endl;
G4cout << "resmass: " << resmass << G4endl;
G4cout << "recoilmom: " << recoilmom << G4endl;
G4cout << "resenerg: " << resenerg << G4endl;
G4cout << "erescm: " << frunAction->GetPrimGener()->Geterescm() << G4endl;
G4cout << "prodm: " << prodm << G4endl;
G4cout << "beammass: " << beammass << G4endl;
G4cout << "e0recoil: " << e0recoil << G4endl;
std::cout << "pppppppppppppppppppppppppppppppppppppp" << std::endl;
 
  gekin = track->GetKineticEnergy()/GeV;
  destep = step->GetTotalEnergyDeposit();
  ngkine = step->GetNumberOfSecondariesInCurrentStep();
  tofg = track->GetGlobalTime()/s;
  sleng = track->GetTrackLength();
  number = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber() + 1;  //OJO Implica que k empieza en 1 -> k - 1
  
  G4LogicalVolume* vol_log = step->GetPreStepPoint()->GetTouchable()->GetVolume()->GetLogicalVolume();
  if(vol_log->GetSensitiveDetector()) 
    { 
     if(vol_log->GetName() == "HSNG")idtype = 1;
     if(vol_log->GetName() == "PMT")idtype = 2;
	 if(vol_log->GetName() == "ENDV")idtype = 3;
     } 
  
  if(step->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetName())
    
  if(track->GetCurrentStepNumber() == 1) {ntmult++;}
  
  if (particle->GetParticleType() == "gamma") itrtyp = 1;
  if (particle->GetParticleType() == "lepton") itrtyp = 2;
  if (particle->GetParticleName() == "neutron") itrtyp = 3;
  if (particle->GetParticleType() == "hadron") itrtyp = 4;
  if (particle->GetParticleType() == "muon") itrtyp = 5;
  if (particle->GetParticleType() == "geantino") itrtyp = 6;
  if (particle->GetParticleType() == "opticalphoton") itrtyp = 7;
  if (particle->GetParticleType() == "nucleus") itrtyp = 8;
  
  G4StepStatus preStatus  = step->GetPreStepPoint()->GetStepStatus();
  G4StepStatus postStatus = step->GetPostStepPoint()->GetStepStatus();

  G4VPhysicalVolume* preVolume = step->GetPreStepPoint()->GetPhysicalVolume();
  G4VPhysicalVolume* postVolume = step->GetPostStepPoint()->GetPhysicalVolume();
  if (preVolume != nullptr && postVolume != nullptr && preVolume == postVolume) inwvol = 0;    
  if (preStatus == fGeomBoundary) inwvol = 1;
  if (postStatus == fGeomBoundary) inwvol = 2;
  if (postStatus == fWorldBoundary ) inwvol = 3;

  if (track->GetCreatorProcess())
      kcase = track->GetCreatorProcess()->GetProcessName();
 
  G4ThreeVector pos = step->GetPostStepPoint()->GetPosition();
  G4double pmag = step->GetPostStepPoint()->GetMomentum().mag();
  G4ThreeVector pdir = step->GetPostStepPoint()->GetMomentumDirection().unit();
  G4ThreeVector pVec = step->GetPostStepPoint()->GetMomentum();

  vect[0] = pos.x()/cm;
  vect[1] = pos.y()/cm;
  vect[2] = pos.z()/cm;
  vect[3] = pdir.x();
  vect[4] = pdir.y();
  vect[5] = pdir.z();
  vect[6] = pmag/GeV; 
 
///////////////////////////////////////////////////////////////////////////////

// Posición actual (PostStepPoint)
G4cout << "Empezando el step" << G4endl;
G4ThreeVector pos1 = track->GetPosition();
G4ThreeVector mom = track->GetMomentum();
G4double E_tot = track->GetTotalEnergy(); 
G4double p_mag = mom.mag(); 

G4cout << "Track Position: x=" << pos1.x()/cm 
       << " cm, y=" << pos1.y()/cm 
       << " cm, z=" << pos1.z()/cm << " cm" << G4endl;

G4cout << "Track Momentum: px=" << mom.x()/GeV 
       << " GeV/c, py=" << mom.y()/GeV 
       << " GeV/c, pz=" << mom.z()/GeV << " GeV/c" 
       << G4endl;
       
G4cout << "Momentum magnitude: " << p_mag/GeV << " GeV/c" << G4endl;

G4cout << "Kinetic Energy: " << track->GetKineticEnergy()/GeV << " GeV" << G4endl;

G4cout << "Total Energy: " << E_tot/GeV << " GeV" << G4endl;

std::cout << "x: " << vect[0] << " cm" << std::endl;
std::cout << "y: " << vect[1] << " cm" << std::endl;
std::cout << "z: " << vect[2] << " cm" << std::endl;
std::cout << "I_px: " << vect[3] << std::endl;
std::cout << "I_py: " << vect[4] << std::endl;
std::cout << "I_pz: " << vect[5] << std::endl;
std::cout << "px: " << vect[3]*vect[6] << " GeV/c" << std::endl;
std::cout << "py: " << vect[4]*vect[6] << " GeV/c" << std::endl;
std::cout << "pz: " << vect[5]*vect[6] << " GeV/c" << std::endl;
std::cout << "p: " << vect[6] << " GeV/c" << std::endl;
///////////////////////////////////////////////////////////////////////////////

  auto touch = step->GetPreStepPoint()->GetTouchable();
  nlevel = touch->GetHistoryDepth();
  std::cout << "nlevel = " << nlevel << std::endl;
  for (int i = 0; i <= nlevel; ++i) 
      {
       std::cout << " Level " << i << ": " << touch->GetVolume(i)->GetName() << std::endl;
       }
  names = step->GetPreStepPoint()->GetTouchable()->GetVolume()->GetName();


 gustep(step);
 
 /*
 std::cout << "TESTING..." << std::endl;
 G4double xm[3], xd[3], xd_endv[3], xdd_endv[3];
 G4String chname = "E1";

  xm[0] = -327.373505;
  xm[1] = -0.573720872;
  xm[2] = 577.177063;
 
  xd[0] = -0.647119582;
  xd[1] = -0.0416779518;
  xd[2] = 44.1201553;
 
 gmtod(xm,xd,1,chname);
 
 std::cout << "xd " << xd[0] << ", " << xd[1] << ", " << xd[2] << std::endl;
 
 
 xd = 
 irot = irot_dipole[k - 1];
 
 gitran(xd, &dx_dipole[0][k - 1], irot, xd);  
 
 
std::cout << "Imprimiendo arreglos: " << std::endl;
std::cout << "Caso 1: " << std::endl;
 
// ----------------- 1D -----------------
std::cout << "irot_dipole: ";
for (int i = 0; i < max_dipole; i++) {
    std::cout << irot_dipole[i] << ", ";
}


std::cout << std::endl << "irot_edipol: ";
for (int i = 0; i < max_edipol; i++) {
    std::cout << irot_edipol[i] << ", ";
}
// ----------------- 2D -----------------
// Dipole
std::cout << std::endl << "jcol_dipole:" << std::endl;
for (int j = 0; j < 2; j++) {
    for (int i = 0; i < max_dipole; i++) {
        std::cout << jcol_dipole[j][i] << " ";
    }
    std::cout << std::endl;
}

std::cout << std::endl << "xcol_dipole:" << std::endl;
for (int j = 0; j < 2; j++) {
    for (int i = 0; i < max_dipole; i++) {
        std::cout << xcol_dipole[j][i] << " ";
    }
    std::cout << std::endl;
}

std::cout << std::endl << "ycol_dipole:" << std::endl;
for (int j = 0; j < 2; j++) {
    for (int i = 0; i < max_dipole; i++) {
        std::cout << ycol_dipole[j][i] << " ";
    }
    std::cout << std::endl;
}

std::cout << std::endl << "dxcol_dipole:" << std::endl;
for (int j = 0; j < 2; j++) {
    for (int i = 0; i < max_dipole; i++) {
        std::cout << dxcol_dipole[j][i] << " ";
    }
    std::cout << std::endl;
}

std::cout << std::endl << "dycol_dipole:" << std::endl;
for (int j = 0; j < 2; j++) {
    for (int i = 0; i < max_dipole; i++) {
        std::cout << dycol_dipole[j][i] << " ";
    }
    std::cout << std::endl;
}

std::cout << std::endl << "dx_dipole:" << std::endl;
for (int j = 0; j < 3; j++) {
    for (int i = 0; i < max_dipole; i++) {
        std::cout << dx_dipole[j][i] << " ";
    }
    std::cout << std::endl;
}

// Multipole
std::cout << std::endl << "jcol_mpole:" << std::endl;
for (int j = 0; j < 2; j++) {
    for (int i = 0; i < max_mpole; i++) {
        std::cout << jcol_mpole[j][i] << " ";
    }
    std::cout << std::endl;
}

std::cout << std::endl << "xcol_mpole:" << std::endl;
for (int j = 0; j < 2; j++) {
    for (int i = 0; i < max_mpole; i++) {
        std::cout << xcol_mpole[j][i] << " ";
    }
    std::cout << std::endl;
}

std::cout << std::endl << "ycol_mpole:" << std::endl;
for (int j = 0; j < 2; j++) {
    for (int i = 0; i < max_mpole; i++) {
        std::cout << ycol_mpole[j][i] << " ";
    }
    std::cout << std::endl;
}

std::cout << std::endl << "dxcol_mpole:" << std::endl;
for (int j = 0; j < 2; j++) {
    for (int i = 0; i < max_mpole; i++) {
        std::cout << dxcol_mpole[j][i] << " ";
    }
    std::cout << std::endl;
}

std::cout << std::endl << "dycol_mpole:" << std::endl;
for (int j = 0; j < 2; j++) {
    for (int i = 0; i < max_mpole; i++) {
        std::cout << dycol_mpole[j][i] << " ";
    }
    std::cout << std::endl;
}

std::cout << std::endl << "dx_mpole:" << std::endl;
for (int j = 0; j < 3; j++) {
    for (int i = 0; i < max_mpole; i++) {
        std::cout << dx_mpole[j][i] << " ";
    }
    std::cout << std::endl;
}

// Sole
std::cout << std::endl << "jcol_sole:" << std::endl;
for (int j = 0; j < 2; j++) {
    for (int i = 0; i < max_sole; i++) {
        std::cout << jcol_sole[j][i] << " ";
    }
    std::cout << std::endl;
}

std::cout << std::endl << "xcol_sole:" << std::endl;
for (int j = 0; j < 2; j++) {
    for (int i = 0; i < max_sole; i++) {
        std::cout << xcol_sole[j][i] << " ";
    }
    std::cout << std::endl;
}

std::cout << std::endl << "ycol_sole:" << std::endl;
for (int j = 0; j < 2; j++) {
    for (int i = 0; i < max_sole; i++) {
        std::cout << ycol_sole[j][i] << " ";
    }
    std::cout << std::endl;
}

std::cout << std::endl << "dxcol_sole:" << std::endl;
for (int j = 0; j < 2; j++) {
    for (int i = 0; i < max_sole; i++) {
        std::cout << dxcol_sole[j][i] << " ";
    }
    std::cout << std::endl;
}

std::cout << std::endl << "dycol_sole:" << std::endl;
for (int j = 0; j < 2; j++) {
    for (int i = 0; i < max_sole; i++) {
        std::cout << dycol_sole[j][i] << " ";
    }
    std::cout << std::endl;
}

std::cout << std::endl << "dx_sole:" << std::endl;
for (int j = 0; j < 3; j++) {
    for (int i = 0; i < max_sole; i++) {
        std::cout << dx_sole[j][i] << " ";
    }
    std::cout << std::endl;
}

// EDipole
std::cout << std::endl << "jcol_edipol:" << std::endl;
for (int j = 0; j < 2; j++) {
    for (int i = 0; i < max_edipol; i++) {
        std::cout << jcol_edipol[j][i] << " ";
    }
    std::cout << std::endl;
}

std::cout << std::endl << "xcol_edipol:" << std::endl;
for (int j = 0; j < 2; j++) {
    for (int i = 0; i < max_edipol; i++) {
        std::cout << xcol_edipol[j][i] << " ";
    }
    std::cout << std::endl;
}

std::cout << std::endl << "ycol_edipol:" << std::endl;
for (int j = 0; j < 2; j++) {
    for (int i = 0; i < max_edipol; i++) {
        std::cout << ycol_edipol[j][i] << " ";
    }
    std::cout << std::endl;
}

std::cout << std::endl << "dxcol_edipol:" << std::endl;
for (int j = 0; j < 2; j++) {
    for (int i = 0; i < max_edipol; i++) {
        std::cout << dxcol_edipol[j][i] << " ";
    }
    std::cout << std::endl;
}

std::cout << std::endl << "dycol_edipol:" << std::endl;
for (int j = 0; j < 2; j++) {
    for (int i = 0; i < max_edipol; i++) {
        std::cout << dycol_edipol[j][i] << " ";
    }
    std::cout << std::endl;
}

std::cout << std::endl << "dx_edipol:" << std::endl;
for (int j = 0; j < 3; j++) {
    for (int i = 0; i < max_edipol; i++) {
        std::cout << dx_edipol[j][i] << " ";
    }
    std::cout << std::endl;
}
*/
 
 std::cout << std::endl << "oooooooooooEND STEPooooooooooo" << std::endl; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}

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
/// \file DRAGON/src/DRAGONActionInitialization.cc
/// \brief Implementation of the DRAGON::DRAGONActionInitialization class

#include "DRAGONActionInitialization.hh"
#include "DRAGONPrimaryGeneratorAction.hh"
#include "DRAGONPhysicsList.hh"
#include "DRAGONDetectorConstruction.hh"
#include "DRAGONRunAction.hh"
#include "DRAGONEventAction.hh"
#include "DRAGONTrackingAction.hh"
#include "DRAGONSteppingAction.hh"
#include "DRAGONStackingAction.hh"

namespace DRAGON
{
	
DRAGONActionInitialization::DRAGONActionInitialization(DRAGONDetectorConstruction* detector, DRAGONPhysicsList* physics) 
:fDetector(detector),fphys(physics)
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
void DRAGONActionInitialization::BuildForMaster() const
{
    auto runAction = new DRAGONRunAction(fDetector, nullptr, fphys, nullptr);
    SetUserAction(runAction);
}

void DRAGONActionInitialization::Build() const
{
    auto fprimary = new DRAGONPrimaryGeneratorAction();
    fprimary->SetDetector(fDetector);
    fprimary->SetPhys(fphys);
    SetUserAction(fprimary);

    auto runAction = new DRAGONRunAction(fDetector, fprimary, fphys, nullptr);
    SetUserAction(runAction);

    auto eventAction = new DRAGONEventAction(runAction);
    SetUserAction(eventAction);

    auto steppingAction = new DRAGONSteppingAction(eventAction, runAction);
    runAction->SetSteppingAction(steppingAction);
    SetUserAction(steppingAction);

    SetUserAction(new DRAGONTrackingAction(eventAction));
    SetUserAction(new DRAGONStackingAction());
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}

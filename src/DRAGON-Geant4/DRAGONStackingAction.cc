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
/// \file DRAGON/src/DRAGONStackingAction.cc
/// \brief Implementation of the DRAGON::DRAGONStackingAction class

#include "G4ClassificationOfNewTrack.hh"      //Geant4
#include "G4Track.hh"
#include "G4VProcess.hh"

#include "DRAGONStackingAction.hh"            //Local


namespace DRAGON
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack DRAGONStackingAction::ClassifyNewTrack(const G4Track* track){

//From gustep_gbox.f
//C.
//C *** Daughter particles that were generated in the current step
//C ***                  are put on the stack
//C.

 G4int pdg = track->GetDefinition()->GetPDGEncoding();

        if (pdg == 12 || pdg == 14 || pdg == 16 || pdg == -12 || pdg == -14 || pdg == -16) 
           {
            return fKill;   
        }

        /*
        auto creatorProc = track->GetCreatorProcess();
        if (creatorProc) {
            std::string procName = creatorProc->GetProcessName();

            if (procName == "Decay") {
                return fKill;  
            }
        }

//From gustep_gbox.f
//C.
//C ***   The charged particle has produced scintillation photons and it
//C ***   is still alive; we put it on the stack and we let the photons be
//C ***   tracked first
//C.

    if (name == "opticalphoton") {
        return fUrgent;
    }

    return fWaiting;
*/


const G4VProcess* creatorProc = track->GetCreatorProcess();
if (creatorProc)
   {
    G4String kcase = creatorProc->GetProcessName();
    if (kcase == "Decay" || kcase == "RESR")
       {
        std::cout << "AAAA9" << std::endl;
        return fUrgent;
        }
}





return fUrgent;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DRAGONStackingAction::NewStage()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DRAGONStackingAction::PrepareNewEvent()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
}


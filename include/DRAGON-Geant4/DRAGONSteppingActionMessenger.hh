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
/// \file DRAGONSteppingActionMessenger.hh
/// \brief Definition of the DRAGONSteppingActionMessenger class
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DRAGONSteppingActionMessenger_h
#define DRAGONSteppingActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class G4UIdirectory;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace DRAGON {

class DRAGONSteppingAction;

/// Stepping Action messenger

class DRAGONSteppingActionMessenger : public G4UImessenger
{
	public:
		DRAGONSteppingActionMessenger(DRAGONSteppingAction*);
		~DRAGONSteppingActionMessenger() override;

		void SetNewValue(G4UIcommand*, G4String) override;

	private:
		DRAGONSteppingAction* fSteppingAction = nullptr;
		G4UIdirectory* fSteppingActionDir = nullptr;
		G4UIcmdWithAString* fMAXSCmd = nullptr;
};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
}

#endif


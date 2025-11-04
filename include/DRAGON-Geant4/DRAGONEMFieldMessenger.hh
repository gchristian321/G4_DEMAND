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
/// \file DRAGONEMFieldMessenger.hh
/// \brief Definition of the DRAGONEMFieldMessenger class
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DRAGONEMFieldMessenger_h
#define DRAGONEMFieldMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace DRAGON {

class DRAGONEMField;

/// EMField messenger

class DRAGONEMFieldMessenger : public G4UImessenger
{
	public:
		DRAGONEMFieldMessenger(DRAGONEMField*);
		~DRAGONEMFieldMessenger() override;

		void SetNewValue(G4UIcommand*, G4String) override;

	private:
		DRAGONEMField* fEMField = nullptr;
		G4UIdirectory* fEMFieldDir = nullptr;
		G4UIcmdWithAString* fSCALCmd = nullptr;
                G4UIcmdWithAnInteger* fStepperCMD = nullptr;
                G4UIcmdWithADouble* fMinStepCMD = nullptr;
                G4UIcmdWithADouble* fDeltaChordCMD = nullptr;
                G4UIcmdWithADouble* fDeltaOneStepCMD = nullptr;
                G4UIcmdWithADouble* fDeltaIntersectionCMD = nullptr;
                G4UIcmdWithADouble* fEpsMinCMD = nullptr;
                G4UIcmdWithADouble* fEpsMaxCMD = nullptr;
};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
}

#endif


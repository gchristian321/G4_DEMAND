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
/// \file DRAGON/include/DRAGONDetectorMessenger.hh
/// \brief Definition of the DRAGON::DRAGONDetectorMessenger class

#ifndef DRAGONDetectorMessenger_h
#define DRAGONDetectorMessenger_h 1

#include "G4UImessenger.hh"              //Geant4

class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcmdWithADouble;
class G4UIcmdWithABool;

namespace DRAGON
{

class DRAGONDetectorConstruction;

/// Messenger class that defines commands for DetectorConstruction.
///

class DRAGONDetectorMessenger: public G4UImessenger
{
  public:
    DRAGONDetectorMessenger(DRAGONDetectorConstruction* );
    ~DRAGONDetectorMessenger() override;

    void SetNewValue(G4UIcommand*, G4String) override;
    
  private:
    DRAGONDetectorConstruction* fDetectorConstruction = nullptr;

    G4UIdirectory* fDirectory = nullptr;
    G4UIdirectory* fDetDirectory = nullptr;    
    G4UIcmdWithAString* fMAXSCmd  = nullptr;
    G4UIcmdWithAnInteger* fTUBECmd = nullptr;
    G4UIcmdWithAnInteger* fTARGCmd = nullptr;
    G4UIcmdWithAnInteger* fMPMTCmd = nullptr;
    G4UIcmdWithAString* fPMTRCmd = nullptr;
    G4UIcmdWithAString* fWALLCmd = nullptr;
    G4UIcmdWithADouble* fBGAPCmd = nullptr;
    G4UIcmdWithAString* fFSIDCmd = nullptr;
    G4UIcmdWithAnInteger* fNDMATCmd = nullptr;
    G4UIcmdWithADouble* fHOLECmd = nullptr;
    G4UIcmdWithADouble* fBLKACmd = nullptr;
    G4UIcmdWithADouble* fREFLCmd = nullptr;
    G4UIcmdWithABool* fOVRLCmd = nullptr;
    G4UIcmdWithAString* fMASKCmd  = nullptr;
    G4UIcmdWithAString* fSHLDCmd = nullptr;

};

}

#endif

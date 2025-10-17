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
// $Id$
// 
/// \file DemandEventAction.hh
/// \brief Definition of the DemandEventAction class

#ifndef DemandEventAction_h
#define DemandEventAction_h 1

#include <vector>

#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"

#include "DemandHit.hh"

#include "globals.hh"

/// Event action class
///
/// In EndOfEventAction(), it prints the accumulated quantities of the energy 
/// deposit and track lengths of charged particles in Absober and Gap layers 
/// stored in the hits collections.

class DemandEventAction : public G4UserEventAction
{
public:
  DemandEventAction();
  virtual ~DemandEventAction();

  virtual void  BeginOfEventAction(const G4Event* event);
  virtual void    EndOfEventAction(const G4Event* event);
	
private:
  // methods
  DemandHitsCollection* GetHitsCollection(
		G4int hcID,	const G4Event* event) const;
  void PrintEventStatistics(
		G4double absoEdep, G4double caloTime, const G4ThreeVector& caloPos) const;
	G4ThreeVector FigureOutMeasuredPosition(const DemandHit& hit) const;

	
  // data members                   
	std::vector<G4int> fDemandHCID;
};
                     
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    

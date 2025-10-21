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
/// \file DemandRunAction.hh
/// \brief Definition of the DemandRunAction class

#ifndef DemandRunAction_h
#define DemandRunAction_h 1

#include <vector>
#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;
class DemandRunMessenger;

/// Run action class
///
/// It accumulates statistic and computes dispersion of the energy deposit 
/// and track lengths of charged particles with use of analysis tools:
/// H1D histograms are created in BeginOfRunAction() for the following 
/// physics quantities:
/// - Edep in absorber
/// - Edep in gap
/// - Track length in absorber
/// - Track length in gap
/// The same values are also saved in the ntuple.
/// The histograms and ntuple are saved in the output file in a format
/// accoring to a selected technology in DemandAnalysis.hh.
///
/// In EndOfRunAction(), the accumulated statistic and computed 
/// dispersion is printed.
///

#include <G4ThreeVector.hh>

class TFile;
class TTree;

struct G3Event
{
	G4ThreeVector fPosition;
	G4ThreeVector fMomentumDirection;
	G4double      fKineticEnergy;
};

class DemandRunAction : public G4UserRunAction
{
  public:
    DemandRunAction();
    virtual ~DemandRunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);
	
	  G4long SetupGeant3Input(const G4String&);
	  void GetGeant3Event(G4long indx, G3Event*) const;

  private:
	  void CleanupGeant3Input();
	
  private:
	  G4int fNumRuns;
	  TFile* fG3File;
  	TTree* fG3Tree;
	  G4float E_n, cost_n, cosp_n, sinp_n, xint, yint, zint;
	  G4float E_rec, cost_r, cosp_r, sinp_r;
	  G4int react, recdet;
	  std::vector<G4long> fEventIndices;

	  DemandRunMessenger* fRunMessenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


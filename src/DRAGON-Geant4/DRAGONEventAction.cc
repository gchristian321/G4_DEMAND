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
/// \file DRAGON/src/DRAGONEventAction.cc
/// \brief Implementation of the DRAGON::DRAGONEventAction class

#include "DRAGONEventAction.hh"
#include "DRAGONRunAction.hh"
#include "G4DigiManager.hh"

#include "G4Event.hh"


namespace DRAGON
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DRAGONEventAction::DRAGONEventAction(DRAGONRunAction* runAction)
: fRunAction(runAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DRAGONEventAction::BeginOfEventAction(const G4Event*)
{
//C.    Initialize all ntuple variables
      E_int = 0.;
      E_rec = 0.;
      cost_r = 99.;
      cosp_r = 99.;
      react = 0;
      Nodec = 0;
      recdet = 0;
      for (int i = 0; i < 15; i++)
          {
           E_g[i] = 0.;
           E_gp[i] = 0.;
           cost_g[i] = 99.;
           cost_gp[i] = 99.;
           phi_g[i] = 99.;
           }
      x_r = 0.;
      y_r = 0.;
      z_r = 0.;
      thet_r = 0.;
      xstop = 0.;
      ystop = 0.;
      zstop = 0.;
      xint = 0.;
      yint = 0.;
      zint = 0.;
      x = 99.;
      y = 99.;
      xp = 99.;
      yp = 99.;
      beamtof = 0.;
      for (int i = 0; i < 10; i++)
          {
           xtest[i] = 99.;
           ytest[i] = 99.;
           etest[i] = 0.;
           }
      dsssdpos = 0.;
      recoil_hit_ENDV = 0;
      num_bgos_hit = 0;
      e_bgos_total = 0;
      num_bgo_first = 0;
      e_bgo_first = 0;
      num_bgo_second = 0;
      e_bgo_second = 0;
      num_bgos_hit_ab = 0;
      num_bgo_first_ab = 0;
      e_bgo_first_ab = 0;
      num_bgo_second_ab = 0;
      e_bgo_second_ab = 0;
      pair_productions = 0;
      gammatof = 0.;
      e0_conv = 0.;
      McpHit = false;

      gutrev();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DRAGONEventAction::EndOfEventAction(const G4Event* event)
{
 G4DigiManager* DM = G4DigiManager::GetDMpointer();
 DM->Digitize("PMT"); 
 DM->Digitize("SCNT"); 
 DM->Digitize("DSSD"); 
 
 guout(event);
}

void DRAGONEventAction::uvinit()
{ 
 react = 0;
 idevt = 0;
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}

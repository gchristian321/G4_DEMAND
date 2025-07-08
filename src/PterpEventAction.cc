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
/// \file PterpEventAction.cc
/// \brief Implementation of the PterpEventAction class

#include "PterpEventAction.hh"
#include "PterpSD.hh"
#include "PterpHit.hh"
#include "PterpAnalysis.hh"
#include "PterpDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace {
const double POSITION_FWHM = 3 * CLHEP::cm;
}

PterpEventAction::PterpEventAction()
 : G4UserEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PterpEventAction::~PterpEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PterpHitsCollection* 
PterpEventAction::GetHitsCollection(G4int hcID,
                                  const G4Event* event) const
{
  auto hitsCollection 
    = static_cast<PterpHitsCollection*>(
        event->GetHCofThisEvent()->GetHC(hcID));
  
  if ( ! hitsCollection ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID; 
    G4Exception("PterpEventAction::GetHitsCollection()",
      "MyCode0003", FatalException, msg);
  }         

  return hitsCollection;
}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PterpEventAction::PrintEventStatistics(
	G4double caloEdep, G4double caloTime, const G4ThreeVector& caloPos)
	const
{
  // print event statistics
  G4cout
     << "   Calorimeter: total energy: " 
     << std::setw(7) << G4BestUnit(caloEdep, "Energy")
		 << "       time: "
		 << std::setw(7) << G4BestUnit(caloTime, "Time")
     << "       position: " 
     << std::setw(7) << G4BestUnit(caloPos.x(), "Position (x), ")
     << std::setw(7) << G4BestUnit(caloPos.y(), "Position (y), ")
     << std::setw(7) << G4BestUnit(caloPos.z(), "Position (z)")
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PterpEventAction::BeginOfEventAction(const G4Event* /*event*/)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PterpEventAction::EndOfEventAction(const G4Event* event)
{  
  // Get hits collections IDs (only once)
	if(fPterpHCID.empty()) {
		G4SDManager* sdManager = G4SDManager::GetSDMpointer();
		G4int hcID = sdManager->GetCollectionID("PterpSensitive/pterpCollection");
		fPterpHCID.push_back(hcID);
	}
	
  // Get hits collections
	G4HCofThisEvent* hce = event->GetHCofThisEvent();
	
	std::vector<PterpHitsCollection*> hPterpHC;
	for(G4int hcid : fPterpHCID) {
		hPterpHC.push_back(static_cast<PterpHitsCollection*>(hce->GetHC(hcid)));
	}

	// Get analysis manager
	auto analysisManager = PterpAnalysis::Instance();
	analysisManager->Clear();

	// Loop hits
	G4double tmin = DBL_MAX;
	PterpHit* firstRealHit = 0;
	for(auto& hitCollection : hPterpHC) {
		for(int j=0; j< hitCollection->entries(); j++) {
			auto hit = (*hitCollection)[j];
			auto sensitiveDetector = dynamic_cast<PterpSD&>(
				*(hit->GetPhysicalVolume()->GetLogicalVolume()->GetSensitiveDetector()));
			G4double thresh = 
				sensitiveDetector.GetThreshold(hit->GetID());

			if(hit->GetEnergyQuenched() > thresh) {
				auto pos = FigureOutMeasuredPosition(*hit);
				analysisManager->AddHit(
					hit->GetEnergyQuenched(),
					hit->GetTime(),
					pos.x(),
					pos.y(),
					pos.z(),
					hit->GetParticleA(),
					hit->GetParticleZ() );
			}
			if(hit->GetTime() < tmin) {
				tmin = hit->GetTime();
				firstRealHit = hit;
			}
		}
	}

	// analyze
	if(firstRealHit) {
		analysisManager->SetFirstInteraction(
			firstRealHit->GetTime(),
			firstRealHit->GetActualPosition().x(),
			firstRealHit->GetActualPosition().y(),
			firstRealHit->GetActualPosition().z());
	}
	analysisManager->Analyze();
	analysisManager->FillGenTree();
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector PterpEventAction::FigureOutMeasuredPosition(const PterpHit& hit) const
{
	auto bar = dynamic_cast<PterpDetectorConstruction::ScintillatorBox*>(
		hit.GetPhysicalVolume()->GetLogicalVolume()->GetSolid());
	if(!bar) {
		throw std::runtime_error(
			"PterpEventAction::FigureOutMeasuredPosition --> No conversion from hit to G4Box");
	}		 

	G4ThreeVector measuredPos = hit.GetPosition(); // centre of detector

	switch(bar->GetReadoutType()) {
	case PterpDetectorConstruction::kCube:
		{
			//keep at centre
			break;
		}
	case PterpDetectorConstruction::kBar:
		{
			const std::array<double, 3> side_lengths = {
				bar->GetXHalfLength() * 2,
				bar->GetYHalfLength() * 2,
				bar->GetZHalfLength() * 2
			};
			auto longSide = std::max_element(side_lengths.begin(),side_lengths.end()) - side_lengths.begin();
			measuredPos[longSide] = hit.GetActualPosition()[longSide] +
				G4RandGauss::shoot(0, POSITION_FWHM / (2*sqrt(2*log(2))));
		break;
		}
	default:
		{
			char buf[4096];
			sprintf(buf, "EventAction::FigureOutMeasuredPosition -->"
							" invalid readout type %i", bar->GetReadoutType());
			throw std::runtime_error(buf);
			break;
		}
	}
	return measuredPos;
}

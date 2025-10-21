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
/// \file DemandRunAction.cc
/// \brief Implementation of the DemandRunAction class

#include "DemandRunAction.hh"
#include "DemandAnalysis.hh"
#include "DemandPrimaryGeneratorAction.hh"
#include "DemandRunMessenger.hh"

#include <sstream>

#include <TFile.h>
#include <TTree.h>

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DemandRunAction::DemandRunAction()
	: G4UserRunAction(),
		fNumRuns(0),
		fG3File(nullptr),
		fG3Tree(nullptr)
{ 
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);     

  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in DemandAnalysis.hh
  auto analysisManager = DemandAnalysis::Instance();
	if(!analysisManager) {
		throw std::logic_error("Couldn't create analysis manager in RunAction!");
	}

	fRunMessenger = new DemandRunMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DemandRunAction::~DemandRunAction()
{
	CleanupGeant3Input();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DemandRunAction::BeginOfRunAction(const G4Run* /*run*/)
{ 
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  
  // Get analysis manager
  auto analysisManager = DemandAnalysis::Instance();

  // Open an output file
  //
	std::stringstream fileName;
	fileName << "Demand_" << fNumRuns++ << ".root";
  analysisManager->OpenFile(fileName.str());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DemandRunAction::EndOfRunAction(const G4Run* run)
{
  // print histogram statistics
  //
  auto analysisManager = DemandAnalysis::Instance();
	
  // save histograms & ntuple
  //
  analysisManager->Write();
 
	G4cout << "Number of events simulated: " <<
		run->GetNumberOfEvent() << " ( " << fEventIndices.size() << " )\n";

	G4cout << "Number of hits above threshold: " <<
		analysisManager->GetEventsAboveThreshold() << G4endl;
	
	G4cout << "Number of events with neutron crossing detector: " <<
		analysisManager->GetEventsCrossingDetector() << G4endl;
	
  analysisManager->CloseFile();
	CleanupGeant3Input();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DemandRunAction::CleanupGeant3Input()
{
	if(fG3File){
		if(fG3Tree){
			fG3Tree->ResetBranchAddresses();
		}
		fG3File->Close();
		delete fG3File;
		fG3File = nullptr;
		fG3Tree = nullptr;
	}
	fEventIndices.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace {
template<typename T> void set_branch_address(TTree* tree,const char* name,T& t)
{
	if( tree->SetBranchAddress(name, &t) != 0 ) {
		throw std::runtime_error(
			Form("ERROR in DemandRunAction::SetupGeant3Input -- "
					 "Bad h1000 tree, no %s branch!", name) );
	}
}}

G4long DemandRunAction::SetupGeant3Input(const G4String& g3fname)
{
	fG3File = TFile::Open(g3fname.c_str());
	if(!fG3File){
		throw std::runtime_error(
			Form("Bad input ROOT file: \"%s\" !",
					 g3fname.c_str()));
	}

	fG3Tree = dynamic_cast<TTree*>(fG3File->Get("h1000"));
	if(!fG3Tree){
		throw std::runtime_error(
			Form("Bad input ROOT file (no h1000 tree): \"%s\" !",
					 g3fname.c_str()));
	}

	set_branch_address(fG3Tree,"E_n",E_n);
	set_branch_address(fG3Tree,"cost_n",cost_n);
	set_branch_address(fG3Tree,"cosp_n",cosp_n);
	set_branch_address(fG3Tree,"sinp_n",sinp_n);

	set_branch_address(fG3Tree,"E_rec",E_rec);
	set_branch_address(fG3Tree,"cost_r",cost_r);
	set_branch_address(fG3Tree,"cosp_r",cosp_r);
	set_branch_address(fG3Tree,"sinp_r",sinp_r);

	set_branch_address(fG3Tree,"xint",xint);
	set_branch_address(fG3Tree,"yint",yint);
	set_branch_address(fG3Tree,"zint",zint);

	set_branch_address(fG3Tree,"react",react);
	set_branch_address(fG3Tree,"recdet",recdet);

	fEventIndices.clear();
	for(G4long entry = 0; entry < fG3Tree->GetEntries(); ++entry){
		fG3Tree->GetEntry ( entry );
		if ( react != 0 && recdet != 0 ) {
			fEventIndices.push_back(entry);
		}
	}

	G4cout << "Number of events in GEANT3 file: " << fEventIndices.size() << G4endl;
	return fEventIndices.size();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DemandRunAction::GetGeant3Event(G4long indx, G3Event* g3evt) const
{
	G4long eventNo;
	try {
		eventNo = fEventIndices.at ( indx );
	}
	catch (std::exception& e){
		G4cerr << "ERROR: bad event index " << indx << ", number of valid events: "
					 << fEventIndices.size() << G4endl;
		throw e;
	}
	
	fG3Tree->GetEntry ( eventNo );
	g3evt->fPosition.set(xint,yint,zint);

	G4double sint_n  = sqrt(1.0 - cost_n*cost_n);
	g3evt->fMomentumDirection.set(
		sint_n * cosp_n,
		sint_n * sinp_n,
		cost_n
		);
	g3evt->fKineticEnergy = E_n;

	const double m_rec = 23.2741642 * GeV; // hard coded - same as GEANT3
	const double p_rec = sqrt(pow(E_rec+m_rec, 2) - m_rec*m_rec);
	const double sint_r = sqrt(1.0 - cost_r*cost_r);
	g3evt->fRecoilLorentzVector.set(
		sint_r * cosp_r, sint_r * sinp_r,	cost_r,	E_rec + m_rec
		);
}


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
/// \file DRAGON/include/DRAGONRunAction.hh
/// \brief Definition of the DRAGON::DRAGONRunAction class

#ifndef DRAGONRunAction_h
#define DRAGONRunAction_h 1

#include "G4UserRunAction.hh"
#include "G4UnitsTable.hh" 
#include "globals.hh"

#include <chrono>

class G4Run;
class G4ParticleDefinition;

namespace DRAGON
{

class DRAGONDetectorConstruction;
class DRAGONRun;
class DRAGONPrimaryGeneratorAction;
class DRAGONHistoManager;
class DRAGONRunActionMessenger;
class DRAGONPhysicsList;
class DRAGONSteppingAction;

class DRAGONRunAction : public G4UserRunAction
{
  public:
    DRAGONRunAction(DRAGONDetectorConstruction* det, DRAGONPrimaryGeneratorAction* prim, DRAGONPhysicsList* phys, DRAGONSteppingAction* step);
	DRAGONRunAction();
	~DRAGONRunAction() override;

    G4Run* GenerateRun() override; 
    void BeginOfRunAction(const G4Run*) override;
    void EndOfRunAction(const G4Run*) override;
   
    void uvinit();
    void uginit();
    void uglast(const G4Run* run);
    G4int GetIrevs() {return irevs;}
    G4int* GetIswit() {return iswit;} 
    G4int Getn_detector(){return n_detector;}
    G4int Getngascell(){return ngascell;}
	G4int Getnend(){return nend;}
	G4int Getnfcm2(){return nfcm2;}
	G4int GetNum_Recoils_Q3(){return Num_Recoils_Q3;}
	G4int GetNum_Recoils_Q8(){return Num_Recoils_Q8;}
	G4int GetNum_BeamPart_ENDV(){return Num_BeamPart_ENDV;}
	G4int Getidevt(){return idevt;}
	G4int Getnreact(){return nreact;}
		
    DRAGONDetectorConstruction* GetDetConst() {return fDetector;}
    DRAGONPrimaryGeneratorAction* GetPrimGener() {return fPrimary;}
    DRAGONPhysicsList* GetPhysList() {return fphys;}
    DRAGONSteppingAction* GetStepAction() {return fstep;}
    
    void SetIrevs(G4int newVal) {irevs = newVal;}
    void SetIswit(const G4String& input);
    void SetSteppingAction(DRAGONSteppingAction* step){fstep = step;}
    void SetPrimaryAction(DRAGONPrimaryGeneratorAction* prim){fPrimary = prim;}
    void Setn_detector(G4int n_detector_){n_detector = n_detector_;}
    void Setngascell(G4int ngascell_){ngascell = ngascell_;}
	void Setnend(G4int nend_){nend = nend_;}
	void Setnfcm2(G4int nfcm2_){nfcm2 = nfcm2_;}
	void SetNum_Recoils_Q3(G4int Num_Recoils_Q3_){Num_Recoils_Q3 = Num_Recoils_Q3_;}
	void SetNum_Recoils_Q8(G4int Num_Recoils_Q8_){Num_Recoils_Q8 = Num_Recoils_Q8_;}
	void SetNum_BeamPart_ENDV(G4int Num_BeamPart_ENDV_){Num_BeamPart_ENDV = Num_BeamPart_ENDV_;}
	void Setidevt(G4int idevt_){idevt = idevt_;}
	void Setnreact_(G4int nreact_){nreact = nreact_;}
				        
    void PrintCustomParticles();
    void PrintIonProperties(G4ParticleDefinition* ion);
    void PrintMaterials();
    void PrintLogicalVolumes();
    void PrintSensitiveDetectors();
	
    DRAGONHistoManager*	GetHistoManager(){return fHistoManager;}

  private:
    friend class DRAGONDetectorConstruction;
    friend class DRAGONPrimaryGeneratorAction;
    friend class DRAGONHistoManager; 
	friend class DRAGONPhysicsList;
       
    class DRAGONDetectorConstruction* fDetector = nullptr;
    class DRAGONPrimaryGeneratorAction* fPrimary = nullptr;
    class DRAGONPhysicsList* fphys = nullptr;
    class DRAGONSteppingAction* fstep = nullptr;
    class DRAGONHistoManager* fHistoManager = nullptr;
    
    DRAGONRun* fRun = nullptr; 
    DRAGONRunActionMessenger* fRunMessenger = nullptr;
    
    G4int idebug;
    G4int iswit[10];
    G4int num = 9;
    G4int irevs;

    G4int ievent, idrun;
	G4int jstop;
                              
    std::chrono::time_point<std::chrono::high_resolution_clock> TIMINT;  //From gctime.inc
  
  public:
         G4int Num_Recoils_Q3, Num_Recoils_Q8, Num_BeamPart_ENDV, 
		       nend, nfcm2, idevt, n_detector;
	 	 G4int ngascell, nreact, ntargexit, nbeamout;  
	 	 

};

}

#endif


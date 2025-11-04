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
/// \file DRAGON/include/DRAGONPrimaryGeneratorAction.hh
/// \brief Definition of the DRAGON::DRAGONPrimaryGeneratorAction class

#ifndef DRAGONPrimaryGeneratorAction_h
#define DRAGONPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"   //Geant4
#include "G4SystemOfUnits.hh"
#include "globals.hh"

#include "DRAGONRunAction.hh"                 //local

#include <vector>                             //std

class G4GeneralParticleSource;
class G4Event;

namespace DRAGON
{

class DRAGONDetectorConstruction;
class DRAGONPrimaryGeneratorActionMessenger;
class DRAGONPhysicsList;

class DRAGONPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    DRAGONPrimaryGeneratorAction();
    ~DRAGONPrimaryGeneratorAction() override;
	
	const G4GeneralParticleSource* GetParticleGun() const {return fGPSParticleGun;}
    void GeneratePrimaries(G4Event*) override;
    
    void uvinit();
	void uhinit();
	void Build250();  
    void Build501();
    void Build523(G4String fFileName);
    G4double HRNDM1(G4int IDD);
    G4double sig(G4double e);
    G4double angdist(G4double angle);
          
    void gukine(G4Event* anEvent, G4GeneralParticleSource* fGPSParticleGun);
    void gukine_full(G4Event* anEvent, G4GeneralParticleSource* fGPSParticleGun);
    void gukine_full_up(G4Event* anEvent, G4GeneralParticleSource* fGPSParticleGun);
    void gukine_full_down(G4Event* anEvent, G4GeneralParticleSource* fGPSParticleGun);
    void gukine_full_hole(G4Event* anEvent, G4GeneralParticleSource* fGPSParticleGun);
    void gukine_full_left(G4Event* anEvent, G4GeneralParticleSource* fGPSParticleGun);
    void gukine_full_right(G4Event* anEvent, G4GeneralParticleSource* fGPSParticleGun);
    void gukine_gbox(G4Event* anEvent, G4GeneralParticleSource* fGPSParticleGun);
    void gukine_mitray(G4Event* anEvent, G4GeneralParticleSource* fGPSParticleGun);

    const G4double Getlm1() const {return lm1;}
    const G4double Getlm2() const {return lm2;}
    const G4double Getbeamo() const {return beamo;}
    const G4double Getbeamenerg() const {return beamenerg;}
	const G4double Geterescm() const {return erescm;}
    const G4double GetEmax() const {return emax;}
    const G4double* GetPkine() const {return pkine;}
    const G4double* GetGkine() const {return gkine;}
    const G4double Geteres() const {return eres;}
    //const G4double GetRefEnerg() const {return refenerg;}
    //const G4double GetRefAtNo() const {return refatno;}
    //const G4double GetRefQ() const {return refq;}
    const G4double* GetOffset() const {return offset;}
	DRAGONRunAction* GetDRAGONRunAction(){return DRAGONRunAction_prim;}	
    
	void SetKINE(const G4String& input, G4int option);
    void SetEGAM(const G4String& input);
    void SetGKIN(const G4String& input);
    void SetRAYFILE(const G4String val) {rayfile = val;}
    //void SetRefEnerg(const G4double newVal) {refenerg = newVal;}
    //void SetRefAtNo(const G4double newVal) {refatno = newVal;}
    //void SetRefQ(const G4double newVal) {refq = newVal;}
    void SetDetector(DRAGONDetectorConstruction* det){fDetector = det;}
    void SetPhys(DRAGONPhysicsList* phys){fphys = phys;}
    void SetRunAction(DRAGONRunAction& RunAction_prim){DRAGONRunAction_prim = &RunAction_prim;}  
	
  private:
    DRAGONRunAction* DRAGONRunAction_prim = nullptr;
    DRAGONDetectorConstruction* fDetector = nullptr; 
    DRAGONPhysicsList* fphys = nullptr;
    G4GeneralParticleSource* fGPSParticleGun = nullptr; 
    DRAGONPrimaryGeneratorActionMessenger* fGPSMessenger = nullptr; 
      
    G4int ikine;         
    G4double pkine[10],  
             mkine,         
             gkine[10];  

    G4double beamenerg,    //OJO Estas variables vinenen todas de DRAGONPhysicsList (no deben ser propiedad de esta clase)
             beammass,
			 beamo,
             sigx,
             amax,
             sigy,
             bmax,
             bunchl,
             emax,
             buncht,
             tofg,         
         	 eres0, 
	         resenerg, 
	         reswidth,
	         eres,
			 erescm,
			 zprod,  //OJO Estas variables vienen de DRAGONPhysicsList y se usan en el .cc (no deberian ser propiedad de esta clase)
             m1,     //Definirlas como locales en los metodos que las usan, sus valores vienen de DRAGONPhysicsList
			 m2, 
	     	 z1, 
			 z2, 
	         er, 
			 gp, 
			 gg, 
			 omg, 
			 ell;
    G4double lm1,
	         lm2,  
			 offset[4],  
	         egamma[10];
	
	G4String rayfile,    
             fFileName;	
    std::vector<G4double> energies,
                          counts,
                          cdf,     
                          sig_energy,
                          sig_values,   
                          angdist_costheta,
                          angdist_values;         
};  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
}
#endif

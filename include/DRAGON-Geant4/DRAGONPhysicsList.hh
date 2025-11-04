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
/// \file hadronic/Hadr01/include/PhysicsList.hh
/// \brief Definition of the PhysicsList class
//
//
//
/////////////////////////////////////////////////////////////////////////
//
// PhysicsList
//
// Created: 31.04.2006 V.Ivanchenko
//
// Modified:
// 04.06.2006 Adoptation of Hadr01 (V.Ivanchenko)
//
////////////////////////////////////////////////////////////////////////
// 

#ifndef DRAGONPhysicsList_h
#define DRAGONPhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "G4Material.hh"
#include "globals.hh"
#include "DRAGONRunAction.hh"

class G4VPhysicsConstructor;

namespace DRAGON
{

class DRAGONPhysicsListMessenger;
class DRAGONDetectorConstruction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DRAGONPhysicsList: public G4VModularPhysicsList
{
public:

  DRAGONPhysicsList(DRAGONDetectorConstruction*);
  ~DRAGONPhysicsList() override;

  void ConstructParticle() override;
  void ConstructProcess() override;
  void SetCuts() override;
  void AddIonGasModels();
  
  void uvinit();
  G4int Getipart() const {return ipart;}
  G4int Getirecoil() const {return irecoil;}
  G4int Getlkine() const {return lkine;}
  G4int Getitckov() const {return itckov;}
  G4double Getfkine(int i) const {return fkine[i];}
  G4double Getbeamenerg() const {return beamenerg;}
  G4double Getbeammom() const {return beammom;}
  G4double Getbeammass() const {return beammass;}
  G4double Getbuncht() const {return buncht;}
  G4double Getprodm() const {return prodm;}
  G4double Geteres0() const {return eres0;}
  G4double Gete0recoil() const {return e0recoil;}
  G4double Getresenerg() const {return resenerg;}
  G4double Getresmass() const {return resmass;}
  G4double Getreswidth() const {return reswidth;}
  G4double Getrecoilmom() const {return recoilmom;}
  G4double Getbunchl() const {return bunchl;}
  G4double Getbeamo() const {return beamo;}
  G4double Getm1() const {return m1;}
  G4double Getm2() const {return m2;}
  G4double Getz1() const {return z1;}
  G4double Getz2() const {return z2;}
  G4double Getzprod() const {return zprod;}
  G4double Geter() const {return er;}
  G4double Getgp() const {return gp;}
  G4double Getgg() const {return gg;}
  G4double Getomg() const {return omg;}
  G4double Getell() const {return ell;}
  G4double Getires() const {return ires;}
  G4double Getatarg() const {return atarg;}
  G4double Getmtarg() const {return mtarg;}
  
  G4double Getamax() const {return amax;}
  G4double Getbmax() const {return bmax;}
  G4double Getemax() const {return emax;}
  G4double Getsigx() const {return sigx;}
  G4double Getsigy() const {return sigy;}
  const G4bool Getalpha() const {return alpha;}

  DRAGONDetectorConstruction* GetDRAGONDet() {return DRAGON_det;}
  DRAGONRunAction* GetDRAGONRunAction(){return DRAGONRunAction_phys;}
   
  void SetFKIN(const G4String& input);
  void SetCUTs(const G4String& input);
  void SetInputCardname(const G4String& input);

  void SetBEAM(G4double ener){beamenerg = ener;}
  void SetTUNE(const G4String& input);
  void SetMTUN(const G4String& input);
  void SetGeom(DRAGONDetectorConstruction* det){DRAGON_det = det;}
  void Setmtarg(G4int mtarg_){mtarg = mtarg_;}
  void SetMtarg(G4Material* Mtarg_){Mtarg = Mtarg_;}
  void Setentdens(G4double entdens_){entdens = entdens_;}
  void Setexitdens(G4double exitdens_){exitdens = exitdens_;}
  void SetRunAction(DRAGONRunAction& RunAction_phys) {DRAGONRunAction_phys = &RunAction_phys;}
  
  void ureact();
  void beaminit();
  void ReadUserReactionNameList();
  void PrintIonProperties(G4ParticleDefinition* ion); //OJO Eliminar
  void DefineCo60copia();
  G4ParticleDefinition* DefineCo60cascade();
  
private:
    G4VPhysicsConstructor* fEmPhysicsList;
    G4VPhysicsConstructor* fDecPhysicsList;
	G4VPhysicsConstructor* fRadDecPhysicsList;   
    G4VPhysicsConstructor* fIonPhysicsList; 
    G4VPhysicsConstructor* fOptPhysicsList;
	G4VPhysicsConstructor* fStepLimiterPhysics;
	G4VPhysicsConstructor* fHadElastPhysics;
	G4VPhysicsConstructor* fHadPhysicsFTFP_BERT;
	
	DRAGONPhysicsListMessenger* fMessenger;
    DRAGONDetectorConstruction* DRAGON_det;
	DRAGONRunAction* DRAGONRunAction_phys;
 
    G4int mtarg;    //OJO Quitar de aqui (uggeom.hh)
    G4Material *Mtarg;
    
    G4String CARDNAME;
    
    G4double beam_mass_excess, recoil_mass_excess;    //From "params.hh"
    G4double part_width, gam_width, spin_stat_fac;
    G4double level[16], life[16], ztarg;
    G4double br[15][10], beamlifetime;
    G4int rstate, md[15][10];
    G4String beamtyp;
    G4String rectyp;

    G4bool alpha;       
    G4double targmass,resmass,resenerg,eres0,reswidth,prodm,atarg;
    G4double zbeam,abeam,zprod; 
    G4double entdens,exitdens;    
    G4int irecoil;
    
    G4double amass;   //OJO Desconozco su origen
    
    G4int ipart;         //From gckine.hh
    G4int lkine;     //From uevent.hh
    G4int iswit[10];     //From gcflag.hh
    G4double fkine[10];
    G4double egamma[10]; //From uevent.hh 
    
    G4double cutgam, cutele;
    
    //From beamcom.hh
    G4double amumev = 0.93149432E+03;
    G4double beamenerg;
    G4double beammass;
    G4double beammom;
    G4double sigx;
    G4double amax;
    G4double sigy;
    G4double bmax;
    G4double bunchl;
    G4double emax;
    G4double buncht;
    G4double tofg;
    G4double refenerg;
    G4double refatno;
    G4double refq;
    G4double offset[4];  
    G4double energscale;
    G4double e0recoil;
    G4double e0beam;
    G4double ex,ey,el; 
    G4double beamvel;
    G4double beamo;
    G4double recoilmom;
    G4double bscale;
    G4double escale;
    G4int ires;
    
    G4double m1, m2, mprod, z1, z2, er, gp, gg, omg, ell;
    
    G4double clight = 29979245800.;  //light velocity in cm s −1;
    
    G4int itckov;   //OJO Implementar Cherenkov

};
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


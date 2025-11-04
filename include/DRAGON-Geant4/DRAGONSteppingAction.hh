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
/// \file DRAGON/include/DRAGONSteppingAction.hh
/// \brief Definition of the DRAGON::DRAGONSteppingAction class

#ifndef DRAGONSteppingAction_h
#define DRAGONSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4Types.hh"
#include "G4SystemOfUnits.hh"
#include "G4String.hh"

class G4LogicalVolume;
class G4Step;
class G4ParticleDefinition;
class G4DynamicParticle;

namespace DRAGON
{

class DRAGONEventAction;
class DRAGONRunAction;
class DRAGONSteppingActionMessenger;


class DRAGONSteppingAction : public G4UserSteppingAction
{
	public:
		DRAGONSteppingAction(DRAGONEventAction* eventAction,DRAGONRunAction* runAction);
		~DRAGONSteppingAction() override = default;

		void UserSteppingAction(const G4Step*) override;
                
                void uvinit();
                void gustep(const G4Step*);
                void gustep_trgt(const G4Step*);
                void gustep_mitray(const G4Step*);
                void gustep_gbox(const G4Step*);
                void PrintDynamicParticleProperties(G4DynamicParticle* dynParticle);
                
                void ghidet();
                void ghipmt(const G4Step*);
                void gureact(const G4Step*);
                void Setipart(G4int ipart_){ipart = ipart_;}
                void Setires(G4int ires_){ires = ires_;}
                void Setirecoil(G4int irecoil_){irecoil = irecoil_;}
                void Setalpha(G4int alpha_){alpha = alpha_;}
                void Setn_detmate(G4int n_detmate_){n_detmate = n_detmate_;}
                void Setlkine(G4int lkine_){lkine = lkine_;}
                void Setitckov(G4int itckov_){itckov = itckov_;}
                void Settargtype(G4int targtype_){targtype = targtype_;}
                void Seteres(G4double eres_){eres = eres_;}
                void Setresmass(G4double resmass_){resmass = resmass_;}
                void Setrecoilmom(G4double recoilmom_){recoilmom = recoilmom_;}
                void Setresenerg(G4double resenerg_){resenerg = resenerg_;}
                void Seterescm(G4double erescm_){erescm = erescm_;}
                void Setprodm(G4double prodm_){prodm = prodm_;}
                void Setbeammass(G4double beammass_){beammass = beammass_;}
                void Sete0recoil(G4double e0recoil_){e0recoil = e0recoil_;}
				void Setlen_max(G4double len_max_){len_max = len_max_;}
				void Setmax_step(G4double max_step_){max_step = max_step_;}
						
		
		G4int Getipart(){return ipart;}
		G4int Getires(){return ires;}
		G4int GetMaxStep() const {return max_step;}
		G4int Getistop() const {return istop;}
		
	private:
		DRAGONEventAction* fEventAction = nullptr;
                DRAGONRunAction* frunAction = nullptr;
		G4LogicalVolume* fScoringVolume = nullptr;
		
		G4int max_step;
		G4double len_max;
		DRAGONSteppingActionMessenger* fSteppingActionMessenger = nullptr;

                G4int inwvol, number, number_old = 0;
                G4int ntmult = 0, ntmult_old = 0;
                G4String names = "", name_old = "", chname_nlevel = "", chname_2level = "", chname_pres = "", chname_old = "";
                G4int nlevel;
				G4int itrtyp;
                G4int ipart, ires, irecoil, lpart;
                G4bool alpha;
                G4String kcase = "";
                G4double vect[7];
                G4int lkine;
                G4int targtype;
                G4double gekin;
                G4double eres;
                G4double destep;
                G4double E_int;
		G4double tofg;
		G4double newm;
		G4double resmass, resenerg, erescm;
		G4double beamtof;
		G4int react;
		G4double xint, yint, zint;
		G4int istop;
                G4int ntargexit,nbeamout;
		G4double prodm;
		G4double e0recoil;
		G4int ngkine0,ngkine;
		G4double sleng;
		G4int* iflgk;
		G4int gkin;
	        //G4String wname, vname;
	        G4int n_detector = 0;
	        G4int pair_productions = 0;
	        G4int itckov, iscnt;
	        G4int n_detmate;
	        
	        G4double tlast, trec;
	        G4double beammass,recoilmom;
	        bool McpHit;
	        G4double xstop, ystop, zstop;
	        G4double xtest[10], ytest[10], etest[10];
	        G4double hits[5];
	        G4int recdet;
	        G4int nevent;
	        G4int dsssdpos;
	        G4int jslit, jstop;   //OJO Quitar de aqui
			
			G4int idtype;
			G4int nstep;
			

//C.            Information on the photon conversion point (i)
//C.           ------------------------------------------------
//C.
//C.    vname(i)   : The full name of the conversion volume 
//C.    ivcopy(2,i): The volume copy number of the conversion volume
//C.    true_conv(3,i): The coordinates of the conversion piont [cm]
//C.
//C.    Index array relating in the event ITRA -> n_gamma	

			static const G4int max_conv = 10;
            G4int ivcopy[2][max_conv];
	        G4double true_conv[3][max_conv];
	        G4String vname[max_conv];
			G4int nconv = 0;

	        static const G4int max_itra = 30;
            G4int index_track_to_gamma[max_itra];  
			

};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


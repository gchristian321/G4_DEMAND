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
/// \file DRAGON/include/DRAGONEventAction.hh
/// \brief Definition of the DRAGON::DRAGONEventAction class

#ifndef DRAGONEventAction_h
#define DRAGONEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

namespace DRAGON
{

class DRAGONRunAction;
class DRAGONEventAction;

/// Event action class

class DRAGONEventAction : public G4UserEventAction
{
  public:
    DRAGONEventAction(DRAGONRunAction* runAction);
    ~DRAGONEventAction() override = default;

    void BeginOfEventAction(const G4Event* event) override;
    void EndOfEventAction(const G4Event* event) override;

    void gutrev();
    void uvinit();
    void guout(const G4Event* event);
    void guout_mitray(const G4Event* event);
    void guout_gbox(const G4Event* event);
    
    void SetRunAction(DRAGONRunAction* runAction){fRunAction = runAction;}
    void Setn_flag(G4int n_flag_){n_flag = n_flag_;}
    void Setpair_productions(G4int pair_productions_){pair_productions = pair_productions_;}
	void SetMcpHit(bool McpHit_){McpHit = McpHit_;}
	void Setrecdet(bool recdet_){recdet = recdet_;}
	void Setjslit(G4int jslit_){jslit = jslit_;}
	void Setjstop(G4int jstop_){jstop = jstop_;}
        void Setedetect(G4double edetect_){edetect = edetect_; }
        void Setntot(G4int ntot_){ntot = ntot_; }	
		
	DRAGONRunAction* GetRunAction(){return fRunAction;}
    G4int Getn_flag(){return n_flag;}
    G4int Getpair_productions(){return pair_productions;}
	bool GetMcpHit(){return McpHit;}
	bool Getrecdet(){return recdet;}
	G4int Getjslit(){return jslit;}
	G4int Getjstop(){return jstop;}
        G4double Getedetect(){return edetect;}
        G4int Getntot(){return ntot; }

  private:
    DRAGONRunAction* fRunAction = nullptr;
  
    G4double E_int, E_rec, E_g[15], E_gp[15], cost_g[15], phi_g[15],
             cost_gp[15], cost_r, cosp_r, x_r, y_r, z_r, thet_r,
             xstop, ystop, zstop, xint, yint, zint,x,y,xp,yp,xtest[10],
             ytest[10],
             etest[10], beamtof;
    G4int Nodec, react, dsssdpos;
      G4int nout = 0;
	const G4String label[3] = {"TA","RG","ET"};
	G4int idevt;

public:
       bool McpHit;
	   G4int recdet;
	   G4int recoil_hit_ENDV;
	   G4int jslit,jstop;
	   G4int n_flag;
	   G4double tot_thrshld = 0.0;   //OJO
	   G4int pair_productions;
	   G4double edetect;
	   G4double z_react;

           G4double e_bgos_total, e_bgo_first, e_bgo_second,
           e_bgo_first_ab, e_bgo_second_ab, gammatof, e0_conv;
           G4int num_bgos_hit, num_bgo_first,
           num_bgo_second, num_bgos_hit_ab,
           num_bgo_first_ab, num_bgo_second_ab;
         
           G4double e_detect;
	   G4int nloss[6]; 
           G4int ievent = 0;
           G4int ntot;
           G4double E_threshold;
           G4double x_mean, y_mean, z_mean;

/*	
C.
C.                Information on the photons
C.                --------------------------
C.
C.    true_e: The true energy of photon [MeV]
C.    true_d: The true photon unit vector [unit-vector-components]
C.
C.    -------------------------------------------------
C.
C.    index_track_to_gamma: ITRA -> n_gamma
C.
C.    ntot: Total number of PMT hits 
C.    ntot: Total number of PMT hits 
C.
C.    ifngr_max: The coordinates of the maximum finger
C.
C.    x_conv, y_conv, z_conv: The coordinates of 'first' conversion [cm]
C.
C.    x_max, y_max, z_max:    The coordinates of max energy deposition [cm]
C.
C.    x_mean, y_mean, z_mean: The coordinates of mean energy deposition [cm]
C.
*/
      static const G4int max_photon = 10;
      
	  G4double true_e[max_photon], true_d[3][max_photon];
 
      G4int ifngr_max;

      G4double x_conv, y_conv, z_conv;
      G4double x_max, y_max, z_max;


};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif



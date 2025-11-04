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
/// \file DRAGON/include/DRAGONDetectorConstruction.hh
/// \brief Definition of the DRAGON::DRAGONDetectorConstruction class

#ifndef DRAGONDetectorConstruction_h
#define DRAGONDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"      //Geant4
#include "G4SystemOfUnits.hh"
#include "G4LogicalVolume.hh"
#include "G4EqMagElectricField.hh"
#include "G4PropagatorInField.hh"
#include "G4TransportationManager.hh"
#include "G4ChordFinder.hh"
#include "G4ClassicalRK4.hh"

#include "DRAGONRunAction.hh"                  //local
#include "DRAGONDetectorMessenger.hh"
#include "DRAGONEMField.hh"


class G4VPhysicalVolume;
class G4UserLimits;

namespace DRAGON
{

class DRAGONPhysicsList;


class DRAGONDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DRAGONDetectorConstruction(DRAGONPhysicsList*);
    ~DRAGONDetectorConstruction() override;
    
    void uvinit();
    void ugeom();
    void udet();
    void udetmitray();
    void ugeo_space();
    void ugeo_defin(); 
    void ugeo_detector();
    void ugeo_finger();
    void ugeo_pmt();
	void ugeo_trgt();
	void ugeo2_trgt();   //OJO no se usa en ugeom_trgt_soltarg
    void ugeo_trgt_small();
    void ugeo_trgt_large();
    void ugeo_trgt_small_up();
    void ugeo_trgt_small_down();
    void ugeo_trgt_small_right();
    void ugeo_trgt_small_left();  
    void ugeo_trgt_small_hole();  
    
    void ugeom_setup();
    void mitray_setup();
    void ugeo_fcup(G4double pos[3]);
    void ugeo_end(G4double pos[3]);
    void ugeo_test(G4int k,G4double pos[3], G4double rot_angles[3]);
    void ugeo_start(G4double pos[3]);
    void ugeo_dipole(G4int k);
    void ugeo_edipol(G4int k);
    void ugeo_mpole(G4int k,G4double rot_angles[14*3]);
    void ugeo_sole(G4int k);
    void ugeo_col(G4double pos[45*3],G4double rot_angles[45*3],G4double data[6],G4String rname);
    void ugeo_mcp(G4double pos[3],G4double data[2],G4String rname,G4int nmcp);
    void ugeo_dssd(G4double pos[3]);
	void neighborhood();
    
    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;
    
    G4double GetTLrms() const {return TLrms;}
    G4double GetRrms() const {return Rrms;}
    G4double Gettargetl() const {return targetl;}
    G4int GetTubeType() const {return tubetype;}
    G4int GetTargType() const {return targtype;}
    G4int GetMTypePMT() const {return mtype_pmt;}
    G4Material* Getmtarg() const {return G4mtarg;}
	G4double GetPMTSize() const {return pmt_size;}
	G4double GetPMTLength() const {return pmt_length;}
	const G4double* GetWall() const{return wall;}
	G4int GetNDetMate() const {return n_detmate;}
	G4double GetAprt() const {return aprt;}
	G4double GetBulkAbs() const {return bulk_absorption;}
	G4double GetPaintAbs() const {return paint_absorption;}
	G4LogicalVolume* GetWRLDLogVolume() {return WRLD_log;}
	G4bool GetCheckOverlaps() {return checkOverlaps;}
    G4double Getentdens() {return entdens;}
    G4double Getexitdens() {return exitdens;}
    G4double Getlen_max() const {return len_max;}
	G4double Getmax_step() const {return max_step;}
	G4int Getntot(const G4String& logicalName);
		
    void SetRunAction(DRAGONRunAction& RunAction_det) {DRAGONRunAction_det = &RunAction_det;}
    void SetTUBE(G4int newVal);
    void SetTARG(G4int newVal) {targtype = newVal;}
    void SetMPMT(G4int newVal) {mtype_pmt = newVal;}
    void SetPMTR(const G4String& input);
    void SetWALL(const G4String& input);
 	void SetBGAP(G4double newValue) {box_width = newValue;}
 	void SetFSID(const G4String& input);
    void SetDMAT(G4int newVal) {n_detmate = newVal;}
    void SetHOLE(G4double newVal) {aprt = newVal;}
 	void SetBLKA(G4double newVal) {bulk_absorption = newVal;}
 	void SetREFL(G4double newVal) {paint_absorption = newVal;}
 	void SetCheckOverlaps(G4bool newBool) {checkOverlaps = newBool;}
    void SetPhys(DRAGONPhysicsList* phys){DRAGON_phys = phys;}
    void SetMaxStepsTrackLength(const G4String& input);
    void SetMASK(const G4String& input);
    void SetSHLD(const G4String& input);
    void PrintVolumesInfo();
 private:
    DRAGONRunAction* DRAGONRunAction_det = nullptr;
    DRAGONPhysicsList* DRAGON_phys = nullptr;
	DRAGONDetectorMessenger* fDetMessenger = nullptr; 
    G4UserLimits* userLimits = nullptr;
	G4LogicalVolume* WRLD_log;
	G4bool checkOverlaps;

	G4int max_step; 
    G4double len_max;
    G4int mtarg, mcent, ment[6], mex[6], mbox;
    G4Material* G4mtarg;
    G4double Rrms, TLrms, 
	         targetl, rent, lent, len1, riren1, rilen1,
			 len2, riren2, rilen2,
			 len3, riren3, rilen3,
             lex1, rilex1, rirex1, lex2, rilex2, rirex2,
             lex3, rilex3, rirex3, lex4, rilex4, rirex4;
    G4double s_finger, z_finger, air_gap, d_air[2], d_mtl, gap;
    G4double zent[6], zex[7];
    G4int  targtype, tubetype;
    G4double shield_end[2], xtnd_block;
    G4int mask[30];     
    G4String MASK;
    G4double hexagon_small_width, hexagon_large_width, depth,
             wall[3], box_width, box_length, aprt, col_length, 
             col_collar_length;
    G4int mtype_pmt;
    G4double pmt_size, pmt_length;
    G4double bulk_absorption, paint_absorption;
     G4int n_detmate;
    G4double entdens,exitdens;
	
	static const G4int max = 30;
	
	
	public:
	    static const G4int max_hexagon = 30;
	    static const G4int Nn = 10;
	    G4int n_fngr[Nn][max] = {0};	    
	    G4double x_fngr[max_hexagon], y_fngr[max_hexagon], z_fngr[max_hexagon];
	    G4int adjacency_matrix[30][30];
    };

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

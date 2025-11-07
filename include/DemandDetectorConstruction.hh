/// \file DemandDetectorConstruction.hh
/// \brief Definition of the DemandDetectorConstruction class

#ifndef DemandDetectorConstruction_h
#define DemandDetectorConstruction_h 1

#include <array>
#include <vector>
#include <string>
#include <map>

#include "G4VUserDetectorConstruction.hh"
#include "G4Box.hh"
#include "globals.hh"
#include "DemandDetectorMessenger.hh"

class G4Material;
class G4VPhysicalVolume;
class G4AssemblyVolume;
class G4LogicalVolume;

class DemandDetectorConstruction: public G4VUserDetectorConstruction {
public:
	static const int kCube = 0;
	static const int kBar  = 1;
	struct Module_t {
		double m_Dist;
		double m_Theta;
		double m_Phi;
		std::array<int,3> m_Voxnum;
		std::array<double,3> m_Voxsize;
		double m_Threshold;
		int  m_ReadoutType;
		bool m_Rotate;
		G4ThreeVector m_Position;
		std::array<int,3> m_HavePosition;
	};

public:
	class ScintillatorBox : public G4Box {
 public:
		ScintillatorBox (const G4String &pName, G4double pX, G4double pY, G4double pZ):
			G4Box(pName,pX,pY,pZ){fReadoutType = kCube;}
		virtual ~ScintillatorBox(){}
		int GetReadoutType() const { return fReadoutType; }
		void SetReadoutType(int type) { fReadoutType = type; }
 private:
		int fReadoutType;
	};
	
public:
  DemandDetectorConstruction();
  virtual ~DemandDetectorConstruction();

	virtual G4VPhysicalVolume* Construct();
	virtual void ConstructSDandField();

	G4AssemblyVolume* CreateAssembly(const Module_t& module);
	void SetModules(const std::vector<Module_t>& modules)
		{ fModules = modules; }
	void AddModule();
	
  G4Material* GetScintillatorMaterial();
	const std::vector<std::string>& GetSensitiveDetectorNames() const
		{ return fSensitiveDetectorNames; }

	const Module_t* GetModule(G4int) const;
	Module_t* GetModuleByOriginalSpecification(size_t);
	const Module_t* GetModuleByOriginalSpecification(size_t) const;

	void SetTargetMaterial(const G4String&);
	void SetTargetThickness(G4double thick) { fTargetThickness = thick; }
	const G4Material* GetTargetMaterial() const { return fTargetMaterial; }
	G4double GetTargetThickness() const { return fTargetThickness; }
	bool GetHaveTarget() const { return fHaveTarget; }
	void SetBGOMask(const G4String& mask){fBGOMask=mask;}
	const G4String& GetBGOMask() const{return fBGOMask;}
	void SetUseDRAGON(bool use){fUseDRAGON=use;}
	bool GetUseDRAGON() const {return fUseDRAGON;}
	
private:
	bool fCheckOverlaps;
	std::vector<Module_t> fModules;
	std::vector<G4AssemblyVolume*> fAssemblies;
	std::vector<std::string> fSensitiveDetectorNames;
	std::map<G4int, const Module_t*> fModuleIndexMap;
	DemandDetectorMessenger* fDetectorMessenger;
	G4Material* fTargetMaterial;
	G4double fTargetThickness;
	bool fHaveTarget;
	G4String fBGOMask;
	bool fUseDRAGON;



	// DRAGON-Geant4 Stuff
private:
	void ugeo_defin();
	void ugeo_detector();
private:
	bool checkOverlaps;
	G4LogicalVolume* WRLD_log;
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




#endif


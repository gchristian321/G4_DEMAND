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
class G4PVPlacement;
namespace CLHEP {
class HepRotation;
}
namespace DRAGON {
class DRAGONDetectorConstruction;
class DRAGONPhysicsList;
}

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
	void SetUseChamber(bool use){fUseChamber=use;}
	bool GetUseChamber() const {return fUseChamber;}
	G4VPhysicalVolume* ConstructWithDragon();
	G4VPhysicalVolume* ConstructWithoutDragon();
	void ConstructNeutronDetectorModules(G4VPhysicalVolume* world);
	void SetDragonDetectorConstruction(DRAGON::DRAGONDetectorConstruction* det) {fDragonDet = det;}
	DRAGON::DRAGONDetectorConstruction* GetDragonDetectorConstruction() const  {return fDragonDet;}
	void SetDragonPhysicsList(DRAGON::DRAGONPhysicsList* phys) {fDragonPhys = phys;}
	DRAGON::DRAGONPhysicsList* GetDragonPhysicsList() const  {return fDragonPhys;}

	void ConstructS2230Detectors(G4VPhysicalVolume*);
	G4LogicalVolume* ConstructS2230CasingAssembly();
	void PlaceCasingAssembly(
		G4LogicalVolume* assembly, G4LogicalVolume* mother, const G4ThreeVector& pos);
	bool GetUseS2230Assembly() const   { return fUseS2230Assembly; }
	void SetUseS2230Assembly(bool use) { fUseS2230Assembly = use;  }
	void SetS2230Threshold(size_t i, double thresh) {
		fS2230Thresholds.at(i) = thresh;
	}
	void SetS2230ThresholdSigma(G4double sig) { fS2230ThresholdSigma = sig; }
	G4double GetS2230ThresholdSigma() const { return fS2230ThresholdSigma; }
	
private:
	bool fCheckOverlaps;
	std::vector<Module_t> fModules;
	std::vector<G4AssemblyVolume*> fAssemblies;
	std::vector<std::string> fSensitiveDetectorNames;
	std::map<G4int, const Module_t*> fModuleIndexMap;
	DemandDetectorMessenger* fDetectorMessenger;
	G4Material* fTargetMaterial;
	G4Material* fMatAir;
	G4Material* fMatSteel;
	G4Material* fMatAl;
	G4Material* fMatVacuum;
	G4Material* fMatMuMetal;
	G4Material* fMatPMTGlass;
	
	G4double fTargetThickness;
	bool fHaveTarget;
	G4String fBGOMask;
	bool fUseDRAGON;
	bool fUseChamber;
	bool fUseS2230Assembly;
	std::vector<double> fS2230Thresholds;
	G4double fS2230ThresholdSigma;
	DRAGON::DRAGONDetectorConstruction* fDragonDet;
	DRAGON::DRAGONPhysicsList* fDragonPhys;

private:
	// class to make components of the outer "skeleton" wireframe
	// not currently used but keep in case of desire to revive some day
	class SkeletonFrame {
 public:
		SkeletonFrame();
		~SkeletonFrame(){}
		G4double GetSFthick()const{return fSFthick;};
		G4LogicalVolume* GetFrameLV()const{return fFrameLV;}
		G4LogicalVolume* GetHLV()const{return fHLV;}
		G4double GetTotalThickness()const{return fTotalThickness;}
		G4double GetHYpos()const{return fHoleYpos;}
 private:
		void ConstructFrame();
		void ConstructH();
 private:
		const G4double fSFthick;
		const G4double fEps;
		const G4ThreeVector fEps3;
		G4double fTotalThickness; // outside->outside
		G4double fHoleYpos;
		G4LogicalVolume* fFrameLV;
		G4LogicalVolume* fHLV;
		G4Material* fMatAl;
	};

};

#endif


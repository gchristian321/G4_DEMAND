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
};




#endif


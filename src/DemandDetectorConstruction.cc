#include <algorithm>
#include <map>

#include "DemandDetectorConstruction.hh"
#include "DemandSD.hh"
// DRAGON
#include "DRAGONDetectorConstruction.hh"
#include "DRAGONPhysicsList.hh"
#include "Materials.hh"

#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
//#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4ExtrudedSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4TwoVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4AssemblyVolume.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4UserLimits.hh"
#include "G4Trd.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4UnionSolid.hh"
#include "CADMesh.hh"
#include "G4IonisParamMat.hh"
#include "G4VisExtent.hh"


// Auto-generated code for skeleton frame parts
#include "skelframe/FrontFrameUpper.inc"
#include "skelframe/FrontFrameLower.inc"
#include "skelframe/BackFrameUpper.inc"
#include "skelframe/BackFrameLower.inc"
#include "skelframe/H1.inc"
#include "skelframe/H2.inc"
#include "skelframe/H3.inc"
#include "skelframe/H4.inc"
#include "skelframe/SideHoles.inc"
#include "skelframe/TopHoles.inc"
#include "skelframe/BottomHoles.inc"
#include "skelframe/VerticalSupport.inc"
#include "skelframe/VerticalSupportTop.inc"
#include "skelframe/HorizontalSupportUpper.inc"
#include "skelframe/HorizontalSupportLower.inc"


using namespace std;
namespace DRAGON {
extern G4double PDE_zpos;
extern G4double PDE_zlen;
}


DemandDetectorConstruction::DemandDetectorConstruction()
	:G4VUserDetectorConstruction(),
	 fCheckOverlaps(false)
{
	fDetectorMessenger = new DemandDetectorMessenger(this);
	fTargetMaterial = nullptr;
	fTargetThickness = 0;
	fHaveTarget = false;
	fBGOMask="";
	fUseDRAGON=false;
	fUseChamber=true;
	fUseS2230Assembly = false;
	for(size_t i=0;i<8;++i){
		fS2230Thresholds.push_back(100*keV);
	}
	fDragonDet=nullptr;
	fDragonPhys=nullptr;

	G4NistManager *nist = G4NistManager::Instance();
  fMatAir   = nist->FindOrBuildMaterial("G4_AIR");
  fMatSteel = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  fMatAl    = nist->FindOrBuildMaterial("G4_Al");
	fMatVacuum = nist->FindOrBuildMaterial("G4_Galactic");
	fMatPMTGlass = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");

	// Mu-metal approximation (80% Ni, 20% Fe) for PMT shield
  auto* elNi = nist->FindOrBuildElement("Ni");
  auto* elFe = nist->FindOrBuildElement("Fe");
  fMatMuMetal = new G4Material("MuMetal_80Ni20Fe", 8.7*g/cm3, 2);
  fMatMuMetal->AddElement(elNi, 0.80);
  fMatMuMetal->AddElement(elFe, 0.20);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DemandDetectorConstruction::~DemandDetectorConstruction()
{
	delete fDetectorMessenger;
}


G4VPhysicalVolume* DemandDetectorConstruction::Construct()
{
	auto worldPV = fUseDRAGON ?
		ConstructWithDragon() :
		ConstructWithoutDragon() ;
	return worldPV;
}

G4VPhysicalVolume* DemandDetectorConstruction::ConstructWithDragon()
{
	// \todo look at tube stuff
	fDragonDet = new DRAGON::DRAGONDetectorConstruction(nullptr);
	fDragonPhys->SetGeom(fDragonDet);
	fDragonDet->SetPhys(fDragonPhys);

	if(!fDragonDet){
		throw std::runtime_error(
			"DemandDetectorConstruction:: UseDRAGON on but fDragonDet not set!");
	}
	fDragonDet->SetTUBE(0);
	fDragonDet->SetPMTR("2.54 2.5");
	fDragonDet->SetHOLE(4.496);
//	fDragonDet->SetTARG(2);
	fDragonDet->SetTARG(0);
	fDragonDet->SetCheckOverlaps(false);
	if(fBGOMask != ""){
		fDragonDet->SetMASK(fBGOMask);
	}
	G4VPhysicalVolume* WRLD_phys = fDragonDet->Construct();

	// neutron detectors
	// place neutron detectors in "DETE" volume, not world
	if (GetUseS2230Assembly()) {
		// use packaged S2230 asselbly w/ detectors inside
		ConstructS2230Detectors(
			G4PhysicalVolumeStore::GetInstance()->GetVolume("DETE")
			);
	}
	else {
		//  ose "old neutron detectors from module commands"
		ConstructNeutronDetectorModules(
			G4PhysicalVolumeStore::GetInstance()->GetVolume("DETE")
			);
		//ConstructNeutronDetectorModules(WRLD_phys);
	}
	
	return WRLD_phys;
}

G4VPhysicalVolume* DemandDetectorConstruction::ConstructWithoutDragon()
{
  G4double worldSizeXY = 10*meter;
  G4double worldSizeZ  = 10*meter;
	   
  //     
  // World
  //
  auto worldS 
    = new G4Box("WRLD",           // its name
                 worldSizeXY/2, worldSizeXY/2, worldSizeZ/2); // its size
                         
  auto worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 fMatVacuum,  // its material
                 "WRLD");         // its name
                                   
  auto world
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume                         
                 "WRLD",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
	
	
	G4VisAttributes * worldAttr = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
  worldAttr->SetVisibility(false);
  worldLV->SetVisAttributes(worldAttr);

	// neutron detectors
	// place neutron detectors in "DETE" volume, not world
	if (GetUseS2230Assembly()) {
		// use packaged S2230 asselbly w/ detectors inside
		ConstructS2230Detectors(world);
	}
	else {
		//  ose "old neutron detectors from module commands"
		ConstructNeutronDetectorModules(world);
	}
	
	//     
  // Target
  //
	if(fTargetMaterial && fTargetThickness >= 1e-5*micrometer) {
		G4double size_XY = 3*cm;
		G4Box* solidTrgt =    
			new G4Box("Target",                    //its name
								0.5*size_XY, 0.5*size_XY, 0.5*fTargetThickness); //its size
      
		G4LogicalVolume* logicTrgt =                         
			new G4LogicalVolume(solidTrgt,             //its solid
													fTargetMaterial,       //its material
													"Target");             //its name

		new G4PVPlacement(0,                       //no rotation
											G4ThreeVector(0,0,0),    // location
											logicTrgt,               //its logical volume
											"Target",                //its name
											worldLV,                 //its mother  volume
											false,                   //no boolean operation
											0,                       //copy number
											fCheckOverlaps);          //overlaps checking
		fHaveTarget = true;
	}	

	if(GetUseChamber())
	{
		//
		// Target Chamber
		//
		G4Box *solidOuterChamber = new G4Box("solidOuterChamber",25.4*mm, 128.65*mm, 85.725*mm);
		G4LogicalVolume *logicOuterChamber = new G4LogicalVolume(solidOuterChamber, fMatAl, "logicOuterChamber");
		G4VPhysicalVolume *physOuterChamber = new G4PVPlacement(0, G4ThreeVector(0.,-96.91*mm,0.), logicOuterChamber, "physOuterChamber", worldLV, false, 0, true);

		G4Box *solidInnerChamber = new G4Box("solidInnerChamber", 22.23*mm, 125.48*mm, 82.555*mm);
		G4LogicalVolume *logicInnerChamber = new G4LogicalVolume(solidInnerChamber, fMatVacuum, "logicInnerChamber");
		G4VPhysicalVolume *physInnerChamber = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), logicInnerChamber, "physInnerChamber", logicOuterChamber, false, 0, true);

		G4Tubs *solidExitHole = new G4Tubs("solidExitHole", 0.0, 4.5*mm, 1.585*mm, 0.0, 360.0);
		G4LogicalVolume *logicExitHole = new G4LogicalVolume(solidExitHole, fMatVacuum, "logicExitHole");
		G4VPhysicalVolume *physExitHole = new G4PVPlacement(0, G4ThreeVector(0., 96.91*mm, 84.14*mm), logicExitHole, "physExitHole", logicOuterChamber, false, 0, true);

		//
		// Beam Pipe Collar
		//
		G4Tubs *solidCollar = new G4Tubs("solidCollar", 4.5*mm, 19.02*mm, 4*mm, 0., 360.);
		G4LogicalVolume *logicCollar = new G4LogicalVolume(solidCollar, fMatAl, "logicCollar");
		G4VPhysicalVolume *physCollar = new G4PVPlacement(0, G4ThreeVector(0.,0.,89.725*mm), logicCollar, "physCollar", worldLV, false, 0, true);

		//
		// Beam Pipe, with collar
		//
		G4Tubs *solidBeamPipe = new G4Tubs("solidBeamPipe", 4.5*mm, 12.68*mm, 21.3875*mm, 0., 360.);
		G4LogicalVolume *logicBeamPipe = new G4LogicalVolume(solidBeamPipe, fMatAl, "logicBeamPipe");
		G4VPhysicalVolume *physBeamPipe = new G4PVPlacement(0, G4ThreeVector(0., 0., 115.1125*mm), logicBeamPipe, "physBeamPipe", worldLV, false, 0, true);

		//
		// Beam pipe, no collar
		//
		//G4Tubs *solidBeamPipe = new G4Tubs("solidBeamPipe", 4.5*mm, 12.68*mm, 25.3875*mm, 0., 360.);
		//G4LogicalVolume *logicBeamPipe = new G4LogicalVolume(solidBeamPipe, fMatAl, "logicBeamPipe");
		//G4VPhysicalVolume *physBeamPipe = new G4PVPlacement(0, G4ThreeVector(0., 0., 111.1125*mm), logicBeamPipe, "physBeamPipe", worldLV, false, 0, true);

		//
		// Test detector
		//
		/*G4Box *solidOGS = new G4Box("solidOGS", 15.0*mm, 15.0*mm, 15.0*mm);
			G4LogicalVolume *logicOGS = new G4LogicalVolume(solidOGS, matAir, "logicOGS");
			G4VPhysicalVolume *physOGS = new G4PVPlacement(0, G4ThreeVector(-30.0*mm, 30.0*mm, 121.5*mm), logicOGS, "physOGS", worldLV, false, 0, true);*/

		//Create simple enclosure box
		/*G4double Enclosure_height = 161*mm;
			G4double Enclosure_width = 43.45*mm;
			G4double Enclosure_length = 41.5*mm;
			G4double Enclosure_thickness = 3*mm;

			//Define the outer enclosure
			G4Box *fullEnclosure = new G4Box("fullEnclosure", Enclosure_width/2, Enclosure_height/2, Enclosure_length/2);
			G4Box *subtractEnclosure = new G4Box("subtractEnclosure", (Enclosure_width-2*Enclosure_thickness)/2, (Enclosure_height-2*Enclosure_thickness)/2, (Enclosure_length-2*Enclosure_thickness)/2);
			G4SubtractionSolid *Enclosure = new G4SubtractionSolid("Enclosure",fullEnclosure,subtractEnclosure);
	
			G4LogicalVolume* logicEnclosure = new G4LogicalVolume(Enclosure, fMatAl, "logicEnclosure");

			new G4PVPlacement(0,G4ThreeVector(0.0,0.0,200*mm),logicEnclosure,"physEnclosure", worldLV, false, 0, true);*/
	}
	
	return world;
}

void DemandDetectorConstruction::ConstructNeutronDetectorModules(G4VPhysicalVolume* world)
{
	for(const auto& module : fModules) {
		G4AssemblyVolume* assembly = CreateAssembly(module);

		G4ThreeVector Tm;
		if(module.m_HavePosition[0] == 1 &&
			 module.m_HavePosition[1] == 1 &&
			 module.m_HavePosition[2] == 1)     {
			Tm = module.m_Position;
		}
		else {
			// use r, theta, phi
			double halfThickness = module.m_Voxnum[2] * (module.m_Voxsize[2]/2);
			double rho = module.m_Dist + halfThickness;
			Tm.set(rho*sin(module.m_Theta)*cos(module.m_Phi),
						 rho*sin(module.m_Theta)*sin(module.m_Phi),
						 rho*cos(module.m_Theta)
				);
		}
			
		G4RotationMatrix Rm;
		if(module.m_Rotate){
			Rm.rotateY(module.m_Theta);
			Rm.rotateZ(module.m_Phi);
		}

		G4Transform3D Trans(Rm,Tm);
		assembly->MakeImprint(
			world->GetLogicalVolume(), Trans);
		fAssemblies.push_back(assembly);
	}
}

void DemandDetectorConstruction::ConstructSDandField() {
	auto demandDetector = new DemandSD("DemandSensitive");	

	if(GetUseS2230Assembly()) {
#if 1
		G4LogicalVolume* DEMAND_scintLV =	G4LogicalVolumeStore::GetInstance()
			->GetVolume("DEMAND_scintLV", true);
		if(!DEMAND_scintLV){throw std::runtime_error("no DEMAND_scintLV found");}

		DEMAND_scintLV->SetSensitiveDetector(demandDetector);
		G4SDManager::GetSDMpointer()->AddNewDetector(demandDetector);

    // std::stringstream sstr; sstr << "CasingPV_" << i;
    // new G4PVPlacement(casingRot,
    //                   casingPos,
    //                   casingAssemblyLV,
    //                   sstr.str().c_str(),
    //                   moduleLV,
    //                   false,
    //                   (int)i,
    //                   true);

		for (auto pv : *G4PhysicalVolumeStore::GetInstance()) {
			if (std::string(pv->GetName()).substr(0,8) == "CasingPV") {
				G4int copyNo = pv->GetCopyNo();
				G4double thresh = 100*keV;
				try { thresh = fS2230Thresholds.at(copyNo); }
				catch(std::exception& e){
					G4cerr << "ERROR: can't set threshold for detno " << copyNo << G4endl;
					throw e;
				}
				demandDetector->AddThreshold(copyNo, thresh);
				G4cout << " S2230 Threshold, det " << copyNo << ", " << thresh/keV << " keVee\n";
			}
		}
#endif
	}
	else{	
		int detno = 1, moduleno = 0;
		for(G4AssemblyVolume* assembly : fAssemblies) {
			unsigned ncubes = assembly->TotalImprintedVolumes ();
			auto itBegin = assembly->GetVolumesIterator();
			auto thisModule = fModules.at(moduleno++);
			G4double thresh = thisModule.m_Threshold;
			
			for(auto it = itBegin; it< itBegin + ncubes; ++it) {
				G4LogicalVolume * cubeLogic = (*it)->GetLogicalVolume();
			
				G4SDManager::GetSDMpointer()->AddNewDetector(demandDetector);
				cubeLogic->SetSensitiveDetector(demandDetector);

				G4int copyNo = (*it)->GetCopyNo();
				demandDetector->AddThreshold(copyNo,thresh);
				fModuleIndexMap.emplace(copyNo, &thisModule);
			
				char name[256];
				sprintf(name,"demandCrystal_%d",detno);	
				fSensitiveDetectorNames.push_back(name);

				++detno;
			}
		}
	}
}

namespace { G4Element* GetElementFromLibrary(const std::string& Name) {
	static G4Element *m_D = nullptr;
	static G4Element *m_T = nullptr;
	static G4Element *m_He3 = nullptr;
	
  if (Name == "D" || Name == "d") {
    if (!m_D) {
			m_D = new G4Element(Name.c_str(), Name.c_str(), 1);
      G4Isotope* isotope
				= new G4Isotope(Name.c_str(), 1, 2, 2.01410178 * g / mole);
      m_D->AddIsotope(isotope, 1);
    }
    return m_D;
  }
	
  else if (Name == "T" || Name == "t") {
    if (!m_T) {
      m_T = new G4Element(Name.c_str(), Name.c_str(), 1);
      G4Isotope* isotope
				= new G4Isotope(Name.c_str(), 1, 3, 3.0160492 * g / mole);
      m_T->AddIsotope(isotope, 1);
    }
    return m_T;
  }

  else if (Name == "He3" || Name == "3He") {
    if (!m_He3) {
      m_He3 = new G4Element(Name.c_str(), Name.c_str(), 1);
      G4Isotope* isotope
				= new G4Isotope(Name.c_str(), 2, 1, 3.0160293 * g / mole);
      m_He3->AddIsotope(isotope, 1);
    }
    return m_He3;
  }

  G4NistManager* man = G4NistManager::Instance();
  return man->FindOrBuildElement(Name.c_str());
} }

G4Material* DemandDetectorConstruction::GetScintillatorMaterial()
{
	static G4Material* material = nullptr;
	if(material == nullptr) {
		G4double density = 1.096 * g / cm3;
		material = new G4Material("Organic-Glass-Scintillator", density, 3);
		material->AddElement(GetElementFromLibrary("H"), 36);
		material->AddElement(GetElementFromLibrary("C"), 42);
    material->AddElement(GetElementFromLibrary("Si"), 1);

		// Set birks constant from fit to published proton quenching
		// data [T.A. Laplace et al 2020 JINST 15 P11020].
		// Data retreived from berkeley lab scintillator library
		// https://scintillator.lbl.gov/organic-glass-quenching-data/
		//
		// For the fit, the scale parameter S is fixed to 1.
		// Fit on GAC OneDrive SMU/projects/DRAGON-alpha_n/S2230/NeutronSinglesAnalysis_sydney/OGS - Quenching.ipynb
		//
		// Note that the Birks quenching is used for heavy ions and protons
		// above 25.5 MeV.  Below that a semi-emperical fit is used for protons.
		material->GetIonisation()->SetBirksConstant(0.07283 * mm/MeV);
	}
	return material;
}


namespace {
double getLocalPosition (int i, int n, double width) {
	double lp = (i-n/2)*width;
	if(n%2 == 0) { lp += width/2; }
	return lp;
} }

G4AssemblyVolume* DemandDetectorConstruction::CreateAssembly
(const Module_t& module)
{
	// Static stuff
	static G4VisAttributes*	VisSquare =
		new G4VisAttributes(G4Colour(0, 1, 1, 0.5));   
	static map<std::array<double,3>, G4LogicalVolume*> squareVolumes;
	auto buildLocalDetector =
		[&]	(){
			auto it = squareVolumes.find(module.m_Voxsize);
			if (it != squareVolumes.end()){
				return it->second;
			}
			auto localBox_ = new ScintillatorBox(
				"Demand_LocalBox_",
				0.5*module.m_Voxsize[0],
				0.5*module.m_Voxsize[1],
				0.5*module.m_Voxsize[2]);
			localBox_->SetReadoutType(module.m_ReadoutType);
			auto squareDetector =	new G4LogicalVolume
				(localBox_,GetScintillatorMaterial(),
				 "logic_Demand_LocalBox_",0,0,0);
			squareDetector->SetVisAttributes(VisSquare);
			squareVolumes.emplace(module.m_Voxsize,squareDetector);
			return squareDetector;
		};

	G4AssemblyVolume* assemblyDetector = new G4AssemblyVolume();
	
	for(int iz=0; iz< module.m_Voxnum[2]; ++iz) {
		double zlocal = getLocalPosition(iz,module.m_Voxnum[2],module.m_Voxsize[2]);
		G4cout << "zlocal ("<<iz<<"): " << zlocal << G4endl;
		// zlocal += ((nz-1)*voxsize[2]/2);

		for(int ix = 0; ix< module.m_Voxnum[0]; ++ix) {
			double xlocal = getLocalPosition(
				ix,module.m_Voxnum[0],module.m_Voxsize[0]);
			
			for(int iy = 0; iy< module.m_Voxnum[1]; ++iy) {
				double ylocal = getLocalPosition(
					iy,module.m_Voxnum[1],module.m_Voxsize[1]);

				G4ThreeVector Ta(xlocal,ylocal,zlocal);
				G4RotationMatrix Ra;
				G4Transform3D Trans_local(Ra,Ta);
				G4LogicalVolume* cubeLV = buildLocalDetector();				
				assemblyDetector->AddPlacedVolume(
					cubeLV,  Trans_local); 
			}
		}
	}
		
	return assemblyDetector;
}

const DemandDetectorConstruction::Module_t*
DemandDetectorConstruction::GetModule(G4int copyNo) const
{
	auto it = fModuleIndexMap.find(copyNo);
	if ( it == fModuleIndexMap.end() ) {
		throw std::invalid_argument(
			"Bad copyNo to DemandDetectorConstruction::GetModule()");
	}
	return it->second;
}

DemandDetectorConstruction::Module_t*
DemandDetectorConstruction::GetModuleByOriginalSpecification(size_t i)
{
	try {
		return &(fModules.at(i));
	}
	catch(std::exception&) {
		char buf[4096];
		sprintf(buf, "Invalid index, %i, to DemandDetectorConstruction::"
						"GetModuleByOriginalSpecification() -- maximum is: %i",
						int(i), int(fModules.size())-1);
		throw std::invalid_argument(buf);
	}
}

const DemandDetectorConstruction::Module_t*
DemandDetectorConstruction::GetModuleByOriginalSpecification(size_t i) const
{
	try {
		return &(fModules.at(i));
	}
	catch(std::exception&) {
		char buf[4096];
		sprintf(buf, "Invalid index, %i, to DemandDetectorConstruction::"
						"GetModuleByOriginalSpecification() -- maximum is: %i",
						int(i), int(fModules.size())-1);
		throw std::invalid_argument(buf);
	}
}

void DemandDetectorConstruction::SetTargetMaterial(const G4String& mat)
{
	static G4Element *D_ = nullptr;
	if(!D_) {
		D_ = new G4Element("D","D",1);
		G4Isotope* isotope = new G4Isotope("D", 1,2,2.01410178*g/mole);
		D_->AddIsotope(isotope,1);
	}
	G4Element* Si_ = G4NistManager::Instance()->FindOrBuildElement(14);
	G4Element* C_  = G4NistManager::Instance()->FindOrBuildElement(6);
	G4Element* He_ = G4NistManager::Instance()->FindOrBuildElement(2);
	G4Element* H_  = G4NistManager::Instance()->FindOrBuildElement(1);

	if(mat == "CD2") {
		static G4Material* material_CD2 = nullptr;
		if(material_CD2 == nullptr) {
			G4double density =1.06*g/cm3;
			material_CD2 = new G4Material("CD2", density,2);
			material_CD2->AddElement(C_,1);
			material_CD2->AddElement(D_,2);
		}
		fTargetMaterial = material_CD2;
	}
	else if(mat == "sH2") { // solid hydrogen
		static G4Material* material_sH2 = nullptr;
		if(material_sH2 == nullptr) {
			G4double density = 0.086*g/cm3;
			material_sH2 = new G4Material("sH2", density, 1);
			material_sH2->AddElement(H_,1);
		}
		fTargetMaterial = material_sH2;
	}
	else if(mat == "sD2") { // solid deuterium
		static G4Material* material_sD2 = nullptr;
		if(material_sD2 == nullptr) {
			G4double density = 0.1967*g/cm3;
			material_sD2 = new G4Material("sD2", density, 1);
			material_sD2->AddElement(D_,1);
		}
		fTargetMaterial = material_sD2;
	}
	else if(mat == "SiHe") { // He-implanted Silicon
		static G4Material* material_SiHe = nullptr;
		if(material_SiHe == nullptr) {
			G4double density = 1.4824*g/cm3;
			material_SiHe = new G4Material("SiHe", density, 2);
			material_SiHe->AddElement(Si_,2);
			material_SiHe->AddElement(He_,1);
		}
		fTargetMaterial = material_SiHe;
	}
	else if(mat == "He") { // He gas
		static G4Material* material_He = nullptr;
		if(material_He == nullptr) {
			G4double pressure = 4.9; // Torr
			G4double density = ((pressure*133.322*4.002602)/(8.314463*298.15*1e6))*g/cm3;
			material_He = new G4Material("He", density, 1);
			material_He->AddElement(He_,1);
		}
		fTargetMaterial = material_He;
	}
	else {
		char buf[4096];
		sprintf(buf, "Invalid material: \"%s\".\n"
						"Valid entries are \"CD2\" or \"sH2\".",
						mat.c_str());
		throw std::invalid_argument(buf);
	}
}


void DemandDetectorConstruction::AddModule()
{
	DemandDetectorConstruction::Module_t module;
	module.m_Dist = 1*m;
	module.m_Theta = 0.;
	module.m_Phi   = 0.;
	module.m_Voxnum = {10,10,1};
	module.m_Voxsize = {20,20,20};
	module.m_Threshold = 100*keV;
	module.m_ReadoutType = kCube;
	module.m_Rotate = true;
	module.m_Position.set(0,0,0);
	module.m_HavePosition = {0,0,0};
	fModules.emplace_back(module);
}

namespace
{
void solid_vis (G4LogicalVolume* lv, G4Color color, double alpha = 1)
{
	color.SetAlpha(alpha);
	auto vis = new G4VisAttributes(color);
	
	vis->SetVisibility(true);
	vis->SetForceSolid(true);      // THIS is the key
// vis->SetForceWireframe(true); // alternative
// vis->SetForceAuxEdgeVisible(true);
	lv->SetVisAttributes(vis);
};
void solid_vis (G4LogicalVolume* lv, double alpha = 1)
{
	solid_vis(lv,G4Color(0.2,0.6,0.9),alpha);
}

G4Colour aluminumGrey(0.80, 0.80, 0.82, 1.0);
G4Colour nearBlack(0.50, 0.50, 0.50, 1.0);
G4Colour scintClear(0.75, 0.85, 1.00, 1.0);

} // namespace


G4LogicalVolume* DemandDetectorConstruction::ConstructS2230CasingAssembly()
{	
	// -------------------------
	// Build ONE casing "assembly" LV that contains:
	// - shell (Al subtraction solid)
	// - inner air cavity (box) as sibling of shell
	// - PMT + collars + scint inside inner air
	// This assembly LV is then placed 8x in moduleLV.
	// -------------------------
	// Outer box full lengths (x,y,z)
	const G4double boxDims[3] = { 160.*mm, 43.45*mm, 40.*mm };
	const G4double wall_thickness = 3.175*mm;

	// Inner (air cavity) full lengths
	G4double boxDimsInner[3];
	for (int i = 0; i < 3; ++i) {
		boxDimsInner[i] = boxDims[i] - 2.0*wall_thickness;
		G4cout << "boxDimsInner["<<i<<"]: " << boxDimsInner[i]/mm << " mm\n";
	}
	 

	// Assembly container is the *outer envelope* (air, invisible)
	auto* casingAssemblySolid = new G4Box("CasingAssemblySolid",
																				0.5*boxDims[0], 0.5*boxDims[1], 0.5*boxDims[2]);
	auto* casingAssemblyLV = new G4LogicalVolume(casingAssemblySolid, fMatAir, "CasingAssemblyLV");
	casingAssemblyLV->SetVisAttributes(G4VisAttributes::GetInvisible());

	// ---- Shell solid (Al): outer minus inner ----
	auto* outerBox = new G4Box("OuterBox",
														 0.5*boxDims[0], 0.5*boxDims[1], 0.5*boxDims[2]);

	// Use a cutter for subtraction ONLY; do not re-use this for placed air.
	auto* innerBoxCut = new G4Box(
		"InnerBoxCut",
		0.5*boxDimsInner[0],
		0.5*boxDimsInner[1],
		0.5*boxDimsInner[2]
		);

	auto* shellSolid_nocutout = new G4SubtractionSolid("AlShell_nocutout", outerBox, innerBoxCut);

	// Cutout at top (a cylindrical notch). This is used only to subtract from shell.
	const G4double cutout_diam = 20*mm;
	auto* cutoutSolid = new G4Tubs("cutoutSolid",
																 0.0,
																 0.5*cutout_diam,
																 0.5*wall_thickness + 1.0*mm, // oversize ok for cutter
																 0.0*deg,
																 360.0*deg);

	auto* rotCutout = new G4RotationMatrix();
	rotCutout->rotateY(90*deg);

	const G4double posCutout = 0.5*boxDimsInner[0] + 0.5*wall_thickness; // along +X

	auto* shellSolid = new G4SubtractionSolid(
		"AlShell", shellSolid_nocutout, cutoutSolid,
		rotCutout, G4ThreeVector(posCutout, 0, 0)
		);

	auto* shellLV = new G4LogicalVolume(shellSolid, fMatAl, "AlShellLV");
	solid_vis(shellLV, aluminumGrey, 0.2);

	// Place shell inside casing assembly
	new G4PVPlacement(nullptr, G4ThreeVector(), shellLV, "AlShellPV",
										casingAssemblyLV, false, 0, true);

	// ---- Inner air cavity (placed) ----
	// IMPORTANT: this must be a sibling of the shell (same mother), NOT a child of shellLV.
	// Shrink slightly to avoid coplanar-face / tolerance overlaps.
	const G4double eps = 0.05*mm;

	auto* innerAirSolid = new G4Box("InnerAirSolid",
																	0.5*boxDimsInner[0] - eps,
																	0.5*boxDimsInner[1] - eps,
																	0.5*boxDimsInner[2] - eps);

	auto* airLV = new G4LogicalVolume(innerAirSolid, fMatAir, "airLV");

	new G4PVPlacement(nullptr, G4ThreeVector(), airLV, "airPV",
										casingAssemblyLV, false, 0, true);

	// -------------------------
	// Internals go inside airLV
	// -------------------------

	// PMT: hollow cylinder (mu-metal)
	const G4double pmt_len = 120*mm;
	const G4double pmt_diam_outer = 31*mm;
	const G4double pmt_diam_inner = 26*mm;
//    const G4double pmt_wall_thick = 0.8*mm; // wrong

	auto* pmtSolid = new G4Tubs("PMTSolid",
															0.5*pmt_diam_inner, // inner radius
															0.5*pmt_diam_outer, // outer radius
															0.5*pmt_len,        // half-length (local Z axis)
															0.0*deg,
															360.0*deg);

	auto* pmtLV = new G4LogicalVolume(pmtSolid, fMatPMTGlass, "PMTLV");
	solid_vis(pmtLV, G4Color(0,0,1), 0.5);

	// rotate cylinder so axis is along +X (local)
	auto* pmtRot = new G4RotationMatrix();
	pmtRot->rotateY(90*deg);

	// place PMT in air: "butt" against +X end of inner cavity
	// Note this is inconsistent with CAD files, which show the PMT
	// arbitrarily placed. In reality they are butted against the
	// far edge of the casing
	double gap = 0.1*mm; // to get rid of overlap warning
	const G4double pmtPos = 0.5*boxDimsInner[0] - 0.5*pmt_len - gap;
	G4cout << "PMT position: " << pmtPos << G4endl;
		
	new G4PVPlacement(pmtRot,
										G4ThreeVector(pmtPos, 0, 0),
										pmtLV,
										"PMTPV",
										airLV,
										false,
										0,
										true);

	// Glass PMT->scint interface
	const G4double glassDiam = 26*mm;
	const G4double glassThickness = 1*mm;
	auto* glassSolid = new G4Tubs("GlassSolid",
																0, // inner radius
																0.5*glassDiam, // outer radius
																0.5*glassThickness, // half-length (local Z axis)
																0.0*deg, 360.0*deg
		);
	auto* glassLV = new G4LogicalVolume(glassSolid, fMatPMTGlass, "GLASSLV");
	solid_vis(pmtLV, G4Color(0,0,1), 0.8);
//		solid_vis(glassLV, aluminumGrey);
	new G4PVPlacement(pmtRot,
										G4ThreeVector(pmtPos - 0.5*pmt_len - 0.1*mm),
										glassLV, "GLASSPV", airLV, false, 0, true
		);

	// Mu Metal shielding
	const G4double muMetalThickness = 0.8*mm;
	const G4double muMetalDiam = pmt_diam_outer + 0.02*mm;
	auto* muMetalSolid = new G4Tubs("MuMetalSolid",
																	0.5*muMetalDiam, // inner R
																	0.5*muMetalDiam + muMetalThickness, // outer R
																	0.5*pmt_len,
																	0*deg, 360*deg
		);
	auto* muMetalLV = new G4LogicalVolume(muMetalSolid, fMatMuMetal, "MU_METAL_LV");
	solid_vis(muMetalLV, nearBlack);
	new G4PVPlacement(pmtRot,
										G4ThreeVector(pmtPos,0,0),
										muMetalLV, "MuMetalPV", airLV, false, 0, true
		);																		

	// PMT collars (Al): square plate with circular hole (aligned to PMT axis)
	// In CAD, collar is 37x37x3.175.  But it actually fits into a notch in the
	// outer casing.  So for GEANT4, make the collar the same size as it's mother
	// airLV in x/y.  This will avoid overlap issues and give the same mechanical
	// geometry
	//
	// auto* innerAirSolid = new G4Box("InnerAirSolid",
	//                             0.5*boxDimsInner[0] - eps, --> long axis
	//                             0.5*boxDimsInner[1] - eps,
	//                             0.5*boxDimsInner[2] - eps);
	// correction orientation is [2,1,thick] for collar
	const G4double collarThick = 3.175*mm;
	const G4double collarDims[3] = {
		boxDimsInner[2]-eps*2, boxDimsInner[1]-eps*2, collarThick
	};

	auto* collarOuter = new G4Box("collarOuter",
																0.5*collarDims[0],
																0.5*collarDims[1],
																0.5*collarDims[2]);

	const G4double collarHoleDiam = muMetalDiam + 2*muMetalThickness; // 31.8*mm;
	auto* collarHole = new G4Tubs("collarHole",
																0.0,
																0.5*collarHoleDiam + 0.01*mm,
																0.5*collarThick + 1.0*mm, // oversize cutter
																0.0*deg,
																360.0*deg);

	auto* collarSolid = new G4SubtractionSolid("collarSolid", collarOuter, collarHole);
	auto* collarLV = new G4LogicalVolume(collarSolid, fMatAl, "collarLV");
	solid_vis(collarLV, G4Color(1,0,0), 1);

	// collar positions (your CAD-based numbers)
	// far end of right collar is 56.225 mm from centre
	// far end of left collar is 25.95 mm from centre
	std::vector<G4double> collar_pos = { //{ (87.6-32)*mm, (166.6-32)*mm };
		-25.95*mm + 0.5*collarThick,
		+56.225*mm - 0.5*collarThick
	}; 
	for (size_t icol = 0; icol < collar_pos.size(); ++icol) {
		std::stringstream sstr; sstr << "collarPV_" << icol;

		new G4PVPlacement(pmtRot,
											G4ThreeVector(collar_pos[icol], 0, 0),
											collarLV,
											sstr.str().c_str(),
											airLV,
											false,
											(int)icol,
											true);
	}

	// OGS scintillator: 30 mm cube
	const G4double scintDims[3] = { 30*mm, 30*mm, 30*mm };
	auto scintBox = new ScintillatorBox(
		"scintBox",
		0.5*scintDims[0],
		0.5*scintDims[1],
		0.5*scintDims[2]);
	scintBox->SetReadoutType(DemandDetectorConstruction::kCube);


	auto* scintLV = new G4LogicalVolume(scintBox, GetScintillatorMaterial(), "DEMAND_scintLV");
	solid_vis(scintLV, scintClear, 1);

	// butt against PMT on -X side
	const G4double scintPos = pmtPos - 0.5*pmt_len - glassThickness - 0.5*scintDims[0];
	G4cout << "ScintPos: " << scintPos << G4endl;
		
	new G4PVPlacement(nullptr,
										G4ThreeVector(scintPos, 0, 0),
										scintLV,
										"DEMAND_scintPV",
										airLV,
										false,
										0,
										true);

	return casingAssemblyLV;
}

void DemandDetectorConstruction::PlaceCasingAssembly(
	G4LogicalVolume* casingAssemblyLV,
	G4LogicalVolume* motherLV,
	const G4ThreeVector& pos)
{
	// Place 8x casings
	const G4double dS = 112*mm;
	const G4double dT = 44.5*mm;

	std::vector<G4double> casing_X = {
		-dS, -dS, 0,  dS,  dS,  dT, 0,  -dT
	};
	std::vector<G4double> casing_Y = {
		0,   dT, dS, dT,  0,  -dS, -dS, -dS
	};
	std::vector<G4double> casing_rot = {
		180*deg, 180*deg, 270*deg, 0*deg, 0*deg, 90*deg, 90*deg, 90*deg
	};
	
	for (size_t i=0; i<casing_rot.size(); ++i) {
		auto* casingRot = new G4RotationMatrix();
		casingRot->rotateZ(casing_rot[i]);
		G4String PVname = "CasingPV_" + std::to_string(i);

		new G4PVPlacement(
			casingRot, pos + G4ThreeVector(
				casing_X[i], casing_Y[i]),
			casingAssemblyLV, PVname, motherLV, false, (int)i, true
			);
	}
}

void DemandDetectorConstruction::ConstructS2230Detectors(G4VPhysicalVolume* mother_phys)
{
	auto motherLV = mother_phys->GetLogicalVolume();

	// Import auto-generated code which creates the outer "skeleton frame"
	// using G4Extruded solids.
	//
	// Code was generated by exporting various CAD parts (in Onshape) as DFX files.
	// Then use a ChatGPT generated python script which parses the DFX files and
	// generates the C++ code to create the part as an extruded solid. The python
	// file can be found in the etc/dfx_to_geant4.py
	//
	// The code for each part is stored in it's own c++ file with the .inc extension.
	// These are located in src/skelframe and #included at the top of this file. Each
	// inc file contains a function Build<name>, which creates the logical volume for
	// part <name>.
	//
	// Here we call all of these functions to build the reqired logical volumes.
	// 
	auto BuildLV = [](G4LogicalVolume* (*f)(), G4Color color) {
		auto LV = f();	solid_vis(LV, G4Color(1,1,1,0.7));		return LV;
	};
	
	auto FrontFrameUpperLV = BuildLV(
		BuildFrontFrameUpper, G4Color(0,0,1));
	auto FrontFrameLowerLV = BuildLV(
		BuildFrontFrameLower, G4Color(0,0,1));
	auto BackFrameUpperLV = BuildLV(
		BuildBackFrameUpper, G4Color(0,0,1));
	auto BackFrameLowerLV = BuildLV(
		BuildBackFrameLower, G4Color(0,0,1));
	auto H1LV = BuildLV(
		BuildH1, G4Color(0,1,0));
	auto H2LV = BuildLV(
		BuildH2, G4Color(1,1,0));
	auto H3LV = BuildLV(
		BuildH3, G4Color(0,1,1));
	auto H4LV = BuildLV(
		BuildH4, G4Color(1,0,1));
	auto SideHolesLV = BuildLV(
		BuildSideHoles, G4Color(0,1,0));
	auto TopHolesLV = BuildLV(
		BuildTopHoles, G4Color(0,1,0));
	auto BottomHolesLV = BuildLV(
		BuildBottomHoles, G4Color(0,1,0));
	auto VerticalSupportLV = BuildLV(
		BuildVerticalSupport, G4Color(1,0,1));
	auto VerticalSupportTopLV = BuildLV(
		BuildVerticalSupportTop, G4Color(1,0,1));
	auto HorizontalSupportUpperLV = BuildLV(
		BuildHorizontalSupportUpper, G4Color(1,0,1));
	auto HorizontalSupportLowerLV = BuildLV(
		BuildHorizontalSupportLower, G4Color(1,0,1));
	
	// Also build the logical volume for the casing assembly.
	// This contains all detectors, PMTs, detector housing, etc.
	auto CasingAssemblyLV = ConstructS2230CasingAssembly();


	// Now place all of the detector parts (casing assembly + skeleton
	// frame) into the mother volume. The z-positions of each part are
	// gotten from CAD measurements, see below.	The entire assembly
	// is placed such that the 1" inner collar (parts H2 and H3) is flush
	// against the "PDE" collar defined in ugeom_gbox.cc.
	// 
  // CAD Measured Distances of Various Parts
	// These measurements can also be found online:
	// https://smuhalifax-my.sharepoint.com/:x:/g/personal/greg_christian_smu_ca/IQCNVfJwkPjdT7PywL-vrzEgAVjHRmvAEBj_J8378qOAJ3A?e=ZoUGjK
	//
	// We use calculated distance from PDE upstream face to center
	// of each part to place everything. This is based off downstream
	// H2/H3 (which define a 1" aperture) being flush against PDE.
	//
	// Note that the position and length of PDE is exposed in ugeom_gbox.cc
	//
  // Part | O --> Far	| O --> Near	| Thickness	| O-->Center	| PDE --> Center
  // Back Frame	23.925	20.75	3.175	22.3375	2.3875
  // Back H1	20.75	17.575	3.175	19.1625	-0.7875
  // Back H2	19.95	17.575	2.375	18.7625	-1.1875
  // Back H3	19.95	17.575	2.375	18.7625	-1.1875
  // Back H4	20.75	17.575	3.175	19.1625	-0.7875
  // Center Assembly	20	-20	40	0	-19.95
  // Front H4	-20.75	-17.575	-3.175	-19.1625	-39.1125
  // Front H3	-19.95	-17.575	-2.375	-18.7625	-38.7125
  // Front H2	-19.95	-17.575	-2.375	-18.7625	-38.7125
  // Front H1	-20.75	-17.575	-3.175	-19.1625	-39.1125
  // Front Frame	-23.925	-20.75	-3.175	-22.3375	-42.2875

	const G4double zPDEFace = DRAGON::PDE_zpos - DRAGON::PDE_zlen;
	
	const G4double zBackFrame  = zPDEFace + +2.3875*mm;
	const G4double zBackH1     = zPDEFace + -0.7875*mm;
	const G4double zBackH2     = zPDEFace + -1.1875*mm;
	const G4double zBackH3     = zPDEFace + -1.1875*mm;
	const G4double zBackH4     = zPDEFace + -0.7875*mm;
	const G4double zCenter     = zPDEFace + -19.95*mm;
	const G4double zFrontH4    = zPDEFace + -39.1125*mm;
	const G4double zFrontH3    = zPDEFace + -38.7125*mm;
	const G4double zFrontH2    = zPDEFace + -38.7125*mm;
	const G4double zFrontH1    = zPDEFace + -39.1125*mm;
	const G4double zFrontFrame = zPDEFace + -42.2875*mm;
	
	// Define some distances which are used for X and Y placements
	// of various frame parts.
	G4double eps_frm = 0.0001*mm;
	G4double frame_thickness = 3.175*mm;

	// Now place everything. First define some variables to set the
	// X and Y positions of each part (these are updated as needed
	// throughout the placement process).
	G4double partX = 0, partY = 0;
	
	// Back frame (upper and lower)
	new G4PVPlacement(
		nullptr,	G4ThreeVector(partX, partY, zBackFrame),
		BackFrameUpperLV, "BackFrameUpper_PV",
		motherLV, false, 0, true);
	
	new G4PVPlacement(
		nullptr,	G4ThreeVector(partX, partY, zBackFrame),
		BackFrameLowerLV, "BackFrameLower_PV",
		motherLV, false, 0, true);
	
	// Back H1
	new G4PVPlacement(
		nullptr, G4ThreeVector(partX,partY,zBackH1),
		H1LV, "BackH1_PV",
		motherLV, false, 0, true);

	// Back H2
	new G4PVPlacement(
		nullptr, G4ThreeVector(partX,partY,zBackH2),
		H2LV, "BackH2_PV",
		motherLV, false, 0, true);

	// Back H3
	new G4PVPlacement(
		nullptr, G4ThreeVector(partX,partY,zBackH3),
		H3LV, "BackH3_PV",
		motherLV,	false, 0, true);
	
	// Back H4
	new G4PVPlacement(
		nullptr, G4ThreeVector(partX,partY,zBackH4),
		H4LV, "BackH4_PV",
		motherLV,	false, 0, true);

	// Detector Housing Assembly
	PlaceCasingAssembly(
		CasingAssemblyLV, motherLV,
		G4ThreeVector(partX,partY,zCenter)
		);
	
	// Side Holes
	auto rotSideHole = new G4RotationMatrix;
	rotSideHole->rotateY(90*deg);

	// left
	partX = -(196.675 - frame_thickness/2 + eps_frm);
	new G4PVPlacement(
		rotSideHole, G4ThreeVector(partX,partY,zCenter),
		SideHolesLV, "SideHoleLeft_PV",
		motherLV,	false, 0, true);

	// right
	partX *= -1;
	new G4PVPlacement(
		rotSideHole, G4ThreeVector(partX,partY,zCenter),
		SideHolesLV, "SideHoleRight_PV",
		motherLV,	false, 0, true);

	// Top Holes
	auto rotTopHole = new G4RotationMatrix;
	rotTopHole->rotateX(90*deg);
	rotTopHole->rotateZ(90*deg);
	partX = 0;
	partY = 196.675*mm - frame_thickness/2 + eps_frm;
	new G4PVPlacement(
		rotTopHole, G4ThreeVector(partX,partY,zCenter),
		TopHolesLV, "TopHole_PV",
		motherLV, false, 0, true);

	// Bottom Holes
	partY = -(196.675*mm - frame_thickness/2 + eps_frm);
	new G4PVPlacement(
		rotTopHole, G4ThreeVector(partX,partY,zCenter),
		BottomHolesLV, "BottomHole_PV",
		motherLV, false, 0, true);

	// Vertical Supports (bottom)
	// Left
	partX = -(69.925*mm - frame_thickness/2 + eps_frm);
	partY = 0;
	new G4PVPlacement(
		rotSideHole, G4ThreeVector(partX,partY,zCenter),
		VerticalSupportLV, "VerticalSupportLeft_PV",
		motherLV, false, 0, true);

	// Right
	partX *= -1;
	new G4PVPlacement(
		rotSideHole, G4ThreeVector(partX,partY,zCenter),
		VerticalSupportLV, "VerticalSupportRight_PV",
		motherLV, false, 0, true);

	// Vertical Supports (top)
	// Left
	partX = -(25.425*mm - frame_thickness/2 + eps_frm);
	new G4PVPlacement(
		rotSideHole, G4ThreeVector(partX,partY,zCenter),
		VerticalSupportTopLV, "VerticalSupportTopLeft_PV",
		motherLV, false, 0, true);

	// Right
	partX *= -1;
	new G4PVPlacement(
		rotSideHole, G4ThreeVector(partX,partY,zCenter),
		VerticalSupportTopLV, "VerticalSupportTopRight_PV",
		motherLV, false, 0, true);

	// Horizontal Supports (upper)
	// Left
	partX = 0;
	partY = 69.925*mm - frame_thickness/2 + eps_frm;
	new G4PVPlacement(
		rotTopHole, G4ThreeVector(partX,partY,zCenter),
		HorizontalSupportUpperLV, "HorizontalSupportUpperLeft_PV",
		motherLV, false, 0, true);

	// Right
	auto rotHSupportRight = new G4RotationMatrix;
	rotHSupportRight->rotateX(90*deg);
	rotHSupportRight->rotateZ(-90*deg);
	new G4PVPlacement(
		rotHSupportRight, G4ThreeVector(partX,partY,zCenter),
		HorizontalSupportUpperLV, "HorizontalSupportUpperRight_PV",
		motherLV, false, 0, true);

	// Horizontal Supports (lower)
	// Left
	partX = 0;
	partY = -(25.425*mm - frame_thickness/2 + eps_frm);
	new G4PVPlacement(
		rotTopHole, G4ThreeVector(partX,partY,zCenter),
		HorizontalSupportLowerLV, "HorizontalSupportLowerLeft_PV",
		motherLV, false, 0, true);

	// Right
	new G4PVPlacement(
		rotHSupportRight, G4ThreeVector(partX,partY,zCenter),
		HorizontalSupportLowerLV, "HorizontalSupportLowerRight_PV",
		motherLV, false, 0, true);
	
	
	// Front H pieces and frame all centered in X, Y
	partX = partY = 0;
	
	// Front H1
	new G4PVPlacement(
		nullptr, G4ThreeVector(partX,partY,zFrontH1),
		H1LV, "FrontH1_PV",
		motherLV, false, 0, true);

	// Front H2
	new G4PVPlacement(
		nullptr, G4ThreeVector(partX,partY,zFrontH2),
		H2LV, "FrontH2_PV",
		motherLV, false, 0, true);

	// Front H3
	new G4PVPlacement(
		nullptr, G4ThreeVector(partX,partY,zFrontH3),
		H3LV, "FrontH3_PV",
		motherLV,	false, 0, true);
	
	// Front H4
	new G4PVPlacement(
		nullptr, G4ThreeVector(partX,partY,zFrontH4),
		H4LV, "FrontH4_PV",
		motherLV,	false, 0, true);

	// Front Frame (upper and lower)
	new G4PVPlacement(
		nullptr,	G4ThreeVector(partX, partY, zFrontFrame),
		FrontFrameUpperLV, "FrontFrameUpper_PV",
		motherLV, false, 0, true);
	
	new G4PVPlacement(
		nullptr,	G4ThreeVector(partX, partY, zFrontFrame),
		FrontFrameLowerLV, "FrontFrameLower_PV",
		motherLV, false, 0, true);	
}



/////////////////////////////
//  class SkeletonFrame    //
/////////////////////////////

namespace{
	G4Box* make_box(const G4String& name, const G4ThreeVector& halfDims){
		return new G4Box(name, halfDims.x(), halfDims.y(), halfDims.z());
	};
}

DemandDetectorConstruction::SkeletonFrame::SkeletonFrame():
	fSFthick(3.175*CLHEP::mm), fEps(0.001*CLHEP::mm), fEps3(fEps,fEps,fEps)
{
	ConstructFrame();
	ConstructH();
	fMatAl = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
}


void DemandDetectorConstruction::SkeletonFrame::ConstructFrame()
{
	// Rectangular Al plate (holes cut out of this)
	const G4ThreeVector frameDims(387*mm, 387*mm, fSFthick);
	auto* sfSolid = new G4Box(
		"SkeletonFrameSolid", 0.5*frameDims[0], 0.5*frameDims[1], 0.5*frameDims[2]
		);
	
	///////////////////////
	// Cutout outer edges//
	///////////////////////
	// Bottom notches
	const G4ThreeVector cutoutDimsBottom(120.25*mm, 164.75*mm, fSFthick);
	auto* cutoutSolidBottom = make_box(
		"OuterCutoutBottom", 0.5*cutoutDimsBottom+fEps3
		);

	const G4ThreeVector cutoutPosBottom(
		+0.5*frameDims[0] - 0.5*cutoutDimsBottom[0],
		-0.5*frameDims[1] + 0.5*cutoutDimsBottom[1],
		+0
		);
	auto* sfWithCutoutBottom = new G4SubtractionSolid(
		"SFWithBottomCutout",
		new G4SubtractionSolid("SFWithBottomCutout1",
													 sfSolid, cutoutSolidBottom,
													 nullptr, cutoutPosBottom
			),
		cutoutSolidBottom, nullptr, G4ThreeVector(
			-cutoutPosBottom.x(), cutoutPosBottom.y(), cutoutPosBottom.z()
			) );

	// Top notches
	const G4ThreeVector cutoutDimsTop(164.75*mm, 120.25*mm, fSFthick);
	auto* cutoutSolidTop = make_box(
		"OuterCutoutTop", 0.5*cutoutDimsTop+fEps3
		);

	const G4ThreeVector cutoutPosTop(
		+0.5*frameDims[0] - 0.5*cutoutDimsTop[0],
		+0.5*frameDims[1] - 0.5*cutoutDimsTop[1],
		+0
		);
	auto* sfWithTopCutout = new G4SubtractionSolid(
		"SFWithTopCutout",
		new G4SubtractionSolid("SFWithTopCutout1",
													 sfWithCutoutBottom, cutoutSolidTop,
													 nullptr, cutoutPosTop
			),
		cutoutSolidTop, nullptr, G4ThreeVector(
			-cutoutPosTop.x(), cutoutPosTop.y(), cutoutPosTop.z()
			) );

	//////////////////////
	//Cutout Inner Holes//
	//////////////////////

	// Left/Right
	const G4ThreeVector cutoutDimsLR(149*mm, 76*mm, fSFthick);
	auto* cutoutSolidLR = make_box(
		"InnerCutoutLR", 0.5*cutoutDimsLR+fEps3
		);

	const G4ThreeVector cutoutPosLR(
		+0.5*frameDims[0] - 0.5*cutoutDimsLR[0] - 6.5*mm,
		0.5*cutoutDimsLR[1] - 15.75*mm,
		+0
		);
	auto* sfWithLRCutout = new G4SubtractionSolid(
		"SFWithLRCutout",
		new G4SubtractionSolid("SFWithTopCutout1",
													 sfWithTopCutout, cutoutSolidLR,
													 nullptr, cutoutPosLR
			),
		cutoutSolidLR, nullptr, G4ThreeVector(
			-cutoutPosLR.x(), cutoutPosLR.y(), cutoutPosLR.z()
			) );

	// Top/bottom
	const G4ThreeVector cutoutDimsInnerTop(31*mm, 149*mm, fSFthick);
	const G4ThreeVector cutoutDimsInnerBottom(120.5*mm, 149*mm, fSFthick);
	auto* cutoutInnerTop = make_box(
		"InnerCutoutTop", 0.5*cutoutDimsInnerTop+fEps3
		);
	auto* cutoutInnerBottom = make_box(
		"InnerCutoutBottom", 0.5*cutoutDimsInnerBottom+fEps3
		);
	const G4ThreeVector cutoutPosInnerTop(
		0, +0.5*frameDims[1] - 0.5*cutoutDimsInnerTop[1] - 6.5*mm, 0
		);
	const G4ThreeVector cutoutPosInnerBottom(
		0, -0.5*frameDims[1] + 0.5*cutoutDimsInnerBottom[1] + 6.5*mm, 0
		);
	auto* sfWithInnerTopCutout = new G4SubtractionSolid(
		"SFWithInnerTopCutout",
		sfWithLRCutout, cutoutInnerTop, nullptr, cutoutPosInnerTop
		);
	auto* sfWithInnerCutouts = new G4SubtractionSolid(
		"SFWithInnerCutouts",
		sfWithInnerTopCutout, cutoutInnerBottom, nullptr, cutoutPosInnerBottom
		);


	// Center hole
	// --> Not complete!
	const G4double centerHoleDiam = 42*mm;
	auto* centerHole = new G4Tubs(
		"centerHole",
		0.0,
		0.5*centerHoleDiam,
		0.5*fSFthick + fEps,
		0, 360*deg
		);
	auto* sfWithCenterHole = new G4SubtractionSolid(
		"SFWithCenterHole", sfWithInnerCutouts, centerHole, nullptr, G4ThreeVector(0,0,0)
		);	
	
	fFrameLV = new G4LogicalVolume(sfWithCenterHole, fMatAl, "SkeletonFrameLV");
	solid_vis(fFrameLV, G4Color(1,0,0), 1);
};

void DemandDetectorConstruction::SkeletonFrame::ConstructH()
{
	// Full plate (cut out from this)
	const G4ThreeVector HDims(132.5*mm, 98.25*mm, fSFthick);
	auto* HSolid = make_box("HSolid",0.5*HDims);

	// cutout top
	const G4ThreeVector HCutoutTopDims(44.5*mm, 34.25*mm, fSFthick+fEps);
	auto* HCutoutTop = make_box("HCutoutTop",HCutoutTopDims/2);
	const G4ThreeVector HCutoutTopPos(0,HDims[1]/2-HCutoutTopDims[1]/2,0);
	auto* HCutoutTopSolid = new G4SubtractionSolid(
		"HWithCutoutTop", HSolid, HCutoutTop, nullptr, HCutoutTopPos
		);

	// cutout center hole
	const double holeCutoutDiam = 25.5*mm;
	fHoleYpos = 32*mm - HDims[1]/2;
	auto* HCutoutHole = new G4Tubs(
		"HCutoutHole", 0, 0.5*holeCutoutDiam, fSFthick+fEps, 0, 360*deg
		);
	auto* HCutoutHoleSolid = new G4SubtractionSolid(
		"HWithCutoutHole", HCutoutTopSolid, HCutoutHole,
		nullptr, G4ThreeVector(0,fHoleYpos,0)
		);

	// left/right cutouts
	const G4ThreeVector HLRCutoutDims(34.25*mm, HDims[1]-6.5*mm,fSFthick+fEps);
	auto* HLRCutout = make_box("HLRCutout",HLRCutoutDims/2+fEps3);
	const G4ThreeVector HLRCutoutPos(
		HDims[0]/2-HLRCutoutDims[0]/2, HDims[1]/2-HLRCutoutDims[1]/2, 0
		);
	auto* HLRCutoutSolid = new G4SubtractionSolid(
		"HLRCutout", new G4SubtractionSolid(
			"HLRCutout1", HCutoutHoleSolid, HLRCutout, nullptr, HLRCutoutPos),
		HLRCutout, nullptr, G4ThreeVector(
			-HLRCutoutPos[0], HLRCutoutPos[1], HLRCutoutPos[2]
			) );
	
	// Logical Volume
	fHLV = new G4LogicalVolume(
		HLRCutoutSolid,fMatAl,"SkeletonFrame_H_LV"
		);
	solid_vis(fHLV,G4Color(1,0,0),1);
}




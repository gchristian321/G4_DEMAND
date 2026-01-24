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

using namespace std;

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
	
	G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;  
  G4double density; 

  // Vacuum
  auto defaultMaterial = new G4Material(
		"Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
		kStateGas, 2.73*kelvin, 3.e-18*pascal);

  G4NistManager *nist = G4NistManager::Instance();
  G4Material *matAir = nist->FindOrBuildMaterial("G4_AIR");
  G4Material *matSteel = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  G4Material *matAl = nist->FindOrBuildMaterial("G4_Al");

  // //OGS
  // G4Material *matOGS = new G4Material("matOGS", 1.096*g/cm3, 3);
  // matOGS->AddElement(nist->FindOrBuildElement("H"),36);
  // matOGS->AddElement(nist->FindOrBuildElement("C"),42);	
  // matOGS->AddElement(nist->FindOrBuildElement("Si"),1);


   
  //     
  // World
  //
  auto worldS 
    = new G4Box("WRLD",           // its name
                 worldSizeXY/2, worldSizeXY/2, worldSizeZ/2); // its size
                         
  auto worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 defaultMaterial,  // its material
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
		G4LogicalVolume *logicOuterChamber = new G4LogicalVolume(solidOuterChamber, matAl, "logicOuterChamber");
		G4VPhysicalVolume *physOuterChamber = new G4PVPlacement(0, G4ThreeVector(0.,-96.91*mm,0.), logicOuterChamber, "physOuterChamber", worldLV, false, 0, true);

		G4Box *solidInnerChamber = new G4Box("solidInnerChamber", 22.23*mm, 125.48*mm, 82.555*mm);
		G4LogicalVolume *logicInnerChamber = new G4LogicalVolume(solidInnerChamber, defaultMaterial, "logicInnerChamber");
		G4VPhysicalVolume *physInnerChamber = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), logicInnerChamber, "physInnerChamber", logicOuterChamber, false, 0, true);

		G4Tubs *solidExitHole = new G4Tubs("solidExitHole", 0.0, 4.5*mm, 1.585*mm, 0.0, 360.0);
		G4LogicalVolume *logicExitHole = new G4LogicalVolume(solidExitHole, defaultMaterial, "logicExitHole");
		G4VPhysicalVolume *physExitHole = new G4PVPlacement(0, G4ThreeVector(0., 96.91*mm, 84.14*mm), logicExitHole, "physExitHole", logicOuterChamber, false, 0, true);

		//
		// Beam Pipe Collar
		//
		G4Tubs *solidCollar = new G4Tubs("solidCollar", 4.5*mm, 19.02*mm, 4*mm, 0., 360.);
		G4LogicalVolume *logicCollar = new G4LogicalVolume(solidCollar, matAl, "logicCollar");
		G4VPhysicalVolume *physCollar = new G4PVPlacement(0, G4ThreeVector(0.,0.,89.725*mm), logicCollar, "physCollar", worldLV, false, 0, true);

		//
		// Beam Pipe, with collar
		//
		G4Tubs *solidBeamPipe = new G4Tubs("solidBeamPipe", 4.5*mm, 12.68*mm, 21.3875*mm, 0., 360.);
		G4LogicalVolume *logicBeamPipe = new G4LogicalVolume(solidBeamPipe, matAl, "logicBeamPipe");
		G4VPhysicalVolume *physBeamPipe = new G4PVPlacement(0, G4ThreeVector(0., 0., 115.1125*mm), logicBeamPipe, "physBeamPipe", worldLV, false, 0, true);

		//
		// Beam pipe, no collar
		//
		//G4Tubs *solidBeamPipe = new G4Tubs("solidBeamPipe", 4.5*mm, 12.68*mm, 25.3875*mm, 0., 360.);
		//G4LogicalVolume *logicBeamPipe = new G4LogicalVolume(solidBeamPipe, matAl, "logicBeamPipe");
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
	
			G4LogicalVolume* logicEnclosure = new G4LogicalVolume(Enclosure, matAl, "logicEnclosure");

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
		// Fit on GAC OneDrive SMU/projects/DRAGON-alpha_n/S2230/NeutronSinglesAnalysis_sydney/OGS - Quenching.ipynb
		material->GetIonisation()->SetBirksConstant(0.104222 * mm/MeV);
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

namespace {
void solid_vis (G4LogicalVolume* lv, G4Color color, double alpha = 1)
{
	color.SetAlpha(alpha);
	auto vis = new G4VisAttributes(color);
	
	vis->SetVisibility(true);
	vis->SetForceSolid(true);        // THIS is the key
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


void DemandDetectorConstruction::ConstructS2230Detectors(G4VPhysicalVolume* mother_phys)
{
  auto* nist  = G4NistManager::Instance();
  G4Material* matAir = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* matAl  = nist->FindOrBuildMaterial("G4_Al");

  // Mu-metal approximation (80% Ni, 20% Fe) for PMT shield
  auto* elNi = nist->FindOrBuildElement("Ni");
  auto* elFe = nist->FindOrBuildElement("Fe");
  auto* muMetal = new G4Material("MuMetal_80Ni20Fe", 8.7*g/cm3, 2);
  muMetal->AddElement(elNi, 0.80);
  muMetal->AddElement(elFe, 0.20);

  // -------------------------
  // Module container (entire assembly)
  // CAD extents ~ 393 x 393 x 47.8 mm^3 (x,y,z)
  // -------------------------
  auto* moduleSolid = new G4Box("ModuleSolid", 0.5*393*mm, 0.5*393*mm, 0.5*47.8*mm);
  auto* moduleLV    = new G4LogicalVolume(moduleSolid, matAir, "S2230_ModuleLV");
  moduleLV->SetVisAttributes(G4VisAttributes::GetInvisible());

  // -------------------------
  // Build ONE casing "assembly" LV that contains:
  // - shell (Al subtraction solid)
  // - inner air cavity (box) as sibling of shell
  // - PMT + collars + scint inside inner air
  // This assembly LV is then placed 8x in moduleLV.
  // -------------------------
  auto ConstructCasingAssemblyLV = [&]() -> G4LogicalVolume*
  {
    // Outer box full lengths (x,y,z)
    const G4double boxDims[3] = { 160.*mm, 43.4*mm, 40.*mm };
    const G4double wall_thickness = 3.2*mm;

    // Inner (air cavity) full lengths
    G4double boxDimsInner[3];
    for (int i = 0; i < 3; ++i) boxDimsInner[i] = boxDims[i] - 2.0*wall_thickness;

    // Assembly container is the *outer envelope* (air, invisible)
    auto* casingAssemblySolid = new G4Box("CasingAssemblySolid",
                                          0.5*boxDims[0], 0.5*boxDims[1], 0.5*boxDims[2]);
    auto* casingAssemblyLV = new G4LogicalVolume(casingAssemblySolid, matAir, "CasingAssemblyLV");
    casingAssemblyLV->SetVisAttributes(G4VisAttributes::GetInvisible());

    // ---- Shell solid (Al): outer minus inner ----
    auto* outerBox = new G4Box("OuterBox",
                               0.5*boxDims[0], 0.5*boxDims[1], 0.5*boxDims[2]);

    // Use a cutter for subtraction ONLY; do not re-use this for placed air.
    auto* innerBoxCut = new G4Box("InnerBoxCut",
                                  0.5*boxDimsInner[0],
                                  0.5*boxDimsInner[1],
                                  0.5*boxDimsInner[2]);

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

    auto* shellSolid = new G4SubtractionSolid("AlShell",
                                              shellSolid_nocutout,
                                              cutoutSolid,
                                              rotCutout,
                                              G4ThreeVector(posCutout, 0, 0));

    auto* shellLV = new G4LogicalVolume(shellSolid, matAl, "AlShellLV");
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

    auto* airLV = new G4LogicalVolume(innerAirSolid, matAir, "airLV");

    new G4PVPlacement(nullptr, G4ThreeVector(), airLV, "airPV",
                      casingAssemblyLV, false, 0, true);

    // -------------------------
    // Internals go inside airLV
    // -------------------------

    // PMT: hollow cylinder (mu-metal)
    const G4double pmt_len = 120*mm;
    const G4double pmt_diam = 31*mm;
    const G4double pmt_wall_thick = 0.8*mm;

    auto* pmtSolid = new G4Tubs("PMTSolid",
                                0.5*pmt_diam - pmt_wall_thick, // inner radius
                                0.5*pmt_diam,                  // outer radius
                                0.5*pmt_len,                   // half-length (local Z axis)
                                0.0*deg,
                                360.0*deg);

    auto* pmtLV = new G4LogicalVolume(pmtSolid, muMetal, "PMTLV");
    solid_vis(pmtLV, nearBlack);

    // rotate cylinder so axis is along +X (local)
    auto* pmtRot = new G4RotationMatrix();
    pmtRot->rotateY(90*deg);

    // place PMT in air: "butt" against +X end of inner cavity
		double gap = 0.1*mm; // to get rid of overlap warning
    const G4double pmtPos = 0.5*boxDimsInner[0] - 0.5*pmt_len - gap;

    new G4PVPlacement(pmtRot,
                      G4ThreeVector(pmtPos, 0, 0),
                      pmtLV,
                      "PMTPV",
                      airLV,
                      false,
                      0,
                      true);

    // PMT collars (Al): square plate with circular hole (aligned to PMT axis)
    const G4double collarDims[3] = { 33.4*mm, 33.4*mm, 3.2*mm };

    auto* collarOuter = new G4Box("collarOuter",
                                  0.5*collarDims[0],
                                  0.5*collarDims[1],
                                  0.5*collarDims[2]);

    auto* collarHole = new G4Tubs("collarHole",
                                  0.0,
                                  0.5*pmt_diam + 0.1*mm,
                                  0.5*collarDims[2] + 1.0*mm, // oversize cutter
                                  0.0*deg,
                                  360.0*deg);

    auto* collarSolid = new G4SubtractionSolid("collarSolid", collarOuter, collarHole);
    auto* collarLV = new G4LogicalVolume(collarSolid, matAl, "collarLV");
    solid_vis(collarLV, aluminumGrey, 0.2);

    // collar positions (your CAD-based numbers)
    std::vector<G4double> collar_pos = { (87.6-32)*mm, (166.6-32)*mm };
    for (size_t icol = 0; icol < collar_pos.size(); ++icol) {
      std::stringstream sstr; sstr << "collarPV_" << icol;

      new G4PVPlacement(pmtRot,
                        G4ThreeVector(collar_pos[icol] - 0.5*boxDims[0], 0, 0),
                        collarLV,
                        sstr.str().c_str(),
                        airLV,
                        false,
                        (int)icol,
                        true);
    }

    // OGS scintillator: 30 mm cube
    const G4double scintDims[3] = { 30*mm, 30*mm, 30*mm };

    // auto* scintBox = new G4Box("scintBox",
    //                            0.5*scintDims[0],
    //                            0.5*scintDims[1],
    //                            0.5*scintDims[2]);
		auto scintBox = new ScintillatorBox(
			"scintBox",
			0.5*scintDims[0],
			0.5*scintDims[1],
			0.5*scintDims[2]);
		scintBox->SetReadoutType(DemandDetectorConstruction::kCube);


    auto* scintLV = new G4LogicalVolume(scintBox, GetScintillatorMaterial(), "DEMAND_scintLV");
    solid_vis(scintLV, scintClear, 1);

    // butt against PMT on -X side
    const G4double scintPos = pmtPos - 0.5*pmt_len - 0.5*scintDims[0];
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
  };

  // Build the casing assembly once
  auto* casingAssemblyLV = ConstructCasingAssemblyLV();

  // -------------------------
  // Place 8 casings into moduleLV
  // -------------------------
  const G4double dS = 112*mm;
  const G4double dT = 44.5*mm;

  std::vector<G4double> casing_X = { -dS, -dS, 0,  dS,  dS,  dT, 0,  -dT };
  std::vector<G4double> casing_Y = {  0,   dT, dS, dT,  0,  -dS, -dS, -dS };
  std::vector<G4double> casing_rot = { 180*deg, 180*deg, 270*deg, 0*deg, 0*deg, 90*deg, 90*deg, 90*deg };

  for (size_t i = 0; i < casing_X.size(); ++i) {
    auto* casingRot = new G4RotationMatrix();
    casingRot->rotateZ(casing_rot[i]);

    auto casingPos = G4ThreeVector(casing_X[i], casing_Y[i], 0);

    std::stringstream sstr; sstr << "CasingPV_" << i;
    new G4PVPlacement(casingRot,
                      casingPos,
                      casingAssemblyLV,
                      sstr.str().c_str(),
                      moduleLV,
                      false,
                      (int)i,
                      true);
  }

  // -------------------------
  // Place module into world (mother_phys)
  // -------------------------
  const G4double module_zpos = 11.312*cm;
  new G4PVPlacement(nullptr,
                    G4ThreeVector(0,0,module_zpos),
                    moduleLV,
                    "S2230_ModulePV",
                    mother_phys->GetLogicalVolume(),
                    false,
                    0,
                    true);

  // NOTE on IDs:
  // - each casing instance has copyNo = i at the CasingPV level
  // - inside each casing, PMT/collars/scint have their own copyNos
  // In a SensitiveDetector you can use touchable history:
  //   touch->GetVolume(0)->GetCopyNo()  (scint copy, etc)
  //   touch->GetVolume(1)->GetCopyNo()  (airPV copy: will be 0 here)
  //   touch->GetVolume(2)->GetCopyNo()  (CasingPV_i copy: i)  <-- useful
	
#if 0
	///
	/// Construct outer wire frame using exported
	/// STL file and CADMesh.
	auto mesh = CADMesh::TessellatedMesh::FromSTL("../TDE3468.stl");

// If the STL numbers are mm (very likely here):
	mesh->SetScale(mm);

// If it turns out the STL numbers are inches, use:
// mesh->SetScale(25.4*mm);

	auto frameSolid = mesh->GetSolid();
	auto frameLV = new G4LogicalVolume(frameSolid, matAl, "FrameLV");

// vis (readable on black bg)
	auto vis = new G4VisAttributes(G4Colour(0.45, 0.45, 0.48, 1.0));
	vis->SetForceSolid(true);
	vis->SetForceAuxEdgeVisible(true);
	frameLV->SetVisAttributes(vis);

// place (start with no rot/offset because STL is centered)
	auto rotFrame = new G4RotationMatrix();
	rotFrame->rotateZ(180*deg);
	new G4PVPlacement(rotFrame,
										G4ThreeVector(0,0,0),
										frameLV, "FramePV",
										//mother_phys->GetLogicalVolume(),
                    moduleLV,
										false, 0, true);
#endif

}



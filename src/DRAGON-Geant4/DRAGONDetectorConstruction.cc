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
/// \file DRAGON/src/DRAGONDetectorConstruction.cc
/// \brief Implementation of the DRAGON::DRAGONDetectorConstruction class


#include "G4RunManager.hh"                   //Geannt4
#include "G4PhysicalVolumeStore.hh"
#include "G4AutoDelete.hh"
#include "G4UserLimits.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4Polyhedra.hh"
#include "G4VSolid.hh"
#include "G4SDManager.hh"

#include "DRAGONDetectorConstruction.hh"     //local
#include "DRAGONPhysicsList.hh"
#include "DRAGONSensitiveDetector.hh"
#include "Materials.hh"


namespace DRAGON
{

DRAGONDetectorConstruction::DRAGONDetectorConstruction(DRAGONPhysicsList* phys)
:DRAGON_phys(phys), checkOverlaps(true), MASK(" "), box_width(5.08 * cm)
{
 uvinit();

 fDetMessenger = new DRAGONDetectorMessenger(this);
 }

DRAGONDetectorConstruction::~DRAGONDetectorConstruction()
{
 delete userLimits;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DRAGONDetectorConstruction::Construct()
{
std::cout << "-GEOMETRIA-" << std::endl;

Materials* materials = Materials::Instance();


materials->atarg = 4;//DRAGON_phys->Getatarg();
materials->bulk_absorption = this->bulk_absorption;
materials->paint_absorption = this->paint_absorption;
materials->ugmate();
materials->ugstmed();  
// DRAGON_phys->Setmtarg(materials->mtarg);
// DRAGON_phys->SetMtarg(materials->Target);

std::cout << "Material: " << materials->Target->GetName() << std::endl;

// DRAGON_phys->Setentdens(materials->entdens);
// DRAGON_phys->Setexitdens(materials->exitdens);

G4cout << *(G4Material::GetMaterialTable()); 

ugeom();
//C     .
//C     .-->   Geometry description
//C     .
//C     .-->   Initialize MITRAY B-field routines by reading RAYTRACE file

//mitray_setup();

G4PhysicalVolumeStore* volumeStore = G4PhysicalVolumeStore::GetInstance();
G4VPhysicalVolume* WRLD_phys = volumeStore->GetVolume("WRLD");
userLimits = new G4UserLimits();
userLimits->SetMaxAllowedStep(max_step);
userLimits->SetUserMaxTrackLength(len_max*cm);
WRLD_phys->GetLogicalVolume()->SetUserLimits(userLimits);

//G4int a = Getntot("HSNG");

return WRLD_phys;

}

void DRAGONDetectorConstruction::uvinit()
{
 //From uvinit.f
 //C.
 //C.    *************** USER DEFAULTS ***************
 //C.
 max_step = 100000;   //! Default maximum number of steps
 len_max = 10000.;    //! Default is 100m track length
 //C.
 n_detmate = 14;
 //C.
 s_finger = 5.08;
 z_finger = 7.62;
 //C.
 d_air[0] = 0.01270;
 d_air[1] = 0.01270;
 d_mtl = 0.00254;
 //C.
 wall[0] = 0.5;
 wall[1] = 0.1;
 wall[2] = 0.1;
 //C.
 //C. gap  = 4.0;
 //C. pinhe= 0.6;
 //C. pinhx = 0.8;
 //C.
 //C. shield_end[0] =  0.0;
 //C. shield_end[1] = 20.0;
 //C.
 mtype_pmt  =  1;
 pmt_size   =  2.5;
 pmt_length = 10.0;
 //C.
 bulk_absorption  = 10.0;
 paint_absorption = 0.20;
 //C.
 xtnd_block  = 1.0;
 
//DRAGON default geometry configuration
d_mtl = 0.063500002;      
d_air[0] =0.035500012; 
d_air[1] =0.317500025; 
z_finger =7.61999998;
pmt_length =2.5;
hexagon_large_width =6.82774401;
hexagon_small_width =5.91299963; 
air_gap =0.127000004;
aprt =0.449600041;
wall[0] =0.1;
wall[1] =0.317500025;
wall[2] =0.497800022;
box_width = 5.07999992;
depth = 8.00099945;
s_finger = 5.588;	
targtype = 0;
}


void DRAGONDetectorConstruction::SetMaxStepsTrackLength(const G4String& input) {
    	std::istringstream iss(input);
	iss >> max_step >> len_max; 
}

void DRAGONDetectorConstruction::SetMASK(const G4String& input)
{
 std::istringstream iss(input);
    	iss >> MASK;
}

void DRAGONDetectorConstruction::SetSHLD(const G4String& input) {
    	std::istringstream iss(input);
    	iss >> shield_end[0] >> shield_end[1]; 

    G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void DRAGONDetectorConstruction::SetPMTR(const G4String& input) 
{
    std::istringstream iss(input);
    iss >> pmt_size >> pmt_length;

    G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void DRAGONDetectorConstruction::SetTUBE(G4int newVal) 
	 {
		 G4cout << "SetTUBE: " << newVal << G4endl;
	  tubetype = newVal;  
	  
	  G4RunManager::GetRunManager()->GeometryHasBeenModified();
	  }

void DRAGONDetectorConstruction::SetWALL(const G4String& input) 
{
    std::istringstream iss(input);
    iss >> wall[0] >> wall[1] >> wall[2];

    G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void DRAGONDetectorConstruction::SetFSID(const G4String& input) 
{
    std::istringstream iss(input);
    iss >> s_finger >> z_finger >> air_gap >> d_air[0] >> d_air[1] >> d_mtl;
	
	 std::cout << "FROM SetFSID" << std::endl;
	 std::cout << "s_finger " << s_finger << std::endl;
	 std::cout << "z_finger " << z_finger << std::endl;
	 std::cout << "air_gap " << air_gap << std::endl;
	 std::cout << "d_air[0] " << d_air[0] << std::endl;
	 std::cout << "d_air[1] " << d_air[1] << std::endl;
	 std::cout << "d_mtl " << d_mtl << std::endl;

    G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DRAGONDetectorConstruction::ConstructSDandField()
     { 
      udet();       //From uginit.f
      udetmitray(); 

      PrintVolumesInfo();
    
      DRAGONEMField* field = new DRAGONEMField();
       
      if(field)
        {
         G4AutoDelete::Register(field);
         G4EqMagElectricField* fEquation = new G4EqMagElectricField(field);
         //G4EqMagElectricField* fEquation = field->GetEquation();
         G4TransportationManager* transportManager = G4TransportationManager::GetTransportationManager();
         G4FieldManager* fFieldMgr = new G4FieldManager();
         G4PropagatorInField* fFieldPropagator = transportManager->GetPropagatorInField();
         fFieldPropagator->SetMaxLoopCount(10000);
         fFieldMgr->SetDetectorField(field);
         G4MagIntegratorStepper* fStepper = new G4ClassicalRK4(fEquation,8);
         //G4MagIntegratorStepper* fStepper = field->GetStepper(); 
         G4ChordFinder* fChordFinder = new G4ChordFinder((G4MagneticField*)this,field->GetfMinStep(),fStepper);
         fChordFinder->SetDeltaChord(field->GetfDeltaChord());
         fFieldMgr->SetAccuraciesWithDeltaOneStep(field->GetfDeltaOneStep());
         fFieldMgr->SetDeltaIntersection(field->GetfDeltaIntersection());
         fFieldPropagator->SetMinimumEpsilonStep(field->GetfEpsMin());
         fFieldPropagator->SetMaximumEpsilonStep(field->GetfEpsMax());
         fFieldMgr->SetChordFinder(fChordFinder);
 
         //G4int numElements = sizeof(field->DRAGONEMelements)/sizeof(field->DRAGONEMelements[0]);
         G4int numElements = (field->DRAGONEMelements).size();
         
         for(G4int i=0; i<numElements;i++)
            {
             G4String name = field->GetAllDRAGONEMFilters()[i].name;
             //G4LogicalVolume* DRAGONEMfilter = G4LogicalVolumeStore::GetInstance()->GetVolume(name);
             G4LogicalVolume* DRAGONEMfilter = G4LogicalVolumeStore::GetInstance()->GetVolume("D1");
             if (DRAGONEMfilter)
                {DRAGONEMfilter->SetFieldManager(fFieldMgr, true);}
             else
                {
                 std::cout << "UNABLE TO READ VOLUMES FOR ELECTRIC OR MAGNETIC FIELD ASSOCIATION" << std::endl;
                 std::cout << "!!! Abort geometry implementation !!!" << std::endl;
                 G4RunManager::GetRunManager()->AbortEvent();
                 }
             } 
         }
      else
         {
          std::cout << "UNABLE TO IMPLEMENT ELECTRIC AND MAGNETIC FIELDS" << std::endl;
          std::cout << "!!! Abort geometry implementation !!!" << std::endl;
          G4RunManager::GetRunManager()->AbortEvent();
          }
     } 


void DRAGONDetectorConstruction::PrintVolumesInfo() {

G4cout << "******************************************************************" << G4endl;
G4cout << "     Summary of Geometry Elements Parameters and Positions        " << G4endl;
G4cout << "******************************************************************" << G4endl;

for (auto lv : *G4LogicalVolumeStore::GetInstance()) {
    G4cout << "\nLogical Volume: " << lv->GetName() << G4endl;
    G4cout << "-------------------------" << G4endl;

    if (lv->GetSensitiveDetector()) 
	   {
        auto sd = static_cast<DRAGONSensitiveDetector*>(lv->GetSensitiveDetector());
        G4cout << "Volumen lógico: " << lv->GetName() << " es sensible." << G4endl;
        G4cout << "Nombre del SD: " << sd->GetName() << G4endl;
        G4cout << "Nombre de la colección de hits: " << sd->GetHitsCollectionName() << G4endl;
        }

    auto material = lv->GetMaterial();
    if (material) {
      G4cout << "  Material: " << material->GetName() << G4endl;
    } else {
      G4cout << "  [No material assigned]" << G4endl;
    }

    G4VSolid* solid = lv->GetSolid();

    if (!solid) {
        G4cout << "  [No solid attached]\n";
        return;
    }

    G4cout << "  Solid type: " << solid->GetEntityType() << G4endl;

    if (auto box = dynamic_cast<G4Box*>(solid)) {
        G4cout << "    X half-length: " << box->GetXHalfLength() / cm << " cm\n"
               << "    Y half-length: " << box->GetYHalfLength() / cm << " cm\n"
               << "    Z half-length: " << box->GetZHalfLength() / cm << " cm\n";
    }
    else if (auto trd = dynamic_cast<G4Trd*>(solid)) {
        G4cout << "    X half-length (−z): " << trd->GetXHalfLength1() / cm << " cm\n"
               << "    X half-length (+z): " << trd->GetXHalfLength2() / cm << " cm\n"
               << "    Y half-length (−z): " << trd->GetYHalfLength1() / cm << " cm\n"
               << "    Y half-length (+z): " << trd->GetYHalfLength2() / cm << " cm\n"
               << "    Z half-length:      " << trd->GetZHalfLength() / cm << " cm\n";
    }
    else if (auto trap = dynamic_cast<G4Trap*>(solid)) {
        G4cout << "    Z half-length: " << trap->GetZHalfLength() / cm << " cm\n"
               << "    Theta: " << trap->GetTheta() / deg << " deg\n"
               << "    Phi:   " << trap->GetPhi() / deg << " deg\n"
               << "    Y1 half: " << trap->GetYHalfLength1() / cm << " cm\n"
               << "    X1 half (−y): " << trap->GetXHalfLength1() / cm << " cm\n"
               << "    X2 half (+y): " << trap->GetXHalfLength2() / cm << " cm\n"
               << "    Alpha1: " << trap->GetTanAlpha1() << " (tan)\n"
               << "    Y2 half: " << trap->GetYHalfLength2() / cm << " cm\n"
               << "    X3 half (−y): " << trap->GetXHalfLength3() / cm << " cm\n"
               << "    X4 half (+y): " << trap->GetXHalfLength4() / cm << " cm\n"
               << "    Alpha2: " << trap->GetTanAlpha2() << " (tan)\n";
    }
    else if (auto tubs = dynamic_cast<G4Tubs*>(solid)) {
        G4cout << "    Inner radius: " << tubs->GetInnerRadius() / cm << " cm\n"
               << "    Outer radius: " << tubs->GetOuterRadius() / cm << " cm\n"
               << "    Z half-length: " << tubs->GetZHalfLength() / cm << " cm\n"
               << "    Start phi: " << tubs->GetStartPhiAngle() / deg << " deg\n"
               << "    Delta phi: " << tubs->GetDeltaPhiAngle() / deg << " deg\n";
    }
    else if (auto cons = dynamic_cast<G4Cons*>(solid)) {
        G4cout << "    Inner radius −Z: " << cons->GetInnerRadiusMinusZ() / cm << " cm\n"
               << "    Outer radius −Z: " << cons->GetOuterRadiusMinusZ() / cm << " cm\n"
               << "    Inner radius +Z: " << cons->GetInnerRadiusPlusZ() / cm << " cm\n"
               << "    Outer radius +Z: " << cons->GetOuterRadiusPlusZ() / cm << " cm\n"
               << "    Z half-length:    " << cons->GetZHalfLength() / cm << " cm\n"
               << "    Start phi: " << cons->GetStartPhiAngle() / deg << " deg\n"
               << "    Delta phi: " << cons->GetDeltaPhiAngle() / deg << " deg\n";
    }

    auto lvStore = G4LogicalVolumeStore::GetInstance();
    for (auto lv : *lvStore) {
        if (lv->GetSolid() == solid) {
            G4cout << "  Attached to logical volume: " << lv->GetName() << G4endl;
            if (lv->GetNoDaughters() != 0)G4cout << "  Daughter volumes: " << G4endl;
            for (int i = 0; i < lv->GetNoDaughters(); ++i) {
                auto phys = lv->GetDaughter(i);
                auto pos = phys->GetTranslation();
                G4cout << "    " << phys->GetName()
                       << " at position: ("
                       << pos.x() / cm << ", "
                       << pos.y() / cm << ", "
                       << pos.z() / cm << ") cm\n";
            }
        }
    }
  }
 G4cout << "******************************************************************" << G4endl << G4endl;
G4cout << "       Summary of Active Sensitive Detectors (SDs)      " << G4endl << G4endl; 
auto sdManager = G4SDManager::GetSDMpointer();
sdManager->ListTree(); 
G4cout << "******************************************************************" << G4endl << G4endl;
}

G4int DRAGONDetectorConstruction::Getntot(const G4String& logicalName)
{
    auto* lvStore = G4LogicalVolumeStore::GetInstance();
    auto* pvStore = G4PhysicalVolumeStore::GetInstance();

    G4LogicalVolume* PMT_log = nullptr;
    for (auto* lv : *lvStore) {
        if (lv->GetName() == logicalName) {
            PMT_log = lv;
            break;
        }
    }

    if (!PMT_log) {
        G4cerr << "No logical volume found with name: " << logicalName << G4endl;
        return 0;
    }

    G4int count = 0;
    for (auto* pv : *pvStore) {
        if (pv->GetLogicalVolume() == PMT_log) {
            count++;
        }
    }

    G4cout << "Logical volume '" << logicalName
           << "' has " << count << " physical copies." << G4endl;

    return count;
}


}

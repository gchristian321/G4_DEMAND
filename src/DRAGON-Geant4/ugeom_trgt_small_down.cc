#include "G4Tubs.hh"                         //Geant4
#include "G4Box.hh" 
#include "G4Cons.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4UserLimits.hh"

#include "DRAGONDetectorConstruction.hh"     //local
#include "Materials.hh"                  

namespace DRAGON
{

void DRAGONDetectorConstruction::ugeo_trgt_small_down()
     { 
      //C     Small Pumping tubes
            
      G4double z, par[10], zcol, shape[3];
      
      Materials* materials = Materials::Instance();
      G4Material* material;
      
      G4LogicalVolumeStore* volumeStore = G4LogicalVolumeStore::GetInstance();
      G4LogicalVolume* DETE_log = volumeStore->GetVolume("DETE");

      //CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      //C     Al entrance tubes
      
      par[0] = len2;
      par[1] = rilen2;
      par[2] = Rrms;
      par[3] = riren2;
      par[4] = Rrms;
      
      G4VSolid* EN2C_solid = new G4Cons("EN2C",par[1]*cm,par[2]*cm,par[3]*cm,par[4]*cm,par[0]*cm,0.0*deg,360.0*deg);
      //TMED->6
      material = materials->Aluminium;
      G4LogicalVolume* EN2C_log = new G4LogicalVolume(EN2C_solid,material,"EN2C");

      par[0] = len3;
      par[1] = rilen3;
      par[2] = Rrms;
      par[3] = riren3;
      par[4] = Rrms;

      G4VSolid* EN3C_solid = new G4Cons("EN3C",par[1]*cm,par[2]*cm,par[3]*cm,par[4]*cm,par[0]*cm,0.0*deg,360.0*deg);
      //TMED->6
      G4LogicalVolume* EN3C_log = new G4LogicalVolume(EN3C_solid,material,"EN3C");

      //CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      //C     Exit apertures (taken from Dave O's Feb 06 drawings)

      par[0] = 3.52/2.;
      par[1] = 1.81/2.;
      par[2] = Rrms;
      par[3] = 1.81/2.;
      par[4] = Rrms;
      G4VSolid* EX2C_solid = new G4Cons("EX2C",par[1]*cm,par[2]*cm,par[3]*cm,par[4]*cm,par[0]*cm,0.0*deg,360.0*deg);
      //TMED->6
      G4LogicalVolume* EX2C_log = new G4LogicalVolume(EX2C_solid,material,"EX2C");

      par[0] = 5.5/2.;
      par[1] = 1.8/2.;
      par[2] = Rrms;
      par[3] = 2.01/2.;
      par[4] = Rrms;
      G4VSolid* EX3C_solid = new G4Cons("EX3C",par[1]*cm,par[2]*cm,par[3]*cm,par[4]*cm,par[0]*cm,0.0*deg,360.0*deg);
      //TMED->6
      G4LogicalVolume* EX3C_log = new G4LogicalVolume(EX3C_solid,material,"EX3C");

      par[0] = 4.16/2.;
      par[1] = 2.67/2.;
      par[2] = Rrms;
      par[3] = 2.67/2.;
      par[4] = Rrms;      
      G4VSolid* EX4C_solid = new G4Cons("EX4C",par[1]*cm,par[2]*cm,par[3]*cm,par[4]*cm,par[0]*cm,0.0*deg,360.0*deg);
      //TMED->6
      G4LogicalVolume* EX4C_log = new G4LogicalVolume(EX4C_solid,material,"EX4C");

      par[0] = 24.9/2.;
      par[1] = 2.71/2.;
      par[2] = Rrms;
      par[3] = 3.67/2.;
      par[4] = Rrms;
      G4VSolid* EX5C_solid = new G4Cons("EX5C",par[1]*cm,par[2]*cm,par[3]*cm,par[4]*cm,par[0]*cm,0.0*deg,360.0*deg);
      //TMED->6
      G4LogicalVolume* EX5C_log = new G4LogicalVolume(EX5C_solid,material,"EX5C");

      par[0] = 22.9/2.;
      par[1] = 5.45/2.;
      par[2] = Rrms;
      par[3] = 6.51/2.;
      par[4] = Rrms;
      G4VSolid* EX6C_solid = new G4Cons("EX6C",par[1]*cm,par[2]*cm,par[3]*cm,par[4]*cm,par[0]*cm,0.0*deg,360.0*deg);
      //TMED->6
      G4LogicalVolume* EX6C_log = new G4LogicalVolume(EX6C_solid,material,"EX6C");

      //CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      //C     Define entrance gas volumes and place the collimators inside
      par[0] = 0.;
      par[1] = Rrms;
      par[2] = len1;
      //C     Position new gas volumes to match new collimator 
      par[2] = 0.5*0.68;

      G4VSolid* EN1G_solid = new G4Tubs("EN1G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
    
      //TMED->28
      material = materials->CentralVacuum;
      G4LogicalVolume* EN1G_log = new G4LogicalVolume(EN1G_solid,material,"EN1G");
      G4UserLimits* EN1GLimits = new G4UserLimits(10.0*cm);  //From ugstmed_trgt.f
      EN1G_log->SetUserLimits(EN1GLimits);
      
      zcol = -(0.5*box_length -wall[1] - col_collar_length + col_length + par[2]);
    
      G4ThreeVector EN1G_pos(0.0*cm,0.0*cm,zcol*cm);
      G4VPhysicalVolume* EN1G_phys = new G4PVPlacement(0,EN1G_pos,EN1G_log,"EN1G",DETE_log,false,0,checkOverlaps);

      z = zent[0];
      
      par[2] = 0.5 * (zent[0] - zent[2] - len1-len2);
      
      G4VSolid* EN2G_solid = new G4Tubs("EN2G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28
      G4LogicalVolume* EN2G_log = new G4LogicalVolume(EN2G_solid,material,"EN2G");
	  	  
      z = .5 *(zent[0] + zent[2]-len1 + len2);
      
      G4ThreeVector EN2G_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* EN2G_phys = new G4PVPlacement(0,EN2G_pos,EN2G_log,"EN2G",DETE_log,false,0,checkOverlaps);  

      //C     Define entrance 3 + 4 gas volumes

      par[0] = 0.;
      par[1] = Rrms;
      par[2] = len2;

      G4VSolid* EN3G_solid = new G4Tubs("EN3G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28
      G4LogicalVolume* EN3G_log = new G4LogicalVolume(EN3G_solid,material,"EN3G");

      G4VPhysicalVolume* EN2C_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0.0*cm),EN2C_log,"EN2C",EN3G_log,false,0,checkOverlaps);  
      z = zent[2];
      G4ThreeVector EN3G_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* EN3G_phys = new G4PVPlacement(0,EN3G_pos,EN3G_log,"EN3G",DETE_log,false,0,checkOverlaps);  
	
      par[2] = 0.5 * (zent[2] - zent[4] - len3 -len2);

      G4VSolid* EN4G_solid = new G4Tubs("EN4G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28
      G4LogicalVolume* EN4G_log = new G4LogicalVolume(EN4G_solid,material,"EN4G");

      z = .5 * (zent[2] + zent[4] +len3 -len2);

      G4ThreeVector EN4G_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* EN4G_phys = new G4PVPlacement(0,EN4G_pos,EN4G_log,"EN4G",DETE_log,false,0,checkOverlaps); 
      
      //C     Define entrance 5 + 6 gas volumes

      par[0] = 0.;
      par[1] = Rrms;
      par[2] = len3;

      G4VSolid* EN5G_solid = new G4Tubs("EN5G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28
      G4LogicalVolume* EN5G_log = new G4LogicalVolume(EN5G_solid,material,"EN5G");

      G4VPhysicalVolume* EN3C_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0.0*cm),EN3C_log,"EN3C",EN5G_log,false,0,checkOverlaps); 

      z = zent[4];

      G4ThreeVector EN5G_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* EN5G_phys = new G4PVPlacement(0,EN5G_pos,EN5G_log,"EN5G",DETE_log,false,0,checkOverlaps);

      par[2] = 0.5 * (TLrms + zent[4] - len3);

      G4VSolid* EN6G_solid = new G4Tubs("EN6G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28
      G4LogicalVolume* EN6G_log = new G4LogicalVolume(EN6G_solid,material,"EN6G");

      z = -TLrms + par[2];

      G4ThreeVector EN6G_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* EN6G_phys = new G4PVPlacement(0,EN6G_pos,EN6G_log,"EN6G",DETE_log,false,0,checkOverlaps);

      //CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      //C     Create exit gas volumes, place collimators within them
      par[0] = 0.;
      par[1] = Rrms;
      par[2] = 6.7/2.;

      G4VSolid* EX2G_solid = new G4Tubs("EX2G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28
      G4LogicalVolume* EX2G_log = new G4LogicalVolume(EX2G_solid,material,"EX2G");
      z = 26.25;
      G4ThreeVector EX2G_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* EX2G_phys = new G4PVPlacement(0,EX2G_pos,EX2G_log,"EX2G",DETE_log,false,0,checkOverlaps);

      par[0] = 0.;
      par[1] = Rrms;
      par[2] = 3.52/2.;

      G4VSolid* EX3G_solid = new G4Tubs("EX3G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28
      G4LogicalVolume* EX3G_log = new G4LogicalVolume(EX3G_solid,material,"EX3G");
      G4VPhysicalVolume* EX2C_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0.0*cm),EX2C_log,"EX2C",EX3G_log,false,0,checkOverlaps);
      z = 31.44;
      G4ThreeVector EX3G_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* EX3G_phys = new G4PVPlacement(0,EX3G_pos,EX3G_log,"EX3G",DETE_log,false,0,checkOverlaps);
      
      par[0] = 0.;
      par[1] = Rrms;
      par[2] = 4.2/2.;

      G4VSolid* EX4G_solid = new G4Tubs("EX4G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28
      G4LogicalVolume* EX4G_log = new G4LogicalVolume(EX4G_solid,material,"EX4G");
      z = 35.3;
      G4ThreeVector EX4G_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* EX4G_phys = new G4PVPlacement(0,EX4G_pos,EX4G_log,"EX4G",DETE_log,false,0,checkOverlaps);

      par[0] = 0.;
      par[1] = Rrms;
      par[2] = 5.5/2.;

      G4VSolid* EX5G_solid = new G4Tubs("EX5G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28
      G4LogicalVolume* EX5G_log = new G4LogicalVolume(EX5G_solid,material,"EX5G");
      G4ThreeVector EX3C_pos(0.0*cm,0.0*cm,0.0*cm);
      G4VPhysicalVolume* EX3C_phys = new G4PVPlacement(0,EX3C_pos,EX3C_log,"EX3C",EX5G_log,false,0,checkOverlaps);
      z = 40.15;
      G4ThreeVector EX5G_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* EX5G_phys = new G4PVPlacement(0,EX5G_pos,EX5G_log,"EX5G",DETE_log,false,0,checkOverlaps);
      
      par[0] = 0.;
      par[1] = Rrms;
      par[2] = 5.24/2.;

      G4VSolid* EX6G_solid = new G4Tubs("EX6G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28
      G4LogicalVolume* EX6G_log = new G4LogicalVolume(EX6G_solid,material,"EX6G");
      z = 45.48;
      G4ThreeVector EX6G_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* EX6G_phys = new G4PVPlacement(0,EX6G_pos,EX6G_log,"EX6G",DETE_log,false,0,checkOverlaps);

      par[0] = 0.;
      par[1] = Rrms;
      par[2] = 4.16/2.;

      G4VSolid* EX7G_solid = new G4Tubs("EX7G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28
      G4LogicalVolume* EX7G_log = new G4LogicalVolume(EX7G_solid,material,"EX7G");
      G4VPhysicalVolume* EX4C_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),EX4C_log,"EX4C",EX7G_log,false,0,checkOverlaps);
      z = 50.22;
      G4ThreeVector EX7G_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* EX7G_phys = new G4PVPlacement(0,EX7G_pos,EX7G_log,"EX7G",DETE_log,false,0,checkOverlaps);
      
      par[0] = 0.;
      par[1] = Rrms;
      par[2] = 7.85/2.;

      G4VSolid* EX8G_solid = new G4Tubs("EX8G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28
      G4LogicalVolume* EX8G_log = new G4LogicalVolume(EX8G_solid,material,"EX8G");
      z = 56.175;
      G4ThreeVector EX8G_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* EX8G_phys = new G4PVPlacement(0,EX8G_pos,EX8G_log,"EX8G",DETE_log,false,0,checkOverlaps);

      par[0] = 0.;
      par[1] = Rrms;
      par[2] = 24.9/2.;

      G4VSolid* EX9G_solid = new G4Tubs("EX9G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28
      G4LogicalVolume* EX9G_log = new G4LogicalVolume(EX9G_solid,material,"EX9G");
      G4VPhysicalVolume* EX5C_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0.0*cm),EX5C_log,"EX5C",EX9G_log,false,0,checkOverlaps);
      z = 72.55;

      //CCC Collimator for alpha test
  	//C 	down


  	//C. **** Do +y border
 
  	shape[0] = 2.5;
  	shape[1] = 2.5;
  	shape[2] = 0.025;
  	G4VSolid* ATC1_solid = new G4Box("ATC1",shape[0]*cm,shape[1]*cm,shape[2]*cm);
  	//TMED->5
  	material = materials->Copper;
  	G4LogicalVolume* ATC1_log = new G4LogicalVolume(ATC1_solid, material,"ATC1");
  	G4VPhysicalVolume* ATC1_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,(-1.2885+0.7+2.5)*cm,8.93048077*cm),ATC1_log,"ATC1",EX9G_log,false,0,checkOverlaps);  
 
  	//C. **** Do -y border

  	shape[0] = 2.5;
  	shape[1] = 2.5;
  	shape[2] = 0.025;
  	G4VSolid* ATC2_solid = new G4Box("ATC2",shape[0]*cm,shape[1]*cm,shape[2]*cm);
  	//TMED->5
  	G4LogicalVolume* ATC2_log = new G4LogicalVolume(ATC2_solid, material,"ATC2");
  	G4VPhysicalVolume* ATC2_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,(-1.2885-0.7-2.5)*cm,8.93048077*cm),ATC2_log,"ATC2",EX9G_log,false,0,checkOverlaps);  
 
  	//C. **** Do +x border
	 
  	shape[0] = 2.5;
  	shape[1] = 2.5;  
  	shape[2] = 0.025;
  	G4VSolid* ATC3_solid = new G4Box("ATC3",shape[0]*cm,shape[1]*cm,shape[2]*cm);
  	//TMED->5
  	G4LogicalVolume* ATC3_log = new G4LogicalVolume(ATC3_solid, material,"ATC3");
  	G4VPhysicalVolume* ATC3_phys = new G4PVPlacement(0,G4ThreeVector((.465+2.5)*cm,0.0*cm,8.93048077*cm),ATC3_log,"ATC3",EX9G_log,false,0,checkOverlaps);  
 

  	//C. **** Do -x border

  	shape[0] = 2.5;
  	shape[1] = 2.5;  
  	shape[2] = 0.025;
  	G4VSolid* ATC4_solid = new G4Box("ATC4",shape[0]*cm,shape[1]*cm,shape[2]*cm);
  	//TMED->5
  	G4LogicalVolume* ATC4_log = new G4LogicalVolume(ATC4_solid, material,"ATC4");
  	G4VPhysicalVolume* ATC4_phys = new G4PVPlacement(0,G4ThreeVector((-0.465-2.5)*cm,0.0*cm,8.93048077*cm),ATC4_log,"ATC4",EX9G_log,false,0,checkOverlaps);  
 
  	//CCCC End collimator

      G4ThreeVector EX9G_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* EX9G_phys = new G4PVPlacement(0,EX9G_pos,EX9G_log,"EX9G",DETE_log,false,0,checkOverlaps);  
     
      par[0] = 0.;
      par[1] = Rrms;
      par[2] = 3.0/2.;
      G4VSolid* EX10_solid = new G4Tubs("EX10",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28
      material = materials->CentralVacuum;
      G4LogicalVolume* EX10_log = new G4LogicalVolume(EX10_solid, material,"EX10");
      EX10_log->SetUserLimits(EN1GLimits);  //From ugstmed_trgt.f
   
      z = 86.5;
      G4ThreeVector EX10_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* EX10_phys = new G4PVPlacement(0,EX10_pos,EX10_log,"EX10",DETE_log,false,0,checkOverlaps);  
 
      //C     Since the following collimator extends into Q1, it is placed
      //C     in the WRLD coordinates with no gas 
      z = 106.05;
      G4ThreeVector EX6C_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* EX6C_phys = new G4PVPlacement(0,EX6C_pos,EX6C_log,"EX6C",WRLD_log,false,0,checkOverlaps);  
 
      //C     Need to define a TEND to change particle charge state
      //c$$$      par(3)=5.
      //c$$$      z = TLrms-par(3) 
      //c$$$      G4VSolid* TEND_solid = new G4Tubs("TEND",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //c$$$      G4LogicalVolume* TEND_log = new G4LogicalVolume(TEND_solid, material,"TEND");
      //c$$$      G4ThreeVector TEND_pos(0.0*cm,0.0*cm,z*cm);
      //c$$$      G4VPhysicalVolume* TEND_phys = new G4PVPlacement(0,TEND_pos,TEND_log,"TEND",DETE_log,false,0,checkOverlaps); 
      
      }
}

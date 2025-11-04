#include "G4PVPlacement.hh"                    //Geant4
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4VSolid.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4UserLimits.hh"

#include "DRAGONDetectorConstruction.hh"       //local
#include "Materials.hh"


namespace DRAGON
{
//C     Large pumping tubes
void DRAGONDetectorConstruction::ugeo_trgt_large() 
     {
     
      G4double z, par[10], zcol;	
      
      Materials* materials = Materials::Instance();
      G4Material* material;
      
      G4LogicalVolumeStore* volumeStore = G4LogicalVolumeStore::GetInstance();
      G4LogicalVolume* DETE_log = volumeStore->GetVolume("DETE"); 
      
      //CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      //C     Al entrance tubes
      //--------------------EN2C----------------------// 
      Rrms = 6.0;
      par[0] = len2;
      par[1] = rilen2;
      par[2] =  Rrms;
      par[3] = riren2;
      par[4] = Rrms;
      
      G4VSolid* EN2C_solid = new G4Cons("EN2C",par[1]*cm,par[2]*cm,par[3]*cm,par[4]*cm,par[0]*cm,0.0*deg,360.0*deg);
      //TMED->6  
      material = materials->Aluminium;
      G4LogicalVolume* EN2C_log = new G4LogicalVolume(EN2C_solid,material,"EN2C");
      //--------------------EN3C----------------------// 
      par[0] = len3;
      par[1] = rilen3;
      par[2] = Rrms;
      par[3] = riren3;
      par[4] = Rrms;
        
      G4VSolid* EN3C_solid = new G4Cons("EN3C",par[1]*cm,par[2]*cm,par[3]*cm,par[4]*cm,par[0]*cm,0.0*deg,360.0*deg);
      //TMED->6  
      G4LogicalVolume* EN3C_log = new G4LogicalVolume(EN3C_solid,material,"EN3C");
   
      //CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      //C     Exit apertures (As Dave O's Feb 06 drawings)
  
      //--------------------EX2C----------------------// 
      par[0] = 12.9/2.;
      par[1] = 2.34/2.;
      par[2] = Rrms;
      par[3] = 2.92/2.;
      par[4] = Rrms;
       
      G4VSolid* EX2C_solid = new G4Cons("EX2C",par[1]*cm,par[2]*cm,par[3]*cm,par[4]*cm,par[0]*cm,0.0*deg,360.0*deg);
      //TMED->6  
      G4LogicalVolume* EX2C_log = new G4LogicalVolume(EX2C_solid,material,"EX2C");
      //--------------------EX3C----------------------// 
      par[0] = 6.03/2.;
      par[1] = 3.78/2.;
      par[2] = Rrms;
      par[3] = 3.78/2.;
      par[4] = Rrms;
      
      G4VSolid* EX3C_solid = new G4Cons("EX3C",par[1]*cm,par[2]*cm,par[3]*cm,par[4]*cm,par[0]*cm,0.0*deg,360.0*deg);
      //TMED->6  
      G4LogicalVolume* EX3C_log = new G4LogicalVolume(EX3C_solid,material,"EX3C");
      //--------------------EX4C----------------------// 
      par[0] = 13.8/2.;
      par[1] = 3.96/2.;
      par[2] = Rrms;
      par[3] = 4.725/2.;
      par[4] = Rrms;
      
      G4VSolid* EX4C_solid = new G4Cons("EX4C",par[1]*cm,par[2]*cm,par[3]*cm,par[4]*cm,par[0]*cm,0.0*deg,360.0*deg);
      //TMED->6  
      G4LogicalVolume* EX4C_log = new G4LogicalVolume(EX4C_solid,material,"EX4C");
      //--------------------EX6C----------------------// 
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
      par[2] =  0.5*0.68;
      //--------------------EN1G----------------------//  
      
      G4VSolid* EN1G_solid = new G4Tubs("EN1G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28
      material = materials->CentralVacuum;
      G4LogicalVolume* EN1G_log = new G4LogicalVolume(EN1G_solid, material,"EN1G");
      G4UserLimits* EN1GLimits = new G4UserLimits(10.0*cm);  //From ugstmed_trgt.f
      EN1G_log->SetUserLimits(EN1GLimits);      
      zcol = -(0.5*box_length - wall[1] - col_collar_length + col_length + par[2]);
      G4ThreeVector EN1G_pos(0.0*cm,0.0*cm,zcol*cm);
      G4PVPlacement* EN1G_phys = new G4PVPlacement(0,EN1G_pos,EN1G_log,"EN1G",DETE_log,false,0,checkOverlaps);
      
      //--------------------EN2G----------------------// 
      z = zent[0];
      
	  par[2] = 0.5 * (zent[0] - zent[2] - len1-len2);
      //TMED->28
      G4VSolid* EN2G_solid = new G4Tubs("EN2G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      G4LogicalVolume* EN2G_log = new G4LogicalVolume(EN2G_solid, material,"EN2G"); 
      
      z = .5 * (zent[0] + zent[2]-len1 + len2);
      
      G4ThreeVector EN2G_pos(0.0*cm,0.0*cm,z*cm);
      G4PVPlacement* EN2G_phys = new G4PVPlacement(0,EN2G_pos,EN2G_log,"EN2G",DETE_log,false,0,checkOverlaps);
  
      //C     Define entrance 3 + 4 gas volumes 
      //--------------------EN3G----------------------// 
      par[0] = 0.;
      par[1] = Rrms;
      par[2] = len2;
        
      G4VSolid* EN3G_solid = new G4Tubs("EN3G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28
      G4LogicalVolume* EN3G_log = new G4LogicalVolume(EN3G_solid, material,"EN3G"); 
      
      G4PVPlacement* EN2C_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),EN2C_log,"EN2C",EN3G_log,false,0,checkOverlaps);
      z = zent[2];
      G4ThreeVector EN3G_pos(0.0*cm,0.0*cm,z*cm);
      G4PVPlacement* EN3G_phys = new G4PVPlacement(0,EN3G_pos,EN3G_log,"EN3G",DETE_log,false,0,checkOverlaps);
      //--------------------EN4G----------------------// 
      par[2] = 0.5 * (zent[2] - zent[4] - len3 -len2);
  
      G4VSolid* EN4G_solid = new G4Tubs("EN4G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28
      G4LogicalVolume* EN4G_log = new G4LogicalVolume(EN4G_solid, material,"EN4G"); 
      
      z = .5 * (zent[2] + zent[4] +len3 -len2);
      G4ThreeVector EN4G_pos(0.0*cm,0.0*cm,z*cm);
      G4PVPlacement* EN4G_phys = new G4PVPlacement(0,EN4G_pos,EN4G_log,"EN4G",DETE_log,false,0,checkOverlaps);
  
      //C     Define entrance 5 + 6 gas volumes 
      //--------------------EN5G----------------------// 
      par[0] = 0.;
      par[1] = Rrms;
      par[2] = len3;
       
      G4VSolid* EN5G_solid = new G4Tubs("EN5G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28
      G4LogicalVolume* EN5G_log = new G4LogicalVolume(EN5G_solid, material,"EN5G"); 
      G4PVPlacement* EN3C_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),EN3C_log,"EN3C",EN5G_log,false,0,checkOverlaps);
      
      z = zent[4];
      
      G4ThreeVector EN5G_pos(0.0*cm,0.0*cm,z*cm);
      G4PVPlacement* EN5G_phys = new G4PVPlacement(0,EN5G_pos,EN5G_log,"EN5G",DETE_log,false,0,checkOverlaps);
  
      //--------------------EN6G----------------------// 
      par[2] = 0.5 * (TLrms + zent[4] - len3);
       
      G4VSolid* EN6G_solid = new G4Tubs("EN6G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28
      G4LogicalVolume* EN6G_log = new G4LogicalVolume(EN6G_solid, material,"EN6G"); 
      
      z = -TLrms + par[2];
      
      G4ThreeVector EN6G_pos(0.0*cm,0.0*cm,z*cm);
      G4PVPlacement* EN6G_phys = new G4PVPlacement(0,EN6G_pos,EN6G_log,"EN6G",DETE_log,false,0,checkOverlaps);
   
      //CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      //C     Create exit gas volumes, place collimators within them
      //c$$$      npar = 3;
      //c$$$      zcol = -(0.5*box_length + 0.945 + 1.026 + 2.54 + 0.624 + 
      //c$$$     &        3.54 + 2.88 + 2.02 + 0.585*2.);
      //c$$$      lex11 = (zex[2] - lex2 + zcol)/2.;
      //c$$$      zex2 = -zcol + lex11;
      //c$$$      par[2] = lex11;
      //c$$$   
      //c$$$
      //c$$$C.
      //c$$$      npar = 3;
      //c$$$      par(1) = 0.;
      //c$$$      par(2) = Rrms;

      //--------------------EX2G----------------------//  
      par[0] = 0.;
      par[1] = Rrms;
      par[2] = 6.97/2.;
       
      G4VSolid* EX2G_solid = new G4Tubs("EX2G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28
      G4LogicalVolume* EX2G_log = new G4LogicalVolume(EX2G_solid, material,"EX2G");
      z = 26.465;
    
      G4ThreeVector EX2G_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* EX2G_phys = new G4PVPlacement(0,EX2G_pos,EX2G_log,"EX2G",DETE_log,false,0,checkOverlaps); 
      
      //--------------------EX3G----------------------//
      par[0] = 0.;
      par[1] = Rrms;
      par[2] = 12.9/2.;
      
      G4VSolid* EX3G_solid = new G4Tubs("EX3G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
       
      G4LogicalVolume* EX3G_log = new G4LogicalVolume(EX3G_solid, material,"EX3G");
      //TMED->28
      G4PVPlacement* EX2C_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),EX2C_log,"EX2C",EX3G_log,false,0,checkOverlaps);
      z = 36.4;
      G4ThreeVector EX3G_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* EX3G_phys = new G4PVPlacement(0,EX3G_pos,EX3G_log,"EX3G",DETE_log,false,0,checkOverlaps);
      
      //--------------------EX4G----------------------//
      par[0] = 0.;
      par[1] = Rrms;
      par[2] = 4.34/2.;
   
      G4VSolid* EX4G_solid = new G4Tubs("EX4G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28 
      G4LogicalVolume* EX4G_log = new G4LogicalVolume(EX4G_solid, material,"EX4G");
      z = 45.02;
      G4ThreeVector EX4G_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* EX4G_phys = new G4PVPlacement(0,EX4G_pos,EX4G_log,"EX4G",DETE_log,false,0,checkOverlaps);
      //--------------------EX5G----------------------//
      par[0] = 0.;
      par[1] = Rrms;
      par[2] = 6.03/2.;

      G4VSolid* EX5G_solid = new G4Tubs("EX5G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28 
      G4LogicalVolume* EX5G_log = new G4LogicalVolume(EX5G_solid, material,"EX5G");
      G4PVPlacement* EX3C_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),EX3C_log,"EX3C",EX5G_log,false,0,checkOverlaps);
      z = 50.205;
      G4ThreeVector EX5G_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* EX5G_phys = new G4PVPlacement(0,EX5G_pos,EX5G_log,"EX5G",DETE_log,false,0,checkOverlaps);
      //--------------------EX6G----------------------//
      par[0] = 0.;
      par[1] = Rrms;
      par[2] =  8.35/2.;

      G4VSolid* EX6G_solid = new G4Tubs("EX6G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28 
      G4LogicalVolume* EX6G_log = new G4LogicalVolume(EX6G_solid, material,"EX6G");
      z = 57.395;
      G4ThreeVector EX6G_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* EX6G_phys = new G4PVPlacement(0,EX6G_pos,EX6G_log,"EX6G",DETE_log,false,0,checkOverlaps);
      
      //--------------------EX7G----------------------//
      par[0] = 0.;
      par[1] = Rrms;
      par[2] = 13.8/2.;

      G4VSolid* EX7G_solid = new G4Tubs("EX7",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28 
      G4LogicalVolume* EX7G_log = new G4LogicalVolume(EX7G_solid, material,"EX7G");
      G4PVPlacement* EX4C_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0.0*cm),EX4C_log,"EX4C",EX7G_log,false,0,checkOverlaps);
      z = 68.47;
      G4ThreeVector EX7G_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* EX7G_phys = new G4PVPlacement(0,EX7G_pos,EX7G_log,"EX7G",DETE_log,false,0,checkOverlaps);
      
      //--------------------EX8G----------------------//
      par[0] = 0.;
      par[1] = Rrms;
      par[2] = 12.63/2.;

      G4VSolid* EX8G_solid = new G4Tubs("EX8G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28 
      G4LogicalVolume* EX8G_log = new G4LogicalVolume(EX8G_solid, material,"EX8G");
      z = 81.685;
      G4ThreeVector EX8G_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* EX8G_phys = new G4PVPlacement(0,EX8G_pos,EX8G_log,"EX8G",DETE_log,false,0,checkOverlaps);
     
      //C     Since the following collimator extends into Q1, it is placed
      //C     in the WRLD coordinates with no gas 
      z = 106.05;
      G4ThreeVector EX6C_pos(0.0*cm,0.0*cm,z*cm);
      G4PVPlacement* EX6C_phys = new G4PVPlacement(0,EX6C_pos,EX6C_log,"EX6C",WRLD_log,false,0,checkOverlaps); 
//c$$$      par[2]=5.;
//c$$$      z = TLrms-par[2];
//c$$$      G4VSolid* TEND_solid = new G4Tubs("TEND",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
//c$$$      G4ThreeVector TEND_pos(0.0*cm,0.0*cm,z*cm);
//c$$$      G4PVPlacement* TEND_phys = new G4PVPlacement(0,TEND_pos,TEND_log,"TEND",DETE_log,false,0,checkOverlaps); 
      }
  
}

#include "G4PVPlacement.hh"                 //Geant4
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4UserLimits.hh"

#include "DRAGONDetectorConstruction.hh"    //local
#include "Materials.hh"

namespace DRAGON
{

void DRAGONDetectorConstruction::ugeo2_trgt() 
     {    
      G4double z, par[10], zcol;
 
      Materials* materials = Materials::Instance();
      G4Material* material;
     
      G4LogicalVolumeStore* volumeStore = G4LogicalVolumeStore::GetInstance();
      G4LogicalVolume* DETE_log = volumeStore->GetVolume("DETE");
      //------------------------------------------------------------------

      //C.--> Define Al Entrance tubes

      //C. Aluminum
      par[0] = len1;
      par[1] = rilen1;
      par[2] =  Rrms;
      par[3] = riren1;
      par[4] = Rrms;

      //      G4VSolid* EN1C_solid = new G4Cons("EN1C",par[1]*cm,par[2]*cm,par[3]*cm,par[4]*cm,par[0]*cm,0.0*deg,360.0*deg);
      //      G4LogicalVolume* EN1C_log = new G4LogicalVolume(EN1C_solid,water_mat,"EN1C");
      
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

      //C.--> Define exit apertures

      par[0] = lex1;
      par[1] = rilex1;
      par[2] = Rrms;
      par[3] = rirex1;
      par[4] = Rrms;

      //      G4VSolid* EX1C_solid = new G4Cons("EX1C",par[1]*cm,par[2]*cm,par[3]*cm,par[4]*cm,par[0]*cm,0.0*deg,360.0*deg);
      //      TMED->6
      //      G4LogicalVolume* EX1C_log = new G4LogicalVolume(EX1C_solid,water_mat,"EX1C");

      par[0] = lex2;
      par[1] = rilex2;
      par[2] = Rrms;
      par[3] = rirex2;
      par[4] = Rrms;

      G4VSolid* EX2C_solid = new G4Cons("EX2C",par[1]*cm,par[2]*cm,par[3]*cm,par[4]*cm,par[0]*cm,0.0*deg,360.0*deg);
      //TMED->6
      G4LogicalVolume* EX2C_log = new G4LogicalVolume(EX2C_solid,material,"EX2C");    

      par[0] = lex3;
      par[1] = rilex3;
      par[2] = Rrms;
      par[3] = rirex3;
      par[4] = Rrms;

      G4VSolid* EX3C_solid = new G4Cons("EX3C",par[1]*cm,par[2]*cm,par[3]*cm,par[4]*cm,par[0]*cm,0.0*deg,360.0*deg);
      //TMED->6
      G4LogicalVolume* EX3C_log = new G4LogicalVolume(EX3C_solid,material,"EX3C"); 

      par[0] = lex4;
      par[1] = rilex4;
      par[2] = Rrms;
      par[3] = rirex4;
      par[4] = Rrms;

      G4VSolid* EX4C_solid = new G4Cons("EX4C",par[1]*cm,par[2]*cm,par[3]*cm,par[4]*cm,par[0]*cm,0.0*deg,360.0*deg);
      //TMED->6
      G4LogicalVolume* EX4C_log = new G4LogicalVolume(EX4C_solid,material,"EX4C"); 

      //C.--> Define entrance 1 + 2 gas volumes

      //C.    nmed = ment[0]; //! Need vacuum for solid carbon target

      par[0] = 0.;
      par[1] = Rrms;
      par[2] = len1;

      //C.--> Position new gas volumes to match new collimator 

      par[2] = 0.5*0.68;
    
      G4VSolid* EN1G_solid = new G4Tubs("EN1G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28
      material = materials->CentralVacuum;      
      G4LogicalVolume* EN1G_log = new G4LogicalVolume(EN1G_solid, material,"EN1G");
      G4UserLimits* EN1GLimits = new G4UserLimits(10.0*cm);  //From ugstmed_trgt.f
      EN1G_log->SetUserLimits(EN1GLimits);
      
      //CORRECTION IMPLEMENTED TO FIX OVERLAP FROM GEANT3 GEOMETRY WITH VOLUME PUJ1
      //MODIFICATION MADE DEC-3-2024 BEN MARLOW
      G4double correction1 = 0.2445;
      
      zcol = -(0.5*box_length -wall[1] - col_collar_length + col_length + par[2]) + correction1;
      G4ThreeVector EN1G_pos(0.0*cm,0.0*cm,zcol*cm);
      G4VPhysicalVolume* EN1G_phys = new G4PVPlacement(0,EN1G_pos,EN1G_log,"EN1G",DETE_log,false,0,checkOverlaps);  

      z = zent[0];

      //C.    nmed = ment[1]; //! Need vacuum for solid carbon target
  
      //CORRECTION IMPLEMENTED TO FIX OVERLAP FROM GEANT3 GEOMETRY WITH VOLUME EN1G
      //MODIFICATION MADE DEC-3-2024 BEN MARLOW
      G4double correction2 = 0.34625;
      
      par[2] = 0.5 * (zent[0] - zent[2] - len1-len2) - correction2;

      G4VSolid* EN2G_solid = new G4Tubs("EN2G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28
      G4LogicalVolume* EN2G_log = new G4LogicalVolume(EN2G_solid, material,"EN2G");
      
      z = .5 * (zent[0] + zent[2]-len1 + len2) - correction2;

      G4ThreeVector EN2G_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* EN2G_phys = new G4PVPlacement(0,EN2G_pos,EN2G_log,"EN2G",DETE_log,false,0,checkOverlaps);  


      //C.--> Define entrance 3 + 4 gas volumes

      //C.    nmed = ment[2]; //! Need vacuum for solid carbon target
 
      par[0] = 0.;
      par[1] = Rrms;
      par[2] = len2;

      G4VSolid* EN3G_solid = new G4Tubs("EN3G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28
      G4LogicalVolume* EN3G_log = new G4LogicalVolume(EN3G_solid, material,"EN3G");
      
      G4VPhysicalVolume* EN2C_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0.0*cm),EN2C_log,"EN2C",EN3G_log,false,0,checkOverlaps);  

      z = zent[2];

      G4ThreeVector EN3G_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* EN3G_phys = new G4PVPlacement(0,EN3G_pos,EN3G_log,"EN3G",DETE_log,false,0,checkOverlaps);  

      //C.    nmed = ment[3]; //! Need vacuum for solid carbon target
   
      par[2] = 0.5 * (zent[2] - zent[4] - len3 -len2);

      G4VSolid* EN4G_solid = new G4Tubs("EN4G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28
      G4LogicalVolume* EN4G_log = new G4LogicalVolume(EN4G_solid, material,"EN4G");
      
      z = .5 * (zent[2] + zent[4] +len3 -len2);

      G4ThreeVector EN4G_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* EN4G_phys = new G4PVPlacement(0,EN4G_pos,EN4G_log,"EN4G",DETE_log,false,0,checkOverlaps);  

      //C.--> Define entrance 5 + 6 gas volumes

      //C.    nmed = ment[4] //! Need vacuum for solid carbon target
     
      par[0] = 0.;
      par[1] = Rrms;
      par[2] = len3;

      G4VSolid* EN5G_solid = new G4Tubs("EN5G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28
      G4LogicalVolume* EN5G_log = new G4LogicalVolume(EN5G_solid, material,"EN5G");
      
      G4VPhysicalVolume* EN3C_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0.0*cm),EN3C_log,"EN5G",EN5G_log,false,0,checkOverlaps);  

      z = zent[4];

      G4ThreeVector EN5G_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* EN5G_phys = new G4PVPlacement(0,EN5G_pos,EN5G_log,"EN5G",DETE_log,false,0,checkOverlaps);  

      //C.    nmed = ment[5]; //! Need vacuum for solid carbon target
     
      par[2] = 0.5 * (TLrms + zent[4] - len3);

      G4VSolid* EN6G_solid = new G4Tubs("EN6G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28
      G4LogicalVolume* EN6G_log = new G4LogicalVolume(EN6G_solid, material,"EN6G");

      z = -TLrms + par[2];
     
      G4ThreeVector EN6G_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* EN6G_phys = new G4PVPlacement(0,EN6G_pos,EN6G_log,"EN6G",DETE_log,false,0,checkOverlaps);  

      //C.--> Define exit 1 + 2  gas volumes

      //C.    nmed = mex[0]; //! Need vacuum for solid carbon target

      par[0] = 0.;
      par[1] = Rrms;
      par[2] = lex1;
      //C     par[2] = 0.5*0.68;
     
      G4VSolid* EX1G_solid = new G4Tubs("EX1G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28
      G4LogicalVolume* EX1G_log = new G4LogicalVolume(EX1G_solid, material,"EX1G");
      z = -zcol;
      G4ThreeVector EX1G_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* EX1G_phys = new G4PVPlacement(0,EX1G_pos,EX1G_log,"EX1G",DETE_log,false,0,checkOverlaps);  

      //C      z = zex[0];

      //C.    nmed = mex[1] //! Need vacuum for solid carbon target

      par[0] = 0.;
      par[1] = Rrms;
      //C     JS Change to par(3) done for new pumping tube system
      //C      par(3) = (zex(3) - zex(1) - lex1 - lex2)/2.
      par[2] = (zex[2] + zcol - lex1 - lex2)/2.;

      G4VSolid* EX2G_solid = new G4Tubs("EX2G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28
      G4LogicalVolume* EX2G_log = new G4LogicalVolume(EX2G_solid, material,"EX2G");
      
      //C     JS Change to z done for new pumping tube system
      //C      z = (zex(1) + zex(3) + lex1 - lex2)/2.
      z = (-zcol + zex[2] + lex1 - lex2)/2.;      

      G4ThreeVector EX2G_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* EX2G_phys = new G4PVPlacement(0,EX2G_pos,EX2G_log,"EX2G",DETE_log,false,0,checkOverlaps);  

      //C.--> Define exit 3 + 4  gas volumes

      //C.    nmed = mex[2] //! Need vacuum for solid carbon target

      par[0] = 0.;
      par[1] = Rrms;
      par[2] = lex2;

      G4VSolid* EX3G_solid = new G4Tubs("EX3G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28
      G4LogicalVolume* EX3G_log = new G4LogicalVolume(EX3G_solid, material,"EX3G");

      G4VPhysicalVolume* EX2C_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0.0*cm),EX2C_log,"EX2C",EX3G_log,false,0,checkOverlaps);  

      z = zex[2];

      G4ThreeVector EX3G_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* EX3G_phys = new G4PVPlacement(0,EX3G_pos,EX3G_log,"EX3G",DETE_log,false,0,checkOverlaps);  

      //C.    nmed = mex[3] //! Need vacuum for solid carbon target

      par[0] = 0.;
      par[1] = Rrms;
      par[2] = (zex[4] - zex[2]- lex2 - lex3)/2.;

      G4VSolid* EX4G_solid = new G4Tubs("EX4G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28
      G4LogicalVolume* EX4G_log = new G4LogicalVolume(EX4G_solid, material,"EX4G");

      z = (zex[2] + zex[4] + lex2 - lex3)/2.;
   
      G4ThreeVector EX4G_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* EX4G_phys = new G4PVPlacement(0,EX4G_pos,EX4G_log,"EX4G",DETE_log,false,0,checkOverlaps);  

      //C.--> Define exit 5 + 6  gas volumes

      //C.    nmed = mex[4] //! Need vacuum for solid carbon target

      par[0] = 0.;
      par[1] = Rrms;
      par[2] = lex3;

      G4VSolid* EX5G_solid = new G4Tubs("EX5G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28
      G4LogicalVolume* EX5G_log = new G4LogicalVolume(EX5G_solid, material,"EX5G");
  
      G4VPhysicalVolume* EX3C_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0.0*cm),EX3C_log,"EX3C",EX5G_log,false,0,checkOverlaps);  

      z = zex[4];

      G4ThreeVector EX5G_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* EX5G_phys = new G4PVPlacement(0,EX5G_pos,EX5G_log,"EX5G",DETE_log,false,0,checkOverlaps);  

      //C.    nmed = mex[5] //! Need vacuum for solid carbon target

      par[0] = 0.;
      par[1] = Rrms;
      par[2] = (zex[6] - zex[4] - lex3- lex4)/2.;

      G4VSolid* EX6G_solid = new G4Tubs("EX6G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //TMED->28
      G4LogicalVolume* EX6G_log = new G4LogicalVolume(EX6G_solid, material,"EX6G");
      
      z = (zex[6] + zex[4] + lex3 - lex4)/2.;

      G4ThreeVector EX6G_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* EX6G_phys = new G4PVPlacement(0,EX6G_pos,EX6G_log,"EX6G",DETE_log,false,0,checkOverlaps);  

      par[0] = 0.;
      par[1] = Rrms;
      par[2] = lex4;
      //C. JS: In new pumping tube system EX7G is inside Quad 1.  This causes serious
      //C. problems with GEANTS interpretation of the magnetic fields.  DH has
      //C. suggested removing EX7G from the simulation
      //C. G4VSolid* EX7G_solid = new G4Tubs("EX7G",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //C. G4LogicalVolume* EX7G_log = new G4LogicalVolume(EX7G_solid, material,"EX7G");
      //C. 
      //C. G4VPhysicalVolume* EX4C_phys = new G4PVPlacement(0,EX4C_pos,EX4C_log,"EX4C",EX7G_log,false,0,checkOverlaps); 

      z = zex[6];
      //C. JS - EX7G removed in new pumping tube system.  See above comment.
      //C     
      //C  G4ThreeVector EX7G_pos(0.0*cm,0.0*cm,z*cm);
      //C  G4VPhysicalVolume* EX7G_phys = new G4PVPlacement(0,EX7G_pos,EX7G_log,"EX7G",DETE_log,false,0,checkOverlaps); 

      //C JS changed par(3) for new pumping tube system
      par[2] = (TLrms+30.0 -zex[6]+lex4)/2.;
      //C      par[2] =  (TLrms -zex[6]+lex4)/2.;
      //C JS changed z for new pumping tube system
      z = TLrms+30.0-par[2];
      //C      z = TLrms-par[2];

      //C      JS has just outright removed TEND from the simulation
      //C      It doesn't work right now that it's inside Q1
      //C      G4VSolid* TEND_solid = new G4Tubs("TEND",par[0]*cm,par[1]*cm,par[2]*cm,0.0*deg,360.0*deg);
      //C      G4LogicalVolume* TEND_log = new G4LogicalVolume(TEND_solid, material,"TEND");
      
      //C      G4ThreeVector TEND_pos(0.0*cm,0.0*cm,z*cm);
      //C      G4VPhysicalVolume* TEND_phys = new G4PVPlacement(0,TEND_pos,TEND_log,"TEND",DETE_log,false,0,checkOverlaps); 
      
      }
}

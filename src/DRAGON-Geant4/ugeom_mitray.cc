#include "G4VSolid.hh"                     //Geant4
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Trap.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4SDManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "DRAGONSensitiveDetector.hh"
#include "G4UserLimits.hh"

#include "DRAGONDetectorConstruction.hh"    //local
#include "Materials.hh"              
#include "geom_dipole.hh"      
#include "geom_edipol.hh"  
#include "geom_mpole.hh"   
#include "geom_sole.hh"           


namespace DRAGON
{
	
void DRAGONDetectorConstruction::ugeo_dipole(G4int k) 
     {
	  //************************************************************************
      //*                                                                      *
      //*                     Define a RAYTRACE dipole magnet                  *
      //*                                                                      *
      //************************************************************************
		 
      G4double r, dr, z11, z22, alpha, beta, phi, theta, rho;
      
      G4double x1, x2, x11, x12, x21, x22;
      G4double dx, dy, dz;
      G4double dxp, dyp, dzp;
      
      G4double shape[11], pos[3];
      
      std::ostringstream kname, lname;
      G4String vvname;
      
      G4VisAttributes* visAttributes = new G4VisAttributes(G4Colour(0.4, 0.4, 1.0, 0.5)); 
      visAttributes->SetVisibility(true); 
      visAttributes->SetForceSolid(true);
      
      G4double degrad = M_PI/180.0;
      G4double raddeg = 180.0/M_PI;
      
      Materials* materials = Materials::Instance();
      G4Material* material;
      
      kname << "D" << k;
      lname << "BN" << k;
  
      std::copy(&pos_dipole[k-1][0], &pos_dipole[k-1][0] + 3, pos);
      gap   =  gap_dipole[k];        
      r     =  std::abs(r_dipole[k]);
      dr    =  dr_dipole[k];
      phi   =  degrad * phi_dipole[k];
      alpha =  degrad * alpha_dipole[k];
      beta  =  degrad * beta_dipole[k];
      z11   =  z11_dipole[k];
      z22   =  z22_dipole[k];

      x21 = r * std::sin(phi/2.0);
      x21 = x21 + r * (1.0 - std::cos(phi/2.0))*std::tan(phi/2.0 - alpha);
      x21 = x21 + z11/std::cos(phi/2.0 - alpha);
      x21 = x21 + dr * std::tan(phi/2.0 - alpha);

      x11 = x21 - (r + dr) * std::tan(phi/2.0 - alpha);

      x22 = r * std::sin(phi/2.0);
      x22 = x22 + r * (1.0 - std::cos(phi/2.0)) * std::tan(phi/2.0 - beta);
      x22 = x22 + z22/std::cos(phi/2.0 - beta);
      x22 = x22 + dr * std::tan(phi/2.0 - beta);

      x12 = x22 - (r + dr) * std::tan(phi/2.0 - beta);

      x1 = (x11 + x12)/2.0;
      x2 = (x21 + x22)/2.0;

      theta = std::atan((x2 - x1 - (r + dr) * std::tan(phi/2.0 - alpha))/(r + dr));

      shape[0] = gap/2.0;
      shape[1] = 0.0;
      shape[2] = 0.0;
      shape[3] = (r + dr)/2.0;
      shape[4] = x1;
      shape[5] = x2;
      shape[6] = raddeg * theta;
      shape[7] = shape[3]; 
      shape[8] = shape[4];
      shape[9] = shape[5];
      shape[10] = shape[6];
    
      //TMED->1
      //material = materials->Vacuum;
      //G4VSolid* kname_solid = new G4Trap(kname.str(), shape[0]*cm, shape[1]*deg, shape[2]*deg, shape[3]*cm, shape[7]*cm, shape[4]*cm, shape[5]*cm, shape[8]*cm, shape[9]*cm, shape[6]*deg, shape[10]*deg);
      G4VSolid* kname_solid = new G4Trap(kname.str(), shape[0]*cm, shape[1]*deg, shape[2]*deg, shape[3]*cm, shape[4]*cm, shape[5]*cm, shape[6]*deg, shape[7]*cm, shape[8]*cm, shape[9]*cm, shape[10]*deg);
      //TMED->2
      material = materials->Vacuum;
      G4LogicalVolume* kname_log = new G4LogicalVolume(kname_solid, material, kname.str());
      kname_log->SetVisAttributes(visAttributes);
     
      dx = -(r + dr)/2.0 * std::tan(theta);
      dx = dx - (x1 - x11);
      dy = -(r + dr)/2.0;
      dz = 0.0;
      
      //CC Vacuum Vessel Componants
      
      //c     Rotation Matrices 
      
      G4RotationMatrix* rotMatrix200 = new G4RotationMatrix();
      rotMatrix200->rotateZ(25*deg);
      G4RotationMatrix* rotMatrix201 = new G4RotationMatrix();
      rotMatrix201->rotateZ(-6*deg);
      G4RotationMatrix* rotMatrix202 = new G4RotationMatrix();
      rotMatrix202->rotateZ(-25*deg);
      G4RotationMatrix* rotMatrix203 = new G4RotationMatrix();
      rotMatrix203->rotateZ(12*deg);
      G4RotationMatrix* rotMatrix204 = new G4RotationMatrix();
      rotMatrix204->rotateZ(-33*deg);
      
      G4RotationMatrix* irot = new G4RotationMatrix();
      irot->rotateX(-90.0*deg);
               
      if (k == 2) 
         {
          // MD2
          shape[0] = r-9.0;
          shape[1] = r-6.0;
          shape[2] = 5.0;
          shape[3] = 90. - raddeg * phi/2.;
          shape[4] = raddeg * phi;   
                    
          vvname = "VV1";
          G4VSolid* VV1_solid = new G4Tubs(vvname, shape[0]*cm, shape[1]*cm, shape[2]*cm, shape[3]*deg, shape[4]*deg);
          //TMED->20
          material = materials->StainlessSteel;          
          G4LogicalVolume* VV1_log = new G4LogicalVolume(VV1_solid, material, vvname);
          new G4PVPlacement(0, G4ThreeVector(dx*cm, dy*cm, dz*cm), VV1_log, vvname, kname_log, false, k, checkOverlaps);
          G4UserLimits* CMBRLimits = new G4UserLimits(10.0*cm);  //From ugstmed_trgt.f
          VV1_log->SetUserLimits(CMBRLimits);
      
          shape[0] = r+6.0;
          shape[1] = r+9.0;
          shape[2] = 5.0;
          shape[3] = 90. - raddeg * phi/2.;
          shape[4] = raddeg * phi;
         
          vvname = "VV2";
          G4VSolid* VV2_solid = new G4Tubs(vvname, shape[0]*cm, shape[1]*cm, shape[2]*cm, shape[3]*deg, shape[4]*deg);
          //TMED->20
          G4LogicalVolume* VV2_log = new G4LogicalVolume(VV2_solid, material, vvname);
          new G4PVPlacement(0, G4ThreeVector(dx*cm, dy*cm, dz*cm), VV2_log, vvname, kname_log, false, k, checkOverlaps);
          VV2_log->SetUserLimits(CMBRLimits);
      
          shape[0] = r-9.0;
          shape[1] = r+9.0;
          shape[2] = 0.375/2.0;
          shape[3] = 90. - raddeg * phi/2.;
          shape[4] = raddeg * phi;

          vvname = "VV3";
          G4VSolid* VV3_solid = new G4Tubs(vvname, shape[0]*cm, shape[1]*cm, shape[2]*cm, shape[3]*deg, shape[4]*deg);
          //TMED->20
          G4LogicalVolume* VV3_log = new G4LogicalVolume(VV3_solid, material, vvname);
          new G4PVPlacement(0, G4ThreeVector(dx*cm, dy*cm, dz*cm + 5.0*cm + (0.375/2.0)*cm), VV3_log, vvname, kname_log, false, k, checkOverlaps);
          VV3_log->SetUserLimits(CMBRLimits);

          shape[0] = r-9.0;
          shape[1] = r+9.0;
          shape[2] = 0.375/2;
          shape[3] = 90. - raddeg * phi/2.;
          shape[4] = raddeg * phi;

          vvname = "VV4";
          G4VSolid* VV4_solid = new G4Tubs(vvname, shape[0]*cm, shape[1]*cm, shape[2]*cm, shape[3]*deg, shape[4]*deg);
          //TMED->20
          G4LogicalVolume* VV4_log = new G4LogicalVolume(VV4_solid, material, vvname);
          new G4PVPlacement(0, G4ThreeVector(dx*cm, dy*cm, dz*cm - 5.0*cm - (0.375/2.0)*cm), VV4_log, vvname, kname_log, false, k, checkOverlaps);  
          VV4_log->SetUserLimits(CMBRLimits);
     
          irot->rotateZ(-(90.0-raddeg*phi)*deg);
          }
      else 
         {
          if (k == 1) 
             {
              // MD1
              
              //Part 1
              dxp = dx - (100.0-6.3/2.0*2.54-0.375/2.0*2.54)*cos(65.*degrad) - 1.75*2.54/2.0*cos(25.0*degrad);
              dyp = dy + (100.0-6.3/2.0*2.54-0.375/2.0*2.54)*sin(65.*degrad) - 1.75*2.54/2.0*sin(25.0*degrad);
              dzp = dz;
              
              shape[0] = 1.75/2.0*2.54;
              shape[1] = 0.375/2.0*2.54;
              shape[2] = 4.0;

              vvname = "VV5";
              G4VSolid* VV5_solid = new G4Box("VV5", shape[0]*cm, shape[1]*cm, shape[2]*cm);
              //TMED->20
              material = materials->StainlessSteel; 
              G4LogicalVolume* VV5_log = new G4LogicalVolume(VV5_solid, material, vvname);
              new G4PVPlacement(rotMatrix200, G4ThreeVector(-dxp*cm, dyp*cm, dzp*cm), VV5_log, vvname, kname_log, false, k, checkOverlaps);
              G4UserLimits* CMBRLimits = new G4UserLimits(10.0*cm);  //From ugstmed_trgt.f
              VV5_log->SetUserLimits(CMBRLimits);

              //Part 2
              dxp = dx - (100.0+6.3/2.0*2.54+0.375/2.0*2.54)*cos(65.*degrad) + 22.8460*cos(25.0*degrad);
              dyp = dy + (100.0+6.3/2.0*2.54+0.375/2.0*2.54)*sin(65.*degrad) + 22.8460*sin(25.0*degrad);
              dzp = dz;
               
              shape[0] = 27.291;
              shape[1] = 0.375/2.0*2.54;
              shape[2] = 4.0;

              vvname = "VV6";
              G4VSolid* VV6_solid = new G4Box("VV6", shape[0]*cm, shape[1]*cm, shape[2]*cm);
              //TMED->20
              G4LogicalVolume* VV6_log = new G4LogicalVolume(VV6_solid, material, vvname);
              new G4PVPlacement(rotMatrix200, G4ThreeVector(-dxp*cm, dyp*cm, dzp*cm), VV6_log, vvname, kname_log, false, k, checkOverlaps);

              //Part 3
              dxp = dx - (100.0+6.3/2.0*2.54+0.375/2.0*2.54)*cos(65.*degrad) + (22.8460+27.291)*cos(25.0*degrad) + 26.515*cos(6.0*degrad);
              dyp = dy + (100.0+6.3/2.0*2.54+0.375/2.0*2.54)*sin(65.*degrad) + (22.8460+27.291)*sin(25.0*degrad) - 26.515*sin(6.0*degrad);
              dzp = dz;
               
              shape[0] = 26.515;
              shape[1] = 0.375/2.0*2.54;
              shape[2] = 4.0;

              vvname = "VV7";
              G4VSolid* VV7_solid = new G4Box("VV7", shape[0]*cm, shape[1]*cm, shape[2]*cm);
              //TMED->20
              G4LogicalVolume* VV7_log = new G4LogicalVolume(VV7_solid, material, vvname);
              new G4PVPlacement(rotMatrix201, G4ThreeVector(-dxp*cm, dyp*cm, dzp*cm), VV7_log, vvname, kname_log, false, k, checkOverlaps);

              //Part 4
              dxp = dx - (100.0+6.3/2.0*2.54+0.375/2.0*2.54)*cos(65.*degrad) + (22.8460+27.291)*cos(25.0*degrad) + 53.03*cos(6.0*degrad) + 6.836935*cos(25.0*degrad);
              dyp = dy + (100.0+6.3/2.0*2.54+0.375/2.0*2.54)*sin(65.*degrad) + (22.8460+27.291)*sin(25.0*degrad) - 53.03*sin(6.0*degrad) - 6.836935*sin(25.0*degrad);
              dzp = dz;
               
              shape[0] = 6.836935;
              shape[1] = 0.375/2.0*2.54;
              shape[2] = 4.0;

              vvname = "VV8";
              G4VSolid* VV8_solid = new G4Box("VV8", shape[0]*cm, shape[1]*cm, shape[2]*cm);
              //TMED->20
              G4LogicalVolume* VV8_log = new G4LogicalVolume(VV8_solid, material, vvname);
              new G4PVPlacement(rotMatrix202, G4ThreeVector(-dxp*cm, dyp*cm, dzp*cm), VV8_log, vvname, kname_log, false, k, checkOverlaps);
               
              //Part 5
              dxp = dx - (100.0-6.3/2.0*2.54-0.375/2.0*2.54)*cos(65.*degrad) + 19.234*cos(12.0*degrad);
              dyp = dy + (100.0-6.3/2.0*2.54-0.375/2.0*2.54)*sin(65.*degrad) + 19.234*sin(12.0*degrad);
              dzp = dz;

              shape[0] = 19.234;
              shape[1] = 0.375/2.0 * 2.54;
              shape[2] = 4.0;

              vvname = "VV9";
              G4VSolid* VV9_solid = new G4Box(vvname, shape[0]*cm, shape[1]*cm, shape[2]*cm);
              //TMED->20
              G4LogicalVolume* VV9_log = new G4LogicalVolume(VV9_solid, material, vvname);
              new G4PVPlacement(rotMatrix203, G4ThreeVector(-dxp*cm, dyp*cm, dzp*cm), VV9_log, vvname, kname_log, false, k, checkOverlaps);
              
              //Part 6
              dxp = dx - (100.0-6.3/2.0*2.54-0.375/2.0*2.54)*cos(65.*degrad) + 2.0*19.234*cos(12.0*degrad) + 20.0495*cos(33.0*degrad);
              dyp = dy + (100.0-6.3/2.0*2.54-0.375/2.0*2.54)*sin(65.*degrad) + 2.0*19.234*sin(12.0*degrad) - 20.0495*sin(33.0*degrad);
              dzp = dz;

              shape[0] = 20.0495;
              shape[1] = 0.375/2.0 * 2.54;
              shape[2] = 4.0;

              vvname = "VV10";
              G4VSolid* VV10_solid = new G4Box(vvname.c_str(), shape[0]*cm, shape[1]*cm, shape[2]*cm);
              //TMED->20
              G4LogicalVolume* VV10_log = new G4LogicalVolume(VV10_solid, material, vvname);
              new G4PVPlacement(rotMatrix204, G4ThreeVector(-dxp*cm, dyp*cm, dzp*cm), VV10_log, vvname, kname_log, false, k, checkOverlaps);

              //Part 7
              dxp = dx - (100.0-6.3/2.0*2.54-0.375/2.0*2.54)*cos(65.*degrad) + 2.0*19.234*cos(12.0*degrad) + 40.099*cos(33.0*degrad) + 6.836935*cos(25.0*degrad);
              dyp = dy + (100.0-6.3/2.0*2.54-0.375/2.0*2.54)*sin(65.*degrad) + 2.0*19.234*sin(12.0*degrad) - 40.099*sin(33.0*degrad) - 6.836935*sin(25.0*degrad);
              dzp = dz;

              shape[0] = 6.836935;
              shape[1] = 0.375/2.0*2.54;
              shape[2] = 4.0;

              vvname = "VV11";
              G4VSolid* VV11_solid = new G4Box(vvname.c_str(), shape[0]*cm, shape[1]*cm, shape[2]*cm);
              //TMED->20
              G4LogicalVolume* VV11_log = new G4LogicalVolume(VV11_solid, material, vvname);
              new G4PVPlacement(rotMatrix202, G4ThreeVector(-dxp*cm, dyp*cm, dzp*cm), VV11_log, vvname, kname_log, false, k, checkOverlaps);

              //Part 8
              dxp = dx - (100.0) * cos(65.0 * degrad) + 20.0 * cos(25.0 * degrad);
              dyp = dy + (100.0) * sin(65.0 * degrad) + 20.0 * sin(25.0 * degrad);
              dzp = dz + 4.0 + 0.375 * 2.54/2.0;

              shape[0] = 20.0 + 1.75 * 2.54;
              shape[1] = (6.3 + 2.0 * 0.375) * 2.54/2.0 + 10.0;
              shape[2] = 0.375/2.0 * 2.54;

              vvname = "VV12";
              G4VSolid* VV12_solid = new G4Box(vvname.c_str(), shape[0]*cm, shape[1]*cm, shape[2]*cm);
              //TMED->20
              G4LogicalVolume* VV12_log = new G4LogicalVolume(VV12_solid, material, vvname);
              new G4PVPlacement(rotMatrix200, G4ThreeVector(-dxp*cm, dyp*cm, dzp*cm), VV12_log, vvname, kname_log, false, k, checkOverlaps);

              //Part 9
              dxp = dx - (100.0) * cos(65.0 * degrad) + 20.0 * cos(25.0 * degrad) + 40.0;
              dyp = dy + (100.0) * sin(65.0 * degrad) + 20.0 * sin(25.0 * degrad) - 10.0;
              dzp = dz + 4.0 + 0.375 * 2.54/2.0;

              shape[0] = 20.0 + 1.75 * 2.54 + 13.0;
              shape[1] = (6.3 + 2.0 * 0.375) * 2.54/2.0 + 10.0 + 32.0;
              shape[2] = 0.375/2.0 * 2.54;

              vvname = "VV13";
              G4VSolid* VV13_solid = new G4Box(vvname.c_str(), shape[0]*cm, shape[1]*cm, shape[2]*cm);
              //TMED->20
              G4LogicalVolume* VV13_log = new G4LogicalVolume(VV13_solid, material, vvname);
              new G4PVPlacement(rotMatrix202, G4ThreeVector(-dxp*cm, dyp*cm, dzp*cm), VV13_log, vvname, kname_log, false, k, checkOverlaps);

              //Part 10
              dxp = dx - (100.0) * cos(65.0 * degrad) + 20.0 * cos(25.0 * degrad);
              dyp = dy + (100.0) * sin(65.0 * degrad) + 20.0 * sin(25.0 * degrad); 
              dzp = dz - 4.0 - 0.375 * 2.54/2.0;

              shape[0] = 20.0+1.75*2.54;
              shape[1] = (6.3+2.0*0.375)*2.54/2.0+10.0;
              shape[2] = 0.375/2.0*2.54;
               
              vvname = "VV14";
              G4VSolid* VV14_solid = new G4Box(vvname.c_str(), shape[0]*cm, shape[1]*cm, shape[2]*cm);
              //TMED->20
              G4LogicalVolume* VV14_log = new G4LogicalVolume(VV14_solid, material, vvname);
              new G4PVPlacement(rotMatrix200, G4ThreeVector(-dxp*cm, dyp*cm, dzp*cm), VV14_log, vvname, kname_log, false, k, checkOverlaps);

              //Part 11
              dxp = dx - (100.0) * cos(65.0 * degrad) + 20.0 * cos(25.0 * degrad) + 40.0;
              dyp = dy + (100.0) * sin(65.0 * degrad) + 20.0 * sin(25.0 * degrad) - 10.0;
              dzp = dz - 4.0 - 0.375 * 2.54/2.0;
               
              shape[0] = 20.0+1.75*2.54+13.0;
              shape[1] = (6.3+2.0*0.375)*2.54/2.0+10.0+32.0;
              shape[2] = 0.375/2.0*2.54;

              vvname = "VV15";
              G4VSolid* VV15_solid = new G4Box(vvname.c_str(), shape[0]*cm, shape[1]*cm, shape[2]*cm);
              //TMED->20
              G4LogicalVolume* VV15_log = new G4LogicalVolume(VV15_solid, material, vvname);
              new G4PVPlacement(rotMatrix202, G4ThreeVector(-dxp*cm, dyp*cm, dzp*cm), VV15_log, vvname, kname_log, false, k, checkOverlaps);
              
              irot->rotateZ((90-raddeg*phi/2.)*deg);  
              }
          else 
              {
               G4cout << "Wow, three MD's, add my vacuum vessel please." << G4endl;
               std::exit(EXIT_FAILURE);
		       }   
          }
 
      //CC End of Vacuum Vessel

      shape[0] = r-dr;
      shape[1] = r+dr;
      shape[2] = gap/2.;
      shape[3] = 90. - raddeg * phi/2.;
      shape[4] = 90. + raddeg * phi/2.;

      //G4VSolid* lname_solid = new G4Tubs(lname.str(), shape[0]*cm, shape[1]*cm, shape[2]*cm, shape[3]*deg, shape[4]*deg);
      //TMED->3
      //material = materials->Vacuum;
      //G4LogicalVolume* lname_log = new G4LogicalVolume(lname_solid, material, lname.str());
      
      //G4VPhysicalVolume* lname_phys = new G4PVPlacement(0,G4ThreeVector(dx*cm, dy*cm, dz*cm),lname_log,lname.str(),kname_log,false,k,checkOverlaps);  

      rho = raddeg * atan((z11/cos(alpha))/r);
      rho = 90.0 - raddeg*phi/2. - rho;
      rho = degrad * rho;

      //C. ***** To position the A-frame at the edge of the entrance fringe field

      //dx = dx - sqrt(r**2+(z11/cos(alpha))**2) * cos(rho);
      //dy = dy + sqrt(r**2+(z11/cos(alpha))**2) * sin(rho);
      //dz = 0.0;

      //C. ***** To position the A-frame at the edge of the entrance EFB

      dx = dx - r * sin(phi/2.);
      dy = dy + r * cos(phi/2.);
      dz = 0.0;

      shape[0] = 2.*(r+dr);
      shape[1] = gap/2.;
      shape[2] = 2.*x2;

      //G4VSolid* mname_solid = new G4Tubs(mname.str(), shape[0]*cm, shape[1]*cm, shape[2]*cm, shape[3]*deg, shape[4]*deg);
      //TMED->3
      //material = materials->Vacuum;
      //G4LogicalVolume* mname_log = new G4LogicalVolume(mname_solid, material, mname.str());
      //G4VPhysicalVolume* mname_phys = new G4PVPlacement(0,G4ThreeVector(dx*cm, dy*cm, dz*cm),mname_log,mname.str(),mname_log,false,k,checkOverlaps);  

      dx_dipole[0][k-1] = dx;
      dx_dipole[1][k-1] = dy;
      dx_dipole[2][k-1] = dz;

      if(k==1)irot_dipole[k-1] = 10;
      if(k==2)irot_dipole[k-1] = 37;
      
//OJO
/*
std::cout << "PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP\n";
std::cout << "D" << k;
std::cout << " dx_dipole " << dx_dipole[0][k] << " " << dx_dipole[1][k] << " " << dx_dipole[2][k] << std::endl;
std::cout << "irot_dipole " << irot_dipole[k] << std::endl;
std::cout << "PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP\n";
*/
      new G4PVPlacement(irot,G4ThreeVector(pos[0]*cm, pos[1]*cm, pos[2]*cm),kname_log,kname.str(),WRLD_log,false,k,checkOverlaps);  
      }


void DRAGONDetectorConstruction::ugeo_edipol(G4int k) 
     {
	  //************************************************************************
      //*                                                                      *
      //*               Define a RAYTRACE electrostatic deflector              *
      //*                                                                      *
      //************************************************************************	 
		
      G4double r, dr, gap, z11, z22, phi, theta, rho;
      
      G4double x1, x2, x11, x12, x21, x22;
      G4double dx, dy, dz;
      
      G4double shape[11], pos[3];
      
      std::ostringstream kname, lname, lnameb;
      
      G4VisAttributes* visAttributes = new G4VisAttributes(G4Colour(1.0, 1., .0, 0.5)); 
      visAttributes->SetVisibility(true); 
      visAttributes->SetForceSolid(true);
      
      G4double degrad = M_PI/180.0;
      G4double raddeg = 180.0/M_PI;
      
      Materials* materials = Materials::Instance();
      G4Material* material;
      
      kname << "E" << k;
      lname << "PN" << k;
      lnameb << "PM" << k;
 
      std::copy(&pos_edipol[k-1][0], &pos_edipol[k-1][0] + 3, pos); 
      gap   =  gap_edipol[k];        
      r     =  std::abs(r_edipol[k]);
      dr    =  dr_edipol[k];
      phi   =  degrad * phi_edipol[k];
      z11   =  z11_edipol[k];
      z22   =  z22_edipol[k];

	  x21 = r*sin(phi/2.);
      x21 = x21 + r*(1.-cos(phi/2.))*tan(phi/2.);
      x21 = x21 + z11/cos(phi/2.);
      x21 = x21 + dr * tan(phi/2.);

      x11 = x21 - (r+dr) * tan(phi/2.);

      x22 = r*sin(phi/2.);
      x22 = x22 + r*(1.-cos(phi/2.))*tan(phi/2.);
      x22 = x22 + z22/cos(phi/2.);
      x22 = x22 + dr * tan(phi/2.);

      x12 = x22 - (r+dr) * tan(phi/2.);

      x1 = (x11 + x12)/2.;
      x2 = (x21 + x22)/2.;

      theta = atan((x2-x1-(r+dr)*tan(phi/2.))/(r+dr));
      
      shape[0] = gap/2.0;
      shape[1] = 0.0;
      shape[2] = 0.0;
      shape[3] = (r + dr)/2.0;
      shape[4] = x1;
      shape[5] = x2;
      shape[6] = raddeg * theta;
      shape[7] = shape[3]; 
      shape[8] = shape[4];
      shape[9] = shape[5];
      shape[10] = shape[6];

      //TMED->1
      //material = materials->Vacuum;
      G4VSolid* kname_solid = new G4Trap(kname.str(), shape[0]*cm, shape[1]*deg, shape[2]*deg, shape[3]*cm, shape[4]*cm, shape[5]*cm, shape[6]*deg, shape[7]*cm, shape[8]*cm, shape[9]*cm, shape[10]*deg);
      //TMED->2
      material = materials->Vacuum;
      G4LogicalVolume* kname_log = new G4LogicalVolume(kname_solid, material, kname.str());
      kname_log->SetVisAttributes(visAttributes);

      dx = -(r + dr)/2.0 * std::tan(theta);
      dx = dx - (x1 - x11);
      dy = -(r + dr)/2.0;
      dz = 0.0;		
   
      shape[0] = r-8.0;
      shape[1] = r-5.0;
      if (k == 1)
         {
          shape[2] = 14.0;
          }
      else
         shape[2] = 15.0;
      
      shape[3] = 90. - raddeg * phi/2.;
      shape[4] = raddeg * phi;	
  
      G4VSolid* lname_solid = new G4Tubs(lname.str(), shape[0]*cm, shape[1]*cm, shape[2]*cm, shape[3]*deg, shape[4]*deg);
      //TMED->20
      material = materials->StainlessSteel;
      G4LogicalVolume* lname_log = new G4LogicalVolume(lname_solid, material, lname.str());
      visAttributes = new G4VisAttributes(G4Colour(0.0, 1., 1.0, 0.5)); 
      lname_log->SetVisAttributes(visAttributes);
      G4VPhysicalVolume* lname_phys = new G4PVPlacement(0,G4ThreeVector(dx*cm, dy*cm, dz*cm),lname_log,lname.str(),kname_log,false,k,checkOverlaps); 
      G4UserLimits* CMBRLimits = new G4UserLimits(10.0*cm);  //From ugstmed_trgt.f
      lname_log->SetUserLimits(CMBRLimits); 
	 
	  shape[0] = r+5.0;
      shape[1] = r+8.0;
      if (k == 1)
         {
          shape[2] = 14.0;
          }
      else
         shape[2] = 15.0;
      
      shape[3] = 90. - raddeg * phi/2.;
      shape[4] = raddeg * phi;
	
	  G4VSolid* lnameb_solid = new G4Tubs(lnameb.str(), shape[0]*cm, shape[1]*cm, shape[2]*cm, shape[3]*deg, shape[4]*deg);
      //TMED->20
      G4LogicalVolume* lnameb_log = new G4LogicalVolume(lnameb_solid, material, lnameb.str());
      lnameb_log->SetVisAttributes(visAttributes);
      
      G4VPhysicalVolume* lnameb_phys = new G4PVPlacement(0,G4ThreeVector(dx*cm, dy*cm, dz*cm),lnameb_log,lnameb.str(),kname_log,false,k,checkOverlaps);  

      rho = raddeg * atan(z11/r);
      rho = 90.0 - raddeg*phi/2. - rho;
      rho = degrad * rho;

      //C. ***** To position the A-frame at the edge of the entrance fringe field

      //dx = dx - sqrt(r**2+z11**2) * cos(rho);
      //dy = dy + sqrt(r**2+z11**2) * sin(rho);
      //dz = 0.0;

      //C. ***** To position the A-frame at the edge of the entrance EFB

      dx = dx - r * sin(phi/2.);
      dy = dy + r * cos(phi/2.);
      dz = 0.0;

      shape[0] = 2.*(r+dr);
      shape[1] = gap/2.;
      shape[2] = 2.*x2;
      
      //G4VSolid* mname_solid = new G4Box(mname.str(),shape[0]*cm,shape[1]*cm,shape[2]*cm);
      //TMED->5
      //material = materials->Copper;
      //G4LogicalVolume* mname_log = new G4LogicalVolume(mname_solid,material,mname.str());
	  
  	  G4RotationMatrix* irot = new G4RotationMatrix();
      irot->rotateX(-90.0*deg);
      
      if (k == 1) irot->rotateZ(30.*deg);     
      else irot->rotateZ(-(90.0-raddeg*phi/2.0)*deg);
  	  //G4VPhysicalVolume * mname_phys = new G4PVPlacement(0,G4ThreeVector(dx*cm, dy*cm, dz*cm),mname_log,lname.str(),kname_log,false,k,checkOverlaps);  
	
	  dx_edipol[0][k-1] = dx;
      dx_edipol[1][k-1] = dy;
      dx_edipol[2][k-1] = dz;
      
      if(k==2)irot_edipol[k-1] = 47;
      if(k==1)irot_edipol[k-1] = 24;

//OJO
/*
std::cout << "PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP\n";
std::cout << "E" << k;
std::cout << " dx_edipol " << dx_edipol[0][k] << " " << dx_edipol[1][k] << " " << dx_edipol[2][k] << std::endl;
std::cout << "irot_edipol " << irot_edipol[k] << std::endl;
 std::cout << "PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP\n";
*/     
      G4VPhysicalVolume* kname_phys = new G4PVPlacement(irot,G4ThreeVector(pos[0]*cm, pos[1]*cm, pos[2]*cm),kname_log,kname.str(),WRLD_log,false,k,checkOverlaps);    
	  }

void DRAGONDetectorConstruction::ugeo_mpole(G4int k,G4double rot_angles[14*3]) 
     {
	  //************************************************************************
      //*                                                                      *
      //*                   Define a RAYTRACE mulipoles magnet                 *
      //*                                                                      *
      //************************************************************************

      static G4int m = 0;
 
      G4double shape[3], pos[3], dx, dy, dz, shape_pip[3];
      
      std::ostringstream kname, pname;
 	  
 	  G4VisAttributes* visAttributes = new G4VisAttributes(G4Colour(0., 1., 1., 0.5)); 
      visAttributes->SetVisibility(true); 
      visAttributes->SetForceSolid(true);
      
	  kname << "Q" << k;
      pname << "BP" << k;
      
      G4double degrad = M_PI/180.0;
      G4double raddeg = 180.0/M_PI;
      
      Materials* materials = Materials::Instance();
      G4Material* material;
	  
	  std::copy(&pos_mpole[k-1][0], &pos_mpole[k-1][0] + 3, pos);
	  
	  shape[0] = 0.0;
      shape[1] = r_mpole[k];
      shape[2] = (efblength_mpole[k]+z11_mpole[k]+z22_mpole[k])/2.;
      if(kname.str() == "Q1")G4cout << "Q1..... " << shape[2] << G4endl;
	
	  //TMED->1
      //material = materials->Vacuum;
      //G4VSolid* kname_solid = new G4Tubs(kname.str(),shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      G4VSolid* kname_solid = new G4Tubs(kname.str(),shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
  
      //TMED->2
      material = materials->Vacuum;
      G4LogicalVolume* kname_log = new G4LogicalVolume(kname_solid,material,kname.str());
   	  kname_log->SetVisAttributes(visAttributes);
      //C. ***** To position the A-frame at the edge of the entrance EFB

      dx = 0.0;
      dy = 0.0;
      dz = (z11_mpole[k]-z22_mpole[k]-efblength_mpole[k])/2.;

      shape[0] = 2.*shape[1];
      shape[1] = 2.*shape[1];
      shape[2] = 2.*shape[2];  
      
      //G4VSolid* mname_solid = new G4Box(mname.str(),shape[0]*cm,shape[1]*cm,shape[2]*cm);
      //material = materials->Vacuum;
      //G4LogicalVolume* mname_log = new G4LogicalVolume(mname_solid,material,mname.str());
      //G4VPhysicalVolume* mname_phys = new G4PVPlacement(0,G4ThreeVector(dx*cm, dy*cm, dz*cm),mname_log,mname.str(),kname_log,false,k,checkOverlaps); 
      
      dx_mpole[0][k-1] = dx;
      dx_mpole[1][k-1] = dy;
      dx_mpole[2][k-1] = dz;
   
//OJO   
/*
std::cout << "PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP\n";
std::cout << "Q" << k;
std::cout << " dx_mpole " << dx_mpole[0][k] << " " << dx_mpole[1][k] << " " << dx_mpole[2][k] << std::endl;
std::cout << "PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP\n";
*/
      G4int a,b,c,d,e,f;
      
      G4RotationMatrix* rotMatrix = new G4RotationMatrix();
      a = m++;
      rotMatrix->rotateZ(rot_angles[a]*deg); 
      b = m++;
      rotMatrix->rotateY(rot_angles[b]*deg); 
      c = m++;
      rotMatrix->rotateX(rot_angles[c]*deg);     
      G4VPhysicalVolume* kname_phys = new G4PVPlacement(rotMatrix,G4ThreeVector(pos[0]*cm, pos[1]*cm, pos[2]*cm),kname_log,kname.str(),WRLD_log,false,k,checkOverlaps);  
       
      //C.--> insert beampipe in same space as multipole
      //C.    (comment this out if not required)
      shape_pip[0] = dxcol_mpole[0][k];                                        //! inner radius
      shape_pip[1] = shape_pip[0]+0.15875;                                     //! plus 1/16th inch
      shape_pip[2] = (efblength_mpole[k]+z11_mpole[k]+z22_mpole[k])/2.;
      
      G4VSolid* pname_solid = new G4Tubs(pname.str(),shape_pip[0]*cm,shape_pip[1]*cm,shape_pip[2]*cm,0.0*deg,360.0*deg);
      //TMED->20
      material = materials->StainlessSteel;
      G4LogicalVolume* pname_log = new G4LogicalVolume(pname_solid,material,pname.str());
      G4VPhysicalVolume* pname_phys = new G4PVPlacement(rotMatrix,G4ThreeVector(pos[0]*cm, pos[1]*cm, pos[2]*cm),pname_log,pname.str(),WRLD_log,false,k,checkOverlaps); 
      G4UserLimits* CMBRLimits = new G4UserLimits(10.0*cm);  //From ugstmed_trgt.f 
      pname_log->SetUserLimits(CMBRLimits); 
      }
      
void DRAGONDetectorConstruction::ugeo_sole(G4int k)
     {
      //************************************************************************
      //*                                                                      *
      //*                    Define a RAYTRACE solenoid magnet                 *
      //*                                                                      *
      //************************************************************************
 
      G4double shape[3], pos[3], dx, dy, dz;

      std::ostringstream kname;
 
	  kname << "S" << k;
      
      G4double degrad = M_PI/180.0;
      G4double raddeg = 180.0/M_PI;
      
      Materials* materials = Materials::Instance();
      G4Material* material;
  
      std::copy(&pos_sole[k-1][0], &pos_sole[k-1][0] + 3, pos);

      shape[0] = 0.0;
      shape[1] = r_sole[k];
      shape[2] = (efblength_sole[k]+z11_sole[k]+z22_sole[k])/2.;

      G4VSolid* kname_solid = new G4Tubs(kname.str(),shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->2
      material = materials->Vacuum;
      G4LogicalVolume* kname_log = new G4LogicalVolume(kname_solid,material,kname.str());

      //C. ***** To position the A-frame at the edge of the entrance EFB

      dx = 0.0;
      dy = 0.0;
      dz = (z11_sole[k]-z22_sole[k]-efblength_sole[k])/2.;

      shape[0] = 2.*shape[1];
      shape[1] = 2.*shape[1];
      shape[2] = 2.*shape[1];
      
      //G4VSolid* mname_solid = new G4Box(mname.str(), shape[0]*cm, shape[1]*cm, shape[2]*cm);
      //G4LogicalVolume* mname_log = new G4LogicalVolume(mname_solid, material, mname.str());
      //TMED->2
      //material = materials->Vacuum;
      //G4LogicalVolume* mname_log = new G4LogicalVolume(mname_solid, material, mname.str());
      //new G4PVPlacement(0, G4ThreeVector(dx*cm,dy*cm,dz*cm), mname_log, mname, kname_log, false, k, checkOverlaps);

      dx_sole[0][k-1] = dx;
      dx_sole[1][k-1] = dy;
      dx_sole[2][k-1] = dz;

//OJO
/*
std::cout << "PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP\n"; 
std::cout << "S" << k;
std::cout << " dx_sole " << dx_sole[0][k] << " " << dx_sole[1][k] << " " << dx_sole[2][k] << std::endl;
std::cout << "PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP\n"; 
*/ 
      new G4PVPlacement(0, G4ThreeVector(dx*cm,dy*cm,dz*cm), kname_log, kname.str(), WRLD_log, false, k, checkOverlaps);
      }


void DRAGONDetectorConstruction::ugeo_col(G4double pos[45*3],G4double rot_angles[45*3],G4double data[6],G4String rname)
     {
	  //************************************************************************
      //*                                                                      *
      //*                        Define a REAL COLLIMATOR                      *
      //*                                                                      *
      //************************************************************************

      static G4int i, m=-3;

      G4double shape[3];
      G4double d[3];
     
      std::ostringstream kname;
      
      G4VisAttributes* visAttributes = new G4VisAttributes(G4Colour(0.5, 0.25, 0., 0.5)); 
      visAttributes->SetVisibility(true); 
      visAttributes->SetForceSolid(true);
      
      Materials* materials = Materials::Instance();
      G4Material* material; 
      
      G4RotationMatrix* rotMatrix = new G4RotationMatrix();
           
      if(true)
        {		 
	     if(data[0] == 0)
	       {
		    i = i + 1;
		    kname << "C" << i;
		    m = m + 3;
        
            rotMatrix->rotateZ(rot_angles[m]*deg);
            rotMatrix->rotateY(rot_angles[m+1]*deg);
            rotMatrix->rotateX(rot_angles[m+2]*deg);
	
	        if(rname == "QSLT") 
	          {
               shape[0] = 20.0*data[3];
               shape[1] = 20.0*data[4];
               } 
            else 
               {
                if(rname == "MSLT") 
                  {
                   shape[0] = 10.0*data[3];
                   shape[1] = 10.0*data[4];
                   } 
                else 
                   {
				    if(rname == "FSLT")	
				   	  {
					   shape[0] = 10.0*data[3];
                       shape[1] = 10.0*data[4];
					   }
				    else
			           {
                        shape[0] = 1.0*data[3];
                        shape[1] = 1.0*data[4];
			            }
                    }
	            }
	        shape[2] = data[5];
	    
            //--------------------C----------------------//
            G4VSolid* kname_solid = new G4Box(kname.str(),shape[0]*cm,shape[1]*cm,shape[2]*cm);
            //TMED->5
            material = materials->Copper;
            G4LogicalVolume* kname_log = new G4LogicalVolume(kname_solid,material,kname.str());
      		kname_log->SetVisAttributes(visAttributes);
	        //C. **** Do +x border
	    
	        d[0] = data[1] + data[3] + shape[0];
            d[1] = data[2];
            d[2] = 0.0;

	        G4VPhysicalVolume* kname_phys = new G4PVPlacement(rotMatrix,G4ThreeVector(pos[m]*cm,pos[m+1]*cm,pos[m+2]*cm),kname_log,kname.str(),WRLD_log,false,0,checkOverlaps);  
 	  
	        //C. **** Do -x border
	    
	        d[0] = data[1] - data[3] - shape[0];

            m = m + 3;

            rotMatrix->rotateZ(rot_angles[m]*deg);
            rotMatrix->rotateY(rot_angles[m+1]*deg);
            rotMatrix->rotateX(rot_angles[m+2]*deg);

 	        kname_phys = new G4PVPlacement(rotMatrix,G4ThreeVector(pos[m]*cm,pos[m+1]*cm,pos[m+2]*cm),kname_log,kname.str(),WRLD_log,false,1,checkOverlaps);  
 	  
 	        i = i + 1;
 	    
 	        kname.str(" ");
	        kname << "C" << i;		 
	     
	        if(rname == "QSLT") 
	          {
               shape[0] = 20.0*data[3];
               shape[1] = 20.0*data[4];
               } 
            else 
               {
                if(rname == "MSLT") 
                  {
                   shape[0] = 10.0*data[3];
                   shape[1] = 10.0*data[4];
                   } 
                else 
                   {
				    if(rname == "FSLT")	
				      { 
					   shape[0] = 10.0*data[3];
                       shape[1] = 10.0*data[4];
					   }
				    else
			           {
                        shape[0] = 1.0*data[3];
                        shape[1] = 1.0*data[4];
			            }
                    }
	            }
	        shape[2] = data[5];

            //--------------------C----------------------//
            kname_solid = new G4Box(kname.str(),shape[0]*cm,shape[1]*cm,shape[2]*cm);
            //TMED->5
            material = materials->Copper;
            kname_log = new G4LogicalVolume(kname_solid,material,kname.str());
			kname_log->SetVisAttributes(visAttributes);
	        //C. **** Do +y border
	    
	        d[0] = data[1];
            d[1] = data[2] + data[4] + shape[1];
            d[2] = 0.0;
	    
            m = m + 3;

            rotMatrix->rotateZ(rot_angles[m]*deg);
            rotMatrix->rotateY(rot_angles[m+1]*deg);
            rotMatrix->rotateX(rot_angles[m+2]*deg);
	   
	        kname_phys = new G4PVPlacement(rotMatrix,G4ThreeVector(pos[m]*cm,pos[m+1]*cm,pos[m+2]*cm),kname_log,kname.str(),WRLD_log,false,0,checkOverlaps);  
 	  
	        //C. **** Do -y border
	    
	        d[1] = data[2] - data[4] - shape[1];

            m = m + 3;

            rotMatrix->rotateZ(rot_angles[m]*deg);
            rotMatrix->rotateY(rot_angles[m+1]*deg);
            rotMatrix->rotateX(rot_angles[m+2]*deg);
        
 	        kname_phys = new G4PVPlacement(rotMatrix,G4ThreeVector(pos[m]*cm,pos[m+1]*cm,pos[m+2]*cm),kname_log,kname.str(),WRLD_log,false,1,checkOverlaps);  
 	        }	 
	     else	 
	        {
		     i = i + 1;
		     kname << "C" << i;

		     shape[0] = data[3];
             shape[1] = data[4];
             shape[2] = data[5];

             //--------------------C----------------------//
		     G4VSolid* kname_solid = new G4Tubs(kname.str(),shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
             //TMED->5
             material = materials->Copper;
             G4LogicalVolume* kname_log = new G4LogicalVolume(kname_solid,material,kname.str());
		     if(i == 30 || i == 31) {
		     	visAttributes = new G4VisAttributes(G4Colour(1.0, 0.6, 0.8, 0.5)); 
		     	kname_log->SetVisAttributes(visAttributes);
		     }
             d[0] = data[1];
             d[1] = data[2];
             d[2] = 0.0;

             m = m + 3;
           
             rotMatrix->rotateZ(rot_angles[m]*deg);
             rotMatrix->rotateY(rot_angles[m+1]*deg);
             rotMatrix->rotateX(rot_angles[m+2]*deg);
           
             G4VPhysicalVolume* kname_phys = new G4PVPlacement(rotMatrix,G4ThreeVector(pos[m]*cm,pos[m+1]*cm,pos[m+2]*cm),kname_log,kname.str(),WRLD_log,false,0,checkOverlaps);  
 	         }
         } 
 	  }

void DRAGONDetectorConstruction::ugeo_end(G4double pos[3])
     {
      //**********************************************************************
      //*                                                                      *
      //*                   Define a RAYTRACE final volume                     *
      //*                                                                      *
      //************************************************************************
	
	  G4int MCP;

      G4double shape[3], dead, det;

	  shape[0] = 2.4;
      shape[1] = 2.4;
      shape[2] = 0.1;
      det = shape[2];
	  
      G4VisAttributes* visAttributes = new G4VisAttributes(G4Colour(0., 0., 1., 0.5)); 
      visAttributes->SetVisibility(true); 
      visAttributes->SetForceSolid(true);
      
      Materials* materials = Materials::Instance();
      G4Material* material;  
	
	  //--------------------ENDV----------------------// 
      G4VSolid* ENDV_solid = new G4Box("ENDV",shape[0]*cm,shape[1]*cm,shape[2]*cm);
      //TMED->19 
      material = materials->Silicon;
      G4LogicalVolume* ENDV_log = new G4LogicalVolume(ENDV_solid, material,"ENDV");
      ENDV_log->SetVisAttributes(visAttributes);
      G4ThreeVector ENDV_pos(pos[0]*cm,pos[1]*cm,pos[2]*cm);
      G4VPhysicalVolume* ENDV_phys = new G4PVPlacement(0,ENDV_pos,ENDV_log,"ENDV",WRLD_log,false,0,checkOverlaps); 
      
      //C     Position deadlayer at surface of volume.
      //C     Note: since GEANT cannot track very accurately through
      //C     thin layers, create a new medium which is modified silicon,
      //C     make it 1000 x thicker but proportionately less dense so that the 
      //C     energy loss calculated is correct.
      dead = 100*3.5E-05; //!! 0.35 micrometres * 100.
      shape[2] = dead/2.;
      //--------------------DEAD----------------------// 
      G4VSolid* DEAD_solid = new G4Box("DEAD",shape[0]*cm,shape[1]*cm,shape[2]*cm);
      //TMED->23
      material = materials->ModifiedSilicon;
      G4LogicalVolume* DEAD_log = new G4LogicalVolume(DEAD_solid, material,"DEAD");
      DEAD_log->SetVisAttributes(visAttributes);
      G4ThreeVector DEAD_pos(pos[0]*cm,pos[1]*cm,pos[2]*cm+(det+shape[2])*cm);
      G4VPhysicalVolume* DEAD_phys = new G4PVPlacement(0,DEAD_pos,DEAD_log,"DEAD",WRLD_log,false,0,checkOverlaps);  
      
      G4UserLimits* DEADLimits = new G4UserLimits(0.00004*cm);   //From ugstmed.f
      DEAD_log->SetUserLimits(DEADLimits);
      
      //C     Add in MCP foil 20 microgram/cm**2, 25.4 cm diameter, 50 cm from 
      //C     DSSSD
      MCP = 0;
      if(MCP == 1)
        {
		 //C     First make foil holder (Aluminium) 10cm x 10 cm x 0.1 cm
		 //--------------------HOLD----------------------// 
		 shape[0] = 5.0;
         shape[1] = 5.0;
         shape[2] = 0.05;
       
         G4VSolid* HOLD_solid = new G4Box("HOLD",shape[0]*cm,shape[1]*cm,shape[2]*cm);
         //TMED->6 
         material = materials->Aluminium;
         G4LogicalVolume* HOLD_log = new G4LogicalVolume(HOLD_solid, material,"HOLD");
         HOLD_log->SetVisAttributes(visAttributes);
         //C     Position a 25.4 mm hole (vacuum) into the holder
         //--------------------HOLE----------------------// 
		 shape[0] = 0.0;
         shape[1] = 1.0;    //!!1.27
         shape[2] = 0.05; 
         
         G4VSolid* HOLE_solid = new G4Tubs("HOLE",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
         //TMED->1 
         material = materials->Vacuum;
         G4LogicalVolume* HOLE_log = new G4LogicalVolume(HOLE_solid, material,"HOLE");
         HOLE_log->SetVisAttributes(visAttributes);
         //G4ThreeVector HOLE_pos(0.0*cm,0.0*cm,0.0*cm);
         //G4VPhysicalVolume* HOLE_phys = new G4PVPlacement(0,HOLE_pos,HOLE_log,"HOLE",HOLD_log,false,0,checkOverlaps);  
         //C     Position the MCP foil into the holder
         //--------------------MCP----------------------// 
		 shape[0] = 0.0;
         shape[1] = 1.27;
         shape[2] = 8.889E-03;  //!! 20 microgram/cm**2 / 2.25 g/cm**3 x 1000
         
         G4VSolid* MCP_solid = new G4Tubs("MCP",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
         //TMED->24  
         material = materials->MCP_Carbon;
         G4LogicalVolume* MCP_log = new G4LogicalVolume(MCP_solid, material,"MCP");
         MCP_log->SetVisAttributes(visAttributes);
         //G4ThreeVector MCP_pos(0.0*cm,0.0*cm,0.0*cm);
         //G4VPhysicalVolume* MCP_phys = new G4PVPlacement(0,MCP_pos,MCP_log,"MCP",HOLE_log,false,0,checkOverlaps);  
         //C     Position MCP at -50cm from DSSSD
         //G4ThreeVector HOLD_pos(pos[0],pos[1],pos[2]+50.*cm);
         //G4LogicalVolume* WRLD_log = volumeStore->GetVolume("WRLD"); 
         //G4VPhysicalVolume* HOLD_phys = new G4PVPlacement(0,HOLD_pos,HOLD_log,"HOLD",WRLD_log_log,false,0,checkOverlaps);  
         G4UserLimits* MCPLimits = new G4UserLimits(0.00004*cm);  //From ugstmed.f
         MCP_log->SetUserLimits(MCPLimits);
         }
      }
      

void DRAGONDetectorConstruction::ugeo_start(G4double pos[3])
     {
	  //************************************************************************
      //*                                                                      *
      //*                   Define a RAYTRACE start volume                     *
      //*                                                                      *
      //************************************************************************
      G4double shape[3]; 
      
      Materials* materials = Materials::Instance();
      G4Material* material; 
    
      shape[0] = 100.;
      shape[1] = 100.;
      shape[2] = 1.0;
      
	  G4VSolid* STRV_solid = new G4Box("STRV",shape[0]*cm,shape[1]*cm,shape[2]*cm);
      //TMED->1
      material = materials->Vacuum;
      G4LogicalVolume* STRV_log = new G4LogicalVolume(STRV_solid,material,"STRV");
      G4ThreeVector STRV_pos(pos[0]*cm,pos[1]*cm,pos[2]*cm);
      G4VPhysicalVolume* STRV_phys = new G4PVPlacement(0,STRV_pos,STRV_log,"STRV",WRLD_log,false,0,checkOverlaps);
	  }


void DRAGONDetectorConstruction::ugeo_dssd(G4double pos[3])
     {
	  //************************************************************************
      //*                                                                      *
      //*                   Define a DSSD detector                             *
      //*                                                                      *
      //************************************************************************ 
      G4double shape[3];
      
      Materials* materials = Materials::Instance();
      G4Material* material;  

      shape[0] = 4.8;
      shape[1] = 4.8;
      shape[2] = 0.1;

	  G4VSolid* DSSD_solid = new G4Box("DSSD",shape[0]*cm,shape[1]*cm,shape[2]*cm);
      //TMED->19
      material = materials->Silicon;
      G4LogicalVolume* DSSD_log = new G4LogicalVolume(DSSD_solid,material,"DSSD");
      G4ThreeVector DSSD_pos(pos[0]*cm,pos[1]*cm,pos[2]*cm);
      G4VPhysicalVolume* DSSD_phys = new G4PVPlacement(0,DSSD_pos,DSSD_log,"DSSD",WRLD_log,false,0,checkOverlaps);
	  }

void DRAGONDetectorConstruction::ugeo_test(G4int k,G4double pos[3], G4double rot_angles[3])
     {
	  //************************************************************************
      //*                                                                      *
      //*                   Define a DSSD test volume (alpha acceptance tests) *
      //*                                                                      *
      //************************************************************************

      G4double shape[3];
      std::ostringstream kname;
      
      Materials* materials = Materials::Instance();
      G4Material* material; 
  	  
  	  G4RotationMatrix* rotMatrix = new G4RotationMatrix();
  	  rotMatrix->rotateZ(rot_angles[0]*deg);
      rotMatrix->rotateY(rot_angles[1]*deg);
      rotMatrix->rotateX(rot_angles[2]*deg);
  	  
	  if(k < 10) 
	    {
		 kname << "TST" << k;
         } 
      else 
         {
			kname << "TS" << k;
          }
	  //--------------------TST_ or TS_----------------------// 
	  shape[0] = 4.8;
      shape[1] = 4.8;
      shape[2] = 0.1;
	  
	  G4VSolid* solid = new G4Box(kname.str(),shape[0]*cm,shape[1]*cm,shape[2]*cm);
      //TMED->1
      material = materials->Vacuum;
      G4LogicalVolume* log = new G4LogicalVolume(solid,material,kname.str());
      G4ThreeVector TST_pos(pos[0]*cm,pos[1]*cm,pos[2]*cm);
      G4VPhysicalVolume* phys = new G4PVPlacement(rotMatrix,TST_pos,log,kname.str(),WRLD_log,false,k,checkOverlaps);  
      }

void DRAGONDetectorConstruction::ugeo_mcp(G4double pos[3],G4double data[2],G4String rname,G4int nmcp)
     {
	  
	  //************************************************************************
      //*                                                                      *
      //*                   Define a MCP detector                              *
      //*                                                                      *
      //************************************************************************

      G4double shape[3];
      
      G4VisAttributes* visAttributes = new G4VisAttributes(G4Colour(0., 1.0, 0., 0.5)); 
      visAttributes->SetVisibility(true); 
      visAttributes->SetForceSolid(true);
      
      Materials* materials = Materials::Instance();
      G4Material* material;

      std::ostringstream kname, lname;
      
      kname << "HLD" << nmcp;
      lname << "HOL" << nmcp;
           
      //C        G4cout << "ugeo_mcp for MCP " << rname << G4endl;
      //C        G4cout << "diameter " << data[0] << G4endl;
      //C        G4cout << "thickness ug/cm^2 " << data[1] << G4endl;
      //C        G4cout << "thickness cm " << data[1]/2250. << G4endl;
      //C        G4cout << "holder name " << kname << G4endl;
      //C        G4cout << "hole name " << lname << G4endl;
      //C     First make foil holder (Aluminium) 10cm x 10 cm x 0.1 cm
      //--------------------HLD----------------------// 
      shape[0] = 5.;
      shape[1] = 5.;
      shape[2] = 0.05;

      //--------------------HLD----------------------//
      G4VSolid* kname_solid = new G4Box(kname.str(),shape[0]*cm,shape[1]*cm,shape[2]*cm);
      //TMED->6
      material = materials->Aluminium;
      G4LogicalVolume* kname_log = new G4LogicalVolume(kname_solid,material,kname.str());
      kname_log->SetVisAttributes(visAttributes);
      //C     Position a hole (vacuum) into the holder
      //--------------------HOL----------------------// 
      shape[0] = 0.;
      shape[1] = data[0];                //!! radius in cm   
      shape[2] = 0.05;
      
      G4VSolid* lname_solid = new G4Tubs(lname.str(),shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->1 
      material = materials->Vacuum;
      G4LogicalVolume* lname_log = new G4LogicalVolume(lname_solid, material,lname.str());
      lname_log->SetVisAttributes(visAttributes);
      G4ThreeVector lname_pos(0.0*cm,0.0*cm,0.0*cm);
      G4VPhysicalVolume* lname_phys = new G4PVPlacement(0,lname_pos,lname_log,lname.str(),kname_log,false,0,checkOverlaps);  
      //C     Position the MCP foil into the holder
      shape[0] = 0.;
      shape[1] = data[0];                //!! radius in cm 
      shape[2] = data[1]/(2250.);     //!!  microgram/cm**2 / 2.25 g/cm**3 x 1000     
      
      G4VSolid* rname_solid = new G4Tubs(rname,shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->24 
      material = materials->MCP_Carbon;
      G4LogicalVolume* rname_log = new G4LogicalVolume(rname_solid, material,rname);
      rname_log->SetVisAttributes(visAttributes);
      G4ThreeVector rname_pos(0.0*cm,0.0*cm,0.0*cm);
      G4VPhysicalVolume* rname_phys = new G4PVPlacement(0,rname_pos,rname_log,rname,lname_log,false,0,checkOverlaps);  
 	  //C     Position MCP at posn
 	  G4ThreeVector kname_pos(pos[0]*cm,pos[1]*cm,pos[2]*cm);
 	  G4VPhysicalVolume* kname_phys = new G4PVPlacement(0,kname_pos,kname_log,kname.str(),WRLD_log,false,0,checkOverlaps);  
      }


void DRAGONDetectorConstruction::ugeo_fcup(G4double pos[3])
     {
      //************************************************************************
      //*                                                                      *
      //*                   Define a Farraday Cup                              *
      //*                                                                      *
      //************************************************************************

      G4double shape[3], p1, p2, p3; 
      
      Materials* materials = Materials::Instance();
      G4Material* material;  

      shape[0] = 0.;
      shape[1] = 2.5;
      shape[2] = 2.;
	  
      p1 = pos[0]; 
      p2 = pos[1];
      p3 = pos[2] + shape[2];  
      //G4cout << "posns " << " " << pos[0] << "," << pos[1] << "," << pos[2] << "," << p1 << "," << p2 << "," << p3 << G4endl;

      G4VSolid* FCUP_solid = new G4Tubs("FCUP",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->19
      material = materials->Silicon;
      G4LogicalVolume* FCUP_log = new G4LogicalVolume(FCUP_solid, material,"FCUP");
      G4ThreeVector FCUP_pos(p1*cm,p2*cm,p3*cm);
      G4RotationMatrix* rotMatrix = new G4RotationMatrix();
      rotMatrix->rotateZ(0.*deg);
      rotMatrix->rotateY(70.*deg);
      rotMatrix->rotateX(0.*deg);
      G4VPhysicalVolume* FCUP_phys = new G4PVPlacement(rotMatrix,FCUP_pos,FCUP_log,"FCUP",WRLD_log,false,0,checkOverlaps);      
      }

void DRAGONDetectorConstruction::udetmitray()
{
//C.
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//C                                                                      C
//C     UDETMITRAY is a user routine which defines certain volumes to be C
//C     detectors and defines what information each detector collects.   C
//C     For convenience the detector type variable (IDTYPE) is used      C
//C     for some volumes:                                                C
//C                                                                      C
//C                       1 - Detector Crystal                           C
//C                       2 - PM Tube                                    C
//C                       3 - Silicon                                    C
//C                       4 - Slits/Collimators                          C
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//C.
//C.   3 - Silicon
//C.
      G4SDManager* DRAGONSDManager = G4SDManager::GetSDMpointer();
      auto DRAGONSDDSSD = new DRAGONSensitiveDetector("DSSD",this);
      DRAGONSDManager->AddNewDetector(DRAGONSDDSSD);
      G4LogicalVolumeStore* volumeStore = G4LogicalVolumeStore::GetInstance();
      G4LogicalVolume* ENDV_log = volumeStore->GetVolume("ENDV");
      
      if (ENDV_log) 
         {ENDV_log->SetSensitiveDetector(DRAGONSDDSSD);} 
      else 
         {G4cerr << "Error: No logical volume named 'ENDV' found!" << G4endl;}
}



}





































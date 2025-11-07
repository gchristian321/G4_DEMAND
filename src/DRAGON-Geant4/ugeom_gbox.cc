//************************************************************************
//*                                                                      *
//*                          Define the Detector                         *
//*                                                                      *
//************************************************************************
 
#include "G4Box.hh"                          //Geant4
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Trd.hh"
#include "G4VSolid.hh"
#include "G4Polyhedra.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"
#include "rescom.hh"    //OJO
#include "G4SDManager.hh"
#include "G4DigiManager.hh"
#include "G4UserLimits.hh"

#include "DRAGONDetectorConstruction.hh"    //local
#include "DRAGONPrimaryGeneratorAction.hh"
#include "DRAGONSensitiveDetector.hh"
#include "DRAGONPhysicsList.hh"
#include "DRAGONDigitizer.hh"
#include "Materials.hh"
           

namespace DRAGON {

void DRAGONDetectorConstruction::ugeo_detector() 
     {
//C *** Local variables
      
      G4double the1, phi1, the2, phi2, the3, phi3;
      
      G4double x, y, z, shape[4];
      G4double bpa_len, bpb_len, bpc_len, bpd_len, bpe_len, 
               bpf_len, bpg_len, bph_len, bpi_len, bpj_len, bpa_in_len, 
               bpj_len2;
      G4double box_height, beam_height;
      
      Materials* materials = Materials::Instance();
      G4Material* material;


//C     std::cout << "ugeom_detector" << std::endl;      
//C.    Define Geometry
      ugeo_defin();

      if(1)//targtype == 1) || DRAGON_phys->Getalpha())
        {
         mcent = 1;
         mbox = 1;
         ment[0] = 1;
         ment[1] = 1;
         ment[2] = 1;
         mex[0] = 1;
         mex[1] = 1;
         mex[2] = 1;
         }
   
      // if(DRAGON_phys->Getalpha()) 
      //   {
      //    mtarg = 1;	
      //    materials->Target = materials->Vacuum;   
		  // }

      hexagon_small_width = s_finger + 2.*d_air[0] + 2.*d_mtl + air_gap;
      hexagon_large_width = 2. * hexagon_small_width / std::sqrt(3.);
      depth = z_finger + d_air[1] + d_mtl;
      col_length = 15.24;
        
//C.--->the x width has an offset of 3.35 to account for two detectors
//C.--->being pulled back by 6.7cm to make room for lead collimators   
      shape[0] = depth + 3.35;                              //! half x width;
      shape[1] = 5.  * hexagon_large_width;                 //! half y width
      shape[2] = 3. * hexagon_small_width + col_length/2.;  //! half z width;
      shape[2] = TLrms;
//G4cout << "-------------- HERE1 ---------------" << G4endl;
			G4cout << "DETE shape[0,1,2],(2),TLrms: "
						 << shape[0] << ", " << shape[1] << ", " << shape[2]<< ", "
						 << 3. * hexagon_small_width + col_length/2. << ", " << TLrms << G4endl;
      
      //--------------------DETE----------------------// 
      G4VSolid* DETE_solid = new G4Box("DETE",shape[0]*cm,shape[1]*cm,shape[2]*cm);
      //TRACKING MEDIA NUMBER
      //TMED->22
      material = materials->CentralVacuum;
      G4LogicalVolume* DETE_log = new G4LogicalVolume(DETE_solid, material,"DETE");
      G4VPhysicalVolume* DETE_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),DETE_log,"DETE",WRLD_log,false,0,checkOverlaps);  
      G4UserLimits* DETELimits = new G4UserLimits(10.0*cm);  //From ugstmed_trgt.f
      DETE_log->SetUserLimits(DETELimits);
    
      box_length = 17.069;
      box_height = 20.0;
      beam_height = - box_height/2. + 3.171;
//G4cout << "-------------- HERE2 ---------------" << G4endl;
     
      //--------------------CMBR----------------------// 
      shape[0] = box_width/2.;
      shape[1] = box_height/2.;
      shape[2] = box_length/2.;
   
      G4VSolid* CMBR_solid = new G4Box("CMBR",shape[0]*cm,shape[1]*cm,shape[2]*cm);
	  std::cout << "CMBR shape[]" << shape[0]/cm << " " << shape[1]/cm  << " " << shape[2]/cm << std::endl;
      //TMED->20
      material = materials->StainlessSteel;    
      G4LogicalVolume* CMBR_log = new G4LogicalVolume(CMBR_solid, material,"CMBR");
      
      G4UserLimits* CMBRLimits = new G4UserLimits(10.0*cm);  //From ugstmed.f
      CMBR_log->SetUserLimits(CMBRLimits);
//G4cout << "-------------- HERE3 ---------------" << G4endl;
      
      //--------------------CMBG----------------------// 
      shape[0] = shape[0] - wall[1];
      shape[1] = shape[1] - wall[1]/2.;
      shape[2] = shape[2] - wall[1];
	 
      G4VSolid* CMBG_solid = new G4Box("CMBG",shape[0]*cm,shape[1]*cm,shape[2]*cm);
      //TMED->35
      material = materials->Baseline;
      G4LogicalVolume* CMBG_log = new G4LogicalVolume(CMBG_solid, material,"CMBG");
	  CMBG_log->SetUserLimits(CMBRLimits);  //From ugstmed_trgt.f
//C.  CMBG needs to be ONLY for solid target sims
      if (targtype == 0)
	  {G4VPhysicalVolume* CMBG_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),CMBG_log,"CMBG",CMBR_log,false,0,checkOverlaps);}
      else
	     {	  
	      if (targtype == 1)
		     {G4VPhysicalVolume* CMBG_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),CMBG_log,"CMBG",CMBR_log,false,0,checkOverlaps);}
	      } 
     
      y = beam_height; 
//c   CR changed to ONLY
      G4ThreeVector CMBR_pos(0.0*cm,y*cm,0.0*cm);
      G4VPhysicalVolume* CMBR_phys = new G4PVPlacement(0,CMBR_pos,CMBR_log,"CMBR",DETE_log,false,0,checkOverlaps);  
//G4cout << "-------------- HERE4 ---------------" << G4endl;
    
//C.---> Collimator hole through outer aluminum  gas cell box
//C. changed material & z -> -z from DG's sim - CR. Need irot?
      //--------------------UHOL----------------------// 
      shape[0] = 0.0;
      shape[1] = 0.4;
      shape[2] = wall[1]/2.;  
  
      G4VSolid* UHOL_solid = new G4Tubs("UHOL",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->29
      material = materials->Entrance1;
      G4LogicalVolume* UHOL_log = new G4LogicalVolume(UHOL_solid, material,"UHOL");
      y = -beam_height;
      z = box_length/2. - wall[1]/2.;
      G4ThreeVector UHOL_pos(0.0*cm,y*cm,-z*cm);
      G4VPhysicalVolume* UHOL_phys = new G4PVPlacement(0,UHOL_pos,UHOL_log,"UHOL",CMBR_log,false,0,checkOverlaps); 
      G4UserLimits* UHOLLimits = new G4UserLimits(10.0*cm);  //From ugstmed_trgt.f
      UHOL_log->SetUserLimits(UHOLLimits);  
          
//C.---> Aluminum Collimator on the inside of the gas cell box  
      //--------------------PUAI----------------------//  
      bpa_in_len = 0.472;
      shape[0] = 0.0;
      shape[1] = 1.905;
      shape[2] = bpa_in_len;
  
      G4VSolid* PUAI_solid = new G4Tubs("PUAI",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->6
      material = materials->Aluminium;
      G4LogicalVolume* PUAI_log = new G4LogicalVolume(PUAI_solid, material,"PUAI");
      y = -beam_height;
      z = box_length/2. - wall[1] - bpa_in_len;
      G4ThreeVector PUAI_pos(0.0*cm,y*cm,-z*cm);
      G4VPhysicalVolume* PUAI_phys = new G4PVPlacement(0,PUAI_pos,PUAI_log,"PUAI",CMBG_log,false,0,checkOverlaps);  
//G4cout << "-------------- HERE5 ---------------" << G4endl;
    
      //--------------------PUBI----------------------//  
      shape[0] = 0.0;
      shape[1] = 0.4;
      shape[2] = bpa_in_len;
  
      G4VSolid* PUBI_solid = new G4Tubs("PUBI",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->29
      material = materials->Entrance1;
      G4LogicalVolume* PUBI_log = new G4LogicalVolume(PUBI_solid, material,"PUBI");
      G4VPhysicalVolume* PUBI_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PUBI_log,"PUBI",PUAI_log,false,0,checkOverlaps);  
      PUBI_log->SetUserLimits(UHOLLimits);  //From ugstmed_trgt.f
 
 //C. ---> collimator end collar detail outside box
      bpa_len = 0.472;
      bpb_len = 1.437;
      bpc_len = 0.321;
      bpd_len = 0.159;
      bpe_len = 2.060;
      bpf_len = 0.499;
      bpg_len = 0.476;
      bph_len = 0.980;
      bpi_len = 0.585;
      bpj_len = 0.350;
//    bpj_len2 = 0.467;
//G4cout << "-------------- HERE6 ---------------" << G4endl;

//--------------------PUA1----------------------// 
      shape[0] = 0.0;
      shape[1] = 1.905;
      shape[2] = bpa_len;                //!end collar detail half thickness
  
      G4VSolid* PUA1_solid = new G4Tubs("PUA1",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->6
      material = materials->Aluminium;
      G4LogicalVolume* PUA1_log = new G4LogicalVolume(PUA1_solid, material,"PUA1");
      z =  box_length/2. + bpa_len;
      G4ThreeVector PUA1_pos(0.0*cm,0.0*cm,-z*cm);
      G4VPhysicalVolume* PUA1_phys = new G4PVPlacement(0,PUA1_pos,PUA1_log,"PUA1",DETE_log,false,0,checkOverlaps);  
      //--------------------PUA2----------------------// 
      shape[0] = 0.0;
      shape[1] = 0.4;
      shape[2] = bpa_len;                //!end collar detail half thickness
  
      G4VSolid* PUA2_solid = new G4Tubs("PUA2",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->29
      material = materials->Entrance1;
      G4LogicalVolume* PUA2_log = new G4LogicalVolume(PUA2_solid, material,"PUA2");
      G4VPhysicalVolume* PUA2_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PUA2_log,"PUA2",PUA1_log,false,0,checkOverlaps);  
      PUA2_log->SetUserLimits(UHOLLimits);  //From ugstmed_trgt.f
   
//C.----> section BPB next to the left of collimator detail outside box
      //--------------------PUB1----------------------// 
      shape[0] = 0.0;
      shape[1] = 1.353;
      shape[2] = bpb_len;                //!end collar detail half thickness
  
      G4VSolid* PUB1_solid = new G4Tubs("PUB1",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->7
      material = materials->Lead;
      G4LogicalVolume* PUB1_log = new G4LogicalVolume(PUB1_solid, material,"PUB1");
      z =  box_length/2. + 2.*bpa_len + bpb_len;
      G4ThreeVector PUB1_pos(0.0*cm,0.0*cm,-z*cm);
      G4VPhysicalVolume* PUB1_phys = new G4PVPlacement(0,PUB1_pos,PUB1_log,"PUB1",DETE_log,false,0,checkOverlaps);  
   
      //--------------------PUB2----------------------// 
      shape[0] = 0.0;
      shape[1] = 1.035;
      shape[2] = bpb_len;                //!end collar detail half thickness
  
      G4VSolid* PUB2_solid = new G4Tubs("PUB2",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->6
      material = materials->Aluminium;
      G4LogicalVolume* PUB2_log = new G4LogicalVolume(PUB2_solid, material,"PUB2");
      G4VPhysicalVolume* PUB2_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PUB2_log,"PUB2",PUB1_log,false,0,checkOverlaps);  
      //--------------------PUB3----------------------// 
      shape[0] = 0.0;
      shape[1] = 0.4;
      shape[2] = bpb_len;                //!end collar detail half thickness
  
      G4VSolid* PUB3_solid = new G4Tubs("PUB3",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->29
      material = materials->Entrance1;
      G4LogicalVolume* PUB3_log = new G4LogicalVolume(PUB3_solid, material,"PUB3");
      G4VPhysicalVolume* PUB3_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PUB3_log,"PUB3",PUB2_log,false,0,checkOverlaps);  
      PUB3_log->SetUserLimits(UHOLLimits);  //From ugstmed_trgt.f
      //--------------------PUC1----------------------// 
      shape[0] = 0.0;
      shape[1] = 1.353;
      shape[2] = bpc_len;                //!end collar detail half thickness
//G4cout << "-------------- HERE7 ---------------" << G4endl;  
      G4VSolid* PUC1_solid = new G4Tubs("PUC1",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->7
      material = materials->Lead;
      G4LogicalVolume* PUC1_log = new G4LogicalVolume(PUC1_solid, material,"PUC1");
      z =  box_length/2. + 2.*bpa_len + 2.*bpb_len + bpc_len;
      G4ThreeVector PUC1_pos(0.0*cm,0.0*cm,-z*cm);
      G4VPhysicalVolume* PUC1_phys = new G4PVPlacement(0,PUC1_pos,PUC1_log,"PUC1",DETE_log,false,0,checkOverlaps);  
      //--------------------PUC2----------------------// 
      shape[0] = 0.0;
      shape[1] = 1.035;
      shape[2] = bpc_len;                //!end collar detail half thickness
  
      G4VSolid* PUC2_solid = new G4Tubs("PUC2",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->6
      material = materials->Aluminium;
      G4LogicalVolume* PUC2_log = new G4LogicalVolume(PUC2_solid, material,"PUC2");
      G4VPhysicalVolume* PUC2_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PUC2_log,"PUC2",PUC1_log,false,0,checkOverlaps);  
      //--------------------PUC3----------------------// 
      shape[0] = 0.0;
      shape[1] = 0.45;
      shape[2] = bpc_len;                //!end collar detail half thickness
  
      G4VSolid* PUC3_solid = new G4Tubs("PUC3",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->30
      material = materials->Entrance2;
      G4LogicalVolume* PUC3_log = new G4LogicalVolume(PUC3_solid, material,"PUC3");
      G4VPhysicalVolume* PUC3_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PUC3_log,"PUC3",PUC2_log,false,0,checkOverlaps);  
      PUC3_log->SetUserLimits(UHOLLimits);  //From ugstmed_trgt.f
     
//C.----> section BPD next in line to left of box
      //--------------------PUD1----------------------// 
      shape[0] = 0.0;
      shape[1] = 2.88;
      shape[2] = bpd_len;                //!end collar detail half thickness
  
      G4VSolid* PUD1_solid = new G4Tubs("PUD1",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->7
      material = materials->Lead;
      G4LogicalVolume* PUD1_log = new G4LogicalVolume(PUD1_solid, material,"PUD1");
      z =  box_length/2. + 2.*bpa_len + 2.*bpb_len + 2.*bpc_len + bpd_len;
      G4ThreeVector PUD1_pos(0.0*cm,0.0*cm,-z*cm);
      G4VPhysicalVolume* PUD1_phys = new G4PVPlacement(0,PUD1_pos,PUD1_log,"PUD1",DETE_log,false,0,checkOverlaps);  
      //--------------------PUD2----------------------// 
      shape[0] = 0.0;
      shape[1] = 1.035;
      shape[2] = bpd_len;                //!end collar detail half thickness
  
      G4VSolid* PUD2_solid = new G4Tubs("PUD2",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->6
      material = materials->Aluminium;
      G4LogicalVolume* PUD2_log = new G4LogicalVolume(PUD2_solid, material,"PUD2");
      G4VPhysicalVolume* PUD2_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PUD2_log,"PUD2",PUD1_log,false,0,checkOverlaps);  
      //--------------------PUD3----------------------// 
      shape[0] = 0.0;
      shape[1] = 0.45;
      shape[2] = bpd_len;                //!end collar detail half thickness
  
      G4VSolid* PUD3_solid = new G4Tubs("PUD3",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->30
      material = materials->Entrance2;
      G4LogicalVolume* PUD3_log = new G4LogicalVolume(PUD3_solid, material,"PUD3");
      G4VPhysicalVolume* PUD3_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PUD3_log,"PUD3",PUD2_log,false,0,checkOverlaps);  
      PUD3_log->SetUserLimits(UHOLLimits);  //From ugstmed_trgt.f
    
//C.----> section BPE next in line to left of box  
      //--------------------PUE1----------------------// 
      shape[0] = 0.0;
      shape[1] = 3.20;
      shape[2] = bpe_len;                //!end collar detail half thickness
  //G4cout << "-------------- HERE8 ---------------" << G4endl;
      G4VSolid* PUE1_solid = new G4Tubs("PUE1",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->7
      material = materials->Lead;
      G4LogicalVolume* PUE1_log = new G4LogicalVolume(PUE1_solid, material,"PUE1");
      z =  box_length/2. + 2.*bpa_len + 2.*bpb_len + 2.*bpc_len + 2.*bpd_len + bpe_len;
      G4ThreeVector PUE1_pos(0.0*cm,0.0*cm,-z*cm);
      G4VPhysicalVolume* PUE1_phys = new G4PVPlacement(0,PUE1_pos,PUE1_log,"PUE1",DETE_log,false,0,checkOverlaps);  
      //--------------------PUE2----------------------// 
      shape[0] = 0.0;
      shape[1] = 2.53;
      shape[2] = bpe_len;                //!end collar detail half thickness
  
      G4VSolid* PUE2_solid = new G4Tubs("PUE2",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->8
      material = materials->Air;
      G4LogicalVolume* PUE2_log = new G4LogicalVolume(PUE2_solid, material,"PUE2");
      G4VPhysicalVolume* PUE2_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PUE2_log,"PUE2",PUE1_log,false,0,checkOverlaps);  
      //--------------------PUE3----------------------// 
      shape[0] = 0.0;
      shape[1] = 2.09;
      shape[2] = bpe_len;                //!end collar detail half thickness
  
      G4VSolid* PUE3_solid = new G4Tubs("PUE3",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->6
      material = materials->Aluminium;
      G4LogicalVolume* PUE3_log = new G4LogicalVolume(PUE3_solid, material,"PUE3");
      G4VPhysicalVolume* PUE3_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PUE3_log,"PUE3",PUE2_log,false,0,checkOverlaps);  
      //--------------------PUE4----------------------// 
      shape[0] = 0.0;
      shape[1] = 0.45;
      shape[2] = bpe_len;                //!end collar detail half thickness
  
      G4VSolid* PUE4_solid = new G4Tubs("PUE4",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->30
      material = materials->Entrance2;
      G4LogicalVolume* PUE4_log = new G4LogicalVolume(PUE4_solid, material,"PUE4");
      G4VPhysicalVolume* PUE4_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PUE4_log,"PUE4",PUE3_log,false,0,checkOverlaps);  
      PUE4_log->SetUserLimits(UHOLLimits);  //From ugstmed_trgt.f
         
 //C.----> section BPF next in line to left of box      
      //--------------------PUF1----------------------// 
      shape[0] = 0.0;
      shape[1] = 3.2;
      shape[2] = bpf_len;                
  
      G4VSolid* PUF1_solid = new G4Tubs("PUF1",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->7
      material = materials->Lead;
      G4LogicalVolume* PUF1_log = new G4LogicalVolume(PUF1_solid, material,"PUF1");
      z =  box_length/2. + 2.*bpa_len + 2.*bpb_len + 2.*bpc_len + 2.*bpd_len + 2.*bpe_len + bpf_len;
      G4ThreeVector PUF1_pos(0.0*cm,0.0*cm,-z*cm);
      G4VPhysicalVolume* PUF1_phys = new G4PVPlacement(0,PUF1_pos,PUF1_log,"PUF1",DETE_log,false,0,checkOverlaps);  
      //--------------------PUF2----------------------// 
      shape[0] = 0.0;
      shape[1] = 2.53;
      shape[2] = bpf_len;               
  
      G4VSolid* PUF2_solid = new G4Tubs("PUF2",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->6
      material = materials->Aluminium;
      G4LogicalVolume* PUF2_log = new G4LogicalVolume(PUF2_solid, material,"PUF2");
      G4VPhysicalVolume* PUF2_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PUF2_log,"PUF2",PUF1_log,false,0,checkOverlaps);  
      //--------------------PUF3----------------------// 
      shape[0] = 0.0;
      shape[1] = 2.09;
      shape[2] = bpf_len;               
  
      G4VSolid* PUF3_solid = new G4Tubs("PUF3",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->6
      G4LogicalVolume* PUF3_log = new G4LogicalVolume(PUF3_solid, material,"PUF3");
      G4VPhysicalVolume* PUF3_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PUF3_log,"PUF3",PUF2_log,false,0,checkOverlaps);  
      //--------------------PUF4----------------------// 
      shape[0] = 0.0;
      shape[1] = 1.25;
      shape[2] = bpf_len;               
  
      G4VSolid* PUF4_solid = new G4Tubs("PUF4",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->31
      material = materials->Entrance3;
      G4LogicalVolume* PUF4_log = new G4LogicalVolume(PUF4_solid, material,"PUF4");
      G4VPhysicalVolume* PUF4_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PUF4_log,"PUF4",PUF3_log,false,0,checkOverlaps);  
      PUF4_log->SetUserLimits(UHOLLimits);  //From ugstmed_trgt.f
      //--------------------PUF5----------------------// 
      shape[0] = 0.0;
      shape[1] = 1.04;
      shape[2] = bpf_len;               
  
      G4VSolid* PUF5_solid = new G4Tubs("PUF5",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->6
      material = materials->Aluminium;
      G4LogicalVolume* PUF5_log = new G4LogicalVolume(PUF5_solid, material,"PUF5");
      G4VPhysicalVolume* PUF5_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PUF5_log,"PUF5",PUF4_log,false,0,checkOverlaps);  
      //--------------------PUF6----------------------// 
      shape[0] = 0.0;
      shape[1] = 0.5;
      shape[2] = bpf_len;               
  
      G4VSolid* PUF6_solid = new G4Tubs("PUF6",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->31
      material = materials->Entrance3;
      G4LogicalVolume* PUF6_log = new G4LogicalVolume(PUF6_solid, material,"PUF6");
      G4VPhysicalVolume* PUF6_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PUF6_log,"PUF6",PUF5_log,false,0,checkOverlaps); 
      PUF6_log->SetUserLimits(UHOLLimits);  //From ugstmed_trgt.f 
      //G4cout << "-------------- HERE9 ---------------" << G4endl;
//C.----> section BPG next in line to left of box
      //--------------------PUG1----------------------// 
      shape[0] = 7.75;
      shape[1] = 11.75;
      shape[2] = bpg_len;               
  
      G4VSolid* PUG1_solid = new G4Box("PUG1",shape[0]*cm,shape[1]*cm,shape[2]*cm);
      //TMED->7
      material = materials->Lead;
      G4LogicalVolume* PUG1_log = new G4LogicalVolume(PUG1_solid, material,"PUG1");
      z =  box_length/2. + 2.*bpa_len + 2.*bpb_len + 2.*bpc_len + 2.*bpd_len + 2.*bpe_len + 2.*bpf_len + bpg_len;
      G4ThreeVector PUG1_pos(0.0*cm,0.0*cm,-z*cm);
      G4VPhysicalVolume* PUG1_phys = new G4PVPlacement(0,PUG1_pos,PUG1_log,"PUG1",DETE_log,false,0,checkOverlaps);  
      //--------------------PUG2----------------------// 
      shape[0] = 0.0;
      shape[1] = 2.53;
      shape[2] = bpg_len;               
  
      G4VSolid* PUG2_solid = new G4Tubs("PUG2",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->6
      material = materials->Aluminium;
      G4LogicalVolume* PUG2_log = new G4LogicalVolume(PUG2_solid, material,"PUG2");
      G4VPhysicalVolume* PUG2_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PUG2_log,"PUG2",PUG1_log,false,0,checkOverlaps);  
      //--------------------PUG3----------------------// 
      shape[0] = 0.0;
      shape[1] = 1.25;
      shape[2] = bpg_len;               
  
      G4VSolid* PUG3_solid = new G4Tubs("PUG3",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->29
      material = materials->Entrance1;
      G4LogicalVolume* PUG3_log = new G4LogicalVolume(PUG3_solid, material,"PUG3");
      G4VPhysicalVolume* PUG3_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PUG3_log,"PUG3",PUG2_log,false,0,checkOverlaps);  
      PUG3_log->SetUserLimits(UHOLLimits);  //From ugstmed_trgt.f
      //--------------------PUG4----------------------// 
      shape[0] = 0.0;
      shape[1] = 1.04;
      shape[2] = bpg_len;               
  
      G4VSolid* PUG4_solid = new G4Tubs("PUG4",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->6
      material = materials->Aluminium;
      G4LogicalVolume* PUG4_log = new G4LogicalVolume(PUG4_solid, material,"PUG4");
      G4VPhysicalVolume* PUG4_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PUG4_log,"PUG4",PUG3_log,false,0,checkOverlaps);  
      //--------------------PUG5----------------------// 
      shape[0] = 0.0;
      shape[1] = 0.5;
      shape[2] = bpg_len;               
  
      G4VSolid* PUG5_solid = new G4Tubs("PUG5",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->31
      material = materials->Entrance3;
      G4LogicalVolume* PUG5_log = new G4LogicalVolume(PUG5_solid, material,"PUG5");
      G4VPhysicalVolume* PUG5_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PUG5_log,"PUG5",PUG4_log,false,0,checkOverlaps);  
      PUG5_log->SetUserLimits(UHOLLimits);  //From ugstmed_trgt.f 
      
//C.----> section BPH next in line to left of box    
      //--------------------PUH1----------------------// 
      shape[0] = 0.0;
      shape[1] = 5.71;
      shape[2] = bph_len;               
  
      G4VSolid* PUH1_solid = new G4Tubs("PUH1",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->6
      material = materials->Aluminium;
      G4LogicalVolume* PUH1_log = new G4LogicalVolume(PUH1_solid, material,"PUH1");
      z =  box_length/2. + 2.*bpa_len + 2.*bpb_len + 2.*bpc_len + 2.*bpd_len + 2.*bpe_len + 2.*bpf_len + 2.*bpg_len + bph_len;
      G4ThreeVector PUH1_pos(0.0*cm,0.0*cm,-z*cm);
      G4VPhysicalVolume* PUH1_phys = new G4PVPlacement(0,PUH1_pos,PUH1_log,"PUH1",DETE_log,false,0,checkOverlaps);   
      //--------------------PUH2----------------------// 
      shape[0] = 0.0;
      shape[1] = 1.25;
      shape[2] = bph_len;               
  
      G4VSolid* PUH2_solid = new G4Tubs("PUH2",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->29
      material = materials->Entrance1;
      G4LogicalVolume* PUH2_log = new G4LogicalVolume(PUH2_solid, material,"PUH2");
      G4VPhysicalVolume* PUH2_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PUH2_log,"PUH2",PUH1_log,false,0,checkOverlaps);  
      PUH2_log->SetUserLimits(UHOLLimits); //From ugstmed_trgt.f
      //--------------------PUH3----------------------// 
      shape[0] = 0.0;
      shape[1] = 1.04;
      shape[2] = bph_len;               
  
      G4VSolid* PUH3_solid = new G4Tubs("PUH3",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->6
      material = materials->Aluminium;
      G4LogicalVolume* PUH3_log = new G4LogicalVolume(PUH3_solid, material,"PUH3");
      G4VPhysicalVolume* PUH3_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PUH3_log,"PUH3",PUH2_log,false,0,checkOverlaps);  
      //--------------------PUH4----------------------// 
      shape[0] = 0.0;
      shape[1] = 0.5;
      shape[2] = bph_len;               
  
      G4VSolid* PUH4_solid = new G4Tubs("PUH4",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->31
      material = materials->Entrance3;
      G4LogicalVolume* PUH4_log = new G4LogicalVolume(PUH4_solid, material,"PUH4");
      G4VPhysicalVolume* PUH4_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PUH4_log,"PUH4",PUH3_log,false,0,checkOverlaps);  
      PUH4_log->SetUserLimits(UHOLLimits);  //From ugstmed_trgt.f 
      
//C.----> section BPI next in line to left of box    
      //--------------------PUI1----------------------// 
      shape[0] = 0.0;
      shape[1] = 2.53;
      shape[2] = bpi_len;               
  
      G4VSolid* PUI1_solid = new G4Tubs("PUI1",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->6
      material = materials->Aluminium;
      G4LogicalVolume* PUI1_log = new G4LogicalVolume(PUI1_solid, material,"PUI1");
      z =  box_length/2. + 2.*bpa_len + 2.*bpb_len + 2.*bpc_len + 2.*bpd_len + 2.*bpe_len + 2.*bpf_len + 2.*bpg_len + 2.*bph_len + bpi_len;
      G4ThreeVector PUI1_pos(0.0*cm,0.0*cm,-z*cm);  
      G4VPhysicalVolume* PUI1_phys = new G4PVPlacement(0,PUI1_pos,PUI1_log,"PUI1",DETE_log,false,0,checkOverlaps);  
      //--------------------PUI2----------------------// 
      shape[0] = 0.0;
      shape[1] = 1.25;
      shape[2] = bpi_len;               
  
      G4VSolid* PUI2_solid = new G4Tubs("PUI2",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->29
      material = materials->Entrance1;
      G4LogicalVolume* PUI2_log = new G4LogicalVolume(PUI2_solid, material,"PUI2");
      G4VPhysicalVolume* PUI2_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PUI2_log,"PUI2",PUI1_log,false,0,checkOverlaps);  
      PUI2_log->SetUserLimits(UHOLLimits); //From ugstmed_trgt.f
      //--------------------PUI3----------------------// 
      shape[0] = 0.0;
      shape[1] = 1.04;
      shape[2] = bpi_len;               
  
      G4VSolid* PUI3_solid = new G4Tubs("PUI3",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->6
      material = materials->Aluminium;
      G4LogicalVolume* PUI3_log = new G4LogicalVolume(PUI3_solid, material,"PUI3");
      G4VPhysicalVolume* PUI3_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PUI3_log,"PUI3",PUI2_log,false,0,checkOverlaps);  
      //--------------------PUI4----------------------// 
      shape[0] = 0.0;
      shape[1] = 0.5;
      shape[2] = bpi_len;               
  
      G4VSolid* PUI4_solid = new G4Tubs("PUI4",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->31
      material = materials->Entrance3;
      G4LogicalVolume* PUI4_log = new G4LogicalVolume(PUI4_solid, material,"PUI4");
      G4VPhysicalVolume* PUI4_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PUI4_log,"PUI4",PUI3_log,false,0,checkOverlaps); 
      PUI4_log->SetUserLimits(UHOLLimits);  //From ugstmed_trgt.f  
                
//C.----> section BPJ next in line to left of box
//C. to fill the gap between the end of the pumping tubes and the         
       //--------------------PUJ1----------------------// 
      shape[0] = 0.0;
      shape[1] = 1.04;
      shape[2] = bpj_len;               
	  
	  z =  box_length/2. + 2.*bpa_len + 2.*bpb_len + 2.*bpc_len + 2.*bpd_len + 2.*bpe_len + 2.*bpf_len + 2.*bpg_len + 2.*bph_len + 2.*bpi_len + bpj_len;
	  G4ThreeVector PUJ1_pos(0.0*cm,0.0*cm,-z*cm); 
	  G4VSolid* PUJ1_solid = new G4Tubs("PUJ1",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->6
      material = materials->Aluminium;
      G4LogicalVolume* PUJ1_log = new G4LogicalVolume(PUJ1_solid, material,"PUJ1");
      
      G4VPhysicalVolume* PUJ1_phys = new G4PVPlacement(0,PUJ1_pos,PUJ1_log,"PUJ1",DETE_log,false,0,checkOverlaps);  
      //--------------------PUJ2----------------------// 
      shape[0] = 0.0;
      shape[1] = 0.5;
      shape[2] = bpj_len;               
  
      G4VSolid* PUJ2_solid = new G4Tubs("PUJ2",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
      //TMED->31
      material = materials->Entrance3;
      G4LogicalVolume* PUJ2_log = new G4LogicalVolume(PUJ2_solid, material,"PUJ2");
      G4VPhysicalVolume* PUJ2_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PUJ2_log,"PUJ2",PUJ1_log,false,0,checkOverlaps); 
      PUJ2_log->SetUserLimits(UHOLLimits);  //From ugstmed_trgt.f  
       
//C.***************************************************************
//C.
//C.    Make the collimator assembly for the gas cell for DOWNSTREAM
//C.
//C.***************************************************************    
      //G4cout << "-------------- HERE6 ---------------" << G4endl;
                              //! opposite side rotating matrix
      
      the1 =  270.*deg;
      phi1 =  180.*deg;
      the2 =   90.*deg;
      phi2 =   90.*deg;
      the3 =  180.*deg;
      phi3 =    0.*deg;
      
 //C.      tubetype = 0     
      
      if(tubetype == 0 || (tubetype > 1 && tubetype < 7))
        {
//C.---> Collimator hole through outer aluminum  gas cell box DOWNSTREAM		  
		 //--------------------DHOL----------------------// 
         shape[0] = 0.0;
         shape[1] = 0.450;
         shape[2] = wall[1]/2.;             
  
         G4VSolid* DHOL_solid = new G4Tubs("DHOL",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
         //TMED->32
         material = materials->Ext1;
         G4LogicalVolume* DHOL_log = new G4LogicalVolume(DHOL_solid, material,"DHOL");
         y = -beam_height;
         z = -(box_length/2. - wall[1]/2.);
         G4ThreeVector DHOL_pos(0.0*cm,y*cm,-z*cm);
         G4VPhysicalVolume* DHOL_phys = new G4PVPlacement(0,DHOL_pos,DHOL_log,"DHOL",CMBR_log,false,0,checkOverlaps);  
         G4UserLimits* DHOLLimits = new G4UserLimits(10.0*cm);  //From ugstmed_trgt.f
         DHOL_log->SetUserLimits(DHOLLimits);
        
//C.---> Aluminum Collimator on the inside of the gas cell box	  
		 //--------------------PDAI----------------------// 
         shape[0]=0.0;
         shape[1]=1.905;
         shape[2]=bpa_in_len;       
 
         G4VSolid* PDAI_solid = new G4Tubs("PDAI",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
         //TMED->6
         material = materials->Aluminium;
         G4LogicalVolume* PDAI_log = new G4LogicalVolume(PDAI_solid, material,"PDAI");
         y = -beam_height;
         z = -(box_length/2. - wall[1] - bpa_in_len);
         G4ThreeVector PDAI_pos(0.0*cm,y*cm,-z*cm);
         G4VPhysicalVolume* PDAI_phys = new G4PVPlacement(0,PDAI_pos,PDAI_log,"PDAI",CMBG_log,false,0,checkOverlaps);  
         //--------------------PDBI----------------------//
         shape[0]=0.0;
         shape[1]=0.450;
         shape[2]=bpa_in_len;       
 
         G4VSolid* PDBI_solid = new G4Tubs("PDBI",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
         //TMED->32
         material = materials->Ext1;
         G4LogicalVolume* PDBI_log = new G4LogicalVolume(PDBI_solid, material,"PDBI");
         G4VPhysicalVolume* PDBI_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PDBI_log,"PDBI",PDAI_log,false,0,checkOverlaps);  
         PDBI_log->SetUserLimits(DHOLLimits);         //From ugstmed_trgt.f
         
         //--------------------PDA1----------------------//
         //---> collimator end collar detail outside box
         shape[0]=0.0;
         shape[1]=1.905;
         shape[2]=bpa_len;        //!end collar detail half thickness
 
         G4VSolid* PDA1_solid = new G4Tubs("PDA1",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
         //TMED->6
         material = materials->Aluminium;
         G4LogicalVolume* PDA1_log = new G4LogicalVolume(PDA1_solid, material,"PDA1");
         z = -(box_length/2. + bpa_len);
         G4ThreeVector PDA1_pos(0.0*cm,0.0*cm,-z*cm);
         G4VPhysicalVolume* PDA1_phys = new G4PVPlacement(0,PDA1_pos,PDA1_log,"PDA1",DETE_log,false,0,checkOverlaps);  
         //--------------------PDA2----------------------//
         shape[0]=0.0;
         shape[1]=0.450;
         shape[2]=bpa_len;        //!end collar detail half thickness
 
         G4VSolid* PDA2_solid = new G4Tubs("PDA2",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
         //TMED->32
         material = materials->Ext1;
         G4LogicalVolume* PDA2_log = new G4LogicalVolume(PDA2_solid, material,"PDA2");
         G4VPhysicalVolume* PDA2_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PDA2_log,"PDA2",PDA1_log,false,0,checkOverlaps); 
         PDA2_log->SetUserLimits(DHOLLimits);         //From ugstmed_trgt.f
         //--------------------PDB1----------------------//
         //----> section PDB next to the RIGHT of collimator detail outside box
         shape[0]=0.0;
         shape[1]=1.035;
         shape[2]=bpb_len;        //!end collar detail half thickness
 
         G4VSolid* PDB1_solid = new G4Tubs("PDB1",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
         //TMED->6
         material = materials->Aluminium;
         G4LogicalVolume* PDB1_log = new G4LogicalVolume(PDB1_solid, material,"PDB1");
         z =  -(box_length/2. + 2.*bpa_len + bpb_len);
         G4ThreeVector PDB1_pos(0.0*cm,0.0*cm,-z*cm);
         G4VPhysicalVolume* PDB1_phys = new G4PVPlacement(0,PDB1_pos,PDB1_log,"PDB1",DETE_log,false,0,checkOverlaps);
         //--------------------PDB2----------------------//
         shape[0]=0.0;
         shape[1]=0.450;
         shape[2]=bpb_len;        //!end collar detail half thickness
 
         G4VSolid* PDB2_solid = new G4Tubs("PDB2",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
         //TMED->32
         material = materials->Ext1;
         G4LogicalVolume* PDB2_log = new G4LogicalVolume(PDB2_solid, material,"PDB2");
         G4VPhysicalVolume* PDB2_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PDB2_log,"PDB2",PDB1_log,false,0,checkOverlaps); 
         PDB2_log->SetUserLimits(DHOLLimits);         //From ugstmed_trgt.f
         //--------------------PDD1----------------------//
         //----> section BPD next in line to left of box
         shape[0]=0.0;
         shape[1]=1.035;
         shape[2]=bpc_len + bpd_len;  //!end collar detail half thickness
 
         G4VSolid* PDD1_solid = new G4Tubs("PDD1",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
         //TMED->6
         material = materials->Aluminium;
         G4LogicalVolume* PDD1_log = new G4LogicalVolume(PDD1_solid, material,"PDD1");
         z = -(box_length/2. + 2.*bpa_len + 2.*bpb_len + bpc_len + bpd_len);
         G4ThreeVector PDD1_pos(0.0*cm,0.0*cm,-z*cm);
         G4VPhysicalVolume* PDD1_phys = new G4PVPlacement(0,PDD1_pos,PDD1_log,"PDD1",DETE_log,false,0,checkOverlaps); 
         //--------------------PDD2----------------------//
         shape[0]=0.0;
         shape[1]=0.520;
         shape[2]=bpc_len + bpd_len;  //!end collar detail half thickness
 
         G4VSolid* PDD2_solid = new G4Tubs("PDD2",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
         //TMED->33
         material = materials->Ext2;
         G4LogicalVolume* PDD2_log = new G4LogicalVolume(PDD2_solid, material,"PDD2");
         G4VPhysicalVolume* PDD2_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PDD2_log,"PDD2",PDD1_log,false,0,checkOverlaps); 
         PDD2_log->SetUserLimits(DHOLLimits);         //From ugstmed_trgt.f
   
         //----> section PDE next in line to RIGHT of box
 
         //--------------------PDE1----------------------//
         shape[0]=0.0;
         shape[1]=2.09;
         shape[2]=bpe_len;  //!end collar detail half thickness
 
         G4VSolid* PDE1_solid = new G4Tubs("PDE1",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
         //TMED->6
         material = materials->Aluminium;
         G4LogicalVolume* PDE1_log = new G4LogicalVolume(PDE1_solid, material,"PDE1");
         z =  -(box_length/2. + 2.*bpa_len + 2.*bpb_len + 2.*bpc_len + 2.*bpd_len + bpe_len);
         G4ThreeVector PDE1_pos(0.0*cm,0.0*cm,-z*cm);
         G4VPhysicalVolume* PDE1_phys = new G4PVPlacement(0,PDE1_pos,PDE1_log,"PDE1",DETE_log,false,0,checkOverlaps);
         //--------------------PDE2----------------------//
         shape[0]=0.0;
         shape[1]=0.520;
         shape[2]=bpe_len;  //!end collar detail half thickness
 
         G4VSolid* PDE2_solid = new G4Tubs("PDE2",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
         //TMED->33
         material = materials->Ext2;
         G4LogicalVolume* PDE2_log = new G4LogicalVolume(PDE2_solid, material,"PDE2");
         G4VPhysicalVolume* PDE2_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PDE2_log,"PDE2",PDE1_log,false,0,checkOverlaps);
         PDE2_log->SetUserLimits(DHOLLimits);         //From ugstmed_trgt.f

//       C.----> section PDF next in line to RIGHT of box
         
         //--------------------PDF1----------------------//
         shape[0]=0.0;
         shape[1]=2.53;
         shape[2]=bpf_len+bpg_len;
 
         G4VSolid* PDF1_tube = new G4Tubs("PDF1_tube",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
         //TMED->6
         material = materials->Aluminium;
         G4LogicalVolume* PDF1_log = new G4LogicalVolume(PDF1_tube, material,"PDF1");
         z = -(box_length/2. + 2.*bpa_len + 2.*bpb_len + 2.*bpc_len + 2.*bpd_len + 2.*bpe_len + bpf_len + bpg_len);
         G4ThreeVector PDF1_pos(0.0*cm,0.0*cm,-z*cm);
         G4VPhysicalVolume* PDF1_phys = new G4PVPlacement(0,PDF1_pos,PDF1_log,"PDF1",DETE_log,false,0,checkOverlaps);
         //--------------------PDF2----------------------//
         shape[0]=0.0;
         shape[1]=1.25;
         shape[2]=bpf_len+bpg_len;
 
         G4VSolid* PDF2_tube = new G4Tubs("PDF2_tube",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
         //TMED->34
         material = materials->Ext3;
         G4LogicalVolume* PDF2_log = new G4LogicalVolume(PDF2_tube, material,"PDF2");
         G4VPhysicalVolume* PDF2_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PDF2_log,"PDF2",PDF1_log,false,0,checkOverlaps);
         PDF2_log->SetUserLimits(DHOLLimits);         //From ugstmed_trgt.f
         //--------------------PDF3----------------------//
         shape[0]=0.0;
         shape[1]=1.04;
         shape[2]=bpf_len+bpg_len;
 
         G4VSolid* PDF3_tube = new G4Tubs("PDF3_tube",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
         //TMED->6
         material = materials->Aluminium;
         G4LogicalVolume* PDF3_log = new G4LogicalVolume(PDF3_tube, material,"PDF3");
         G4VPhysicalVolume* PDF3_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PDF3_log,"PDF3",PDF2_log,false,0,checkOverlaps);
         //--------------------PDF4----------------------//
         shape[0]=0.0;
         shape[1]=0.591;
         shape[2]=bpf_len+bpg_len;
 
         G4VSolid* PDF4_tube = new G4Tubs("PDF4_tube",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
         //TMED->34
         material = materials->Ext3;
         G4LogicalVolume* PDF4_log = new G4LogicalVolume(PDF4_tube, material,"PDF4");
         G4VPhysicalVolume* PDF4_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PDF4_log,"PDF4",PDF3_log,false,0,checkOverlaps);
         PDF4_log->SetUserLimits(DHOLLimits);         //From ugstmed_trgt.f
 
//       C.----> section BPH next in line to left of box

         //--------------------PDH1----------------------//
         shape[0]=0.0;
         shape[1]=5.71;
         shape[2]=bph_len;
 
         G4VSolid* PDH1_solid = new G4Tubs("PDH1",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
         //TMED->6
         material = materials->Aluminium;
         G4LogicalVolume* PDH1_log = new G4LogicalVolume(PDH1_solid, material,"PDH1");
         z =  -(box_length/2. + 2.*bpa_len + 2.*bpb_len + 2.*bpc_len + 2.*bpd_len + 2.*bpe_len + 2.*bpf_len + 2.*bpg_len + bph_len);
         G4ThreeVector PDH1_pos(0.0*cm,0.0*cm,-z*cm);
         G4VPhysicalVolume* PDH1_phys = new G4PVPlacement(0,PDH1_pos,PDH1_log,"PDH1",DETE_log,false,0,checkOverlaps);
         //--------------------PDH2----------------------//
         shape[0]=0.0;
         shape[1]=1.25;
         shape[2]=bph_len;
 
         G4VSolid* PDH2_solid = new G4Tubs("PDH2",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
         //TMED->34
         material = materials->Ext3;
         G4LogicalVolume* PDH2_log = new G4LogicalVolume(PDH2_solid, material,"PDH2");
         G4VPhysicalVolume* PDH2_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PDH2_log,"PDH2",PDH1_log,false,0,checkOverlaps);
         PDH2_log->SetUserLimits(DHOLLimits);         //From ugstmed_trgt.f
         //--------------------PDH3----------------------//
         shape[0]=0.0;
         shape[1]=1.04;
         shape[2]=bph_len;
 
         G4VSolid* PDH3_solid = new G4Tubs("PDH3",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
         //TMED->6
         material = materials->Aluminium;
         G4LogicalVolume* PDH3_log = new G4LogicalVolume(PDH3_solid, material,"PDH3");
         G4VPhysicalVolume* PDH3_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PDH3_log,"PDH3",PDH2_log,false,0,checkOverlaps);
         //--------------------PDH4----------------------//
         shape[0]=0.0;
         shape[1]=0.591;
         shape[2]=bph_len;
 
         G4VSolid* PDH4_solid = new G4Tubs("PDH4",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
         //TMED->34
         material = materials->Ext3;
         G4LogicalVolume* PDH4_log = new G4LogicalVolume(PDH4_solid, material,"PDH4");
         G4VPhysicalVolume* PDH4_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PDH4_log,"PDH4",PDH3_log,false,0,checkOverlaps);
         PDH4_log->SetUserLimits(DHOLLimits);         //From ugstmed_trgt.f

//       C.----> section BPI next in line to left of box
         //--------------------PDI1----------------------//
         shape[0]=0.0;
         shape[1]=2.53;
         shape[2]=bpi_len;
 
         G4VSolid* PDI1_solid = new G4Tubs("PDI1",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
         //TMED->6
         material = materials->Aluminium;
         G4LogicalVolume* PDI1_log = new G4LogicalVolume(PDI1_solid, material,"PDI1");
         z =  -(box_length/2. + 2.*bpa_len + 2.*bpb_len + 2.*bpc_len + 2.*bpd_len + 2.*bpe_len + 2.*bpf_len + 2.*bpg_len + 2.*bph_len + bpi_len);
         G4ThreeVector PDI1_pos(0.0*cm,0.0*cm,-z*cm);
         G4VPhysicalVolume* PDI1_phys = new G4PVPlacement(0,PDI1_pos,PDI1_log,"PDI1",DETE_log,false,0,checkOverlaps);
         //--------------------PDI2----------------------//
         shape[0]=0.0;
         shape[1]=1.25;
         shape[2]=bpi_len;
 
         G4VSolid* PDI2_solid = new G4Tubs("PDI2",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
         //TMED->34
         material = materials->Ext3;
         G4LogicalVolume* PDI2_log = new G4LogicalVolume(PDI2_solid, material,"PDI2");
         G4VPhysicalVolume* PDI2_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PDI2_log,"PDI2",PDI1_log,false,0,checkOverlaps);
         PDI2_log->SetUserLimits(DHOLLimits);         //From ugstmed_trgt.f
         //--------------------PDI3----------------------//
         shape[0]=0.0;
         shape[1]=1.04;
         shape[2]=bpi_len;
 
         G4VSolid* PDI3_solid = new G4Tubs("PDI3",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
         //TMED->6
         material = materials->Aluminium;
         G4LogicalVolume* PDI3_log = new G4LogicalVolume(PDI3_solid, material,"PDI3");
         G4VPhysicalVolume* PDI3_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PDI3_log,"PDI3",PDI2_log,false,0,checkOverlaps);
         //--------------------PDI4----------------------//
         shape[0]=0.0;
         shape[1]=0.591;
         shape[2]=bpi_len;
 
         G4VSolid* PDI4_solid = new G4Tubs("PDI4",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
         //TMED->34
         material = materials->Ext3;
         G4LogicalVolume* PDI4_log = new G4LogicalVolume(PDI4_solid, material,"PDI4");
         G4VPhysicalVolume* PDI4_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PDI4_log,"PDI4",PDI3_log,false,0,checkOverlaps); 
         PDI4_log->SetUserLimits(DHOLLimits);         //From ugstmed_trgt.f
//--------------------PDJ1----------------------//
//C.
//C.----> section BPJ next in line to left of box
//C. to fill the gap between the end of the pumping tubes and the
//c      shape[0] = 0.0;
//c      shape[1] = 1.04;
//c      shape[2] = bpj_len2;
//c      z =  box_length/2. + 2.*bpa_len + 2.*bpb_len + 2.*bpc_len + 2.*bpd_len + 2.*bpe_len + 2.*bpf_len + 2.*bpg_len + 2.*bph_len + 2.*bpi_len + bpj_len2;
//c      G4VSolid* PDJ1_solid = new G4Tubs("PDJ1",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
//c      TMED->29
//c      material = materials->Aluminium;
//c      G4LogicalVolume* PDJ1_log = new G4LogicalVolume(PDJ1_solid, material,"PDJ1");
//c      G4VPhysicalVolume* PDJ1_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,z*cm),PDJ1_log,"PDJ1",DETE_log,false,0,checkOverlaps);
//--------------------PDJ2----------------------//
//c      shape[0] = 0.0;
//c      shape[1] = 0.591;
//c      shape[2] = bpj_len2;
//c      G4VSolid* PDJ2_solid = new G4Tubs("PDJ2",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
//c      TMED->29
//c      material = materials->Entrance1;
//c      G4LogicalVolume* PDJ2_log = new G4LogicalVolume(PDJ2_solid, material,"PDJ2");
//c      G4VPhysicalVolume* PDJ2_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0.0*cm),PDJ2_log,"PDJ2",PDJ1_log,false,0,checkOverlaps);
  	     }
      else
         {
		  if (tubetype == 1)
		     {
              //C. Collar on right side of box wall
              //c         col_collar_length = 0.41
              //--------------------CLRD----------------------//
              shape[0]=0.0;
              shape[1]=1.91;
              shape[2]=0.41/2.;
 
              G4VSolid* CLRD_solid = new G4Tubs("CLRD",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
              //TMED->6
              material = materials->Aluminium;
              G4LogicalVolume* CLRD_log = new G4LogicalVolume(CLRD_solid, material,"CLRD");
              y = -beam_height;
              z = box_length/2. - wall[1] - 0.41/2.;
              //--------------------CL1G----------------------//
              shape[1]=0.642;
              G4VSolid* CL1G_solid = new G4Tubs("CL1G",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
              //TMED->32
              material = materials->Ext1;
              G4LogicalVolume* CL1G_log = new G4LogicalVolume(CL1G_solid, material,"CL1G");
              G4ThreeVector CL1G_pos(0.0*cm,0.0*cm,0.0*cm);
              G4VPhysicalVolume* CL1G_phys = new G4PVPlacement(0,CL1G_pos,CL1G_log,"CL1G",CLRD_log,false,0,checkOverlaps); 
  	          G4ThreeVector CLRD_pos(0.0*cm,y*cm,z*cm);
  	          G4VPhysicalVolume* CLRD_phys = new G4PVPlacement(0,CLRD_pos,CLRD_log,"CLRD",CMBG_log,false,0,checkOverlaps); 
  	          G4UserLimits* CL1GLimits = new G4UserLimits(10.0 * cm);  //From ugstmed_trgt.f
              CL1G_log->SetUserLimits(CL1GLimits);
//c           shape[1] = 1.905;
//c           shape[2] = 0.214;
//c           G4VSolid* CLRB_solid = new G4Tubs("CLRB",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
//c           y = -beam_height;
//c           z = box_length/2. - wall[1] - 0.517/2. - 0.517/2. - 0.214;
//c           shape[3] = 0.8474;
//c           G4VSolid* CT1G_solid = new G4Tubs("CT1G",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
//c           G4ThreeVector CT1G_pos(0.0*cm,0.0*cm,0.0*cm);
//c           G4VPhysicalVolume* CT1G_phys = new G4PVPlacement(0,CT1G_pos,CT1G_log,"CT1G",CLRB_log,false,0,checkOverlaps); 
//c           G4ThreeVector CLRB_pos(0.0*cm,y*cm,z*cm);
//c           G4VPhysicalVolume* CLRB_phys = new G4PVPlacement(0,CLRB_pos,CLRB_log,"CLRB",CMBG_log,false,0,checkOverlaps);  	     
			  //C. Cut through box wall
		      //--------------------APTU----------------------//
              shape[0]=0.0;
              shape[1]=0.642;
              shape[2]= wall[1]/2.;
 
              G4VSolid* APTU_solid = new G4Tubs("APTU",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
              //TMED->32
              material = materials->Ext1;
              G4LogicalVolume* APTU_log = new G4LogicalVolume(APTU_solid, material,"APTU");
              y = -beam_height;
              z = box_length/2. - wall[1]/2.;
              G4ThreeVector APTU_pos(0.0*cm,y*cm,z*cm);
              G4VPhysicalVolume* APTU_phys = new G4PVPlacement(0,APTU_pos,APTU_log,"APTU",CMBR_log,false,0,checkOverlaps);
              APTU_log->SetUserLimits(CL1GLimits);   //From ugstmed_trgt.f
			
			  //C. Collar on right side of box
		      //--------------------CLLD----------------------//
              shape[0]=0.0;
              shape[1]=1.91;
              shape[2]=0.945/2.;
 
              G4VSolid* CLLD_solid = new G4Tubs("CLLD",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
              //TMED->6
              material = materials->Aluminium;
              G4LogicalVolume* CLLD_log = new G4LogicalVolume(CLLD_solid, material,"CLLD");
              z = box_length/2. + 0.945/2.;
              //--------------------CL2G----------------------//
              shape[1]=0.642;
              G4VSolid* CL2G_solid = new G4Tubs("CL2G",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
              //TMED->32
              material = materials->Ext1;
              G4LogicalVolume* CL2G_log = new G4LogicalVolume(CL2G_solid, material,"CL2G");
              G4VPhysicalVolume* CL2G_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),CL2G_log,"CL2G",CLLD_log,false,0,checkOverlaps);
		      G4ThreeVector CLLD_pos(0.0*cm,0.0*cm,z*cm);
              G4VPhysicalVolume* CLLD_phys = new G4PVPlacement(0,CLLD_pos,CLLD_log,"CLLD",DETE_log,false,0,checkOverlaps);
              CL2G_log->SetUserLimits(CL1GLimits);   //From ugstmed_trgt.f
              //C. First part of collimator tube			
              //--------------------CM1T----------------------//
              shape[0]=0.0;
              shape[1]=1.27;
              shape[2]=1.026/2.;
 
              G4VSolid* CM1T_solid = new G4Tubs("CM1T",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
              //TMED->6
              material = materials->Aluminium;
              G4LogicalVolume* CM1T_log = new G4LogicalVolume(CM1T_solid, material,"CM1T");
              z = box_length/2. + 0.945 + 1.026/2.;
              //--------------------CM1G----------------------//
              shape[1]=0.642;
              G4VSolid* CM1G_solid = new G4Tubs("CM1G",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
              //TMED->32
              material = materials->Ext1;
              G4LogicalVolume* CM1G_log = new G4LogicalVolume(CM1G_solid, material,"CM1G");
              G4VPhysicalVolume* CM1G_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),CM1G_log,"CM1G",CM1T_log,false,0,checkOverlaps);	
              G4ThreeVector CM1T_pos(0.0*cm,0.0*cm,z*cm);
              G4VPhysicalVolume* CM1T_phys = new G4PVPlacement(0,CM1T_pos,CM1T_log,"CM1T",DETE_log,false,0,checkOverlaps);
              CM1G_log->SetUserLimits(CL1GLimits);   //From ugstmed_trgt.f			
			  //C. Second part of collimator tube
              //--------------------CM2T----------------------//
              shape[0]=0.0;
              shape[1]=1.27;
              shape[2]=2.54/2.;
 
              G4VSolid* CM2T_solid = new G4Tubs("CM2T",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
              //TMED->6
              material = materials->Aluminium;
              G4LogicalVolume* CM2T_log = new G4LogicalVolume(CM2T_solid, material,"CM2T");
              z = box_length/2. + 0.945 + 1.026 + 2.54/2.;
              //--------------------CM2G----------------------//
              shape[1]=0.715;
              G4VSolid* CM2G_solid = new G4Tubs("CM2G",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
              //TMED->32
              material = materials->Ext1;
              G4LogicalVolume* CM2G_log = new G4LogicalVolume(CM2G_solid, material,"CM2G");
              G4VPhysicalVolume* CM2G_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),CM2G_log,"CM2G",CM2T_log,false,0,checkOverlaps);	
              G4ThreeVector CM2T_pos(0.0*cm,0.0*cm,z*cm);
              G4VPhysicalVolume* CM2T_phys = new G4PVPlacement(0,CM2T_pos,CM2T_log,"CM2T",DETE_log,false,0,checkOverlaps);
              CM2G_log->SetUserLimits(CL1GLimits);   //From ugstmed_trgt.f				
			
              //C. Tube up to first part of sheilding assembly	
              
              
              //--------------------CM3T----------------------//
              shape[0]=0.0;
              shape[1]=1.27;
              shape[2]=0.624/2.;
 
              G4VSolid* CM3T_solid = new G4Tubs("CM3T",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
              //TMED->6
              material = materials->Aluminium;
              G4LogicalVolume* CM3T_log = new G4LogicalVolume(CM3T_solid, material,"CM3T");
              z = box_length/2. + 0.945 + 1.026 + 2.54/2. + 0.624/2.;
              //--------------------CM3G----------------------//
              shape[1]=0.95;
              G4VSolid* CM3G_solid = new G4Tubs("CM3G",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
              //TMED->33
              material = materials->Ext2;
              G4LogicalVolume* CM3G_log = new G4LogicalVolume(CM3G_solid, material,"CM3G");
              G4VPhysicalVolume* CM3G_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),CM3G_log,"CM3G",CM3T_log,false,0,checkOverlaps);	
              G4ThreeVector CM3T_pos(0.0*cm,0.0*cm,z*cm);
              G4VPhysicalVolume* CM3T_phys = new G4PVPlacement(0,CM3T_pos,CM3T_log,"CM3T",DETE_log,false,0,checkOverlaps);	
              CM3G_log->SetUserLimits(CL1GLimits);   //From ugstmed_trgt.f	
			  //--------------------PDE1----------------------//
              shape[0]=0.0;
              shape[1]=2.09;
              shape[2]=3.54/2.;                     
 
              G4VSolid* PDE1_solid = new G4Tubs("PDE1",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
              //TMED->6
              material = materials->Aluminium;
              G4LogicalVolume* PDE1_log = new G4LogicalVolume(PDE1_solid, material,"PDE1");
              z = -(box_length/2. + 0.945 + 1.026 + 2.54 + 0.624 + 3.54/2.);
              G4ThreeVector PDE1_pos(0.0*cm,0.0*cm,-z*cm);
              G4VPhysicalVolume* PDE1_phys = new G4PVPlacement(0,PDE1_pos,PDE1_log,"PDE1",DETE_log,false,0,checkOverlaps);		
              //--------------------PDE2----------------------//
              shape[0]=0.0;
              shape[1]=0.95;
              shape[2]=3.54/2.;                     
 
              G4VSolid* PDE2_solid = new G4Tubs("PDE2",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
              //TMED->33
              material = materials->Ext2;
              G4LogicalVolume* PDE2_log = new G4LogicalVolume(PDE2_solid, material,"PDE2");
              G4VPhysicalVolume* PDE2_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PDE2_log,"PDE2",PDE1_log,false,0,checkOverlaps);	
              PDE2_log->SetUserLimits(CL1GLimits);   //From ugstmed_trgt.f	
			
	          //C.----> section PDF next in line to RIGHT of box		
	          //--------------------PDF1----------------------//
              shape[0]=0.0;
              shape[1]=2.53;
              shape[2]=2.88/2.;                     
 
              G4VSolid* PDF1_solid = new G4Tubs("PDF1",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
              //TMED->6
              material = materials->Aluminium;
              G4LogicalVolume* PDF1_log = new G4LogicalVolume(PDF1_solid, material,"PDF1");
              z = -(box_length/2. + 0.945 + 1.026 + 2.54 + 0.624 + 3.54 + 2.88/2.);
              G4ThreeVector PDF1_pos(0.0*cm,0.0*cm,-z*cm);
              G4VPhysicalVolume* PDF1_phys = new G4PVPlacement(0,PDF1_pos,PDF1_log,"PDF1",DETE_log,false,0,checkOverlaps);

              //--------------------PDF2----------------------//
              shape[0]=0.0;
              shape[1]=1.48;
              shape[2]=2.88/2.;                    
 
              G4VSolid* PDF2_solid = new G4Tubs("PDF2_solid",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
              //TMED->34
              material = materials->Ext3;
              G4LogicalVolume* PDF2_log = new G4LogicalVolume(PDF2_solid, material,"PDF2");
              G4VPhysicalVolume* PDF2_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PDF2_log,"PDF2",PDF1_log,false,0,checkOverlaps);
              PDF2_log->SetUserLimits(CL1GLimits);   //From ugstmed_trgt.f	
			
              //--------------------PDF3----------------------//
              shape[0]=0.0;
              shape[1]=1.27;
              shape[2]=2.88/2.;                      
 
              G4VSolid* PDF3_solid = new G4Tubs("PDF3",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
              //TMED->6
              material = materials->Aluminium;
              G4LogicalVolume* PDF3_log = new G4LogicalVolume(PDF3_solid, material,"PDF3");
              G4VPhysicalVolume* PDF3_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PDF3_log,"PDF3",PDF2_log,false,0,checkOverlaps);
              //--------------------PDF4----------------------//
              shape[0]=0.0;
              shape[1]=0.95;
              shape[2]=2.88/2.;                      
 
              G4VSolid* PDF4_solid = new G4Tubs("PDF4",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
              //TMED->34
              material = materials->Ext3;
              G4LogicalVolume* PDF4_log = new G4LogicalVolume(PDF4_solid, material,"PDF4");
              G4VPhysicalVolume* PDF4_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PDF4_log,"PDF4",PDF3_log,false,0,checkOverlaps);
              PDF4_log->SetUserLimits(CL1GLimits);   //From ugstmed_trgt.f
									
              //C.----> section BPH next in line to left of box			
	          //--------------------PDH1----------------------//
              shape[0]=0.0;
              shape[1]=5.71;
              shape[2]=2.02/2.;                      
 
              G4VSolid* PDH1_solid = new G4Tubs("PDH1",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
              //c changed to Alum.
              //TMED->6
              material = materials->Aluminium;
              G4LogicalVolume* PDH1_log = new G4LogicalVolume(PDH1_solid, material,"PDH1");
              z = -(box_length/2. + 0.945 + 1.026 + 2.54 + 0.624 + 3.54 + 2.88 + 2.02/2.);
			  G4ThreeVector PDH1_pos(0.0*cm,0.0*cm,-z*cm);
              G4VPhysicalVolume* PDH1_phys = new G4PVPlacement(0,PDH1_pos,PDH1_log,"PDH1",DETE_log,false,0,checkOverlaps);	
              //--------------------PDH2----------------------//
              shape[0]=0.0;
              shape[1]=1.48;
              shape[2]=2.02/2.;                      
 
              G4VSolid* PDH2_solid = new G4Tubs("PDH2",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
              //c changed to Vac.
              //TMED->6
              material = materials->Ext3;
              G4LogicalVolume* PDH2_log = new G4LogicalVolume(PDH2_solid, material,"PDH2");
              G4VPhysicalVolume* PDH2_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PDH2_log,"PDH2",PDH1_log,false,0,checkOverlaps);
              PDH2_log->SetUserLimits(CL1GLimits);   //From ugstmed_trgt.f	
               //--------------------PDH3----------------------//
              shape[0]=0.0;
              shape[1]=1.27;
              shape[2]=2.02/2.;                      
 
              G4VSolid* PDH3_solid = new G4Tubs("PDH3",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
              //c changed to Vac.
              //TMED->6
              material = materials->Aluminium;
              G4LogicalVolume* PDH3_log = new G4LogicalVolume(PDH3_solid, material,"PDH3");
              G4VPhysicalVolume* PDH3_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PDH3_log,"PDH3",PDH2_log,false,0,checkOverlaps);	
              //--------------------PDH4----------------------//
              shape[0]=0.0;
              shape[1]=0.95;
              shape[2]=2.02/2.;                      
 
              G4VSolid* PDH4_solid = new G4Tubs("PDH4",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
              //c changed to Vac.
              //TMED->34
              material = materials->Ext3;
              G4LogicalVolume* PDH4_log = new G4LogicalVolume(PDH4_solid, material,"PDH4");
              G4VPhysicalVolume* PDH4_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PDH4_log,"PDH4",PDH3_log,false,0,checkOverlaps);
              PDH4_log->SetUserLimits(CL1GLimits);   //From ugstmed_trgt.f
		
		      //C.----> section BPI next in line to left of box	  
		      //--------------------PDI1----------------------//
		      shape[0]=0.0;
              shape[1]=2.53;
              shape[2]=bpi_len;                   
 
              G4VSolid* PDI1_solid = new G4Tubs("PDI1",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
              //c changed to Alum.
              //TMED->6
              material = materials->Aluminium;
              G4LogicalVolume* PDI1_log = new G4LogicalVolume(PDI1_solid, material,"PDI1");
              z = -(box_length/2. + 0.945 + 1.026 + 2.54 + 0.624 + 3.54 + 2.88 + 2.02 + bpi_len);
              G4ThreeVector PDI1_pos(0.0*cm,0.0*cm,-z*cm);
              G4VPhysicalVolume* PDI1_phys = new G4PVPlacement(0,PDI1_pos,PDI1_log,"PDI1",DETE_log,false,0,checkOverlaps);
              //--------------------PDI2----------------------//
		      shape[0]=0.0;
              shape[1]=1.48;
              shape[2]=bpi_len;                   
 
              G4VSolid* PDI2_solid = new G4Tubs("PDI2",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
              //c changed to Vac.
              //TMED->34
              material = materials->Ext3;
              G4LogicalVolume* PDI2_log = new G4LogicalVolume(PDI2_solid, material,"PDI2");
              G4VPhysicalVolume* PDI2_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PDI1_log,"PDI1",PDI1_log,false,0,checkOverlaps);
              PDI2_log->SetUserLimits(CL1GLimits);   //From ugstmed_trgt.f
              //--------------------PDI3----------------------//
		      shape[0]=0.0;
              shape[1]=1.27;
              shape[2]=bpi_len;                   
 
              G4VSolid* PDI3_solid = new G4Tubs("PDI3",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
              //c changed to Alum.
              //TMED->6
              material = materials->Aluminium;
              G4LogicalVolume* PDI3_log = new G4LogicalVolume(PDI3_solid, material,"PDI3");
              G4VPhysicalVolume* PDI3_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PDI3_log,"PDI3",PDI2_log,false,0,checkOverlaps);
              //--------------------PDI4----------------------//
		      shape[0]=0.0;
              shape[1]=0.95;
              shape[2]=bpi_len;                   
 
              G4VSolid* PDI4_solid = new G4Tubs("PDI4",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
              //c changed to Vac.
              //TMED->34
              material = materials->Ext3;
              G4LogicalVolume* PDI4_log = new G4LogicalVolume(PDI4_solid,material,"PDI4");
              G4VPhysicalVolume* PDI4_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),PDI4_log,"PDI4",PDI3_log,false,0,checkOverlaps);
              PDI4_log->SetUserLimits(CL1GLimits);   //From ugstmed_trgt.f

//C.
//C.----> section BPJ next in line to left of box
//C. to fill the gap between the end of the pumping tubes and the
//c      shape[0] = 0.0;
//c      shape[1] = 1.04;
//c      shape[2] = bpj_len2;
//c      z =  box_length/2. + 2.*bpa_len + 2.*bpb_len + 2.*bpc_len + 2.*bpd_len + 2.*bpe_len + 2.*bpf_len + 2.*bpg_len + 2.*bph_len + 2.*bpi_len + bpj_len2;
//c      G4VSolid* PDJ1_solid = new G4Tubs("PDJ1",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
//c      G4VPhysicalVolume* PDJ1_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,z*cm),PDJ1_log,"PDJ1",DETE_log,false,0,checkOverlaps);
//c      shape[0] = 0.0;
//c      shape[1] = 0.71755;
//c      shape[2] = bpj_len2;
//c      G4VSolid* PDJ2_solid = new G4Tubs("PDJ2",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
//c      G4VPhysicalVolume* PDJ2_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0.0*cm),PDJ2_log,"PDJ2",PDJ1_log,false,0,checkOverlaps);
			  }
	      }

      
//C.***************************************************************
//C.
//C.    Make the inner gas cell assembly
//C.
//C.***************************************************************      

//C.    Inner gas cell assembly only needed for gas targets
      if (targtype == 0)
         {   
//C.---> Aluminum outer shell of the cell      
          
          //--------------------CELL----------------------//
          shape[0] = 2.115;		//! bottom of trap auto calculate
          shape[1] = 6.759;		//! top of trap dimension
          shape[2] = 1.905;	    //! width of trap
          shape[3] = 4.208;		//! height of trap     
      
                                    //! orient the trapezoid box correctly
        
          the1 =  180.*deg;		    //! theta x
          phi1 =    0.*deg;		    //! phi x
          the2 =  270.*deg;		    //! theta y
          phi2 =  180.*deg;		    //! phi y
          the3 =   90.*deg;		    //! theta z
          phi3 =   90.*deg;		    //! phi z
          
          G4RotationMatrix* irot_box = new G4RotationMatrix;
          irot_box->rotateY(90*deg);
          irot_box->rotateX(90*deg);
      
          G4VSolid* CELL_solid = new G4Trd("CELL",shape[0]*cm,shape[1]*cm,shape[2]*cm,shape[2]*cm,shape[3]*cm);
          //TMED->6
          material = materials->Aluminium;
          G4LogicalVolume* CELL_log = new G4LogicalVolume(CELL_solid, material,"CELL");
          //C.--->move inner cell so that apertures are centred on beam line
          //C.--->0.976 offset is height from box top to inner cell top
          y = .5 * box_height - shape[3] - 0.976;
          G4ThreeVector CELL_pos(0.0*cm,y*cm,0.0*cm);
          G4VPhysicalVolume* CELL_phys = new G4PVPlacement(irot_box,CELL_pos,CELL_log,"CELL",CMBG_log,false,0,checkOverlaps);
			  
          //C.---> Vacuum inner shell of the cell            
          //--------------------CELG----------------------//
 		  shape[0]=1.798;
          shape[1]=6.094;
          shape[2]=1.588;    
          shape[3]=3.891;
         
          G4VSolid* CELG_solid = new G4Trd("CELG",shape[0]*cm,shape[1]*cm,shape[2]*cm,shape[2]*cm,shape[3]*cm);
          //TMED->27
		  material = materials->Target;
          G4LogicalVolume* CELG_log = new G4LogicalVolume(CELG_solid, material,"CELG");
          G4VPhysicalVolume* CELG_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0*cm),CELG_log,"CELG",CELL_log,false,0,checkOverlaps);
		
          //C.---> Entrance collimator to the trapezoid target  
          //--------------------EAPG----------------------// 
          shape[0]=0.0;
          if (tubetype == 0 || (tubetype > 1 || tubetype < 7))
             {shape[1] = 0.3;}
          else
             {
			  if(tubetype == 1)
			    {
				 //cc       MT updates gas target entrance collimator 21 Oct 2003.
                 shape[1] = 0.2;
			     }
			  }
          shape[2] = 0.5;
          G4VSolid* EAPG_solid = new G4Tubs("EAPG",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
                  
                     	          //! orient the collimators correctly

          the1 =  180.*deg;		//! theta x
          phi1 =    0.*deg;		//! phi x
          the2 =   90.*deg;		//! theta y
          phi2 =   90.*deg;		//! phi y
          the3 =   90.*deg;		//! theta z
          phi3 =    0.*deg;		//! phi z
          
          G4RotationMatrix* irot_col = new G4RotationMatrix;
          irot_col->rotateY(90*deg);
          z = 2.008;
          x = 5.315;
          G4ThreeVector EAPG_pos(x*cm,0.0*cm,z*cm);
		  G4LogicalVolume* EAPG_log = new G4LogicalVolume(EAPG_solid,material,"EAPG");   
          G4VPhysicalVolume* EAPG_phys = new G4PVPlacement(irot_col,EAPG_pos,EAPG_log,"EAPG",CELL_log,false,0,checkOverlaps);
          
          //C.---> Exit collimator to the trapezoid target
          //--------------------XAPG----------------------//
          shape[0]=0.0;
          if (tubetype == 0 || (tubetype > 1 || tubetype < 7))
             {shape[1] = 0.4;}
          else
             {
			  if(tubetype == 1)
			    {
                 //cc       MT updates gas target exit collimator 21 Oct 2003.
                 shape[1] = 0.5;
			     }
			  }
          shape[2] = 0.5;
     
          G4VSolid* XAPG_solid = new G4Tubs("XAPG",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
          //TMED->27
          material = materials->Target;
          G4LogicalVolume* XAPG_log = new G4LogicalVolume(XAPG_solid,material,"XAPG");
          z = 2.008;
          x = -5.315;
          G4ThreeVector XAPG_pos(x*cm,0.0*cm,z*cm);
          G4VPhysicalVolume* XAPG_phys = new G4PVPlacement(irot_col,XAPG_pos,XAPG_log,"XAPG",CELL_log,false,0,checkOverlaps);
					
          }
//C
//C.***************************************************************
//C.
//C.    Make the solid target disc
//C.
//C.***************************************************************     
//C.      Only need solid target disc for, ah, solid target work  
          if (targtype == 1) {
              y = -1.0 * beam_height;
			  
//C       entdens is length from targ entrance to targ center
//C       exitdens is length from targ exit to targ center
//C       Rather than defining a new variable for solid target half lengths (and
//C       be forced to change several other files) JS has used the same variables
//C       despite the poor name for solid targets.
//C
//C       The width and density are exaggerated 100x
          entdens = 1.111e-3;
          exitdens = 1.111e-3;
          shape[0] = 0.000;           //! inner disc radius
          shape[1] = 1.200;           //! outer disc radius
          shape[2] = entdens;            //! half width of disc
//C.      create rotation matrix to aline z-axis of disk with x-axis of CELG
//C.      If the inner gas cell assembly is ever used for solid targets
          
          //--------------------CTAR----------------------//
          G4VSolid* CTAR_solid = new G4Tubs("CTAR",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);
           //TMED->21
          material = materials->Target_Carbon;
          G4LogicalVolume* CTAR_log = new G4LogicalVolume(CTAR_solid,material,"CTAR");
          //G4ThreeVector CTAR_pos(0.0*cm,0.0*cm,z*cm);
          //G4VPhysicalVolume* CTAR_phys = new G4PVPlacement(irot_col,CTAR_pos,CTAR_log,"CTAR",CELG_log,false,0,checkOverlaps);
          G4ThreeVector CTAR_pos(0.0*cm,y*cm,0.0*cm);
          G4VPhysicalVolume* CTAR_phys = new G4PVPlacement(0,CTAR_pos,CTAR_log,"CTAR",CMBG_log,false,0,checkOverlaps);
          
          G4UserLimits* CTARLimits = new G4UserLimits(0.00004*cm);  //From ugstmed.f
          CTAR_log->SetUserLimits(CTARLimits);
          }

								//G4cout << "-------------- HERE8 ---------------" << "targtype:"<<targtype << G4endl;

	 }

void DRAGONDetectorConstruction::ugeo_finger()
     {
      //************************************************************************
      //*                                                                      *
      //*                  Define the BSO finger modules                       *
      //*                                                                      *
      //************************************************************************

	  //C *** Local variables

      G4int i, j, k, n, jm, icopy;
      G4double x, y, z, shape[10];
      G4double the1, phi1, the2, phi2, the3, phi3;
      
      Materials* materials = Materials::Instance();
      G4Material* material;

      G4LogicalVolumeStore* volumeStore = G4LogicalVolumeStore::GetInstance();
      G4LogicalVolume* DETE_log = volumeStore->GetVolume("DETE");
   /*
     std::cout << "d_mtl " << d_mtl << std::endl;
     std::cout << "d_air[0] " << d_air[0] << std::endl;
	  std::cout << "d_air[1] " << d_air[1] << std::endl;
	 std::cout << "z_finger " << z_finger << std::endl;
	 std::cout << "pmt_length " << pmt_length << std::endl;
	 std::cout << "hexagon_large_width " << hexagon_large_width << std::endl;
	 std::cout << "hexagon_small_width " << hexagon_small_width << std::endl;   
	 std::cout << "air_gap " << air_gap << std::endl;
	 std::cout << "aprt " << aprt << std::endl;
	 std::cout << "wall " << wall[0] << " " << wall[1] << " " << wall[2] << std::endl; 
	 std::cout << "box_width " << box_width << std::endl;
	 std::cout << "depth " << depth << std::endl;                                                       //OJO depth toma valores incorrectos en el codigo*/
	 
	 //C std::cout << "ugeo_finger " << std::endl;   
   
      //C.*** Make the housing for scintillator and PMT
      //--------------------HSNG----------------------//
      shape[0] =   0.0*deg;
      shape[1] = 360.0*deg;
      shape[2] =   6;
      shape[3] =   2;
      shape[4] = -((d_mtl+d_air[1]+z_finger+pmt_length)/2.)*cm;
      shape[5] =   0.0*cm;
      shape[6] = (hexagon_small_width/2.)*cm;
      shape[7] =  ((d_mtl+d_air[1]+z_finger+pmt_length)/2.)*cm;
      shape[8] =   0.0*cm;
      shape[9] = shape[6];

    /*  std::cout << "HSNG " << std::endl;
      std::cout << "shape[0] " << shape[0]/deg << std::endl;
      std::cout << "shape[1] " << shape[1]/deg << std::endl;
      std::cout << "shape[2] " << shape[2] << std::endl;
      std::cout << "shape[3] " << shape[3] << std::endl;
      std::cout << "shape[4] " << shape[4] << std::endl;
      std::cout << "shape[5] " << shape[5] << std::endl;
      std::cout << "shape[6] " << shape[6] << std::endl;
      std::cout << "shape[7] " << shape[7] << std::endl;
      std::cout << "shape[8] " << shape[8] << std::endl;
      std::cout << "shape[9] " << shape[9] << std::endl;*/
      
      G4double zPlane[]={shape[4],shape[7]};
      G4double rInner[]={shape[5],shape[8]};
      G4double rOuter[]={shape[6],shape[6]};

      G4VSolid* HSNG_solid = new G4Polyhedra("HSNG",shape[0],shape[1],shape[2],shape[3],zPlane,rInner,rOuter);
      //TMED->8
      material = materials->Air;
      G4LogicalVolume* HSNG_log = new G4LogicalVolume(HSNG_solid,material,"HSNG");
    
      G4VisAttributes* visAttributes = new G4VisAttributes(G4Colour(1.0, 0., 0., 0.)); 
      visAttributes->SetVisibility(true); 
      visAttributes->SetForceSolid(true);
      HSNG_log->SetVisAttributes(visAttributes);
      
      //     C *** Make the fingers ***
      //--------------------FNGR----------------------//
      shape[4] = shape[4] + (pmt_length/2.)*cm;
      shape[6] = shape[6] - (air_gap/2.)*cm;
      shape[7] = shape[7] - (pmt_length/2.)*cm;
      shape[9] = shape[9] - (air_gap/2.)*cm;
      
      zPlane[0]=shape[4];
      zPlane[1]=shape[7];
      rOuter[0]=shape[6];
      rOuter[1]=shape[6];

     /*
      std::cout << "FNGR " << std::endl;
      std::cout << "shape[0] " << shape[0]/deg << std::endl;
      std::cout << "shape[1] " << shape[1]/deg << std::endl;
      std::cout << "shape[2] " << shape[2] << std::endl;
      std::cout << "shape[3] " << shape[3] << std::endl;
      std::cout << "shape[4] " << shape[4] << std::endl;
      std::cout << "shape[5] " << shape[5] << std::endl;
      std::cout << "shape[6] " << shape[6] << std::endl;
      std::cout << "shape[7] " << shape[7] << std::endl;
      std::cout << "shape[8] " << shape[8] << std::endl;
      std::cout << "shape[9] " << shape[9] << std::endl;
      std::cout << "pos 0 0 " << z << std::endl;  
      */
     
      G4VSolid* FNGR_solid = new G4Polyhedra("FNGR",shape[0],shape[1],shape[2],shape[3],zPlane,rInner,rOuter);
      //TMED->6
      material = materials->Aluminium;
      G4LogicalVolume* FNGR_log = new G4LogicalVolume(FNGR_solid,material,"FNGR");
      visAttributes = new G4VisAttributes(G4Colour(0.1, 0.1, 0.1, 1.)); 
      FNGR_log->SetVisAttributes(visAttributes);
      z = -pmt_length/2.;
      G4ThreeVector FNGR_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* FNGR_phys = new G4PVPlacement(0,FNGR_pos,FNGR_log,"FNGR",HSNG_log,false,0,checkOverlaps);  

      //--------------------MGOR----------------------//
      shape[6] = shape[6] - d_mtl*cm;
      shape[9] = shape[9] - d_mtl*cm;
      shape[4] = shape[4] + (d_mtl/2.)*cm;
      shape[7] = shape[7] - (d_mtl/2.)*cm;
   
      G4VSolid* MGOR_solid = new G4Polyhedra("MGOR",shape[0],shape[1],shape[2],shape[3],zPlane,rInner,rOuter);
      //TMED->16
      material = materials->MgO;    //OJO    OJO
      G4LogicalVolume* MGOR_log = new G4LogicalVolume(MGOR_solid,material,"MGOR");
      z = d_mtl/2.;
      G4ThreeVector MGOR_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* MGOR_phys = new G4PVPlacement(0,MGOR_pos,MGOR_log,"MGOR",FNGR_log,false,0,checkOverlaps);

      /*
      std::cout << "MGOR " << std::endl;
      std::cout << "shape[0] " << shape[0]/deg << std::endl;
      std::cout << "shape[1] " << shape[1]/deg << std::endl;
      std::cout << "shape[2] " << shape[2] << std::endl;
      std::cout << "shape[3] " << shape[3] << std::endl;
      std::cout << "shape[4] " << shape[4] << std::endl;
      std::cout << "shape[5] " << shape[5] << std::endl;
      std::cout << "shape[6] " << shape[6] << std::endl;
      std::cout << "shape[7] " << shape[7] << std::endl;
      std::cout << "shape[8] " << shape[8] << std::endl;
      std::cout << "shape[9] " << shape[9] << std::endl;
      std::cout << "pos 0 0 " << z << std::endl;
      */

      //--------------------SCNT----------------------//
      shape[6] = shape[6] - d_air[0]*cm;
      shape[9] = shape[9] - d_air[0]*cm;
      shape[4] = shape[4] + (d_air[1]/2.)*cm;
      shape[7] = shape[7] - (d_air[1]/2.)*cm;
      
      G4VSolid* SCNT_solid = new G4Polyhedra("SCNT",shape[0],shape[1],shape[2],shape[3],zPlane,rInner,rOuter);
      //TMED->14
      material = materials->BGO;
	  //material = materials->GetMatByIndex(n_detmate);	  
	  
      G4LogicalVolume* SCNT_log = new G4LogicalVolume(SCNT_solid,material,"SCNT");
      z = d_air[1]/2.;
      G4ThreeVector SCNT_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* SCNT_phys = new G4PVPlacement(0,SCNT_pos,SCNT_log,"SCNT",MGOR_log,false,0,checkOverlaps);

      /*
      std::cout << "SCNT " << std::endl;
      std::cout << "shape[0] " << shape[0]/deg << std::endl;
      std::cout << "shape[1] " << shape[1]/deg << std::endl;
      std::cout << "shape[2] " << shape[2] << std::endl;
      std::cout << "shape[3] " << shape[3] << std::endl;
      std::cout << "shape[4] " << shape[4] << std::endl;
      std::cout << "shape[5] " << shape[5] << std::endl;
      std::cout << "shape[6] " << shape[6] << std::endl;
      std::cout << "shape[7] " << shape[7] << std::endl;
      std::cout << "shape[8] " << shape[8] << std::endl;
      std::cout << "shape[9] " << shape[9] << std::endl;
      std::cout << "pos 0 0 " << z << std::endl;
      */
   
      //C *** Make two position rotation matrices ***
      G4RotationMatrix* irot_front = new G4RotationMatrix;
      irot_front->rotateZ(-90*deg);
      irot_front->rotateX(-90*deg);
      G4RotationMatrix* irot_back = new G4RotationMatrix;
      irot_back->rotateZ(90*deg);
      irot_back->rotateX(-90*deg);
      G4RotationMatrix* irot_back_inv = new G4RotationMatrix;
      irot_back_inv->rotateZ(-90*deg);
      irot_back_inv->rotateX(-90*deg);
      G4ThreeVector HSNG_pos(0.0*cm,0.0*cm,0.0*cm);
      
      //C.--> Fill square DETE with hexagonal housings HSNG
      
      icopy = 0;
      G4cout << G4endl << "--------- BGO No. X, Y, Z: ---------" << G4endl;
      
      k = 1;
      //C. Detector 1.
      x =   d_mtl/2. + d_air[1]/2. - pmt_length/2. - 7.2;
      y = - aprt - wall[2] - (1./2.) * hexagon_large_width-0.6;
      z =   (5./2.)*hexagon_small_width;
      
      x_fngr[icopy]  = 0.0;
      y_fngr[icopy]  =  y;
      z_fngr[icopy]  =  z;
      
      G4VPhysicalVolume* HSNG_phys;
      if(mask[icopy] == 1)
			{HSNG_phys = new G4PVPlacement(irot_front,G4ThreeVector(-x*cm,y*cm,-z*cm),HSNG_log,"HSNG",DETE_log,false,icopy,checkOverlaps);
      
				G4cout << icopy << "\t" << -x << "\t" << y << "\t" << -z << G4endl;
			}
			icopy++;
			
      k = 2;
      //C. Detector 2.
      x =   d_mtl/2. + d_air[1]/2. - pmt_length/2.;
      y = - aprt - wall[2] - (5./4.) * hexagon_large_width-0.6;
      z =   2.*hexagon_small_width;
      
      x_fngr[icopy]  = 0.0;
      y_fngr[icopy]  =  y;
      z_fngr[icopy]  =  z;
   
      if(mask[icopy] == 1)
			{HSNG_phys = new G4PVPlacement(irot_front,G4ThreeVector(-x*cm,y*cm,-z*cm),HSNG_log,"HSNG",DETE_log,false,icopy,checkOverlaps);
      G4cout << icopy << "\t" << -x << "\t" << y << "\t" << -z << G4endl;	
			}
			icopy++;
      
      //C. Detector 3.
      x =   d_mtl/2. + d_air[1]/2. - pmt_length/2. - 6.7;
      y =   aprt + wall[2] + (1./2.)* hexagon_large_width+0.6;
      z =   2.*hexagon_small_width;
           
      x_fngr[icopy]  = 0.0;
      y_fngr[icopy]  =  y;
      z_fngr[icopy]  =  z;

      if(mask[icopy] == 1)
			{HSNG_phys = new G4PVPlacement(irot_front,G4ThreeVector(-x*cm,y*cm,-z*cm),HSNG_log,"HSNG",DETE_log,false,icopy,checkOverlaps);
        
      G4cout << icopy << "\t" << -x << "\t" << y << "\t" << -z << G4endl;
			}
			icopy++;
			
      k = 3;
      //C. Detector 4.
      x =   d_mtl/2. + d_air[1]/2. - pmt_length/2.;
      y =   aprt + wall[2] + (5./4.)* hexagon_large_width+0.6;
      z =   (3./2.)*hexagon_small_width;
         
      x_fngr[icopy]  = 0.0;
      y_fngr[icopy]  =  y;
      z_fngr[icopy]  =  z;
   
      if(mask[icopy] == 1)
			{HSNG_phys = new G4PVPlacement(irot_front,G4ThreeVector(-x*cm,y*cm,-z*cm),HSNG_log,"HSNG",DETE_log,false,icopy,checkOverlaps);
        
				G4cout << icopy << "\t" << -x << "\t" << y << "\t" << -z << G4endl;
			}
			icopy++;
			
      k = 5;
      //C. Detector 5.
      y =  (9./8.) * hexagon_large_width;
      z =  hexagon_small_width/2.;
        
      x_fngr[icopy]  = 0.0;
      y_fngr[icopy]  =  y;
      z_fngr[icopy]  =  z;
   
      if(mask[icopy] == 1)
			{HSNG_phys = new G4PVPlacement(irot_front,G4ThreeVector(-x*cm,y*cm,-z*cm),HSNG_log,"HSNG",DETE_log,false,icopy,checkOverlaps);
        
      G4cout << icopy << "\t" << -x << "\t" << y << "\t" << -z << G4endl;
			}
						icopy++;
						
      k = 7;
      //C. Detector 6.
      y =  (9./8.) * hexagon_large_width;
      z = -hexagon_small_width/2.;
           
      x_fngr[icopy]  = 0.0;
      y_fngr[icopy]  =  y;
      z_fngr[icopy]  =  z;
   
      if(mask[icopy] == 1)
			{HSNG_phys = new G4PVPlacement(irot_front,G4ThreeVector(-x*cm,y*cm,-z*cm),HSNG_log,"HSNG",DETE_log,false,icopy,checkOverlaps);
        
      G4cout << icopy << "\t" << -x << "\t" << y << "\t" << -z << G4endl;
			}
						icopy++;

      k = 9;
      //C. Detector 7.
      y =   aprt + wall[2] + (5./4.)* hexagon_large_width+0.6;
      z = - (3./2.) * hexagon_small_width;
      
      x_fngr[icopy]  = 0.0;
      y_fngr[icopy]  =  y;
      z_fngr[icopy]  =  z;
   
      if(mask[icopy] == 1)
			{HSNG_phys = new G4PVPlacement(irot_front,G4ThreeVector(-x*cm,y*cm,-z*cm),HSNG_log,"HSNG",DETE_log,false,icopy,checkOverlaps);
        
      G4cout << icopy << "\t" << -x << "\t" << y << "\t" << -z << G4endl;
			}
			icopy++;
      
      k = 10;
      //C. Detector 8.
      y = - aprt - wall[2] - (5./4.) * hexagon_large_width-0.6;
      z = - 2.* hexagon_small_width;
          
      x_fngr[icopy]  = 0.0;
      y_fngr[icopy]  =  y;
      z_fngr[icopy]  =  z;
   
      if(mask[icopy] == 1)
			{HSNG_phys = new G4PVPlacement(irot_front,G4ThreeVector(-x*cm,y*cm,-z*cm),HSNG_log,"HSNG",DETE_log,false,icopy,checkOverlaps);
        
      G4cout << icopy << "\t" << -x << "\t" << y << "\t" << -z << G4endl;
			}
						icopy++;

      //C. Detector 9.
      x =   d_mtl/2. + d_air[1]/2. - pmt_length/2.;
      y =   aprt + wall[2] + (1./2.)* hexagon_large_width+0.6;
      z = - 2.*hexagon_small_width;
      
      x_fngr[icopy]  = 0.0;
      y_fngr[icopy]  =  y;
      z_fngr[icopy]  =  z;
   
      if(mask[icopy] == 1)
			{HSNG_phys = new G4PVPlacement(irot_front,G4ThreeVector(-x*cm,y*cm,-z*cm),HSNG_log,"HSNG",DETE_log,false,icopy,checkOverlaps);
        
      G4cout << icopy << "\t" << -x << "\t" << y << "\t" << -z << G4endl;
			}
						icopy++;
						
      k = 11;
      //C. Detector 10.
      x =   d_mtl/2. + d_air[1]/2. - pmt_length/2.;
      y = - aprt - wall[2] - (1./2.) * hexagon_large_width-0.6;
      z = - (5./2.)*hexagon_small_width;
      
      
      x_fngr[icopy]  = 0.0;
      y_fngr[icopy]  =  y;
      z_fngr[icopy]  =  z;
   
      if(mask[icopy] == 1)
			{HSNG_phys = new G4PVPlacement(irot_front,G4ThreeVector(-x*cm,y*cm,-z*cm),HSNG_log,"HSNG",DETE_log,false,icopy,checkOverlaps);
        
      G4cout << icopy << "\t" << -x << "\t" << y << "\t" << -z << G4endl;
			}
			icopy++;
     
      //C. Detectors 11 through 30.
      for (k = 9; k >= 3; k -= 1) 
          {
		   z = -3./2.*hexagon_small_width + (k-3)*hexagon_small_width/2.;
		   
		   jm = 2;
		   if(k == 3 || k == 5 || k == 7 || k == 9)
		     jm = 1;
		 
		   for (j = 1; j <= jm; ++j) 
		       { 
		        x = box_width/2. + (depth + pmt_length)/2.;
		   
		        if(jm == 1)
		           {y = - (3./8.) * hexagon_large_width;}
		        else
		            y = - (9./8.) * hexagon_large_width;
		        
		        y = y + (j-1) * (3./2.) * hexagon_large_width;
		         
		 

                x_fngr[icopy]  = x - pmt_length/2.;
                y_fngr[icopy]  = y;
                
                if(mask[icopy] == 1)
								{HSNG_phys = new G4PVPlacement(irot_back,G4ThreeVector(-x*cm,y*cm,-z*cm),HSNG_log,"HSNG",DETE_log,false,icopy,checkOverlaps);

   	            G4cout << icopy << "\t" << -x << "\t" << y << "\t" << -z << G4endl;
								}
								icopy++;
		        x = -x;
		         
	

                x_fngr[icopy]  = x + pmt_length/2.;
                y_fngr[icopy]  = y;
                z_fngr[icopy]  = z;
		         
								if(mask[icopy] == 1)
								{HSNG_phys = new G4PVPlacement(irot_back_inv,G4ThreeVector(-x*cm,y*cm,-z*cm),HSNG_log,"HSNG",DETE_log,false,icopy,checkOverlaps);
                G4cout << icopy << "\t" << -x << "\t" << y << "\t" << -z << G4endl;
								}
								icopy++;
	            }
	       } 
  
   std::cout << "RAKATARAKATA" << std::endl;
   std::cout << "x_fngr[]" << " " << "y_fngr[]" << " " << "z_fngr[]" << std::endl;
   for (int i = 0; i < max_hexagon; ++i) 
		 std::cout << i << " " << x_fngr[i] << " " << y_fngr[i] << " " << z_fngr[i] << std::endl;
   std::cout << "n_fngr[][]" << std::endl;
      for (int i = 0; i < Nn; ++i) 
	      {
           for (int j = 0; j < max_hexagon; ++j) 
		       {
				 std::cout << n_fngr[i][j] << std::endl;
		        }
           }
		   		   
      if (1)//DRAGONRunAction_det->GetIswit()[3] == 1) 
	     {
          for (int i = 0; i < max_hexagon; ++i) 
		      {std::cout << i << " " << x_fngr[i] << " " << y_fngr[i] << " " << z_fngr[i] << std::endl;} //OJO no estoy seguro a que archivo va esto (mitray.dat)
          }
      // std::cout << "5" << std::endl;
      for (int i = 0; i < Nn; ++i) 
	      {
           for (int j = 0; j < max_hexagon; ++j) 
		       {
				// std::cout << n_fngr[i][j] << std::endl;
		        }
           }
      
	  neighborhood();

		std::cout << "n_fngr[][]" << max_hexagon<<" " <<Nn<<std::endl;
      for (int i = 0; i < max_hexagon; ++i) 
	      {
           for (int j = 0; j < Nn; ++j) 
		       {
						 std::cout << n_fngr[j][i] << " fngr[i,j] "<<i<<","<<j<<std::endl;
		        }
				}
	  
      std::cout << "end ugeo_finger" << std::endl;
//C.
//C. Filling nearest neighour adding back adjacency matrix
//C.
      k = 0;
      n = 0;

      for (int i = 0; i < 30; ++i) 
	      {
           k = k + 1;
           if (k == 2) k = k - 2;

           for (int j = 0; j < 30; ++j) 
		       {
                n = n + 1;
                if (n == 2) n = n - 2;

                if (i > 9 && j > 9 && n != k) 
				   {;} 
				else 
				   {
                    if (std::sqrt(std::pow(x_fngr[i] - x_fngr[j],2) + std::pow(y_fngr[i] - y_fngr[j],2) + std::pow(z_fngr[i] - z_fngr[j],2)) < 10.0) 
					   {adjacency_matrix[i][j] = 1;} 
				    else 
					   {adjacency_matrix[i][j] = 0;}
                    }
                }
           }
	  adjacency_matrix[4][2]=1;
      adjacency_matrix[2][4]=1;
      adjacency_matrix[8][5]=1;
      adjacency_matrix[5][8]=1;
	  
	   std::cout << "adjacency_matrix[][] " << std::endl;
	    for (int i = 0; i < 30; ++i) 
		     {
			  for (int j = 0; j < 30; ++j) 
			  {std::cout << adjacency_matrix[i][j]  << std::endl;}
		      }
	  
	//  std::cout << "RAKATARAKATA" << std::endl;
      }

void DRAGONDetectorConstruction::ugeo_pmt()
     {
	  //************************************************************************
      //*                                                                      *
      //*                     Define the PMT plate and PMTs                    *
      //*                                                                      *
      //************************************************************************

      G4double z, shape[3];
      
      Materials* materials = Materials::Instance();
      G4Material* material;

	  G4cout << "ugeo_pmt" << G4endl;
	  G4VSolid* PMT_solid;
	  if(mtype_pmt == 1)
	    {
		 shape[0] =  0.0;
         shape[1] =  pmt_size;
         shape[2] =  pmt_length/2.;
         
         PMT_solid = new G4Tubs("PMT",shape[0]*cm,shape[1]*cm,shape[2]*cm,0.0*deg,360.0*deg);  
         }
      else
         {
		  shape[0] =  pmt_size;
          shape[1] =  pmt_size;
          shape[2] =  pmt_length/2.;
         
          PMT_solid = new G4Box("PMT",shape[0]*cm,shape[1]*cm,shape[2]*cm);
          }
      
	  G4VisAttributes* visAttributes = new G4VisAttributes(G4Colour(0.6, 0., 0.3, 1.)); 
      visAttributes->SetVisibility(true); 
      visAttributes->SetForceSolid(true);
      
     
      //TMED->17
	  material = materials->Glass;   
      G4LogicalVolume* PMT_log = new G4LogicalVolume(PMT_solid, material,"PMT");
      PMT_log->SetVisAttributes(visAttributes);
	  G4LogicalVolumeStore* volumeStore = G4LogicalVolumeStore::GetInstance();
      G4LogicalVolume* HSNG_log = volumeStore->GetVolume("HSNG");
      z = depth/2.;
	  G4ThreeVector PMT_pos(0.0*cm,0.0*cm,z*cm);
      G4VPhysicalVolume* PMT_phys = new G4PVPlacement(0,PMT_pos,PMT_log,"PMT",HSNG_log,false,0,checkOverlaps);	
      
     }

void DRAGONDetectorConstruction::udet()
{
//C.
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//C                                                                      C
//C     UDET is a user routine which defines certain volumes to be       C
//C     detectors and defines what information each detector collects.   C
//C     For convenience the detector type variable (IDTYPE) is used      C
//C     for some volumes:                                                C
//C                                                                      C
//C                       1 - Detector Crystal                           C
//C                       2 - PM Tube                                    C
//C                       3 - Silicon                                    C
//C                                                                      C
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//C.

//C.
//************************************************************************
//*                                                                      *
//*                          Detector Scintillator               *
//*                                                                      *
//************************************************************************
//C.
//C.   1 - Detector Crystal 
//C.
      G4SDManager* DRAGONSDManager = G4SDManager::GetSDMpointer();
      auto DRAGONSDSCNT = new DRAGONSensitiveDetector("SCNT",this);
      DRAGONSDManager->AddNewDetector(DRAGONSDSCNT);
      G4LogicalVolumeStore* volumeStore = G4LogicalVolumeStore::GetInstance();
      G4LogicalVolume* HSNG_log = volumeStore->GetVolume("HSNG");
      
      if (HSNG_log) 
         {HSNG_log->SetSensitiveDetector(DRAGONSDSCNT);} 
      else 
         {G4cerr << "Error: No logical volume named 'HSNG' found!" << G4endl;}

//C.
//************************************************************************
//*                                                                      *
//*                                 PM Tube                              *
//*                                                                      *
//************************************************************************
//C.   
//C.   2 - PM Tube
//C.     
      auto DRAGONSDPMT = new DRAGONSensitiveDetector("PMT",this);
      DRAGONSDManager->AddNewDetector(DRAGONSDPMT);
      G4LogicalVolume* PMT_log = volumeStore->GetVolume("PMT");

      if (PMT_log) 
         {PMT_log->SetSensitiveDetector(DRAGONSDPMT);} 
      else 
         {G4cerr << "Error: No logical volume named 'PMT' found!" << G4endl;}

      /*
      auto DRAGONSDCTAR = new DRAGONSensitiveDetector("CTAR");
      DRAGONSDManager->AddNewDetector(DRAGONSDCTAR);
      G4LogicalVolume* CTAR_log = volumeStore->GetVolume("PMT");
      if (CTAR_log) 
         {CTAR_log->SetSensitiveDetector(DRAGONSDCTAR);} 
      else 
         {G4cerr << "Error: No logical volume named 'CTAR' found!" << G4endl;}

      auto DRAGONSDENDV = new DRAGONSensitiveDetector("ENDV");
      DRAGONSDManager->AddNewDetector(DRAGONSDENDV);
      G4LogicalVolume* ENDV_log = volumeStore->GetVolume("PMT");
      if (ENDV_log) 
         {ENDV_log->SetSensitiveDetector(DRAGONSDENDV);} 
      else 
         {G4cerr << "Error: No logical volume named 'ENDV' found!" << G4endl;}  

      auto DRAGONSDDEAD = new DRAGONSensitiveDetector("DEAD");
      DRAGONSDManager->AddNewDetector(DRAGONSDDEAD);
      G4LogicalVolume* DEAD_log = volumeStore->GetVolume("PMT");
      if (DEAD_log) 
         {DEAD_log->SetSensitiveDetector(DRAGONSDDEAD);} 
      else 
         {G4cerr << "Error: No logical volume named 'DEAD' found!" << G4endl;}
         
      auto DRAGONSDMCP = new DRAGONSensitiveDetector("MCP");
      DRAGONSDManager->AddNewDetector(DRAGONSDMCP);
      G4LogicalVolume* MCP_log = volumeStore->GetVolume("PMT");
      if (MCP_log) 
         {MCP_log->SetSensitiveDetector(DRAGONSDMCP);} 
      else 
         {G4cerr << "Error: No logical volume named 'HSNG' found!" << G4endl;}
           
      auto DRAGONSDFCUP = new DRAGONSensitiveDetector("FCUP");
      DRAGONSDManager->AddNewDetector(DRAGONSDFCUP);
      G4LogicalVolume* FCUP_log = volumeStore->GetVolume("PMT");
      if (FCUP_log) 
         {FCUP_log->SetSensitiveDetector(DRAGONSDFCUP);} 
      else 
         {G4cerr << "Error: No logical volume named 'FCUP' found!" << G4endl;}  */
     }

}
























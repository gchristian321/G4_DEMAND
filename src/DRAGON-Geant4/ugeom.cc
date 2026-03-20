#include "G4Box.hh"                         //Geant4
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalVolumeStore.hh"

#include "DRAGONDetectorConstruction.hh"    //local 
#include "Materials.hh"    

#include <fstream>                           //std


namespace DRAGON {

void DRAGONDetectorConstruction::ugeom() 
     {
//C.
//C.    ******************************************************************
//C.    *                                                                *
//C.    *       User routine to define the geometry of the set-up        *
//C.    *                                                                *
//C.    ******************************************************************
//C.

 for (G4int i = 0; i< max_hexagon; i++)
      mask[i] = 1;
 
 if(MASK != " ")
   {
    std::ifstream maskFile(MASK);
    if (!maskFile.is_open()) 
       {
        G4cerr << " No MASKING file opened for simulation! " << MASK << G4endl;
        return;
        }
    int index, value;
    for (int i = 0; i < max_hexagon && maskFile >> index >> value; ++i) 
        {mask[i] = value;}
    maskFile.close();
    }

 G4cout << "MASK: ";
 for(int i=0;i<max_hexagon;++i){
	 G4cout<<i<<": "<<mask[i];
	 if(i<max_hexagon-1)G4cout<<", ";
	 else G4cout<<G4endl;
 }
	    
//C.======================================================================
//C.
//C.    Space 
//C.
//C.======================================================================

ugeo_space();

//C.======================================================================
//C.
//C.    Detector
//C.
//C.======================================================================

ugeo_detector();

//C.======================================================================
//C.
//C.    Target
//C.
//C.======================================================================  
//C.

//tubetype = 2;
std::cout << "YUYTRTRTR "<< tubetype << std::endl;

if(tubetype == 1) 
  ugeo_trgt_large();
else if(tubetype == 0) 
       ugeo_trgt_small();
else if(tubetype == 2) 
       ugeo_trgt_small_up();
else if(tubetype == 3) 
       ugeo_trgt_small_down();
else if(tubetype == 4) 
       ugeo_trgt_small_left();
else if(tubetype == 5)                   //OJO problema
       ugeo_trgt_small_right();
else if(tubetype == 6)  
       ugeo_trgt_small_hole();

//C.======================================================================
//C.
//C.    BSO Fingers
//C.
//C.======================================================================

ugeo_finger();

//************************************************************************
//*                                                                      *
//*                     Define the PMT plate and PMTs                    *
//*                                                                      *
//************************************************************************
        
ugeo_pmt();
}


void DRAGONDetectorConstruction::ugeo_space()
     {
//C.
//************************************************************************
//*                                                                      *
//*                       Define the mother spaces                       *
//*                                                                      *
//************************************************************************
//C.     
//C.                              VACUUM SPACE 
//                              ****************
#if 1
      G4double shape[3];
      shape[0]=1500.;     //! square box space
      shape[1]=shape[0];
      shape[2]=shape[0];
      
      G4Material* material = Materials::Instance()->Vacuum;

      auto WRLD_solid = new G4Box("WRLD",shape[0]*cm,shape[1]*cm,shape[2]*cm); 
      //TMED->1
      WRLD_log = new G4LogicalVolume(WRLD_solid,material,"WRLD"); 
      auto WRLD_phys = new G4PVPlacement(nullptr,G4ThreeVector(),WRLD_log,"WRLD",nullptr,false,0,checkOverlaps); 

#else
			
			WRLD_log = G4LogicalVolumeStore::GetInstance()->GetVolume("WRLD");
			
#endif

		 }

}

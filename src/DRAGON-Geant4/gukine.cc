#include "DRAGONPrimaryGeneratorAction.hh"  //Geant4

#include "DRAGONDetectorConstruction.hh"    //local
#include "DRAGONPhysicsList.hh"           


namespace DRAGON
{
    
void DRAGONPrimaryGeneratorAction::gukine(G4Event* anEvent, G4GeneralParticleSource* fGPSParticleGun)
{
    //************************************************************************
    //*                                                                      *
    //*             GEANT4 user routine to generate Kinematics               *
    //*                        for primary tracks                            *
    //*                                                                      *
    //************************************************************************

G4cout << "YYYYYYYY" << G4endl;
G4cout << "ikine: " << ikine << G4endl;
G4cout << "mkine: " << mkine << G4endl;
G4cout << "lkine: " << fphys->Getlkine() << G4endl;
G4cout << "tubetype: " << fDetector->GetTubeType() << G4endl;
G4cout << "YYYYYYYY" << G4endl;

    if (ikine != 0)
       {gukine_mitray(anEvent, fGPSParticleGun);}
    else if (mkine != 0)
            {gukine_gbox(anEvent, fGPSParticleGun);}
    else if (fphys->Getlkine() != 0)
            {
             if (fDetector->GetTubeType() == 0)
                {gukine_full(anEvent, fGPSParticleGun);  std::cout << "000" << std::endl;}
             else if (fDetector->GetTubeType() == 1)
                     {gukine_full(anEvent, fGPSParticleGun);  std::cout << "111" << std::endl;}
             else if (fDetector->GetTubeType() == 2)
                     {gukine_full_up(anEvent,fGPSParticleGun);  std::cout << "222" << std::endl;}
             else if (fDetector->GetTubeType() == 3)
                     {gukine_full_down(anEvent,fGPSParticleGun);   std::cout << "333" << std::endl;}
             else if (fDetector->GetTubeType() == 4)
                     {gukine_full_left(anEvent,fGPSParticleGun);   std::cout << "444" << std::endl;}
             else if (fDetector->GetTubeType() == 5)
                     {gukine_full_right(anEvent,fGPSParticleGun);    std::cout << "555" << std::endl;}
             else if (fDetector->GetTubeType() == 6)
                     {gukine_full_hole(anEvent,fGPSParticleGun);   std::cout << "666" << std::endl;}
	     }

}

}

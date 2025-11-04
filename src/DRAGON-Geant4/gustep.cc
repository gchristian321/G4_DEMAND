#include "G4Step.hh"                    //Geant4
#include "G4TouchableHandle.hh"       
#include "G4VPhysicalVolume.hh"     
#include "G4EventManager.hh"  

#include "DRAGONSteppingAction.hh"      //local
#include "DRAGONRunAction.hh"
#include "DRAGONPhysicsList.hh"

//C.
//************************************************************************
//*                                                                      *
//*     GEANT3 user routine called at the end of each tracking step      *
//*                                                                      *
//************************************************************************
//C.

namespace DRAGON {


void DRAGONSteppingAction::gustep(const G4Step* step)
{
  chname_2level = step->GetPreStepPoint()->GetTouchable()->GetVolume(nlevel-1)->GetName();
  chname_pres = step->GetPreStepPoint()->GetTouchable()->GetVolume(nlevel-2)->GetName();
  
  names = chname_pres;

  std::cout << "oooooooooooNEW STEPooooooooooo" << std::endl;
  std::cout << "ipart: " << ipart << std::endl;
  std::cout << "itrtyp: " << itrtyp << std::endl;
  std::cout << "inwvol: " << inwvol << std::endl;
  std::cout << "chcase: " << kcase << std::endl;
  std::cout << "vect " << vect[0] << ", " << vect[1] << ", " << vect[2] << ", " << vect[3] << ", " << vect[4] << ", " << vect[5] << ", " << vect[6] << std::endl;
  std::cout << "vect " << vect[0] << ", " << vect[1] << ", " << vect[2] << ", " << vect[3]*vect[6] << ", " << vect[4]*vect[6] << ", " << vect[5]*vect[6] << ", " << vect[6] << std::endl;
  std::cout << "chname_nlevel " << chname_nlevel << std::endl;
  std::cout << "chname_pres " << chname_pres << std::endl;
  std::cout << "chname_2level " << chname_2level << std::endl;

  
  if (ipart >= 80)
     {
      if (chname_2level == "DETE")
         {
          std::cout << "------gustep_trgt 1-----" << std::endl;
          gustep_trgt(step);
          std::cout << "------gustep_trgt 2-----" << std::endl;
          }
      else
         {
          {
           std::cout << "------gustep_mitray 1-----" << std::endl;
           gustep_mitray(step);
           std::cout << "------gustep_mitray 2-----" << std::endl;
           }
          }
      }
  else
     {
//    step->GetTrack()->SetTrackStatus(fStopAndKill);
      
      std::cout << "------gustep_gbox 1-----" << std::endl;
      //gustep_gbox(step);
      std::cout << "------gustep_gbox 2-----" << std::endl;
      }
	  
  }

}



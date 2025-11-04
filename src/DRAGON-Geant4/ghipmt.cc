#include "G4AnalysisManager.hh"     //Geant4
#include "G4Track.hh"

#include "DRAGONSteppingAction.hh"  //Local
#include "DRAGONRunAction.hh"
#include "DRAGONDetectorConstruction.hh"

#include "global_variables.hh"

namespace DRAGON {


void DRAGONSteppingAction::ghipmt(const G4Step* step)
{
	auto gv = GlobalVariables::GetInstance();

//C.
//C *** Description: Give actions when the optical photon enters PMT
//C.
      G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
      G4Track* track = step->GetTrack();
  
      G4double hits[5];
      
      if(idtype == 2)
        {
         hits[0] = frunAction->GetDetConst()->x_fngr[number];   //OJO nlevel-1
         hits[1] = frunAction->GetDetConst()->y_fngr[number];   //OJO nlevel-1
         hits[2] = frunAction->GetDetConst()->z_fngr[number];   //OJO nlevel-1
         hits[3] = tofg;
         hits[4] = 1.0;

         analysisManager->FillH1(gv->IDMap[81],1.*nstep+0.5);
         analysisManager->FillH1(gv->IDMap[82],sleng); 

        istop = 3;
        track->SetTrackStatus(fStopButAlive);
        }

}

}



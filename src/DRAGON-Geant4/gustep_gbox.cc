#include "G4Step.hh"                    //Geant4
#include "G4VProcess.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4AnalysisManager.hh"

#include "DRAGONSteppingAction.hh"      //local
#include "DRAGONDetectorConstruction.hh"
#include "DRAGONRunAction.hh"
#include "DRAGONEventAction.hh"
#include "Materials.hh"

#include "geant3functions.hh"   

namespace DRAGON {


void DRAGONSteppingAction::gustep_gbox(const G4Step* step)
     {
	  auto gv = GlobalVariables::GetInstance();	 
	  
	  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
		 
      G4int in_new_vol, itr;
	  
	  G4int ubuf1, ubuf4, ubuf5;
	  G4double ubuf6;
	  G4Material* ubuf2;
	  G4String ubuf3;
	  
	  G4double vert[3];
	  
	  G4double xm[3], xd[3];
      
      G4String chcase, kcase;
     
      G4Track* track = step->GetTrack();
       
      Materials* mates = Materials::Instance();
      const G4Material* numed_ = step->GetPreStepPoint()->GetMaterial();
      G4Material* n_detmate_ = mates->GetMatByIndex(n_detmate); 
    
      chname_nlevel = names;       
//C     Stop particles, other than beam and recoils, that try enter spectrometer
//C     Some charged particles (ex. positrons) in B fields confused GEANT

      if (ipart < 80 && chname_nlevel == "WRLD")
         {
          istop = 2;
          track->SetTrackStatus(fKillTrackAndSecondaries);  
          }
//C.
//C.    Because INWVOL = 1 can mean either that a new volume has been entered
//C.    or that a new track has been started, define a new variable IN_NEW_VOL
//C.    which specifically indicates a new volume.
//C.
      in_new_vol = 0;
      if (inwvol == 1)
         {
          if(name_old != names || number_old != number)
            {
             if(ntmult == ntmult_old)
                in_new_vol = 1;
             std::cout << "Entrando a nuevo volumen " << std::endl;
             std::cout << "in_new_vol " << in_new_vol << std::endl;
             }
          }
//C.
//C.    Find if any of original photon hit the sleeve of photon detector module
//C.
      if (in_new_vol == 1 && track->GetParentID() == 0)
         {
          if(chname_nlevel == "FNGR" || chname_nlevel == "SCNT")
            {
             if(track->GetDefinition()->GetParticleName() == "gamma")
               {
                if(fEventAction->n_flag == 0)
                  {frunAction->n_detector = frunAction->n_detector + 1;}
                fEventAction->n_flag = fEventAction->n_flag + 1;
                }    
             }
          }
 
//C.
//C.    Stop shower particles if they enter (not simulated) PMT volume
//C.
      if(itrtyp == 1 || itrtyp == 2)
        {
         if(chname_nlevel == "PMT")
           {
            istop = 2;
            track->SetTrackStatus(fKillTrackAndSecondaries);  
            goto _200_;
            }
        }
//C.
//C.    Deal with scintillation photons when they are elsewhere
//C.
      if(in_new_vol == 1 && itrtyp == 7)
        {
        if(chname_nlevel == "PMT")
          {
           ghipmt(step);
           goto _200_;
           }
        else if(chname_nlevel != "SCNT" && chname_nlevel != "MGOR")
               {
                istop = 3;
                track->SetTrackStatus(fStopButAlive);
                goto _200_;
                }
         }

//C.
//C *** Store hits
//C.
      std::cout << "numed: " << numed_->GetName() << std::endl;
      std::cout << "n_detmate: " << n_detmate_->GetName() << std::endl;
      if(numed_ == n_detmate_ && itrtyp != 7)ghidet();
//C.
//C *** Daughter particles that were generated in the current step
//C ***                  are put on the stack
//C.

      if(ngkine > 0)
        {
//C.
        const G4VProcess* creatorProcess = track->GetCreatorProcess();
        if (creatorProcess)G4String kcase = creatorProcess->GetProcessName();

        chcase = kcase;
		itr = 0;
        for (int i = 0 ; i < ngkine; i++) 
            {
             if (ipart == 1) 
                {
                 if (chcase == "GammaConversion") 
                     fEventAction->pair_productions = fEventAction->pair_productions + 1;
                 if (chcase == "GammaConversion" || chcase == "compt" || chcase == "phot") 
                    {
                     if (track->GetParentID() == 0) 
                        {
						 track->SetTrackStatus(fAlive);
                         
						 //From gustep_gbox.f
						 ///////////////////////////////////////////////////////////////////
						 if(nlevel > 2)
						   {
							G4String wname = track->GetVolume()->GetName();
                            G4String vname = track->GetVolume()->GetMotherLogical()->GetName();
						    if(wname == "DETE")
						      {
  						       ubuf1 = track->GetParticleDefinition()->GetPDGEncoding();
						       ubuf2 = track->GetMaterial();
							   ubuf3 = step->GetPreStepPoint()->GetTouchable()->GetVolume(nlevel)->GetName();
						       if(vname == "H")
							     {
							      ubuf4  = track->GetVolume()->GetCopyNo();
                                  ubuf5 = track->GetVolume()->GetCopyNo();
                                  ucopy(vect,xm,3);
							      gmtod(xm,xd,1,chname_nlevel);
						 		  }
							   ubuf6 = xd[3];
							   
					           if(nconv + 1 < max_conv)
							     { 
							      nconv = nconv + 1;
							      //vname[nconv - 1] = track->GetVolume()->GetName();
							      ivcopy[0][nconv - 1] = ubuf4;
                                  ivcopy[1][nconv - 1] = ubuf5;
								  vert[0] = track->GetPosition().x();
								  vert[1] = track->GetPosition().y();
								  vert[2] = track->GetPosition().z();								  
                                  //ucopy(vert, true_conv[0][nconv - 1], 2);
                                  true_conv[2][nconv - 1] = ubuf6;
                                  analysisManager->FillH1(gv->IDMap[23], 1.0*ivcopy[1][nconv - 1]+0.5);
                			      }
							   
							   int it = 1;
							   itr = it;
							   G4int nt = ngkine;
							   G4int ntbeam = step->GetTrack()->GetTrackID();
                               do 
                                 {
                                  if(itr <= max_itra) 
                                     index_track_to_gamma[itr] = index_track_to_gamma[ntbeam];
                                  it++;
                                  }
                                while(it < nt); 
								
							   }
						    }
							///////////////////////////////////////////////////////////////////
						 }
                     }
                 }
			}
         }

//C.

//C.
//C *** Scintillation photons generated?
//C.

  _200_: 
  
      ngkine = 0;
//C.
      name_old   = names;
      number_old = number;
      ntmult_old = ntmult;
//C.
      chname_old = name_old;
      }
}



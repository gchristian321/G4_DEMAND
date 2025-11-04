#include "G4AnalysisManager.hh"   //Geant4
#include "G4Event.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

#include "DRAGONEventAction.hh"   //Local
#include "DRAGONRunAction.hh"
#include "DRAGONSteppingAction.hh"
#include "DRAGONDetectorConstruction.hh"

#include "global_variables.hh"
#include "geant3functions.hh"


namespace DRAGON {


void DRAGONEventAction::guout(const G4Event* event)
{
//C.
//C.    ******************************************************************
//C.    *                                                                *
//C.    *                                                                *
//C.    *       User routine called at the end of each event             *
//C.    *                                                                *
//C.    *                                                                *
//C.    ******************************************************************
//C.


  guout_mitray(event);

  guout_gbox(event);
}



void DRAGONEventAction::guout_mitray(const G4Event* event)
{
	/*
	auto gv = GlobalVariables::GetInstance();
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
 
      G4double theta, phi, pbeam, dir0[3];

      G4int i;
      G4int nt, it, ipart;

      G4double pvert[3], vert[3], amass;
	  
	  G4int nvertx = event->GetNumberOfPrimaryVertex();
      i = 0;

       do 
	       {
            G4PrimaryVertex* vertex = event->GetPrimaryVertex(i);
			G4ThreeVector vertr = vertex->GetPosition(); 
			G4int nt = vertex->GetNumberOfParticle();
			vert[0] = vertex->GetPosition().x();
			vert[1] = vertex->GetPosition().y();
			vert[2] = vertex->GetPosition().z();
            
            it = 0;
            do 
			   {
                G4PrimaryParticle* particle = vertex->GetPrimary(it);
                ipart = particle->GetPDGcode();
				pvert[0] = particle->GetMomentum().x();
				pvert[1] = particle->GetMomentum().y();
				pvert[2] = particle->GetMomentum().z();
                 if (ipart == 61) 
				    { 
                     ucopy(pvert,dir0,3);
                     return;  
                      }
                it++;
                } while (it < nt);

             i++;
            } while (i < nvertx);
	  
     vunit(dir0,3);

     theta = 0.0; 
     if(dir0[0] != 0.0 || dir0[2] != 0.0)
         theta = 1000.0*std::atan2(dir0[0], dir0[2]);
     phi = 1000.0*std::asin(dir0[1]);

     if(fRunAction->GetStepAction()->Getistop() == 100)
       {
        nout  = nout  + 1;
        idevt = idevt + 1;

        pbeam = std::sqrt(pvert[0]*pvert[0] + pvert[1]*pvert[1] + pvert[2]*pvert[2]);
		
		G4ParticleDefinition* particle = G4ParticleTable::GetParticleTable()->FindParticle(ipart);
		amass = particle->GetPDGMass(); 

        std::cout << vert[0] << " " << theta << " " << vert[1] << " " << phi << " "
                  << pbeam/1.E3 << " " << amass << " " << 1.0 << " " << 1.0 << " "
                  << vert[2]/1.E2 << " " << 0.0 << " " << 0.0 << " " << ievent << " "
                  << label[0] << " " << label[1] << " " << label[2] << std::endl;   
        }

     if(fRunAction->GetStepAction()->Getistop()== 100 || fRunAction->GetStepAction()->Getistop() == 200)
       {
        if(jslit > 0) 
          { 
           const long* seeds = CLHEP::HepRandom::getTheSeeds();
           std::cout << " SLIT! " << ievent << " " << seeds[0] << " " << seeds[1] << std::endl;
           }
        
        analysisManager->FillH1(gv->IDMap[5], vert[0]);
        analysisManager->FillH1(gv->IDMap[6], vert[1]);
        analysisManager->FillH1(gv->IDMap[7], theta);
        analysisManager->FillH1(gv->IDMap[8], phi);

        pbeam = std::sqrt(pvert[0]*pvert[0] + pvert[1]*pvert[1] + pvert[2]*pvert[2]);
        pbeam = 100.0*(pbeam/0.25855443-1.0);

        analysisManager->FillH1(gv->IDMap[10], pbeam);

        analysisManager->FillH2(gv->IDMap[103], vert[0], vert[1]);
        analysisManager->FillH2(gv->IDMap[104], theta, phi);
        analysisManager->FillH2(gv->IDMap[107], vert[0], theta);
        analysisManager->FillH2(gv->IDMap[108], vert[1], phi);
        }
	*/	
}


void DRAGONEventAction::guout_gbox(const G4Event* event)
{
	/*
	auto gv = GlobalVariables::GetInstance();
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

      G4int i, ntbeam, nubuf;
      G4int nt, it, ipart;
 
      G4double vert[3], pvert[3], tofg, ubuf[6];

      G4double egamma, gamma[3];
    
      G4double dist, theta;
//C.      
//C.                Information on the photons
//C.                --------------------------
//C.
//C.    true_e: The true energy of photon [MeV]
//C.    true_d: The true polar angle of photon [deg]
//C.
//C.            Information on the photon conversion point (i)
//C.            ----------------------------------------------
//C.
      G4int n_gamma;
//C.
//C.    Index array relating in the event ITRA -> n_gamma
//C.    -------------------------------------------------
//C.
//C.    index_track_to_gamma: ITRA -> n_gamma
//C.
      G4int imax, jmax, jord;

      n_gamma = 0;

      vzero(index_track_to_gamma,max_itra);
//C.
//C.    vname(i)   : The full name of the conversion volume 
//C.    ivcopy(2,i): The volume copy number of the conversion volume
//C.    true_conv(3,i): The coordinates of the conversion piont [cm]
//C.
      egamma = 0.0;

      vfill(true_e,max_photon,-999.0);  
      vfill(&true_d[0][0],3*max_photon,-999.0);

       G4int nvertx = event->GetNumberOfPrimaryVertex();
       i = 0;

       do 
	       {
            G4PrimaryVertex* vertex = event->GetPrimaryVertex(i);
            G4int nt = vertex->GetNumberOfParticle();
            it = 0;
			itr = it;

            do 
			   {
                G4PrimaryParticle* particle = vertex->GetPrimary(it);
                ipart = particle->GetPDGcode();
				if (ipart == 1) 
				    { 
				     pvert[0] = particle->GetMomentum().x();
				     pvert[1] = particle->GetMomentum().y();
				     pvert[2] = particle->GetMomentum().z();
                     ucopy(pvert,gamma,3);
					 egamma = Vmod(gamma,3);
					 vunit(gamma,3);
					 if(n_gamma + 1 <= max_photon)
					   {
                        n_gamma = n_gamma + 1;
                        true_e[n_gamma - 1] = 1000. * egamma;
                        true_d[2][n_gamma - 1] = raddeg * std::acos(gamma[2]);
                        if(itr <= max_itra)index_track_to_gamma[itr] = n_gamma;
					    analysisManager->FillH1(gv->IDMap[21], true_e[n_gamma]);
                        analysisManager->FillH1(gv->IDMap[22], true_d[2][n_gamma]);
                        analysisManager->FillH2(gv->IDMap[24], tofg*1.e9, z_react);
					    std::cout << tofg*1.e9 << ", " << z_react << std::endl;   
                        ucopy(gamma, true_d, 3, n_gamma);
					    }
					  return;
                      }
                it++;
				itr = it;
                } while (it < nt);

		if (track->GetDefinition() == G4Gamma::GammaDefinition()) 
		   {
            const G4VProcess* creatorProcess = track->GetCreatorProcess();
            if (creatorProcess) 
			   {
                G4String processName = creatorProcess->GetProcessName();
                
				if (processName == "InternalConversion") 
				   {
                    const G4StepPoint* preStepPoint = step->GetPreStepPoint();
                    const G4TouchableHandle& touchable = preStepPoint->GetTouchableHandle();
                    const G4VPhysicalVolume* volume = touchable->GetVolume();

                    G4String volumeName = volume->GetName();
                    G4int copyNumber = touchable->GetCopyNumber();

                    G4cout << ">>> Internal conversion photon detected in volume: "
                           << volumeName << " (copy " << copyNumber << ")"
                          << " at position " << preStepPoint->GetPosition()/mm << " mm"
                          << G4endl;

                    // Si quieres almacenar esta info globalmente:
                    // GlobalVariables::GetInstance()->AddICPhoton(volumeName, copyNumber);
                    }
                }
            }

        if(nubuf == 6 && int(ubuf[0]) == 1)
          {  
           if(nconv + 1 <= max_conv)
             {
              nconv++;
              uhtoc(iq[jvolum+int(ubuf[2])], 4, vname[nconv], 4);
              ivcopy[0][nconv] = int(ubuf[3]);
              ivcopy[1][nconv] = int(ubuf[4]);
              ucopy(vert, true_conv[0][nconv], 2);
              true_conv[2][nconv] = ubuf[5];
              analysisManager->FillH1(gv->IDMap[23], 1.0*ivcopy[1][nconv]+0.5);
              }

           int it = 1;
           do 
            {
             itr = q[jv+7+it];
             if(itr <= max_itra) 
               index_track_to_gamma[itr] = index_track_to_gamma[ntbeam];
             it++;
             }
           while(it < nt);
           }

       i++;
    } while(i < nvertx);

//C.
//C      if(ntot == 0)goto _999_;
//C      if(edetect <= tot_thrshld)goto _999_;
//C.

 z_conv = -999.0;
 it = 0;
 do 
  {
   if(z_conv == -999.0)
     {
      x_conv = true_conv[0][it];
      y_conv = true_conv[1][it];
      z_conv = true_conv[2][it];
      }
    it++;
    } 
  while(it < nconv);

  if(x_conv != -999.0 && y_conv != -999.0 && z_conv != -999.0)
    {
     analysisManager->FillH1(gv->IDMap[46], z_conv);

     dist = std::sqrt((x_conv-x_max)*(x_conv-x_max) + (y_conv-y_max)*(y_conv-y_max));
     analysisManager->FillH1(gv->IDMap[48], dist);
     dist = std::sqrt((x_conv-x_max)*(x_conv-x_max) + (y_conv-y_max)*(y_conv-y_max) + (z_conv-z_max)*(z_conv-z_max));
     analysisManager->FillH1(gv->IDMap[49], dist);
     analysisManager->FillH2(gv->IDMap[121], x_conv, y_conv);
     }

  analysisManager->FillH1(gv->IDMap[47], z_max);
  dist = std::sqrt((x_mean-x_max)*(x_mean-x_max) + (y_mean-y_max)*(y_mean-y_max));
  analysisManager->FillH1(gv->IDMap[50], dist);

  analysisManager->FillH1(gv->IDMap[60], x_mean);
  analysisManager->FillH1(gv->IDMap[61], y_mean);
  analysisManager->FillH1(gv->IDMap[62], z_mean);

  analysisManager->FillH2(gv->IDMap[122], x_max, y_max);
  analysisManager->FillH2(gv->IDMap[123], x_mean, y_mean);
  analysisManager->FillH2(gv->IDMap[124], fRunAction->GetDetConst()->x_fngr[ifngr_max], fRunAction->GetDetConst()->y_fngr[ifngr_max]);
  analysisManager->FillH2(gv->IDMap[125], fRunAction->GetDetConst()->z_fngr[ifngr_max], fRunAction->GetDetConst()->x_fngr[ifngr_max]);
  analysisManager->FillH2(gv->IDMap[126], fRunAction->GetDetConst()->z_fngr[ifngr_max], fRunAction->GetDetConst()->y_fngr[ifngr_max]);

  it = 0;
  do
   {
    if(nclu >= it)
      {
       if(it <= 3)
         {
          analysisManager->FillH1(gv->IDMap[93+it], true_e[it]-eclu[it]);
          theta = raddeg * std::acos(DRAGON::vdotn(true_d, dir_clu, it));
	      analysisManager->FillH1(gv->IDMap[96+it], theta);
          }
       }
    it++;
    } while(it < n_gamma);

  _999_:;
  */
}



}



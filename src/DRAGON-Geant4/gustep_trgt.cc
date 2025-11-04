#include "G4Step.hh"                    //Geant4
#include "G4AnalysisManager.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4Track.hh"
#include "G4EventManager.hh"
#include "G4IonTable.hh"
#include "G4VProcess.hh"

#include "DRAGONSteppingAction.hh"      //local
#include "DRAGONEventAction.hh"
#include "DRAGONRunAction.hh"
#include "DRAGONPrimaryGeneratorAction.hh"

#include "global_variables.hh"

namespace DRAGON {

void DRAGONSteppingAction::gustep_trgt(const G4Step* step)
     {
      G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	  auto gv = GlobalVariables::GetInstance();
	  
      size_t idx_pos = 0, idy_pos = 0;
      int idx = 0, idy = 0;

      G4double fstep, vectl[7] = {}, tlast; 

      G4double xtarg, ytarg, xprime, yprime, deltap;

      G4String chcase;
      G4int in_new_vol;
         
      std::cout << "Desde gustep_trgt.f" << std::endl;  
      std::cout << "itrtyp " << itrtyp << std::endl;  
	  
      if(itrtyp != 8) return;
     
      chname_nlevel = names;
      
      in_new_vol = 0;
        
      std::cout << "chname_nlevel " << chname_nlevel << std::endl; 
      std::cout << "name_old " << name_old << std::endl; 
      std::cout << "number " << number << std::endl; 
      std::cout << "number_old " << number_old << std::endl; 
      std::cout << "ntmult " << ntmult << std::endl; 
      std::cout << "ntmult_old " << ntmult_old << std::endl; 
      std::cout << "in_new_vol " << in_new_vol << std::endl; 
      std::cout << "inwvol " << inwvol << std::endl;
                       
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
	
//C Increment counters 
        std::cout << "ipart " << ipart << std::endl;
        std::cout << "irecoil " << irecoil << std::endl;
        std::cout << "ngascell " << fEventAction->GetRunAction()->Getngascell() << std::endl;
        if(chname_nlevel == "EAPG" && ipart == 80) 
//cc mt        if(chname_nlevel == "EAPG" && ipart == irecoil) 
           fEventAction->GetRunAction()->ngascell = frunAction->ngascell + 1;
        std::cout << "ngascell " << fEventAction->GetRunAction()->ngascell << std::endl;
//C.
//C
//C *** Change beam to recoil charge state
//C

      tlast = 0.;
      std::cout << "alpha " << alpha << std::endl;
      if (alpha);
      else if(chname_nlevel == "EX2G" && ipart == 80)   
             {
              std::cout << "AAAA1" << std::endl;
              G4Track* beamTrack = step->GetTrack();
              const G4ParticleDefinition* recoilDef = G4ParticleTable::GetParticleTable()->FindParticle(irecoil);
              auto beamDyn = const_cast<G4DynamicParticle*>(beamTrack->GetDynamicParticle());
              beamDyn->SetCharge(recoilDef->GetPDGCharge());
              }
       
//C.-------------------Check for resonance crossing-----------------------
//C.
//C. Check that it's in gas volume or in solid target
//C      print *,'<<<', lkine,ipart,eres,gekin,destep
//C     JS adds CMBG condition for solid target, as no CELG
      idx_pos = chname_nlevel.find('P');
      idx = (idx_pos != G4String::npos) ? static_cast<int>(idx_pos) + 1 : 0;

      std::cout << "lkine " << lkine << std::endl;
      std::cout << "targtype " << targtype << std::endl;
      std::cout << "index(chname_nlevel,P) " << idx << std::endl;

      if((lkine > 0) && (ipart == 80) &&    
      ((targtype == 0 && idx == 1)                        //!pumping tubes
      || (targtype == 0 && chname_nlevel == "UHOL")       //! listing all other volumes containing gas
      || (targtype == 0 && chname_nlevel == "CMBG")
      || (targtype == 0 && chname_nlevel == "EAPG")
      || (targtype == 0 && chname_nlevel == "CELG")
      || (targtype == 0 && chname_nlevel == "XAPG")
      || (targtype == 0 && chname_nlevel == "DHOL")
      || (targtype == 1 && chname_nlevel != "CMBG")))
     {
      std::cout << "Dentro de volumenes con gas" << std::endl;
     
      eres = frunAction->GetPrimGener()->Geteres();
     
      std::cout << "gekin " << gekin << std::endl;
      std::cout << "eres " << eres << std::endl;
      std::cout << "vect " << vect[0] << ", " << vect[1] << ", " << vect[2] << ", " << vect[3] << ", " << vect[4] << ", " << vect[5] << ", " << vect[6] << std::endl;
      std::cout << "vectl " << vectl[0] << ", " << vectl[1] << ", " << vectl[2] << ", " << vectl[3] << ", " << vectl[4] << ", " << vectl[5] << ", " << vectl[6] << std::endl;
      std::cout << "destep " << destep << std::endl;

      if(gekin > eres )
          {
           std::cout << "AAAA2" << std::endl;
           for (int i = 0; i < 7; i++) 
              vectl[i] = vect[i];
           std::cout << "vect " << vect[0] << ", " << vect[1] << ", " << vect[2] << ", " << vect[3] << ", " << vect[4] << ", " << vect[5] << ", " << vect[6] << std::endl;
        std::cout << "vectl " << vectl[0] << ", " << vectl[1] << ", " << vectl[2] << ", " << vectl[3] << ", " << vectl[4] << ", " << vectl[5] << ", " << vectl[6] << std::endl;
           }
//C        else if (gekin+destep > eres) 
        else 
//C           goto _888_
            {
             std::cout << "AAAA3" << std::endl;
             std::cout << "destep " << destep << std::endl;
             if (destep > 0.0) 
                {
                 std::cout << "AAAA4" << std::endl;
                 
                 fstep = (eres-gekin)/destep;
                              
                 std::cout << "fstep " << fstep << std::endl;
                 int i = 0;
                 do {
                     vect[i] = vect[i] - fstep * (vect[i] - vectl[i]);
                     i++;
                     } while (i < 7);
                 std::cout << "vect " << vect[0] << ", " << vect[1] << ", " << vect[2] << ", " << vect[3] << ", " << vect[4] << ", " << vect[5] << ", " << vect[6] << std::endl;
                 std::cout << "vectl " << vectl[0] << ", " << vectl[1] << ", " << vectl[2] << ", " << vectl[3] << ", " << vectl[4] << ", " << vectl[5] << ", " << vectl[6] << std::endl;
                }  

//C.
//C.
//C.--> Set E_int ntuple
          std::cout << "AAAA5" << std::endl;
          E_int = eres*1000.;
          std::cout << "E_int " << E_int << std::endl;

//C.--> Reaction is occurring, so set the value of beamtof to TOFG
          beamtof = tofg;  
          std::cout << "beamtof " << beamtof << std::endl;
          std::cout << "tofg " << tofg << std::endl;       
//C.
//C. ***    change mass of resonating particle                
//C.
          erescm = frunAction->GetPrimGener()->Geterescm();
          newm = resmass + (erescm-resenerg)/1000.;
//C          std::cout << "<<<<, " << erescm << ", " << resenerg << ", " << resmass << ", " <<  
//C              resmass+(erescm-resenerg)/1000. << std::endl;
             
             std::cout << "Eres (CM) and newmass " << erescm << ", " << newm << std::endl;  
             std::cout << "resmass " << resmass << std::endl;
             std::cout << "resenerg " << resenerg << std::endl;
             std::cout << "newm " << newm << std::endl;                          
   
           std::cout << "Calling gureact" << std::endl;
          std::cout << "vect " << vect[0] << ", " << vect[1] << ", " << vect[2] << ", " << vect[3] << ", " << vect[4] << ", " << vect[5] << ", " << vect[6] << std::endl;
          std::cout << "vect " << vect[0] << ", " << vect[1] << ", " << vect[2] << ", " << vect[3]*vect[6] << ", " << vect[4]*vect[6] << ", " << vect[5]*vect[6] << ", " << vect[6] << std::endl;
          gureact(step);
          std::cout << "vect " << vect[0] << ", " << vect[1] << ", " << vect[2] << ", " << vect[3] << ", " << vect[4] << ", " << vect[5] << ", " << vect[6] << std::endl;
          std::cout << "vect " << vect[0] << ", " << vect[1] << ", " << vect[2] << ", " << vect[3]*vect[6] << ", " << vect[4]*vect[6] << ", " << vect[5]*vect[6] << ", " << vect[6] << std::endl;
          ipart = 80;

          std::cout << "nreact " <<  frunAction->nreact << std::endl;
          react = 1;
          frunAction->nreact =  frunAction->nreact +1;

//C.          if (vect[2] < -5 || vect[2] > 5) 
//C.         sytd::cout << "Suspicious! "  << ", " << vect << ", " << gekin << ", " << fstep << ", " << vectl << std::endl;
//C.
  
          std::cout << "nreact " << frunAction->nreact << std::endl;
          std::cout << "vect " << vect[0] << ", " << vect[1] << ", " << vect[2] << ", " << vect[3] << ", " << vect[4] << ", " << vect[5] << ", " << vect[6] << std::endl;
          std::cout << "std::sqrt(vect[0]*vect[0]+vect(1)*vect[1] " << std::sqrt(vect[0]*vect[0]+vect[1]*vect[1]) << std::endl;
         
          analysisManager->FillH1(gv->IDMap[205],vect[2]);
          analysisManager->FillH1(gv->IDMap[520],vect[0]);
          analysisManager->FillH1(gv->IDMap[521],vect[1]);
          analysisManager->FillH2(gv->IDMap[213],vect[2],std::sqrt(vect[0]*vect[0]+vect[1]*vect[1]));

          xint = vect[0];
          yint = vect[1];
          zint = vect[2];
          
          std::cout << "xint " << xint << std::endl;
          std::cout << "yint " << yint << std::endl;
          std::cout << "zint " << zint << std::endl;
          
//C.          if (zint > 10 || zint < -10) std::cout << chname_nlevel << std::endl;
//C.
          std::cout << "istop " << istop << std::endl;
          istop = 1;
          step->GetTrack()->SetTrackStatus(fStopAndKill);
          std::cout << "istop " << istop << std::endl;
//C.
          goto _999_; 
       }   
//C.
    }
 
 _888_:
//C.
//C.-------------------Fill collimator histograms-----------------------
//C.
//C.--> If particle has stopped then...
//C.
//C.
//C. Pumping tube volumes all end with C

      idy_pos = chname_nlevel.find('C');
      idy = (idy_pos != G4String::npos) ? static_cast<int>(idy_pos) + 1 : 0;

	  std::cout << "Fill collimator histograms" << std::endl; 
	  std::cout << "istop " << istop << std::endl;
	  std::cout << "vect " << vect[0] << ", " << vect[1] << ", " << vect[2] << ", " << vect[3] << ", " << vect[4] << ", " << vect[5] << ", " << vect[6] << std::endl;
	  std::cout << "prodm " << prodm << std::endl;
	  std::cout << "ntargexit " << ntargexit << std::endl;
	  std::cout << "nbeamout " << nbeamout << std::endl;
	  std::cout << "tlast " << tlast << std::endl;
	  std::cout << "chname_nlevel " << chname_nlevel << std::endl;
	  std::cout << "index C " << idy << std::endl;
	  std::cout << "ipart " << ipart << std::endl;
	  std::cout << "irecoil " << irecoil << std::endl;
      
      if(idx == 4 && ipart > irecoil)
        {
         std::cout << "AAAA8" << std::endl;
         tlast=1000.*(std::sqrt(prodm*prodm+vect[6]*vect[6])-prodm);
         std::cout << "tlast " << tlast << std::endl;
         }
      
      if(istop == 2)
      {
//C.
//C.--> If its a recoil....
//C.
       std::cout << "AAAA6" << std::endl;
        if(ipart == irecoil)
        {
         std::cout << "It is a recoil" << std::endl;
         analysisManager->FillH1(gv->IDMap[201],vect[2]);
         std::cout << vect[0] << ", " << vect[1] << ", " << vect[2] << ", " << tlast << std::endl;
//C         std::cout << "recoil stopped trgt" << vect(1) << ", " << vect(2) << ", " << vect(3) << ", " << tlast << std::endl;
//C         std::cout << sleng  <<", volume," << chname_nlevel;
         analysisManager->FillH2(gv->IDMap[211],vect[2],std::sqrt(vect[0]*vect[0]+vect[1]*vect[1]));
         analysisManager->FillH2(gv->IDMap[212],vect[2],std::sqrt(vect[0]*vect[0]+vect[1]*vect[1]));
//C.
//C.--> If its in either target collimator...
//C.
//C        JS puts in CMBG condition for slid target, as no CELL
         if((targtype == 0 && chname_nlevel == "CELL")
      || (targtype == 1 && chname_nlevel == "CMBG") ) 
         {
          std::cout << "its in either target collimator " << std::endl;
          std::cout << "vect " << vect[0] << ", " << vect[1] << ", " << vect[2] << ", " << vect[3] << ", " << vect[4] << ", " << vect[5] << ", " << vect[6] << ", " << std::endl;
           if(vect[2] < 0.0)
            analysisManager->FillH1(gv->IDMap[202],std::sqrt(vect[0]*vect[0]+vect[1]*vect[1]));
          else
             analysisManager->FillH1(gv->IDMap[203],std::sqrt(vect[0]*vect[0]+vect[1]*vect[1]));
          }
       }
//C.
//C.
//C.--> If particle is beam...
//C.
        if(ipart == 80)
          {
           std::cout << "particle is beam ... " << std::endl;
//C.
//C.--> If it is either target collimator...
//C.
//C        JS puts in CMBG condition for solid target, as no CELL
         if((targtype == 0 && chname_nlevel == "CELL")
      || (targtype == 1 && chname_nlevel == "CMBG") ) 
          analysisManager->FillH1(gv->IDMap[206],vect[2]);
        }   
      istop = 1;
      step->GetTrack()->SetTrackStatus(fStopAndKill);
      
    }
    
          if(inwvol == 1 && chname_nlevel == "TEND") 
      {
      std::cout << "AAAA7" << std::endl;
      std::cout << "ntargexit " << frunAction->ntargexit << std::endl;
      std::cout << "nbeamout " << nbeamout << std::endl;
      
      if (ipart == irecoil) frunAction->ntargexit = frunAction->ntargexit +1;
      if (ipart == 80) frunAction->nbeamout = frunAction->nbeamout +1;
      

        analysisManager->FillH1(gv->IDMap[204],1.E6*tofg);
        analysisManager->FillH1(gv->IDMap[207],1.E3*gekin);
        analysisManager->FillH2(gv->IDMap[211],vect[2],std::sqrt(vect[0]*vect[0]+vect[1]*vect[1]));
        analysisManager->FillH2(gv->IDMap[214],vect[0],vect[1]);
        analysisManager->FillH2(gv->IDMap[215],vect[0],vect[3]);
        analysisManager->FillH2(gv->IDMap[216],vect[1],vect[4]);

        xtarg = vect[0]-vect[2]/vect[5]*vect[3];
        ytarg = vect[1]-vect[2]/vect[5]*vect[4];

        analysisManager->FillH2(gv->IDMap[217],xtarg,vect[3]);
        analysisManager->FillH2(gv->IDMap[218],ytarg,vect[4]);

        xprime = 1000.*atan(vect[3]/vect[5]);
        yprime = 1000.*asin(vect[4]);
        deltap = 100.*(vect[6]/std::sqrt(e0recoil*(e0recoil+2*prodm))-1.);
        
        std::cout << "xtarg " << xtarg << std::endl;
        std::cout << "ytarg " << ytarg << std::endl;
        std::cout << "xprime " << xprime << std::endl;
        std::cout << "yprime " << yprime << std::endl;
        std::cout << "deltap " << deltap << std::endl;
        std::cout << "e0recoil " << e0recoil << std::endl;
//C.
//C.  old ray output
//C        std::cout << xtarg << ", " << xprime << ", " << ytarg << ", " << yprime << ", " << dummy << ", " << dummy << ", " << deltap << std::endl;
//C.
//CCC        istop = 1;
//           step->GetTrack()->SetTrackStatus(fStopAndKill);
//C.
      }

      if(ipart == irecoil && lkine > 0 && lpart != irecoil)
       {
//C.   Initial position and momentum of recoil part
          std::cout << "Initial position and momentum of recoil part " << std::endl;
          analysisManager->FillH1(gv->IDMap[5],vect[0]);
          analysisManager->FillH1(gv->IDMap[6],vect[1]);
          analysisManager->FillH1(gv->IDMap[7],vect[3]*1000.);
          analysisManager->FillH1(gv->IDMap[8],vect[4]*1000.);
          analysisManager->FillH1(gv->IDMap[10],deltap);
      }
      lpart = ipart;
      std::cout << "lpart " << lpart << std::endl;
//C
  _999_:
//C.
//C.      
      std::cout << "ngkine " << ngkine << std::endl;
      if(ngkine > 0)
        std::cout << "AAAA8" << std::endl;

      std::cout << "chcase " << chcase << std::endl;

      std::cout << "istop " << istop << std::endl;
      if(alpha)
         ;
      else if(istop != 0)
      {
        std::cout << " Whats stopping me??? (in trgt)" << std::endl;
        std::cout << " istop: " << istop << ", Volume: " << chname_nlevel << ",  " << sleng/cm << ", " <<
                 " ipart: " << ipart << std::endl;
      }

      ngkine = 0;

      name_old   = names;
      number_old = number;
      ntmult_old = ntmult;
      
      std::cout << "name_old " << names << std::endl;
      std::cout << "number_old " << number_old << std::endl;
      std::cout << "ntmult " << ntmult << std::endl;
      }

}



#include "G4Step.hh"                    //Geant4
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"
#include "G4VProcess.hh"
#include "G4IonTable.hh"

#include "DRAGONSteppingAction.hh"      //local
#include "DRAGONDetectorConstruction.hh"
#include "DRAGONRunAction.hh"
#include "DRAGONEventAction.hh"
#include "Materials.hh"

#include "geom_dipole.hh"
#include "geom_edipol.hh"
#include "geom_mpole.hh"
#include "geom_sole.hh"

#include "geant3functions.hh"            

#include "global_variables.hh"

namespace DRAGON {


void DRAGONSteppingAction::gustep_mitray(const G4Step* step)
     {
	  auto gv = GlobalVariables::GetInstance();	 
		 
      G4int i, j, k, irot, kstop, ihit;
	  G4double radius, dr, theta, xlo, xhi, trec, hits[5];

      G4int in_new_vol;	  
      G4String chcase, kdname, kkdname;
      G4Material* chtmed;
      
      G4double xm[3], xd[3], xd_endv[3], xdd_endv[3];
      
      G4double tlast, Zrec, a0, a1, aa, b0, b1, bb, phd, tphd;
      G4double tof, x_mcp, y_mcp, z_mcp, disp, tflight;
      
      G4double clight = 29979245800.;
	  kstop = 0;
     
      G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
     
//	  const G4ParticleDefinition* recoilDef;
//      G4DynamicParticle* beamDyn;
 
//      const G4VProcess* creatorProcess = nullptr;
      const G4ParticleDefinition* beamDef = nullptr;
      
      G4Track* track = step->GetTrack();
      
//C.
//C *** Because INWVOL = 1 can mean either that a new volume has been 
//C *** entered or that a new track has been started, define a new 
//C *** variable IN_NEW_VOL which specifically indicates a new volume.
//C.

      Materials* mates = Materials::Instance();
      const G4Material* numed_ = step->GetPreStepPoint()->GetMaterial();
      G4Material* natmed = mates->GetMatByIndex(n_detmate);
	  
	  chtmed = natmed;
	  chname_nlevel = names;
      
      std::cout << "chtmed " << chtmed << std::endl;
      std::cout << "chname_nlevel " << chname_nlevel << std::endl;
      std::cout << "inwvol " << inwvol << std::endl;
      std::cout << "name_old " << name_old << std::endl;
      std::cout << "number " << number << std::endl;
      std::cout << "number_old " << number_old << std::endl;
      std::cout << "ntmult " << ntmult << std::endl;
      std::cout << "ntmult_old " << ntmult_old << std::endl;
      std::cout << "ipart " << ipart << std::endl;
      std::cout << "irecoil " << irecoil << std::endl;
      std::cout << "nend " << fEventAction->GetRunAction()->Getnend() << std::endl;
      std::cout << "nfcm2 " << fEventAction->GetRunAction()->Getnfcm2() << std::endl;
      std::cout << "Num_Recoils_Q3 " << fEventAction->GetRunAction()->Num_Recoils_Q3 << std::endl;
	  std::cout << "Num_Recoils_Q8 " << fEventAction->GetRunAction()->Num_Recoils_Q8 << std::endl;
      std::cout << "Num_BeamPart_ENDV " << fEventAction->GetRunAction()->Num_BeamPart_ENDV << std::endl;
      
      in_new_vol = 0;
      if(inwvol == 1)
        {
         if(name_old != names || number_old != number)
           {
            if(ntmult != ntmult_old)
              { 
               in_new_vol = 1;
//cc mt            if(ntmult == ntmult_old)in_new_vol = 1;
               if(chname_nlevel == "ENDV" && ipart == irecoil)
			      fEventAction->GetRunAction()->nend = fEventAction->GetRunAction()->nend + 1;
				  
//cc mt           if(chname_nlevel == "C7" && ipart == irecoil) 
               if(chname_nlevel == "C7" && ipart == 80) 
                  fEventAction->GetRunAction()->nfcm2 = fEventAction->GetRunAction()->nfcm2 + 1;  
//C  MT adds counters for the number of recoils that make it
//C     to Q3 (Sext 1) and to Q8 (Quad 6).
               if(chname_nlevel == "Q3" && ipart == irecoil)
                  fEventAction->GetRunAction()->Num_Recoils_Q3 = fEventAction->GetRunAction()->Num_Recoils_Q3 + 1; 
               if(chname_nlevel == "Q8" && ipart == irecoil)
                  fEventAction->GetRunAction()->Num_Recoils_Q8 = fEventAction->GetRunAction()->Num_Recoils_Q8 + 1; 
//C  MT adds counter for the number of beam particles that reach
//C     the end detector.
               if(chname_nlevel == "ENDV" && ipart == 80)
                  fEventAction->GetRunAction()->Num_BeamPart_ENDV = fEventAction->GetRunAction()->Num_BeamPart_ENDV + 1; 
               }
            }
        }
     
//C.      CALL uhtoc(kcase,4,chcase,4)
//C.
  //   const G4VProcess* creatorProc = track->GetCreatorProcess();
  //   if (creatorProc)
  //      {
  //       kcase = creatorProc->GetProcessName();
  //       chcase = kcase;
  //       }
 
      std::cout << "sleng " << sleng << std::endl;
      std::cout << "len_max " << len_max << std::endl;
      
      if(sleng > len_max)
        {
         istop = 6;
		 track->SetTrackStatus(fStopButAlive);
         goto _999_;
        }

//C
//C *** Change beam particle charge state to the same as the recoil
//C  
//         if(ipart == 80)
//		   {
//            recoilDef = G4ParticleTable::GetParticleTable()->FindParticle(ipart);
//            std::cout << recoilDef->GetPDGCharge() << ", " << fphys->Getfkine(1) << std::endl;
//            beamDyn = const_cast<G4DynamicParticle*>(track->GetDynamicParticle());
//            beamDyn->SetCharge(fphys->Getfkine(1));
//		    }
//C
//C *** Calculate recoil (or beam) kinetic energy
//C

      std::cout << "prodm " << prodm << std::endl;
      std::cout << "vect " << vect[0] << ", " << vect[1] << ", " << vect[2] << ", " << vect[3] << ", " << vect[4] << ", " << vect[5] << ", " << vect[6] << std::endl;
      
      if(ipart == irecoil)
        { 
         if(in_new_vol == 1)
           {
            tlast = 1000.*(std::sqrt(prodm*prodm+vect[6]*vect[6])-prodm);
            std::cout << "tlast " << tlast << std::endl;
            }
         trec = std::sqrt( prodm*prodm + vect[6]*vect[6] ) - prodm;
         std::cout << "trec " << trec << std::endl;
         }
      else if(ipart == 80)
	         {trec = std::sqrt( beammass*beammass + vect[6]*vect[6] ) - beammass;
		      std::cout << "trec " << trec << std::endl;
		      std::cout << "beammass " << beammass << std::endl;
		 }
      
      trec = trec*1000.;
      std::cout << "trec " << trec << std::endl;
   
//C
//C *** MCP hit to determine start of TAC
//C     If the current volume is 'MCP0' and we are leaving
//C     that volume, then for the recoil and the beam both
//C     calculate the time of flight using the particle's 
//C     mass and energy. Record the current coordinates.
//C     Note: trec [GeV], mass [GeV/c**2] -> units of 
//C     1/c for time of flight - so have to correct by 'clight'. 
      std::cout << "McpHit " << fEventAction->McpHit << std::endl;
      if(chname_nlevel == "MCP0" && inwvol == 2)
        {
         if(ipart == irecoil)
           {
            tof = std::sqrt( 0.5*prodm/(trec/1000.) );
            tof = tof/clight;
            std::cout << "tof " << tof << std::endl;
            }
         else if(ipart == 80)
                {
                 tof = std::sqrt( 0.5*beammass/(trec/1000.) );
                 tof = tof/clight;
                 std::cout << "tof " << tof << std::endl;
                 }
         fEventAction->McpHit = true;
         x_mcp = vect[0];
         y_mcp = vect[1];
         z_mcp = vect[2];
         std::cout << "vect " << vect[0] << ", " << vect[1] << ", " << vect[2] << ", " << vect[3] << ", " << vect[4] << ", " << vect[5] << ", " << vect[6] << std::endl;
         std::cout << "x_mcp " << x_mcp << std::endl;
         std::cout << "y_mcp " << y_mcp << std::endl;
         std::cout << "z_mcp " << z_mcp << std::endl;
         }

      std::cout << "McpHit " << fEventAction->McpHit << std::endl;   
//C *** If the current volume is 'MCP1' and MCP0 was hit, calculate distance travelled, time-of-flight
//C     between MCPs. 
      if(chname_nlevel == "MCP1" && inwvol == 2)
        {
         if(fEventAction->McpHit)
           {
            disp = std::sqrt( (vect[0]-x_mcp)*(vect[0]-x_mcp) + (vect[1]-y_mcp)*(vect[1]-y_mcp) + (vect[2]-z_mcp)*(vect[2]-z_mcp));
            std::cout << "disp " << disp << std::endl;
            if(tof != 0.0)
              {
               tflight = tof * disp;
               tflight = tflight * 1.E+09; //!! into nanoseconds
               std::cout << "tflight " << tflight << std::endl;
//C              std::cout << "ToF ," << tflight << std::endl;
               analysisManager->FillH1(gv->IDMap[518],tflight);
              }
           }
        }   

//C          
//C
//C *** If particle is escaped beam from ED1
//C
      if(ipart == 80 && chname_nlevel == "D1" && inwvol == 2) 
        {
         analysisManager->FillH1(gv->IDMap[503],1.0);
         std::cout << "Dentro de D1 " << std::endl;
         }
//C
//C *** If recoils are stopped
//C
      std::cout << "istop " << istop << std::endl;
      if(ipart == irecoil && istop != 0 && chname_nlevel != "ENDV") 
        {
         analysisManager->FillH1(gv->IDMap[504],vect[2],vect[0]);
         xstop = vect[0];
         ystop = vect[1];
         zstop = vect[2];
         std::cout << "xstop " << xstop << std::endl;
         std::cout << "ystop " << ystop << std::endl;
         std::cout << "zstop " << zstop << std::endl;
//C       std::cout << "recoil disappeared! " << vect[0] << ", " << vect[1] << ", " << vect[3] << std::endl;
          std::cout << vect[0] << ", " << vect[1] << ", " << vect[2] << ", " << tlast << std::endl;
//C       std::cout << "recoil stopped mit, " << vect[0] << ", " << vect[1] << ", " <<  vect[2] << ", " << tlast << std::endl;
//C       std::cout << sleng << ",  volume  ," << chname_nlevel << std::endl;
         }

//C     alpha acceptance tests - plot x-y coords at Q1 fringe field start
      std::cout << "alpha " << alpha << std::endl;
      if(alpha && ipart == irecoil && chname_nlevel == "Q1" && in_new_vol == 1)
        {
         analysisManager->FillH2(gv->IDMap[522],vect[0],vect[1]);
         dsssdpos = 0;
         xtest[dsssdpos] = vect[0];
         ytest[dsssdpos] = vect[1];
         etest[dsssdpos] = trec;
         std::cout << "xtest[dsssdpos] " << xtest[dsssdpos] << std::endl;
         std::cout << "ytest[dsssdpos] " << ytest[dsssdpos] << std::endl;
         std::cout << "etest[dsssdpos] " << etest[dsssdpos] << std::endl;
         std::cout << vect[0] << ", " << vect[1] << ", " << vect[2] << ", " << trec << std::endl;
         }
      if(alpha && ipart == irecoil && chname_nlevel == "TST1" && in_new_vol == 1)
        {
         dsssdpos = 1;
         xtest[dsssdpos] = (vect[0]+68.5685349)*std::cos(50.*3.1415926/180.) + (vect[2]-359.253174)*std::sin(50.*3.1415926/180.);
         ytest[dsssdpos] = vect[1];
         etest[dsssdpos] = trec;
         std::cout << "xtest[dsssdpos] " << xtest[dsssdpos] << std::endl;
         std::cout << "ytest[dsssdpos] " << ytest[dsssdpos] << std::endl;
         std::cout << "etest[dsssdpos] " << etest[dsssdpos] << std::endl;
         std::cout << (vect[0]+68.5685349)*std::cos(50.*3.1415926/180.) + (vect[2]-359.253174)*std::sin(50.*3.1415926/180.) << ", " << vect[1] << ", " << trec << std::endl;
         }
      if(alpha && ipart == irecoil && chname_nlevel == "TST2" && in_new_vol == 1)
        {
         dsssdpos = 2;
         xtest[dsssdpos] = (vect[0]-(-509.712341+8.1*std::cos(20*3.1415926/180.)))*std::cos(70.*3.1415926/180.) + (vect[2]-(661.428772-8.1*std::sin(20.*3.145926/180.)))*std::sin(70.*3.1415926/180.);
         ytest[dsssdpos] = vect[1];
         etest[dsssdpos] = trec;
         std::cout << "xtest[dsssdpos] " << xtest[dsssdpos] << std::endl;
         std::cout << "ytest[dsssdpos] " << ytest[dsssdpos] << std::endl;
         std::cout << "etest[dsssdpos] " << etest[dsssdpos] << std::endl;
         std::cout << (vect[0]-(-509.712341+8.1*std::cos(20*3.1415926/180.)))*std::cos(70.*3.1415926/180.) + (vect[2]-(661.428772-8.1*std::sin(20.*3.145926/180.)))*std::sin(70.*3.1415926/180.) << ", " << vect[1] << ", " << trec << std::endl;
         }
      if(alpha && ipart == irecoil && chname_nlevel == "TST3" && in_new_vol == 1)
        {
         dsssdpos = 3;
         xtest[dsssdpos] = (vect[0]+1024.617632);
         ytest[dsssdpos] = vect[1];
         etest[dsssdpos] = trec;
         std::cout << "xtest[dsssdpos] " << xtest[dsssdpos] << std::endl;
         std::cout << "ytest[dsssdpos] " << ytest[dsssdpos] << std::endl;
         std::cout << "etest[dsssdpos] " << etest[dsssdpos] << std::endl;
         std::cout << (vect[0]+1024.617632) << ", " << vect[1] << ", " << trec << std::endl;
         }

//C.
//C *** If particle is in ENDV volume
//C.
      if(chname_nlevel == "ENDV" && in_new_vol == 1 )
        {
//C. JS adjusts flag recoil_hit_ENDV to 1 when recoil reaches end detector
         if(ipart == irecoil) 
           {
            fEventAction->recoil_hit_ENDV = 1;
            fEventAction->recdet = 1;
            analysisManager->FillH1(gv->IDMap[510],vect[2]);
            }
    
        std::cout << "recoil_hit_ENDV " << fEventAction->recoil_hit_ENDV << std::endl;
        std::cout << "recdet " << fEventAction->recdet << std::endl;
        std::cout << "etest[dsssdpos] " << etest[dsssdpos] << std::endl;
        std::cout << "vect " << vect[0] << ", " << vect[1] << ", " << vect[2] << ", " << vect[3] << ", " << vect[4] << ", " << vect[5] << ", " << vect[6] << std::endl;
        std::cout << "xm " << xm[0] << ", " << xm[1] << ", " << xm[2] << std::endl;
        std::cout << "xd_endv " << xd_endv[0] << ", " << xd_endv[1] << ", " << xd_endv[2] << std::endl;
        std::cout << "xdd_endv " << xdd_endv[0] << ", " << xdd_endv[1] << ", " << xdd_endv[2] << std::endl;  
//C.
        ucopy(vect,xm,3);    //CALL ucopy(vect(1),xm(1),3) 
        gmtod(xm,xd_endv,1,chname_nlevel);
        std::cout << "xm " << xm[0] << ", " << xm[1] << ", " << xm[2] << std::endl;
        std::cout << "xd_endv " << xd_endv[0] << ", " << xd_endv[1] << ", " << xd_endv[2] << std::endl;
        ucopy(vect,xm,3);    //CALL ucopy(vect(4),xm(1),3)
        std::cout << "vect " << vect[0] << ", " << vect[1] << ", " << vect[2] << ", " << vect[3] << ", " << vect[4] << ", " << vect[5] << ", " << vect[6] << std::endl;
        std::cout << "xm " << xm[0] << ", " << xm[1] << ", " << xm[2] << std::endl;
        gmtod(xm,xd_endv,1,chname_nlevel);
        std::cout << "xm " << xm[0] << ", " << xm[1] << ", " << xm[2] << std::endl;
        std::cout << "xd_endv " << xd_endv[0] << ", " << xd_endv[1] << ", " << xd_endv[2] << std::endl;
//C.
        G4int ievent = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID() + 1;
        
        std::cout << "iswit(8) " << frunAction->GetIswit()[7] << std::endl;
        std::cout << "nevent " << nevent << std::endl;
        std::cout << "ievent " << ievent << std::endl;
        std::cout << "xstop " << xstop << std::endl;
        std::cout << "xstop " << xstop << std::endl;
        
        if(frunAction->GetIswit()[7] == 1)nevent = ievent;
        std::cout << "nevent " << nevent << std::endl;
//C.
        if(xdd_endv[0] != 0.0 || xdd_endv[2] != 0.0)
          {
           xdd_endv[0] = 1000.*std::atan2(xdd_endv[0],xdd_endv[2]);
           std::cout << "xdd_endv(1) " << xdd_endv[0] << std::endl;
           }
        else
          xdd_endv[0] = 0.0;
        
        xdd_endv[1] = 1000.*std::asin(xdd_endv[1]);
        std::cout << "xdd_endv(1) " << xdd_endv[0] << std::endl;
//C.
        hits[0] = vect[0];
        hits[1] = vect[1];
        hits[2] = vect[2];
        hits[2] = 0.;
        hits[4] = trec;

        std::cout << "hits " << hits[0] << ", " << hits[1] << ", " << hits[2] << ", " << hits[3] << ", " << hits[4] << ", " << hits[5] << std::endl;

//C     Calculate polar angle of recoil in lab
//C     note: vect(6) is equal to px/P, i.e is the UNIT vector component
//C     of the momentum in the z-direction. In order to calculate the 
//C     polar angle of the recoil in the plave of the end detector, we 
//C     take into account the fact that the normal vector to the plane 
//C     of the detector, in the direction of ions travelling along the 
//C     optical axis, is exactly antiparallel to the z-axis (though 
//C     obviously not collinear). Therefore, the polar angle is given 
//C     by the inverse cosine of the negative unit z-momentum:
        theta = std::acos(-vect[5]);

//C     Extract the charge (atomic number) of the incident ion
//C     in order to calculate the Pulse Height Defect
        
		beamDef = track->GetDynamicParticle()->GetParticleDefinition();
        Zrec = beamDef->GetAtomicNumber();
        std::cout << "Zrec " << Zrec << std::endl;

//C     Calculate the Pulse Height Defect
//C     
        a0 = 0.804;
        a1 = 1.13E-04;
        b0 = -0.462;
        b1 = -1.625;
        aa = a0 + a1*Zrec*Zrec;
        bb = b0 + b1/Zrec;
        phd = std::pow(10.0, bb) * std::pow(trec, aa);

        tphd = trec - phd;
        std::cout << "tphd " << tphd << std::endl;
   
//C.        CALL gsahit(iset,idet,itra,numbv,hits,ihit)
//C.
        radius = std::sqrt(xd_endv[0]*xd_endv[0]+xd_endv[1]*xd_endv[1]);
        std::cout << "radius " << radius << std::endl;
        std::cout << "xd_endv " << xd_endv[0] << ", " << xd_endv[1] << ", " << xd_endv[2] << std::endl;
//C        dr = std::sqrt(1.-vect[5]*vect[5]);
//C.
//C        theta = 0.0;
//C        if(dr != 0.0 || vect[5] != 0.0)
//C          theta = 1000.*std::atan2(dr,vect[5]);
//C        

//C
//C *** If MCP was hit, calculate distance travelled, time-of-flight
//C     between MCP and DSSSD. 
        if(fEventAction->McpHit)
          {
           disp = std::sqrt( (vect[0]-x_mcp)*(vect[0]-x_mcp) + (vect[1]-y_mcp)*(vect[1]-y_mcp) + (vect[2]-z_mcp)*(vect[2]-z_mcp) );
           std::cout << "disp " << disp << std::endl;
           std::cout << "tof " << tof << std::endl;
           if(tof != 0.0)
             {
              tflight = tof * disp;
              tflight = tflight * 1.E+09;  //!! into nanoseconds
              analysisManager->FillH1(gv->IDMap[517],tflight);
              std::cout << "tflight " << tflight << std::endl;
              }
           }
//C.
//C.      DSSSD hit-pattern
        G4int nstrip = 16;
        G4double pitch = 0.3;
         int i = 0;
       do 
         {
          double xlo = -(static_cast<double>(nstrip) / 2.0) * pitch + static_cast<double>(i - 1) * pitch;
          double xhi = -(static_cast<double>(nstrip) / 2.0) * pitch + static_cast<double>(i) * pitch;
          std::cout << "xlo " << xlo << std::endl;
          std::cout << "xhi " << xhi << std::endl;

          if (xd_endv[0] >= xlo && xd_endv[0] < xhi) 
             { std::cout << "i " << static_cast<float>(i) << std::endl; 
               analysisManager->FillH1(gv->IDMap[401],static_cast<float>(i));}
          if (xd_endv[1] >= xlo && xd_endv[1] < xhi) 
             {std::cout << "i " << static_cast<float>(i); 
              analysisManager->FillH1(gv->IDMap[402],static_cast<float>(i));} 
          i++;
          } while (i < nstrip);
//C.
        analysisManager->FillH1(gv->IDMap[11],xd_endv[0]);
        analysisManager->FillH1(gv->IDMap[12],xd_endv[1]);
        analysisManager->FillH1(gv->IDMap[13],xdd_endv[0]);
        analysisManager->FillH1(gv->IDMap[14],xdd_endv[1]);
        analysisManager->FillH1(gv->IDMap[15],trec);
        analysisManager->FillH1(gv->IDMap[516],tphd);
        analysisManager->FillH1(gv->IDMap[19],gekin*1000.);

        std::cout << "xd_endv(1) " << xd_endv[0] << std::endl;
        std::cout << "xd_endv(2) " << xd_endv[1] << std::endl;
        std::cout << "xdd_endv(1) " << xdd_endv[0] << std::endl;
        std::cout << "xdd_endv(2) " << xdd_endv[1] << std::endl;
        std::cout << "trec " << trec << std::endl;
        std::cout << "tphd " << tphd << std::endl;
        std::cout << "gekin " << gekin << std::endl;

//C.
//c MT histograms
//c
        analysisManager->FillH2(gv->IDMap[310],xd_endv[0],trec);
        analysisManager->FillH1(gv->IDMap[311],trec,xd_endv[0]);

        analysisManager->FillH1(gv->IDMap[111],xd_endv[0],xd_endv[1]);
        analysisManager->FillH1(gv->IDMap[112],xdd_endv[0],xdd_endv[1]);
        analysisManager->FillH1(gv->IDMap[113],xd_endv[0],xdd_endv[0]);
        analysisManager->FillH1(gv->IDMap[114],xd_endv[1],xdd_endv[1]);
//C.
        analysisManager->FillH1(gv->IDMap[115],radius,theta);
        std::cout << "radius " << radius << std::endl;
        std::cout << "theta " << theta << std::endl;

        std::cout << "iswit(7) " << frunAction->GetIswit()[6] << std::endl;   
        if(frunAction->GetIswit()[6] == 2)
          {
           kstop = 1;
           istop = 100;
		   track->SetTrackStatus(fStopButAlive);
           }
        else if(frunAction->GetIswit()[6] == 3)
        {
          std::cout << "idevt " << fEventAction->GetRunAction()->Getidevt() << std::endl;
		  fEventAction->GetRunAction()->idevt = fEventAction->GetRunAction()->idevt + 1;  //OJO fijate que se llama doble a RunAction y en las otras no
          std::cout << "idevt " << fEventAction->GetRunAction()->Getidevt() << std::endl;

          std::cout << xd_endv[0] << ", " << xdd_endv[0] << ", " << xd_endv[1] << ", " << xdd_endv[1] << ", " << 100.*(vect[6]/recoilmom*1000.-1.0) << ", " << xd[2] << std::endl;
          kstop = 1;
          istop = 200;
		  track->SetTrackStatus(fStopButAlive);
          }
        else if(frunAction->GetIswit()[6] == 1)
        {
//C.
//c$$$          std::cout << " x_final: " << xd_endv[0] << " y_final: " << xd_endv[1] << " theta_final: " << xdd_endv[0] << " phi_final: " << xdd_endv[1] << std::endl;
          kstop = 1;
          istop = 200;
		  track->SetTrackStatus(fStopButAlive);
        }

        goto _999_;
        }
      else if(chtmed->GetName() == "G4_Cu")
          {
           std::cout << "chtmed " << chtmed << std::endl;
//C.
//C *** If particle is in COPPER (jaws, slits)
//C.
        jslit = 1;

      if(frunAction->GetIswit()[6] == 1)
        {
         std::cout << "sleng " << sleng << std::endl;
         analysisManager->FillH1(gv->IDMap[16],sleng);
         }

      if(ipart == irecoil)
         analysisManager->FillH1(gv->IDMap[300],sleng);
      if(ipart == 80)
         analysisManager->FillH1(gv->IDMap[301],sleng);

         kstop = 1;
          istop = 5;
		  track->SetTrackStatus(fStopButAlive);
          goto _999_;
        }
 

      if(chname_nlevel == "ENDV")goto _1111_;
      if(chname_nlevel == "STRV")goto _1111_;
      if(chname_nlevel == "DEAD")goto _1111_;
//c
//c
//c  MT adds histograms
//c
//c
//c      if(ipart == irecoil && istop != 0)
//c         analysisManager->FillH1(gv->IDMap[300],sleng);
//c      if(ipart == 80 && istop != 0)
//c         analysisManager->FillH1(gv->IDMap[301],sleng);
//c
//c
//c

//C.
//C *** Check collimators in all RAYTRACE elements
//C.
  
   if (in_new_vol == 1 || inwvol == 2) 
      {
       if (in_new_vol == 1)j = 1;
       if (inwvol == 2)j = 2;
       k = number;
       kdname = names;
	   kdname = kdname.substr(0, 1);
	   std::cout << "k " << k << std::endl;
       std::cout << "kkdname " << names << std::endl;
	   std::cout << "kdname " << kdname << std::endl;
       std::cout << "vect " << vect[0] << ", " << vect[1] << ", " << vect[2] << ", " << vect[3] << ", " << vect[4] << ", " << vect[5] << ", " << vect[6] << std::endl;
       std::cout << "xm " << xm[0] << ", " << xm[1] << ", " << xm[2] << std::endl;
       ucopy(vect,xm,3); 
       std::cout << "vect " << vect[0] << ", " << vect[1] << ", " << vect[2] << ", " << vect[3] << ", " << vect[4] << ", " << vect[5] << ", " << vect[6] << std::endl;
       std::cout << "xm " << xm[0] << ", " << xm[1] << ", " << xm[2] << std::endl;
       std::cout << "xd " << xd[0] << ", " << xd[1] << ", " << xd[2] << std::endl;
       gmtod(xm,xd,1,names);
       std::cout << "xd " << xd[0] << ", " << xd[1] << ", " << xd[2] << std::endl;

       if (kdname == "D") 
          {
           irot = irot_dipole[k];
		   std::cout << "D-----" << std::endl;
           std::cout << "xd " << xd[0] << ", " << xd[1] << ", " << xd[2] << std::endl;
		   std::cout << "irot " << irot << std::endl;
		   std::cout << "dx_dipole(1,k) " << dx_dipole[0][k] << std::endl; 
           gitran(xd, &dx_dipole[0][k - 1], irot, xd);  
           std::cout << "xd " << xd[0] << ", " << xd[1] << ", " << xd[2] << std::endl;
           if (in_new_vol == 1) 
              {
               if (jcol_dipole[j - 1][k] == 1) 
                  {
                   std::cout << "jcol_dipole(j,k) " << jcol_dipole[j - 1][k] << std::endl; 
                   std::cout << "xcol_dipole(j,k) " << xcol_dipole[j - 1][k] << std::endl;
                   std::cout << "dxcol_dipole(j,k)] " << dxcol_dipole[j - 1][k] << std::endl;
                   std::cout << "ycol_dipole(j,k) " << ycol_dipole[j - 1][k] << std::endl;
				   std::cout << "dycol_dipole(j,k) " << dycol_dipole[j - 1][k] << std::endl;
                   if (std::pow(xd[0] - xcol_dipole[j - 1][k], 2) / std::pow(dxcol_dipole[j - 1][k], 2) +
                       std::pow(xd[1] - ycol_dipole[j - 1][k], 2) / std::pow(dycol_dipole[j - 1][k], 2) > 1.0) 
                      {
                       istop = 3;
                       track->SetTrackStatus(fStopButAlive);
                       }
                  } 
               else 
                  {
                   if (std::abs(xd[1] - ycol_dipole[j - 1][k]) > dycol_dipole[j - 1][k]) 
                      {
                       istop = 3;
                       track->SetTrackStatus(fStopButAlive);
                       }
                   }
               } 
           else 
              {
               if (std::abs(xd[1] - ycol_dipole[j - 1][k]) > dycol_dipole[j - 1][k]) 
                  {
                   istop = 3;
                   track->SetTrackStatus(fStopButAlive);
                   }
               }
          } else if (kdname == "Q") 
                    {
                     irot = 0;
					 std::cout << "Q-----" << std::endl;
                     std::cout << "xd " << xd[0] << ", " << xd[1] << ", " << xd[2] << std::endl;
		             std::cout << "irot " << irot << std::endl;
		             std::cout << "dx_mpole(1,k) " << dx_mpole[0][k] << std::endl; 
                     gitran(xd, &dx_mpole[0][k - 1], 0, xd);
					 std::cout << "xd " << xd[0] << ", " << xd[1] << ", " << xd[2] << std::endl;
                     if (jcol_mpole[j - 1][k] == 1) 
                        {
                         std::cout << "jcol_mpole(j,k) " << jcol_mpole[j - 1][k] << std::endl; 
                         std::cout << "xcol_mpole(j,k) " << xcol_mpole[j - 1][k] << std::endl;
                         std::cout << "dxcol_mpole(j,k) " << dxcol_mpole[j - 1][k] << std::endl;
                         std::cout << "ycol_mpole(j,k) " << ycol_mpole[j - 1][k] << std::endl;
				         std::cout << "dycol_mpole(j,k) " << dycol_mpole[j - 1][k] << std::endl;
                         if (std::pow(xd[0] - xcol_mpole[j - 1][k], 2) / std::pow(dxcol_mpole[j - 1][k], 2) +
                             std::pow(xd[1] - ycol_mpole[j - 1][k], 2) / std::pow(dycol_mpole[j - 1][k], 2) > 1.0) 
                            {
                             istop = 3;
                             track->SetTrackStatus(fStopButAlive);
                             }
                        } 
                     else 
                        {
                         if (std::abs(xd[0] - xcol_mpole[j - 1][k]) > dxcol_mpole[j - 1][k]) 
                            {
                             istop = 3;
                             track->SetTrackStatus(fStopButAlive);
                             }
                         if (std::abs(xd[1] - ycol_mpole[j - 1][k]) > dycol_mpole[j - 1][k]) 
                            {
                             istop = 3;
                             track->SetTrackStatus(fStopButAlive);
                             }
                         }
                     } else if (kdname == "S") 
                               {
                                irot = 0;
								std::cout << "S-----" << std::endl;
                                std::cout << "xd " << xd[0] << ", " << xd[1] << ", " << xd[2] << std::endl;
								std::cout << "irot " << irot << std::endl;
                                std::cout << "dx_sole(1,k) " << dx_sole[0][k] << std::endl; 
                                gitran(xd, &dx_sole[0][k], 0, xd);
                                std::cout << "xd " << xd[0] << ", " << xd[1] << ", " << xd[2] << std::endl;
                                if (jcol_sole[j - 1][k] == 1) 
                                   {
									std::cout << "jcol_sole(j,k) " << jcol_sole[j - 1][k] << std::endl; 
                                    std::cout << "xcol_sole(j,k) " << xcol_sole[j - 1][k] << std::endl;
                                    std::cout << "dxcol_sole(j,k)] " << dxcol_sole[j - 1][k] << std::endl;
                                    std::cout << "ycol_sole(j,k) " << ycol_sole[j - 1][k] << std::endl;
				                    std::cout << "dycol_sole(j,k) " << dycol_sole[j - 1][k] << std::endl;
                                    if (std::pow(xd[0] - xcol_sole[j - 1][k], 2) / std::pow(dxcol_sole[j - 1][k], 2) +
                                        std::pow(xd[1] - ycol_sole[j - 1][k], 2) / std::pow(dycol_sole[j - 1][k], 2) > 1.0) 
                                       {
                                        istop = 3;
                                        track->SetTrackStatus(fStopButAlive);
                                        }
                                    } else 
                                         {
                                          if (std::abs(xd[0] - xcol_sole[j - 1][k]) > dxcol_sole[j - 1][k]) 
                                             {
                                              istop = 3;
                                              track->SetTrackStatus(fStopButAlive);
                                              }
                                          if (std::abs(xd[1] - ycol_sole[j - 1][k]) > dycol_sole[j - 1][k]) 
                                             {
                                              istop = 3;
                                              track->SetTrackStatus(fStopButAlive);
                                              }
                                          }
                                std::cout << "xd " << xd[0] << ", " << xd[1] << ", " << xd[2] << std::endl;
                                std::cout << "xm " << xm[0] << ", " << xm[1] << ", " << xm[2] << std::endl;
                                gmtod(xm,xd,1,names);
                                std::cout << "xd " << xd[0] << ", " << xd[1] << ", " << xd[2] << std::endl;
                                std::cout << "xm " << xm[0] << ", " << xm[1] << ", " << xm[2] << std::endl;
                                }
      }
   
     if(in_new_vol == 1 || inwvol == 2)
       {
        if(kdname == "E")
          {
           std::cout << "S-----" << std::endl;
           irot = irot_edipol[k]; 
           std::cout << "xd " << xd[0] << ", " << xd[1] << ", " << xd[2] << std::endl;
		   std::cout << "irot " << irot << std::endl;
           std::cout << "dx_edipol(1,k) " << dx_edipol[0][k] << std::endl; 
           gitran(xd,&dx_edipol[1][k],irot,xd);
           std::cout << "xd " << xd[0] << ", " << xd[1] << ", " << xd[2] << std::endl;
           
           if(in_new_vol == 1)
             {
			  std::cout << "jcol_edipol(j,k) " << jcol_edipol[j - 1][k] << std::endl; 
              std::cout << "xcol_edipol(j,k) " << xcol_edipol[j - 1][k] << std::endl;
              std::cout << "dxcol_edipol(j,k) " << dxcol_edipol[j - 1][k] << std::endl;
              std::cout << "ycol_edipol(j,k) " << ycol_edipol[j - 1][k] << std::endl;
			  std::cout << "dycol_edipol(j,k) " << dycol_edipol[j - 1][k] << std::endl;
              if(jcol_edipol[j - 1][k] == 1)
                {
                 if(std::pow(xd[0]-xcol_edipol[j - 1][k],2)/std::pow(dxcol_edipol[j - 1][k],2) + std::pow(xd[1]-ycol_edipol[j - 1][k],2)/std::pow(dycol_edipol[j - 1][k],2) > 1.0)
                   istop = 3;
                   track->SetTrackStatus(fStopButAlive);
                 }
              else
                 {
                  if(std::abs(xd[0]-xcol_edipol[j - 1][k]) > dxcol_edipol[j - 1][k])
                    {
                     istop = 3;
                     track->SetTrackStatus(fStopButAlive);
                     }
                  if(std::abs(xd[1]-ycol_edipol[j - 1][k]) > dycol_edipol[j - 1][k])
                     {
                      istop = 3;
                      track->SetTrackStatus(fStopButAlive);
                      }
                  }
              }
          else
             {
              if(abs(xd[1]-ycol_edipol[j - 1][k]) > dycol_edipol[j - 1][k])
                 {
                  istop = 3;
                  track->SetTrackStatus(fStopButAlive);
                  }
              }
          }
       }
 
   
      if(istop == 3)
        {
         kstop = 1;
         std::cout << "sleng " << sleng << std::endl;
         analysisManager->FillH1(gv->IDMap[16],sleng);
      if(ipart == irecoil)
         analysisManager->FillH1(gv->IDMap[300],sleng);
      if(ipart == 80)
         analysisManager->FillH1(gv->IDMap[301],sleng);
      
         std::cout << vect[0] << ", " << vect[1] << ", " << vect[2] << ", " << tlast << std::endl;
//C         std::cout <<  "stopped mit" << vect[0] << ", " << vect[1] << ", " << vect[2] << ", " << tlast << std::endl;
//C         std::cout << sleng << ",  volume  " << chname_nlevel << std::endl;
        analysisManager->FillH1(gv->IDMap[17],sleng,std::sqrt(xd[0]*xd[0]+xd[1]*xd[1]));
        }
          
_1111_:
//C.
//C *** Daughter particles that were generated in the current step
//C ***                  are put on the stack
//C.
      std::cout << "ngkine " << ngkine << std::endl;
     
_999_: 

      if(fEventAction->jstop != 0)
      {
        istop = 1;
        step->GetTrack()->SetTrackStatus(fStopAndKill);
        kstop = 1;
        analysisManager->FillH1(gv->IDMap[16],sleng);
//c
      if(ipart == irecoil)
          analysisManager->FillH1(gv->IDMap[300],sleng);
      if(ipart == 80)
          analysisManager->FillH1(gv->IDMap[301],sleng);
//c
        std::cout << " *** Problem!!! *** " << std::endl;
      }
//C      if(kstop == 0 && istop != 0)
      if(alpha);
      else if(istop != 0 && frunAction->GetIswit()[0] == 1)
             {
              std::cout << " Whats stopping me??? (in mit)" << std::endl;
              std::cout << "istop: " << istop << ", Volume: " << chname_nlevel << ",  " << sleng << ", ipart: " << ipart << std::endl;
              }
      fEventAction->jstop = 0;
//C.
      ngkine = 0;
//C.
      name_old   = names;
      number_old = number;
      ntmult_old = ntmult;
      std::cout << "name_old " << name_old << std::endl;
      std::cout << "number_old " << number_old << std::endl;
      std::cout << "ntmult_old " << ntmult_old << std::endl;
      }	  
	  
}



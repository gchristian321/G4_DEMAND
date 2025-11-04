#include "Randomize.hh"                      //Geant4
#include "G4IonTable.hh"
#include "G4GeneralParticleSource.hh"

#include "DRAGONPrimaryGeneratorAction.hh"   //local
#include "DRAGONDetectorConstruction.hh"
#include "DRAGONHistoManager.hh"
#include "DRAGONPhysicsList.hh"


#include "global_variables.hh"

namespace DRAGON
{

void DRAGONPrimaryGeneratorAction::gukine_full(G4Event* anEvent, G4GeneralParticleSource* fGPSParticleGun)
     {
 //C======================================================================C
 //C      Initial beam and reaction information passed via card           C
 //C      FKIN - LKINE FKINE(10)                                          C
 //C                                                                      C
 //C      LKINE reaction number                                           C
 //C                                                                      C
 //C      ( 1) 13N(p,g)14O                                                C
 //C      ( 2) 15O(a,g)19Ne                                               C
 //C      ( 3) 25Al(p,g)26Si                                              C
 //C      ( 4) 17F(p,g)18Ne                                               C
 //C      ( 5) 18F(p,g)19Ne                                               C
 //C      ( 6) 19Ne(p,g)20Na                                              C
 //C      ( 7) 20Na(p,g)21Mg                                              C
 //C      ( 8) 21Na(p,g)22Mg                                              C
 //C      ( 9) 23Mg(p,g)24Al                                              C
 //C      (10) 26mAl(p,g)27                                               C
 //C      (11) 7Be(p,g)8B                                                 C
 //C                                                                      C
 //C======================================================================C
      auto gv = GlobalVariables::Instance;

      G4int ip;
      G4double cosp_r,cost_r, x_r, y_r, z_r;  
      G4double vertex[3], plab[3], plabu, beammom;

      G4double rndm[2], phi, thet, mag, r, inittheta;

      G4double beame, beamv;
      G4double beamx, beamy, beamz, beama, beamb, beamt, tdrift;
      G4double cosx, cosy, cosz;
	  
	  G4double clight = 29979245800.;

//C==
//C     Emittances from Laxtal area divided by 4 pi
//C     delx, dely are 2 sigma values for spot size
//C===
      G4double betagamma;

//C=== 
//C     Pencil beam
//C===
//C      G4double ex0 = 0., ey0 = 0., el0 = 0.;
//C      G4double betagamma0 = 0.01851;
//C      G4double delx = 0., dely = 0.;
      
      G4double  xx, yy, x, xp, y, yp;	  
//C.
//C.  Initialize all ntuple variables

     std::cout << "gukine_full.f" << std::endl;
     
//C.
//C.--> if beam, reset its charge to default
//C. 

     if(fphys->Getipart() == 80)
       {fGPSParticleGun->GetCurrentSource()->SetParticleCharge(fphys->Getfkine(0));}
   
      if(fphys->Getalpha()) 
        {
	     std::cout << "alpha: " << fphys->Getalpha() << std::endl;
		 std::cout << "YUYUYUYUYUYUYUYUYUY" << std::endl;
         G4int irecoil = 80;  	     
         //c         rndm[0] = G4RandGauss::shoot(0,1.0); 
         //c         beame = fphys->Getbeamenerg()*(1. + rndm[0]*0.01);
         Build523(fFileName);
         beame = HRNDM1(523)/1000.;
		 beammass = fphys->Getbeammass();
         beammom = std::sqrt(beame*(beame+2000.*beammass))*0.001;
         //std::cout << beammom << ", " << alpha << ", " << beame << ", " << beammass << std::endl; 
         ip = 80;             
         rndm[0] = G4RandGauss::shoot(0,1.0);
         rndm[1] = G4RandGauss::shoot(0,1.0);  
         rndm[0] = G4UniformRand();   
         rndm[1] = G4UniformRand(); 
         //!even distribution over a 2mm radius circle on axis
         r = std::sqrt(rndm[0]*(0.04));
         inittheta = rndm[1]*2.*3.141592654;
         vertex[0] = r*std::cos(inittheta); //!rndm[0]*0.0851 ! FWHM of 2mm
         vertex[1] = r*std::sin(inittheta); //!rndm[1]*0.0851
         vertex[0] = vertex[0] + offset[0];
         vertex[1] = vertex[1] + offset[1];
         vertex[2] = 0. + offset[2];
         rndm[0] = G4UniformRand();  
         rndm[1] = G4UniformRand();
         phi = 2*3.141592654*rndm[0];
         thet = std::acos(1.0-rndm[1]*0.0003125); //! between 0 and 25 mrad;
         //c         thet = 0.;
         cost_r = std::cos(thet);
         cosp_r = std::cos(phi);
         //! std::cout << phi << ", " << thet;
         //!thet = 0;
         x_r = std::sin(thet) * std::cos(phi);
         y_r = std::sin(thet) * std::sin(phi);
         z_r = std::cos(thet);
         cosx = std::sin(thet) * std::cos(phi);
         cosy = std::sin(thet) * std::sin(phi);
         cosz = std::cos(thet);
         mag = std::sqrt(cosx*cosx + cosy*cosy + cosz*cosz);
         cosx = cosx/mag;
         cosy = cosy/mag;
         cosz = cosz/mag;
      
         plab[0] = beammom*cosx;
         plab[1] = beammom*cosy;
         plab[2] = beammom*cosz;
         	
    	 std::cout << "beamenerg " << fphys->Getbeamenerg() << std::endl;
		 std::cout << "beame " << beame << std::endl;
		 std::cout << "beammass " << beammass << std::endl;
		 std::cout << "beammom " << beammom << std::endl;
    	 std::cout << "inittheta " << inittheta << std::endl;
		 std::cout << "vertex(1) " << vertex[0] << std::endl;
		 std::cout << "vertex(2) " << vertex[1] << std::endl;
		 std::cout << "vertex(3) " << vertex[2] << std::endl;
		 std::cout << "plab(1) " << plab[0] << std::endl;
		 std::cout << "plab(2) " << plab[1] << std::endl;
		 std::cout << "plab(3) " << plab[2] << std::endl;
		 std::cout << "YUYUYUYUYUYUYUYUYUY" << std::endl;
         
    	 G4ParticleDefinition* ion_ = G4ParticleTable::GetParticleTable()->FindParticle(ip);
         G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(ion_->GetAtomicNumber(), ion_->GetAtomicMass());
         fGPSParticleGun->SetParticleDefinition(ion);
         fGPSParticleGun->SetParticleCharge(fphys->Getfkine(0));
         fGPSParticleGun->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(vertex[0]*cm,vertex[1]*cm,vertex[2]*cm));                              //OJO D1
         fGPSParticleGun->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(plab[0],plab[1],plab[2]));
         fGPSParticleGun->GetCurrentSource()->GetEneDist()->SetMonoEnergy(beame*MeV);
         fGPSParticleGun->GeneratePrimaryVertex(anEvent);
		 }
      else 
        { 
//C.--> fill upright long. ellipse at the beam waist

          //rndm[0] = G4UniformRand();  
          //rndm[1] = G4UniformRand();
          
          rndm[0] = G4RandGauss::shoot(0, 1);
          rndm[1] = G4RandGauss::shoot(0, 1);
                  
          beamenerg = fphys->Getbeamenerg();
          emax = fphys->Getemax();
          beame = beamenerg + rndm[0]*emax;
// std::cout << beamenerg << ", " << emax << ", " << beame;
          buncht = fphys->Getbuncht();
          beamt = rndm[1]*buncht;
          
//C.--> ERES in MeV

          rndm[0] = G4UniformRand();  
          rndm[1] = G4UniformRand();
 
          G4cout << "FUEEEGGGGOOOOO" << G4endl;
          G4cout << "rndm[0]: " << rndm[0] << G4endl;
          G4cout << "rndm[1]: " << rndm[1] << G4endl;
          G4cout << "beamenerg: " << beamenerg << G4endl;
          G4cout << "emax: " << emax << G4endl;
          G4cout << "beame: " << beame << G4endl;
          G4cout << "buncht: " << buncht << G4endl;          
          
          if (fphys->Getlkine() != 19)
             {
             eres0 = fphys->Geteres0();
             resenerg = fphys->Getresenerg();
             reswidth = fphys->Getreswidth();
             eres = eres0*(resenerg + rndm[0]*reswidth);
             }

          G4cout << "eres0: " << eres0 << G4endl;
          G4cout << "resenerg: " << resenerg << G4endl;
          G4cout << "reswidth: " << reswidth << G4endl;
          G4cout << "eres: " << eres << G4endl;        

	  if (fphys->Getlkine() >= 19)
	     {
		  Build501();
	      do{
	         xx = HRNDM1(501);   
                 //C.    xx = 0.01; //!CR adds way to let beam never react, uncomment if needed
             eres0 = fphys->Geteres0();
             eres = eres0*xx;
             erescm = xx;                          
		 }while (eres > beame * 0.001);
              
		   G4cout << "xx: " << xx << G4endl;
		   G4cout << "erescm: " << erescm << G4endl; 
		   G4cout << "eres: " << eres << G4endl; 
	       G4cout << "FUEEEGGGGOOOOO" << G4endl;
	
	  //c      std::cout << "Erec (CM)" << ", " << erescm;
	        beammass = fphys->Getbeammass();
          //c      yy = xx**2/(2.*(resmass-beammass-xx/1000.));
          //c      std::cout << eres << ", " << beame*0.001 << ", " << xx << ", " << eres0;
		  //CALL hfill(502,xx,0.,1.0)
          //G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
          //analysisManager->FillH1(gv->IDMap[502],xx);		  
          } 

      //!     beamv = clight*std::sqrt(2.*beame/abeam/amumev);
      //!     betagamma = beamv/std::sqrt(clight**2+beamv*beamv);
      beammom = std::sqrt(beame*(beame+2000.*beammass));  //!in Mev/c        
      beamv = clight*beammom*.001/beammass;
      betagamma = beamv/clight/std::sqrt(1.0 - (beamv/clight)*(beamv/clight));  

      //C.--> fill upright trans. ellipses at the beam waist

      rndm[0] = G4RandGauss::shoot(0,1.0);
      rndm[1] = G4RandGauss::shoot(0,1.0);  
      
      amax = fphys->Getamax();
      sigx = fphys->Getsigx();
      beamx = rndm[0]*sigx;
      beama = rndm[1]*amax;
      x = beamx;
      xp = beama;

      rndm[0] = G4RandGauss::shoot(0,1.0);
      rndm[1] = G4RandGauss::shoot(0,1.0); 

      bmax = fphys->Getbmax();
      sigy = fphys->Getsigy();
      
      beamy = rndm[0]*sigy;
      beamb = rndm[1]*bmax;
      y = beamy;
      yp = beamb;

	  std::cout << "beammom " << beammom << std::endl;
	  std::cout << "beamv " << beamv << std::endl;
   	  std::cout << "betagamma " << betagamma << std::endl;
	  std::cout << "beammass " << beammass << std::endl;
	  std::cout << "sigx " << sigx << std::endl;
	  std::cout << "sigy " << sigy << std::endl;
	  std::cout << "amax " << amax << std::endl;
	  std::cout << "bmax " << bmax << std::endl;
	  std::cout << "beamx " << beamx << std::endl;
	  std::cout << "beamy " << beamy << std::endl;
	  std::cout << "beama " << beama << std::endl;
	  std::cout << "beamb " << beamb << std::endl;
	  std::cout << "offset[0] " << offset[0] << std::endl;
	  std::cout << "offset[1] " << offset[1] << std::endl;
	  std::cout << "offset[2] " << offset[2] << std::endl;
	  std::cout << "offset[3] " << offset[3] << std::endl;

      //c      std::cout << "amax & bmax: " << ", " << amax << "," << bmax

      //C.--> add positional offset from mistuned beam
      
      beamx = beamx + offset[0];
      beamy = beamy + offset[1];

      //C.--> distribute particles over a 2 sig range about beam waist at z=0.0

      rndm[0] = G4RandGauss::shoot(0,1.0);
      rndm[1] = G4RandGauss::shoot(0,1.0); 

      bunchl = fphys->Getbunchl();
      beamz = bunchl*rndm[0];

      std::cout << "bunchl " << bunchl << std::endl;

      //C.--> add beam direction axis components
     
      beama = beama + offset[2];
      beamb = beamb + offset[3];
      //c      beama = 0.
      //c      beamb = 0.
      
      //C.--> calculate direction cosines

      cosz = 1./sqrt(1+beama*beama+beamb*beamb);
      cosx = cosz*beama;
      cosy = cosz*beamb;
      //C.      std::cout << beama << ", " << beamb << ", " << cosz << ", " << cosx << ", " << cosy;

      //C.--> DEFINE PARTICLE MOMENTUM (GEV/c ) AND TYPE

      plabu = beammom*0.001;

      plab[0] = plabu*cosx;
      plab[1] = plabu*cosy;
      plab[2] = plabu*cosz;
	  
	  std::cout << "cosx " << cosx << std::endl;
	  std::cout << "cosy " << cosy << std::endl;
	  std::cout << "cosz " << cosz << std::endl;
	  std::cout << "plabu " << plabu << std::endl;
	  std::cout << "plab[0] " << plab[0] << std::endl;
	  std::cout << "plab[1] " << plab[1] << std::endl;
	  std::cout << "plab[2] " << plab[2] << std::endl;
      
      //C.--> DEFINE PARTICLE ORIGIN (VERTEX)

      //C.--> drift back half the target length to obtain entrance coords (x,y,z)
   
      tdrift = (fDetector->GetTLrms()*.999)/beamv;   //OJO
	  
	  std::cout << "beamv " << beamv << std::endl;
	  std::cout << "TLrms " << fDetector->GetTLrms() << std::endl;
	  std::cout << "tdrift " << tdrift << std::endl;

      vertex[0] = beamx - beamv*cosx*tdrift;
      vertex[1] = beamy - beamv*cosy*tdrift;
      vertex[2] = beamz - beamv*cosz*tdrift;
	  
      std::cout << "vertex[0] " << vertex[0] << std::endl;
      std::cout << "vertex[1] " << vertex[1] << std::endl;
      std::cout << "vertex[2] " << vertex[2] << std::endl;
      
      ip = 80;
                      
      tofg = beamt;      

      G4ParticleDefinition* ion_ = G4ParticleTable::GetParticleTable()->FindParticle(ip);
      G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(ion_->GetAtomicNumber(), ion_->GetAtomicMass());
      fGPSParticleGun->SetParticleDefinition(ion);
      fGPSParticleGun->SetParticleCharge(fphys->Getfkine(0));
      fGPSParticleGun->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(vertex[0]*cm,vertex[1]*cm,vertex[2]*cm));                              
      fGPSParticleGun->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(plab[0],plab[1],plab[2]));
      fGPSParticleGun->GetCurrentSource()->GetEneDist()->SetMonoEnergy(beame*MeV);
      fGPSParticleGun->GeneratePrimaryVertex(anEvent);
	
      G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
      analysisManager->FillH1(gv->IDMap[1],vertex[0]);
      analysisManager->FillH1(gv->IDMap[2],vertex[1]);  
      analysisManager->FillH1(gv->IDMap[3],cosx*1000.); 
      analysisManager->FillH1(gv->IDMap[4],cosy*1000.); 
      analysisManager->FillH1(gv->IDMap[9],plabu*1000.); 
      }    
   }
}
  /*   
      fGPSParticleGun->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(vertex[0],vertex[1],vertex[2]));
      //fGPSParticleGun->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(-3.28836578*m,-0.00754969895*m,5.75433411*m));  //OJO E1
      //fGPSParticleGun->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(-0.11,0.,0.1));                    //OJO E1
      //fGPSParticleGun->GetCurrentSource()->GetEneDist()->SetMonoEnergy(25.0775137*MeV);                                                //OJO E2
      //fGPSParticleGun->GeneratePrimaryVertex(anEvent);
      //fGPSParticleGun->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(-9.7*m,0*m,5.0*m));                             //OJO E2
      //fGPSParticleGun->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(-0.5,0.,-0.7));                    //OJO E2
      //fGPSParticleGun->GetCurrentSource()->GetEneDist()->SetMonoEnergy(8.5775137*MeV);                                                 //OJO E2
      //fGPSParticleGun->GeneratePrimaryVertex(anEvent);
      //fGPSParticleGun->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(-6.8*m,0*m,7.23*m));                            //OJO D2
      //fGPSParticleGun->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(-0.5,0.,0.17));                    //OJO D2
      //fGPSParticleGun->GetCurrentSource()->GetEneDist()->SetMonoEnergy(550.0*MeV);                                                     //OJO D2
      //fGPSParticleGun->GeneratePrimaryVertex(anEvent);
      //fGPSParticleGun->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(0.*m,0.*m,2.2*m));                              //OJO D1
      //fGPSParticleGun->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));                      //OJO D1
      //fGPSParticleGun->GetCurrentSource()->GetEneDist()->SetMonoEnergy(100.0*MeV);                                                       //OJO D1
      //fGPSParticleGun->GeneratePrimaryVertex(anEvent);
  */

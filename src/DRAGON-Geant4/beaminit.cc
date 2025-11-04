#include "G4ParticleTable.hh"         //Local
#include "G4SystemOfUnits.hh"        
#include "G4EmCalculator.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4LossTableManager.hh"
#include "G4Types.hh"
#include "G4ios.hh"
#include "G4VEmProcess.hh"             
#include "G4VEmModel.hh"             
#include "G4ParticleDefinition.hh"   
#include "G4IonTable.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Proton.hh"

#include "DRAGONPhysicsList.hh"       //Geant4

#include <cmath>



namespace DRAGON {


void DRAGONPhysicsList::beaminit() 
     {  
      std::cout << "uginit.f" << std::endl;     
      
      G4double etot;
      G4double recoilenerg,refmom;  //!In MeV
      G4double gamma,totmass,eint,excit,ereccm,toten,momm,betacm,
               gamcm,
               erec, trec, treco, eloss, e0rec;
//C
//C===
//C     Emittances from Laxtal area divided by 4 pi
//C     delx, dely are 2 sigma values for spot size
//C===
      G4double ex0 = 2.7e-4, ey0 = 2.7e-4, el0= 5.0e-12;
      G4double betagamma, betagamma0 = 0.01851;
      G4double delx = 0.25, dely = 0.25;
//C
      G4int partid,imate;
      G4Material* iMate;
      G4double dedx, beamm;
//C.
      G4double etaref, etatune;  

//C
//C If the resonance is being placed in the center of the target...
//C ex. (/DRAGON/gps/beam == 0)
      if (beamenerg == 0.)      
         { 
//C First calculate nominal beam and recoil product energies at
//C center of target        
          beammass = G4ParticleTable::GetParticleTable()->FindParticle(80)->GetPDGMass()/1000.; 
          resmass = G4ParticleTable::GetParticleTable()->FindParticle(81)->GetPDGMass()/1000.; 
          
		  std::cout << " *** " << ", " << targmass << ", " << resmass << ", " << beammass << ", " << resenerg << ", " << prodm << std::endl;
          eres0 = (prodm/targmass)*.001;  //!! non-relativistic transformation from lab to CM 
          std::cout << prodm << ", " << targmass << ", " << prodm/targmass << ", " << eres0 << std::endl;
          e0beam = (resmass-beammass-targmass)*(resmass+beammass+targmass)/(2.0*targmass);
          std::cout << e0beam << std::endl;
          etot = e0beam+beammass;
          e0recoil=(std::pow(resmass,2.)+std::pow(prodm,2.))*(etot+targmass)/(2*std::pow(resmass,2.0)) - prodm;   //!ErecoilCM * gamma = value for 90deg CM gamma
//!  This is the place any angular distribution equations would come in.
//!  Lorentz boosts need to be made for each angle in 3-d kinematics. 
//!
//!  Now corrections for energy loss in the target gas
          mtarg = (atarg == 12.) ? 62 : 27; 
          imate = mtarg; 
          iMate = Mtarg;
          partid = 80;

          if(Mtarg)
{G4cout << "SI- " << G4endl;}
else
{G4cout << "NO- " << G4endl;}


          G4cout << "imate " << imate << G4endl;
          G4cout << "targmass " << targmass << G4endl;
          G4cout << "resmass " << resmass << G4endl;
          G4cout << "beammass " << beammass << G4endl;
          G4cout << "resenerg " << resenerg << G4endl;
          G4cout << "prodm " << prodm << G4endl;
          G4cout << "eres0 " << eres0 << G4endl;
          G4cout << "etot " << etot << G4endl;
          G4cout << "e0recoil " << e0recoil << G4endl;
          G4cout << "imate " << imate << G4endl;
          G4cout << "Material: " << Mtarg->GetName() << G4endl;  
          G4cout << *Mtarg << G4endl; 
          G4cout << "e0beam " << e0beam << G4endl;    
          G4cout << "partid: " << partid << G4endl; 
          G4cout << "Particle name: " << G4ParticleTable::GetParticleTable()->FindParticle(partid)->GetParticleName() << G4endl; 
//C.
//C.
//C     Divide target into 100 slices, add up energyloss to get beam energy
//c        eloss = 0.
//c        beamenerg = e0beam;
//c        for (int i = 0; i<100; i++)
//c            {
//c             beamenerg = beamenerg + eloss
                G4EmCalculator emCal;
                G4ParticleDefinition* DRAGONion = G4ParticleTable::GetParticleTable()->FindParticle(partid);  
                G4int Z = DRAGONion->GetAtomicNumber();
                G4int A = DRAGONion->GetAtomicMass();
                G4double excit = 0.0;
                
                if (auto ionDef = dynamic_cast<G4Ions*>(DRAGONion)) 
                    excit = ionDef->GetExcitationEnergy();
              
                G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z, A, excit);
                G4double dedx = emCal.ComputeTotalDEDX(e0beam, ion, iMate);  

  //#######################################################################
  /*G4EmCalculator calc_;
  double energy_ = 8.0*MeV;
  G4Material* material_ = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ge");
  G4ParticleDefinition* ion_ = G4IonTable::GetIonTable()->GetIon(3, 7);
  G4double dedx_ = calc_.ComputeTotalDEDX(
  energy_,
  ion_,
  material_
  );

  std::cout << *(material_) << std::endl;    
  std::cout << "Particle: " << ion_->GetParticleName() << std::endl;
  std::cout << "Material: " << material_->GetName() << std::endl; 
  std::cout << "Kinetic Energy: " << energy_ << " MeV" << std::endl;
  std::cout << "dE/dx [G4EMCalculator]: " << dedx_ << " MeV/mm" << std::endl;
  std::cout << "---------------------------------" << std::endl; */
  //#######################################################################

  //#######################################################################
  G4EmCalculator calc_;
  double energy_ = 8.*MeV;
  G4Material* material_ = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ge");
  G4ParticleDefinition* ion_ = G4IonTable::GetIonTable()->GetIon(3, 7);
  G4double dedx_ = calc_.ComputeTotalDEDX(
  energy_,
  ion_,
  material_
  );

  std::cout << *(material_) << std::endl;    
  std::cout << "Particle: " << ion_->GetParticleName() << std::endl;
  std::cout << "Material: " << material_->GetName() << std::endl; 
  std::cout << "Kinetic Energy: " << energy_ << " MeV" << std::endl;
  std::cout << "dE/dx [G4EMCalculator]: " << dedx_ << " MeV/mm" << std::endl;
  std::cout << "---------------------------------" << std::endl; 
  //#######################################################################


                G4cout << "4Ion: " << ion->GetParticleName() << " Energia: " << e0beam << " Material: " << iMate->GetName() << " dE/dx: " << dedx << G4endl;

//c             eloss = (entdens/100.)*dedx * 0.001;
//c             }
//c        beamenerg = beamenerg*1000.;        
//C.==
//C.=     entdens - The gas thickness to target center
//C.=     beamenerg, dedx in MeV/cm , e0beam in GeV
//C.==
        
//c        if(beamenerg != 0.)
//c          {continue;}
//c        else {
	
dedx = 0.00679321261;	//OJO  Problema aqui

          beamenerg = e0beam*1000. + dedx * entdens/cm;
//c       }
	  	  std::cout << "entdens: " << entdens/cm << std::endl;
	  	  G4cout << "exitdens " << exitdens/cm << G4endl;  
          std::cout << "dedx: " << dedx << std::endl;
		  std::cout << "beamenerg: " << beamenerg << std::endl;
//c        
          beammom=std::sqrt(beamenerg*(beamenerg+2000.*beammass));  //!in Mev/c        
          beamvel = clight*beammom*.001/beammass;
          gamma = 1.0/sqrt(1.0 - std::pow((beamvel/clight),2.0));
          betagamma = beamvel/clight*gamma;  //!sign wrong on original?

	  	  std::cout << "beammom: " << beammom << std::endl;
          std::cout << "beamvel: " << beamvel << std::endl;
		  std::cout << "gamma: " << gamma << std::endl;
		  std::cout << "betagamma: " << betagamma << std::endl;
//C.
          std::cout << "+++++++++++++++++BEAM AND TARGET+++++++++++++++" << std::endl;
          std::cout << "Beam energy       " << beamenerg << " MeV" << ", " << "Momentum " << beammom << " MeV/c" << std::endl;
          std::cout << "Gas half thickness " << entdens/cm << " cm" << std::endl;
          std::cout << "dE/dx in target   " << dedx << " MeV/cm" << std::endl;
          std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
          beamo = beamenerg - dedx*(entdens/cm)*2.;
          std::cout << "Beam Mean Energy leaving target " << beamo <<  " MeV" << std::endl;

//C.
//C.-->   Calculate initial beam distribution parameters
//C.
          sigx = delx/2.;
          sigy = dely/2.;
//C.
//C.-->   assumed Gaussian buncht = 1 sigma
//C.
          buncht = 1.E-9;
          bunchl = buncht*beamvel;
//C.
//C.-->   scale ex and ey
//C.
          ex = ex0*betagamma0/betagamma;
          ey = ey0*betagamma0/betagamma;
          el = el0*betagamma0/betagamma;
//C.
         if(sigx == 0.)
           {
            amax = 0.;
            bmax = 0.;
            emax = 0.;
            }
         else
            {
             amax = ex/sigx;
             bmax = ey/sigy;
             emax = el/buncht;
             }

//!  Now energy loss for the recoil
         partid = irecoil;                   
//C.
//C.
         e0rec = e0recoil;
		 
         DRAGONion = G4ParticleTable::GetParticleTable()->FindParticle(partid);  
         Z = DRAGONion->GetAtomicNumber();
         A = DRAGONion->GetAtomicMass();
         excit = 0.0;
                
         if (auto ionDef = dynamic_cast<G4Ions*>(DRAGONion)) 
             excit = ionDef->GetExcitationEnergy();
              
         ion = G4IonTable::GetIonTable()->GetIon(Z, A, excit);
         dedx = emCal.ComputeTotalDEDX(e0recoil, ion, iMate);    
		 dedx = 0.0129431868;	                                 //OJO  Problema aqui
		 
         e0recoil = e0recoil -0.001*dedx*(exitdens/cm);
         recoilenerg = e0recoil*1000.;    //!in MeV
         recoilenerg = recoilenerg*(1.+energscale);
         recoilmom = std::sqrt(recoilenerg*(recoilenerg+2000.*prodm));
 
		G4cout << "irecoil " << irecoil << G4endl;         
		G4cout << "exitdens " << exitdens/cm << G4endl;       
		G4cout << "energscale " << energscale << G4endl;  
		G4cout << "refenerg " << refenerg << G4endl;   
		G4cout << "refatno " << refatno << G4endl;               
		G4cout << "amumev " << amumev << G4endl;   
		G4cout << "refq " << refq << G4endl; 
 
         std::cout << "++++++++++++++++RECOIL+++++++++++++++++++++++" << std::endl;
         std::cout << "Recoil Mean Energy from reaction " << e0rec*1000. << " MeV" << std::endl;
         std::cout << "Gas half thickness " << exitdens/cm << " cm" << std::endl;
         std::cout << "dE/dx in target   " << dedx << " MeV/cm" << std::endl;
         std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
         std::cout << "Recoil Mean Energy leaving target " << e0recoil*1000. << " MeV" << ", " << "Momentum " << recoilmom << " MeV/c" << std::endl;
         std::cout << "Recoil energy tuned to " << recoilenerg << " MeV" << std::endl;

//!
//!  Determine scaling parameters for this reaction c/w the reference tune
//!
         etaref=refenerg/(2* std::pow(refatno*amumev,2.0));
         etatune=recoilenerg/(2*std::pow(prodm*1000.,2.0));
         refmom = std::sqrt(refenerg*(refenerg+2*refatno*amumev));
         bscale = recoilmom/fkine[1]/ (refmom/refq);
         escale = recoilenerg*refq/(refenerg*fkine[1])*(1+etatune)/(1+2*etatune)*(1+2*etaref)/(1+etaref);
//C. ----- Trick to make beam particles get through if needed 
//C. ----- (rescale Electric Dipoles)
//C.       escale = escale*prodm/beammass
//C. ----- 
         std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
         std::cout << "Magnetic element scale factor " << bscale << std::endl;
         std::cout << "Electric element scale factor " << escale << std::endl;
         std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl; 
         }       
      else
         {
          if(alpha)
            {
             std::cout << "--------- alpha source ----------" << std::endl;
             beammass = G4ParticleTable::GetParticleTable()->FindParticle(80)->GetPDGMass()/1000.; 
             beammom = std::sqrt(beamenerg*(beamenerg+2000.*beammass));
              
              prodm = beammass;
              recoilenerg = beamenerg*(1.+energscale);
              recoilmom = std::sqrt(recoilenerg*(recoilenerg+2000.*prodm));
              etaref=refenerg/(2*std::pow(refatno*amumev,2.0));
              etatune=recoilenerg/(2*std::pow(prodm*1000.,2.0));
              refmom = sqrt(refenerg*(refenerg+2*refatno*amumev));
              bscale = recoilmom/fkine[1]/ (refmom/refq);
              escale = recoilenerg*refq/(refenerg*fkine[1])*(1+etatune)/(1+2*etatune)*(1+2*etaref)/(1+etaref);
              std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
              std::cout << "Magnetic element scale factor " << bscale << std::endl;
              std::cout << "Electric element scale factor " << escale << std::endl;
              std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
              std::cout << "beaminit.f" << std::endl;
              
			  std::cout << "beammass " << beammass << std::endl;
			  std::cout << "beamenerg " << beamenerg << std::endl;
			  std::cout << "prodm " << prodm << std::endl;
			  std::cout << "energscale " << energscale << std::endl;
			  std::cout << "recoilmom " << recoilmom << std::endl;
			  std::cout << "refatno " << refatno << std::endl;
			  std::cout << "recoilenerg " << recoilenerg << std::endl;
			  std::cout << "refmom " << refmom << std::endl;
			  std::cout << "refatno " << refatno << std::endl;
			  std::cout << "refenerg " << refenerg << std::endl;			  
			  std::cout << "refq " << refq << std::endl;
			  std::cout << "etatune " << etatune << std::endl;
			  std::cout << "etaref " << etaref << std::endl;
			  std::cout << "fkine(1) " << fkine[0] << std::endl;
			  std::cout << "fkine(2) " << fkine[1] << std::endl;
			  }
           else
               {
                std:: cout << "********** Using BEAM FFCARD **********" << std::endl;
                beammass = G4ParticleTable::GetParticleTable()->FindParticle(80)->GetPDGMass()/1000.; 
                resmass = G4ParticleTable::GetParticleTable()->FindParticle(81)->GetPDGMass()/1000.; 
                eres0 = prodm/targmass*0.001;
                mtarg = (atarg == 12.) ? 62 : 27;
                imate = mtarg;
		iMate = Mtarg;
                partid = 80;
                G4EmCalculator emCal;
                G4IonTable* ionTable = G4ParticleTable::GetParticleTable()->GetIonTable();
                //G4double dedx = emCal.ComputeTotalDEDX(e0beam, G4ParticleTable::GetParticleTable()->FindParticle(partid), imate);    //OJO
                beammom=sqrt(beamenerg*(beamenerg+2000.*beammass));
                beamvel = clight*beammom*.001/beammass;
                gamma = 1.0/sqrt(1.0-std::pow(beamvel/clight,2.0));
                betagamma = beamvel/clight*gamma;
                beamo = beamenerg - dedx*entdens*2. *cm;
                beamm = beamenerg - dedx*entdens *cm;
//C.
                std::cout << "+++++++++++++++BEAM AND TARGET+++++++++++++++++" << std::endl;
                std::cout << "Beam energy    " << beamenerg << " MeV" << ", " << "Momentum " << beammom << " MeV/c" << std::endl;
                std::cout << "Gas half thickness " << entdens/cm << " cm" << std::endl;
                std::cout << "dE/dx in target  "<< dedx << " MeV/cm" << std::endl;
                std::cout << "Beam energy at target exit " << beamo << " MeV" << std::endl;
                std::cout << "Beam energy at target centre " << beamm << " MeV" << std::endl;
                std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;             
//C.
//C.       Initial beam distribution parameters
//C.       
                sigx = delx/2.;
                sigy = dely/2.;
//C.
//C.-->   assumed Gaussian buncht = 1 sigma
//C.
                buncht = 1.E-9;
                bunchl = buncht*beamvel;
//C.
//C.-->   scale ex and ey
//C.
                ex = ex0*betagamma0/betagamma;
                ey = ey0*betagamma0/betagamma;
                el = el0*betagamma0/betagamma;
//C.
                if (sigx == 0.)
                   {
                    bmax = 0.;
                    emax = 0.;
                    }
                else
                   {
                    amax = ex/sigx;
                    bmax = ey/sigy;
                    emax = el/buncht;
                    }               
//C.
//C.--> Recoil energy at center of target
//C.
                totmass = std::sqrt( std::pow((beammass+targmass),2.0) + 2.*targmass*(beamm/1000.));
                eint = totmass - beammass - targmass;
                excit = (beammass+targmass-prodm) + eint;
                ereccm = sqrt( std::pow(prodm,2.0) + std::pow(excit,2.0) );
                toten = (beamm/1000.) + beammass;
                momm = sqrt( std::pow(toten,2.0) - std::pow(beammass,2.0) );
                betacm = momm/(toten+targmass);
                gamcm = (toten+targmass)/sqrt(std::pow((toten+targmass),2.0)-std::pow(momm,2.0));
                erec = gamcm*ereccm;
                trec = erec - prodm;
//c        std::cout << totmass << ", " << eint << ", " excit << ", " << ereccm << std::endl;
//c        std::cout << momm << ", " << betacm << ", " << gamcm << ", " << erec << ", " << trec << std::endl;
                partid = irecoil;
                //G4double dedx = emCal.ComputeTotalDEDX(e0beam, G4ParticleTable::GetParticleTable()->FindParticle(partid), imate);     //OJO
                treco = trec - 0.001*dedx*exitdens *cm;
                treco = treco*(1.+energscale);
                recoilmom = 1000*(std::sqrt(std::pow((treco + prodm),2.0) - std::pow(prodm,2.0)));
                std::cout << "++++++++++++++++RECOIL+++++++++++++++++++++++" << std::endl;
                std::cout << "Recoil Mean Energy at target centre (90 deg gamma)" << trec*1000. << " MeV" << std::endl;
                std::cout << "Gas half thickness" << exitdens/cm << " cm" << std::endl;
                std::cout << "dE/dx in target   " << dedx << " MeV/cm" << std::endl;
                std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
                std::cout << "Recoil Mean Energy leaving target" << treco*1000. << " MeV" << ", " << "Momentum" << recoilmom << " MeV/c" << std::endl;
                recoilenerg = treco*1000.;
//!
//!  Determine scaling parameters for this reaction c/w the reference tune
//!
                etaref=refenerg/(2*std::pow((refatno*amumev),2.0));
                etatune=recoilenerg/(2*std::pow((prodm*1000.),2.0));
                refmom = sqrt(refenerg*(refenerg+2*refatno*amumev));
                bscale = recoilmom/fkine[1]/ (refmom/refq);
                escale = recoilenerg*refq/(refenerg*fkine[1])*(1+etatune)/(1+2*etatune)*(1+2*etaref)/(1+etaref);
//C.----- Trick to make beam particles get through if needed 
//C.  ----- (rescale Electric Dipoles)
//C.        escale = escale*prodm/beammass
//C. -----                
              std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
              std::cout << "Magnetic element scale factor " << bscale << std::endl;
              std::cout << "Electric element scale factor " << escale << std::endl;
              std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;  
              }     
          }
      std::cout << "beaminit.f" << std::endl;
     } 
      
}

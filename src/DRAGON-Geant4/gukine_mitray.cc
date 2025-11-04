#include "DRAGONPrimaryGeneratorAction.hh"      //Geant4

#include "Randomize.hh"
#include "DRAGONHistoManager.hh"
#include "G4GeneralParticleSource.hh"
#include "DRAGONPhysicsList.hh"

#include "global_variables.hh"

namespace DRAGON
{

void DRAGONPrimaryGeneratorAction::gukine_mitray(G4Event* anEvent, G4GeneralParticleSource* fGPSParticleGun)
{
 //************************************************************************
 //*                                                                      *
 //*  For iswit[3] == 0 and iswit[4] == 0                                 *
 //*  -----------------------------------                                 *
 //*  KINE card:  IKINE    : GEANT particle type                          *
 //*              PKINE 1  : particle origin x position (cm)              *
 //*                    2  : particle origin y position (cm)              *
 //*                    3  : particle origin z position (cm)              *
 //*                    4  : particle central momentum value (MeV/c)      *
 //*                    5  : RAYTRACE theta [mr]                          *
 //*                    6  : RAYTRACE phi   [mr]                          *
 //*                                                                      *
 //*  For iswit[3] == 0 and iswit[4] == 1                                 *
 //*  -----------------------------------                                 *
 //*  KINE card:  IKINE    : GEANT particle type                          *
 //*              PKINE 1  : particle origin x position (cm)              *
 //*                    2  : particle origin y position (cm)              *
 //*                    3  : particle origin z position (cm)              *
 //*                    4  : particle central momentum value (MeV/c)      *
 //*                    5  : width of x origin [1cm]                      *
 //*                    6  : width of y origin [1cm]                      *
 //*                    7  : horizontal emittance [100mr]                 *
 //*                    8  : vertical   emittance [100mr]                 *
 //*                    9  : 1 +- deltaP/P                                *
 //*                                                                      *
 //************************************************************************
auto gv = GlobalVariables::GetInstance();

    G4int ip;

    G4double vertex[3], vrtx[3], plab[3], pbeam, scale;
    G4double x_max, y_max, rndm[2];
   
    G4int nout = 0, isotok = 0;
    G4double bmass, tmp[3];
    G4String label[3];
    
    G4double theta, phi;
    
    G4double xd[3], xdd[3];

    if (DRAGONRunAction_prim->GetIswit()[3] == 0 && DRAGONRunAction_prim->GetIswit()[4] == 0)
       {
        ip = ikine;
        
        vertex[0] = pkine[0];
	    vertex[1] = pkine[1];
		vertex[2] = pkine[2];
		
		pbeam = pkine[3]/1.E3;              //! GEANT wants GeV/c
		
		plab[0] = std::sin(pkine[4]*1.E-3)*std::cos(pkine[5]*1.E-3);
		plab[1] =                          std::sin(pkine[5]*1.E-3);
		plab[2] = std::cos(pkine[4]*1.E-3)*std::cos(pkine[5]*1.E-3);

        if (DRAGONRunAction_prim->GetIrevs() == 1) plab[2] = -plab[2];
        }
    else
       {
		if (DRAGONRunAction_prim->GetIswit()[3] >= 1 && DRAGONRunAction_prim->GetIswit()[4] == 0)
           {
			std::ifstream infile;
            if (nout == 0)
               {
                 infile.open((rayfile.c_str()));
                
                if (!infile.is_open()) 
                   {
                    G4cerr << " No RAYFILE file opened for simulation! " << G4endl;
                    std::exit(1);
                    } 
                }
                
            while (infile >> pkine[0] >> pkine[4] >> pkine[1] >> pkine[5] >> pkine[3]
                   >> bmass >> tmp[0] >> tmp[1] >> pkine[2] >> tmp[2] >> tmp[3]
                   >> nout >> label[0] >> label[1] >> label[2]) 
                  {isotok = isotok + 1;} 
                
           if(DRAGONRunAction_prim->GetIswit()[3] > 1 && (isotok < DRAGONRunAction_prim->GetIswit()[3])) ;
           
           pkine[2] = 1.E2 * pkine[2];
 
           ip = 61;
  
           vertex[0] = pkine[0];
		   vertex[1] = pkine[1];
		   vertex[2] = pkine[2];
    
           pkine[3] = 1.E3 * pkine[3];   //! read in GeV/c -> MeV/c (as above)
    
           pbeam = pkine[3] / 1.E3;      //! GEANT wants GEV/c
    
           plab[0] = std::sin(pkine[4]*1.E-3)*std::cos(pkine[5]*1.E-3);
           plab[1] =                          std::sin(pkine[5]*1.E-3);
           plab[2] = std::cos(pkine[4]*1.e-3)*std::cos(pkine[5]*1.E-3);
    
           if(DRAGONRunAction_prim->GetIrevs() == 1)
             {
			  plab[0]=-plab[0]; 
			  plab[1]=-plab[1]; 
			  plab[2]=-plab[2];   
		      }
           }
        else
           {
			if(DRAGONRunAction_prim->GetIswit()[3] == 0 && DRAGONRunAction_prim->GetIswit()[4] == 1)
			  {
			  ip = ikine;
			  
				vertex[0] = pkine[0];
				vertex[1] = pkine[1];
				vertex[2] = pkine[2];
				
			    rndm[0] = G4RandGauss::shoot(0,1.0);
                rndm[1] = G4RandGauss::shoot(0,1.0); 
			     
                vertex[0] = vertex[0] + (2.0*rndm[0]-1.0) * pkine[4]/2.;
                vertex[1] = vertex[1] + (2.0*rndm[1]-1.0) * pkine[5]/2.;
        
                rndm[0] = G4RandGauss::shoot(0,1.0);
                pbeam = (1.0+(pkine[8]-1.0)*(2.0*rndm[0]-1.0)) * pkine[3];
                pbeam = pbeam / 1.E3;              //! GEANT wants GeV/c
 
                //C *** Pick a uniform direction for the particle into a rectangular aperture
                
                x_max = std::tan(pkine[6]*1.E-3/2.);
                y_max = std::tan(pkine[7]*1.E-3/2.);
 
                rndm[0] = G4RandGauss::shoot(0,1.0);
                rndm[1] = G4RandGauss::shoot(0,1.0); 

                plab[0] = 2.*x_max*rndm[0] - x_max;
                plab[1] = 2.*y_max*rndm[1] - y_max;
                plab[2] = 1.0;

                if(DRAGONRunAction_prim->GetIrevs() == 1)plab[2] = -plab[2];
  			    }
  		    }
        }
        
    if(plab[0] || plab[1] || plab[2]) 
      {       
       plab[0] = plab[0]/(std::sqrt(plab[0]*plab[0]+plab[1]*plab[1]+plab[2]*plab[2]));   
       plab[1] = plab[1]/(std::sqrt(plab[0]*plab[0]+plab[1]*plab[1]+plab[2]*plab[2]));   
       plab[2] = plab[2]/(std::sqrt(plab[0]*plab[0]+plab[1]*plab[1]+plab[2]*plab[2]));        
       }

    xd[0] = vertex[0];
	xd[1] = vertex[1];
	xd[2] = vertex[2];
    	
    xdd[0] = plab[0];
	xdd[1] = plab[1];
	xdd[2] = plab[2];
    
    //CALL gdtom(xd,vertex,1)
    //CALL gdtom(xdd, plab,2)
    
    plab[0] = plab[0]*pbeam;
    plab[1] = plab[1]*pbeam;
    plab[2] = plab[2]*pbeam;
    
    beamenerg = std::sqrt(fphys->Getbeammom()*fphys->Getbeammom() + fphys->Getbeammass()*fphys->Getbeammass()) - fphys->Getbeammass();
    fGPSParticleGun->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(vertex[0],vertex[1],vertex[2]));
    fGPSParticleGun->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.0));
    fGPSParticleGun->GetCurrentSource()->GetEneDist()->SetMonoEnergy(fphys->Getbeamenerg());
    fGPSParticleGun->GeneratePrimaryVertex(anEvent);
    
    if(DRAGONRunAction_prim->GetIswit()[3] >= 1)
      {
       scale = -vertex[2]/plab[2];
       vrtx[2] = vertex[2]+plab[2]*scale; 
       }
    else
       {
 	    vrtx[0] = vertex[0];
 	    vrtx[1] = vertex[1];
 	    vrtx[2] = vertex[2];
	    }        
    
    G4AnalysisManager* analysis = G4AnalysisManager::Instance();  
    analysis->FillH1(gv->IDMap[1], vrtx[0]);
    analysis->FillH1(gv->IDMap[2], vrtx[1]);
    
    analysis->FillH1(gv->IDMap[101], vrtx[0],vrtx[1]);
    
    if(plab[0] || plab[1] || plab[2]) 
      {       
       plab[0] = plab[0]/(std::sqrt(plab[0]*plab[0]+plab[1]*plab[1]+plab[2]*plab[2]));   
       plab[1] = plab[1]/(std::sqrt(plab[0]*plab[0]+plab[1]*plab[1]+plab[2]*plab[2]));   
       plab[2] = plab[2]/(std::sqrt(plab[0]*plab[0]+plab[1]*plab[1]+plab[2]*plab[2]));        
       }
    
    theta = 0.0;
    if((plab[0] != 0.0) || (plab[2] != 0.0))
      theta = 1000.*std::atan2(plab[0],plab[2]);
    phi = 1000.*std::asin(plab[1]);
      
    analysis->FillH1(gv->IDMap[3], theta); 
    analysis->FillH1(gv->IDMap[4], phi); 
    
    analysis->FillH1(gv->IDMap[102], theta, phi); 
    analysis->FillH1(gv->IDMap[105], vrtx[0], theta); 
    analysis->FillH1(gv->IDMap[106], vrtx[1], phi); 
    
    pbeam = 100.*(pbeam/0.25855443-1.0);
    
    analysis->FillH1(gv->IDMap[9], pbeam);    
    }
    
    }

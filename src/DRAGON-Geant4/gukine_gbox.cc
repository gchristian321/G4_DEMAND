#include "DRAGONPrimaryGeneratorAction.hh"      //Geant4
#include "DRAGONDetectorConstruction.hh"  

#include "Randomize.hh"
#include "DRAGONHistoManager.hh"
#include "G4GeneralParticleSource.hh"
#include "DRAGONPhysicsList.hh"

#include "gukine_gbox.hh"                      //local

#include "rescom.hh"
#include "beamcom.hh"


namespace DRAGON
{

void DRAGONPrimaryGeneratorAction::gukine_gbox(G4Event* anEvent, G4GeneralParticleSource* fGPSParticleGun)
{

/***********************************************************************
*   	                                                               *
*        (generation of isotropic photon beam along line source)       *
*        -------------------------------------------------------       *
*                                                                      *
*  KINE card: MKINE    : number of photons at initial vertex           *
*             GKINE 0  : x of photon origin distribution [cm]          *
*                   1  : y of photon origin distribution [cm]          *
*                   2  : z of photon origin distribution [cm]          *
*                   3  : half length of photon origin x-dimension [cm] *
*                   4  : half length of photon origin y-dimension [cm] *
*                   5  : half length of photon origin z-dimension [cm] *
*             GKINE 6  : photon energy [MeV]                           *
*                   7  : theta [degree]                                *
*                   8  : phi [degree]                                  *
*                   9  : emittance angle [degree]                      *
*                                                                      *
*	         Adapted into Geant4 by Ben Marlow Oct-15-2024             *
***********************************************************************/
	auto gv = GlobalVariables::GetInstance();
	
	G4double beammom;
	
	//From beamcom.hh
	G4double beamenerg;
	
	//For now locally define the above KINE card, needs change when we make commands OJO
	//G4double gkine[10] = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 3.0, 0.0, 0.0, 0.0};
    
    G4double degrad = 180.0/3.141592654;
        G4int ip;
	G4int ngamma;
	G4double pbeam, vertex[3], plab[3];
	G4double theta, ctheta, stheta, phi, cphi, sphi;
	G4double dir[3], costh, sinth, cosph, sinph;
	G4bool rotate;
	G4double rndm[3];
	
    ip = 1;				//! Original particle is a photon
	
//C.
//C. --> First deal with the photon origin
//C.
    vertex[0] = gkine[0];
    vertex[1] = gkine[1];
    vertex[2] = gkine[2]; 
    
    
    for(G4int i=0; i<3; i++) {
    	rndm[i] = G4UniformRand();
    }
    
    vertex[0] = vertex[0] + gkine[3] * (2.0*rndm[0]-1.0);
    vertex[1] = vertex[1] + gkine[4] * (2.0*rndm[1]-1.0);
    vertex[2] = vertex[2] + gkine[5] * (2.0*rndm[2]-1.0);
    
    fGPSParticleGun->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(vertex[0],vertex[1],vertex[2]));
    
    ngamma = 0;
    
    if(gkine[7] == 0.0 && gkine[8] == 0.0 && gkine[9] == 0.0) {
    	rndm[0] = G4UniformRand();
    	rndm[1] = G4UniformRand();
    	
    	ctheta = 2.0 * rndm[0] - 1.0;
    	stheta = std::sqrt((1.0 - ctheta) * (1.0 + ctheta));
    	
    	phi = 2.0 * 3.141592654 * rndm[1];
    	cphi = std::cos(phi);
    	sphi = std::sin(phi);
    	
    	plab[0] = stheta * cphi;
    	plab[1] = stheta * sphi;
    	plab[2] = ctheta;
    }
    
    else {
//C.
//C. -->  Establish cone direction
//C.
    	if(gkine[7] != 0.0) {
    		theta = gkine[7];
    		phi = gkine[8];
    		
    		dir[0] = std::sin(degrad * theta) * std::cos(degrad * phi);
          	dir[1] = std::sin(degrad * theta) * std::sin(degrad * phi);
          	dir[2] = std::cos(degrad * theta);
    	}
    	else {
    		dir[0] = 0.0;
    		dir[1] = 0.0;
    		dir[2] = 1.0;	
    	}	
    	
    	//--> Establish direction in a cone
    	if(gkine[9] != 0) {
    		rndm[0] = G4UniformRand();
    		rndm[1] = G4UniformRand();
    		
    		ctheta = 1. - rndm[0] * (1. - std::cos(degrad * gkine[9]/2.));
			stheta = std::sqrt((1. - ctheta) * (1. + ctheta));

			phi  = 2.0 * 3.141592654 * rndm[1];
			cphi = std::cos(phi);
			sphi = std::sin(phi);

			plab[0] =  stheta * cphi;
			plab[1] =  stheta * sphi;
			plab[2] =  ctheta;
    	}
    	else {
    		plab[0] = 0.0;
    		plab[1] = 0.0;
    		plab[2] = 1.0;
    	}
//C.
//C. -->  Rotate cone into direction
//C.
    	G4double dux = dir[0];  // x-component
		G4double duy = dir[1];  // y-component
		G4double duz = dir[2];  // z-component
		G4double dsith2, dsith, dnorm;
		
		rotate = true;

		if (std::abs(duz) >= 0.85) {
		    dsith2 = dux * dux + duy * duy;

		    if (dsith2 > 0.0) {
		        costh = std::copysign(std::sqrt(1.0 - dsith2), duz);  // cos(θ)
		        dsith = std::sqrt(dsith2);  // sin(θ)
		        sinth = dsith;
		        cosph = dux / dsith;  // cos(φ)
		        sinph = duy / dsith;  // sin(φ)
		    } else if (duz > 0.0) {
		        rotate = false;
		        costh = 1.0;
		        sinth = 0.0;
		        cosph = 1.0;
		        sinph = 0.0;
		    } else {
		        costh = -1.0;
		        sinth = 0.0;
		        cosph = 1.0;
		        sinph = 0.0;
		    }
		} else {
		    costh = duz;  // cos(θ)
		    dsith = std::sqrt((1.0 + duz) * (1.0 - duz));  // sin(θ)
		    sinth = dsith;
		    dnorm = 1.0 / std::sqrt(dux * dux + duy * duy);
		    cosph = dux * dnorm;  // cos(φ)
		    sinph = duy * dnorm;  // sin(φ)
   		}
    	
    	if(rotate) {
    		G4double d0 = dir[0];
    		G4double d1 = dir[1];
    		G4double d2 = dir[2];
    		dir[0]=d0*costh*cosph - d1*sinph + d2*sinth*cosph;
      		dir[1]=d0*costh*sinph + d1*cosph + d2*sinth*sinph;
      		dir[2]=-d0*sinth                 + d2*costh;
    	}
    	
    	ngamma = ngamma+1;
    	
    	pbeam = gkine[6]; //Energy in MeV
    	if(ngamma <= 10) {
    		//"egamma" is an array of length 10 that is set in the ffcards and
    		//   contains the gamma energies
    		//Filled with dummy values for now, will need to change when we implement commands OJO
    		if(egamma[ngamma-1] > 0.0) pbeam = egamma[ngamma-1];
    	}
    	
    	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    	G4ParticleDefinition* particle = particleTable->FindParticle("gamma");
  		fGPSParticleGun->SetParticleDefinition(particle);
		
		//Geant4 needs energy, so convert from momentum relativistically
      	beamenerg = std::sqrt(beammom*beammom + beammass*beammass) - beammass;
      
    	fGPSParticleGun->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(plab[0],plab[1],plab[2]));
    	fGPSParticleGun->GetCurrentSource()->GetEneDist()->SetMonoEnergy(pbeam*MeV);
    	fGPSParticleGun->GeneratePrimaryVertex(anEvent);
    }
    
}    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
}

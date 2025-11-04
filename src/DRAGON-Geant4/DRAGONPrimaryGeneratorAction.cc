//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file DRAGON/src/DRAGONPrimaryGeneratorAction.cc
/// \brief Implementation of the DRAGON::DRAGONPrimaryGeneratorAction class
#include "G4GeneralParticleSource.hh"                    //Geant4
#include "G4SPSEneDistribution.hh"
#include "G4SPSPosDistribution.hh"
#include "G4SPSAngDistribution.hh"
#include "G4RunManager.hh"
#include "G4IonTable.hh"
#include "G4DecayTable.hh"
#include "G4VDecayChannel.hh"
#include "G4ParticleDefinition.hh"
#include "G4UnitsTable.hh"
#include "G4Decay.hh"

#include "DRAGONPrimaryGeneratorAction.hh"               //local
#include "DRAGONDetectorConstruction.hh"
#include "DRAGONPhysicsList.hh"
#include "DRAGONPrimaryGeneratorActionMessenger.hh"
#include "geant3functions.hh"


namespace DRAGON
{

DRAGONPrimaryGeneratorAction::DRAGONPrimaryGeneratorAction() 
:emax(0.),buncht(0.),egamma{0}, offset{0.}, rayfile("rayfile.dat")
{
  uvinit();
//! uses a 0.6% detector efficiency convolution with 2 exponentials and a fermi function to fit the observed spectrum
fFileName = "hist6p2e";
//C *** Define the beam energy and tune scale for non-resonant reactions  //From ugffgo.f
  beamenerg = 0;   //!Default for resonance  case
  
  fGPSParticleGun = new G4GeneralParticleSource();
  fGPSMessenger = new DRAGONPrimaryGeneratorActionMessenger(this); 
     
  fGPSParticleGun->GetCurrentSource()->SetNumberOfParticles(1);
     
  //fGPSParticleGun->GetCurrentSource()->SetParticleTime(0.);
  G4SPSPosDistribution *posDist = fGPSParticleGun->GetCurrentSource()->GetPosDist();
  posDist->SetPosDisShape("Point");
  
  G4SPSEneDistribution *eneDist = fGPSParticleGun->GetCurrentSource()->GetEneDist();
  eneDist->SetEnergyDisType("Mono"); 
  
  
  /*
  fGPSParticleGun = new G4GeneralParticleSource();
  fGPSMessenger   = new DRAGONPrimaryGeneratorActionMessenger(this); 
     
  fGPSParticleGun->SetNumberOfParticles(1);

    // Posición del haz (origen)
    G4SPSPosDistribution* posDist = fGPSParticleGun->GetCurrentSource()->GetPosDist();
    posDist->SetPosDisShape("Point");
    posDist->SetCentreCoords(G4ThreeVector(0.,0.,0.));

    // Energía y dirección
    G4SPSEneDistribution* eneDist = fGPSParticleGun->GetCurrentSource()->GetEneDist();
    eneDist->SetEnergyDisType("Mono"); 
    eneDist->SetMonoEnergy(1.0*MeV);  // Energía cinética de 1 MeV

    fGPSParticleGun->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
*/

 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



void DRAGONPrimaryGeneratorAction::uvinit()
{
 // From uvinit.f
 //C.
 //C.    *************** USER DEFAULTS ***************
 //C.
 //C.
 ikine = 0;                    //! No DRAGON test particle
 //C.
 pkine[0] =   0.0;             //! Particle origin x [cm]
 pkine[1] =   0.0;             //! Particle origin y [cm]
 pkine[2] =   0.0;             //! Particle origin z [cm]
 pkine[3] = 258.55443;         //! Particle momentum [MeV/c]
 pkine[4] =   0.0;             //! horizontal emittance [mr]
 pkine[5] =   0.0;             //! vertical   emittance [mr]
 pkine[6] =   0.0;             //! RAYTRACE theta [mr]
 pkine[7] =   0.0;             //! RAYTRACE phi   [mr]
 pkine[8] =   1.0;             //! 1 +- deltaP/P
  
 //C.
 mkine     =   0;           //! number of photons at initial vertex
 gkine[0] =   0.0;          //! x of photon origin distribution [cm]
 gkine[1] =   0.0;          //! y of   "      "         "       [cm]
 gkine[2] =   0.0;          //! z of   "      "         "       [cm]
 gkine[3] =   0.0;          //! length of photon origin x-dimension [cm]
 gkine[4] =   0.0;          //! length of photon origin y-dimension [cm]
 gkine[5] =   0.0;          //! length of photon origin z-dimension [cm]
 gkine[6] =   4.03;         //! photon energy [MeV]
 gkine[7] =   0.0;          //! theta [degree]
 gkine[8] =   0.0;          //! phi [degree]
 gkine[9] =   0.0;          //! emittance [mrad] 
 //C.

 vzero(egamma,10);
}

DRAGONPrimaryGeneratorAction::~DRAGONPrimaryGeneratorAction()
{
  delete fGPSParticleGun;
  delete fGPSMessenger; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DRAGONPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
   gukine(anEvent, fGPSParticleGun);

/*
    G4ParticleDefinition* ion = G4ParticleTable::GetParticleTable()->FindParticle("Co60copia");
    if (!ion) {
        G4Exception("DRAGONPrimaryGeneratorAction", "Co60copiaNotFound",
                    FatalException, "Co60copia no encontrada en la tabla de partículas.");
        return;
    }

    fGPSParticleGun->GetCurrentSource()->SetParticleDefinition(ion);

    fGPSParticleGun->GeneratePrimaryVertex(anEvent); 
*/
	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void DRAGONPrimaryGeneratorAction::SetEGAM(const G4String& input) {
    std::istringstream iss(input);
    G4double value;
    int count = 0;

    while (iss >> value && count < 10) {
        egamma[count++] = value;
    }     
}


void DRAGONPrimaryGeneratorAction::SetKINE(const G4String& input, G4int option) {
    //Function for setting values from "KINE" command
    if(option == 0) {
    	std::istringstream iss(input);
		iss >> ikine;
		for (G4int i = 0; i < 10; ++i) {
		    if (!(iss >> pkine[i])) {
		        std::cout << "Invalid input string in method SetKINE. Expected 10 numbers.";
		    }
    	}
    }
    
    else if(option == 1) {
    	std::istringstream iss(input);
		iss >> ikine;
		for (G4int i = 0; i < 10; ++i) {
		    if (!(iss >> pkine[i])) {
		        std::cout << "Invalid input string in method SetKINE. Expected 10 numbers.";
		    }
    	}
    }
}

void DRAGONPrimaryGeneratorAction::SetGKIN(const G4String& input) {
    //Function for setting values from "GKIN" command
    std::istringstream iss(input);
    iss >> mkine;
    for (G4int i = 0; i < 10; ++i) {
        if (!(iss >> gkine[i])) {
            std::cout << "Invalid input string in method SetGKIN. Expected 11 numbers.";
        }
    }
}

void DRAGONPrimaryGeneratorAction::Build523(G4String fFileName) 
     {
      std::ifstream file(fFileName+".txt");   
      if (!file.is_open()) 
         { 
    	  std::cout << "Cannot open input file: " << fFileName << ".txt" << std::endl;
      	  G4RunManager::GetRunManager()->AbortRun();
          }
      
      G4double energy, count;
      G4double total_count = 0.0;
      
      while (file >> energy >> count) 
            {
    	     energies.push_back(energy);
             counts.push_back(count);
             total_count += count;
             cdf.push_back(total_count);
	     }
      for (auto& value : cdf) 
          {value /= total_count;}
  
      file.close();
      }

    //user defined angular distribution
void DRAGONPrimaryGeneratorAction::Build250() {
	//C.    user defined angular distribution
  G4double num_bins = 1000.;
  G4double min = -1.;
  G4double max = 1.;
  
  G4double bin_width = (max-min)/num_bins;
  G4double costheta;
  
   for (int i = 1; i <= num_bins; ++i) 
       {
		costheta = min + (i - 0.5) * bin_width;
        angdist_costheta.push_back(costheta);
        angdist_values.push_back(angdist(costheta));
        }
}

G4double DRAGONPrimaryGeneratorAction::angdist(G4double costheta){
   //	c---- 67---- gamma angular distribution
   G4double pi = 3.14159265358979;

   //C  A uniform angular distribution for gammas
   G4double angdist = 1;
   //CC  A dipole angular distribution for gammas
   //C    angdist = (3./(8.*pi))*(1.-costheta*costheta);
   //CC  A quad. angular distribution for gammas
   //C   angdist = (15./(8.*pi))*(1.-costheta*costheta)*costheta*costheta;

   return angdist;
   }

void DRAGONPrimaryGeneratorAction::Build501() {
    m1 = fphys->Getm1();
    m2 = fphys->Getm2();
    zprod = fphys->Getzprod();
    er = fphys->Geter();
    gp = fphys->Getgp();
    gg = fphys->Getgg();
    omg = fphys->Getomg();
    ell = fphys->Getell();
    beamo = fphys->Getbeamo();
    emax = fphys->Getemax();
    beamenerg = fphys->Getbeamenerg();
    
    lm1 = (1.-0.005)*beamo*m2/(m1+m2);
    lm2 = (1.+3.*(emax/beamenerg))*beamenerg*m2/(m1+m2);
  
  //creation of a 1-dim id and filling with FUNC
  G4double num_bins = 1000.;
  G4double bin_width = (lm2-lm1)/num_bins;
  G4double e;
  
   for (int i = 1; i <= num_bins; ++i) 
       {
        e = lm1 + (i - 0.5) * bin_width;
        energies.push_back(e);
        sig_values.push_back(sig(e));
        }
  }


G4double DRAGONPrimaryGeneratorAction::sig(G4double e) 
{
 //c----67---- (p,g) cross-section

    G4int q = 1;  //OJO Este es el Q de la reaccion
  
    G4double cv = 931.494;
    
    G4double pi = 3.141592654;
    
    G4double rm = (m1/cv) * (m2/cv)/((m1/cv)+(m2/cv));
    
    //total width
    G4double gt = gg + gp;
    
    //Sommerfeld parameter (on resonance)
    G4double etar = 0.1575*z1*z2*std::sqrt(rm/er);
    
    //wavenumber (on resonance)
    G4double wnr = 0.219*std::sqrt(rm*er);
    
    //BW cross-section (on resonance)
    G4double sigr = (omg*pi/(wnr*wnr))*gp*gg/((0.5*gt)*(0.5*gt));
    
    sigr = sigr/100.; //into barns
    
    //S-Factor (on resonance)
    G4double srr = sigr*er*std::exp(2.0*pi*etar);
    
    //-----67----- energy dependent part
    
    G4double eta = 0.1575 * z1 * z2 *std::sqrt(rm/e);
    G4double wn = 0.219 *std::sqrt(rm*e);
    
    //Equations 10 P. Decrock et al, Phys. Rev. C 48 (1993), 2057
    G4double gge = gg * std::pow((q + e), (2.0 * ell + 1.0)) / std::pow((q + er), (2.0 * ell + 1.0));
    //G4double gpe = gp * std::exp(-std::sqrt(eg / e)) / std::exp(-std::sqrt(eg / er));
    G4double gpe = gp*std::exp(-2.0*pi*(eta-etar));
   
    G4double gte = gge+gpe;
    
    G4double d1 = std::exp(2.0 * pi * eta) * gpe * gge * std::pow(0.5 * gt, 2);
    G4double d2 = std::exp(2.0 * pi * etar) * gp * gg * (std::pow(e - er, 2) + std::pow(0.5 * gte, 2));
    G4double d3 = d1 / d2;
    G4double sr = srr * d3;
    
    G4double sig = (omg * pi / std::pow(wn, 2)) * gpe * gge / (std::pow(e - er, 2) + std::pow(0.5 * gte, 2));
    sig /= 100.0;
    G4double sn2 = sig * e * std::exp(2.0 * pi * eta);

    //!      Equation 15 P. Decrock et al, Phys. Rev. C 48 (1993), 2057
    
    // snr = c*c *s * 3.44e-4 * std::exp(-0.605 *e)
    G4double snr = 0;
    
    //std::cout << eta << ", " << wn << ", " << gge << ", " << gpe << ", " << gte << ", " << sr;
    
    //convert into keV.barns
    
    snr = snr*1000.0;
    sr = sr*1000.0;
    
    //interference phase factor
    
    G4double del = std::atan(gte/(2.*(e-er)));
    if(del > 0.) del = pi-del;
    
    //addition of resonant and non-resonant parts
    
    G4double se = sr + snr + 2.*std::sqrt(sr*snr)*std::cos(del); //keV.b
    sig = se *exp(-2.0*pi*eta)/(e*1000.0); //barns
   
return sig;
}


G4double DRAGONPrimaryGeneratorAction::HRNDM1(G4int IDD) {
 	if(IDD == 523 || IDD == 501) {
 	G4double rndnum = G4UniformRand(); 
    auto it = std::lower_bound(cdf.begin(), cdf.end(), rndnum);
    size_t index = std::distance(cdf.begin(), it);

    if (index >= energies.size()) return energies.back();
    	return energies[index];
    }

 return 0.;
 }
  



}




























#include "DRAGONSensitiveDetector.hh"
#include "DRAGONDetectorConstruction.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"  
#include "G4ios.hh"
#include "DRAGONHit.hh"

#include "Materials.hh"

#include "geant3functions.hh" 

namespace DRAGON
{ 

DRAGONSensitiveDetector::DRAGONSensitiveDetector(const G4String& name,DRAGONDetectorConstruction* det)
    : G4VSensitiveDetector(name),
      fHitsCollection(nullptr),DRAGON_det(det),
      fHitCollectionID(-1) {
    collectionName.insert(name+"HitsCollection");
}

DRAGONSensitiveDetector::~DRAGONSensitiveDetector() {}

void DRAGONSensitiveDetector::Initialize(G4HCofThisEvent* hce) {
   if (!hce) {
        G4cerr << "DRAGONSensitiveDetector::Initialize: null G4HCofThisEvent pointer" << G4endl;
        return;
    }
    fHitsCollection = new G4THitsCollection<DRAGONHit>(SensitiveDetectorName, collectionName[0]);

    if (fHitCollectionID < 0) {
        fHitCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    }

    hce->AddHitsCollection(fHitCollectionID, fHitsCollection);
}

G4bool DRAGONSensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory* /*history*/) {
	   in_new_vol = 0;
       if (step->GetPreStepPoint()->GetStepStatus() == fGeomBoundary)
		   in_new_vol = 1;
	   G4StepStatus preStatus  = step->GetPreStepPoint()->GetStepStatus();
       G4StepStatus postStatus = step->GetPostStepPoint()->GetStepStatus();
	   G4VPhysicalVolume* preVolume = step->GetPreStepPoint()->GetPhysicalVolume();
       G4VPhysicalVolume* postVolume = step->GetPostStepPoint()->GetPhysicalVolume();
       if (preVolume != nullptr && postVolume != nullptr && preVolume == postVolume) inwvol = 0;    
       if (preStatus == fGeomBoundary) inwvol = 1;
       if (postStatus == fGeomBoundary) inwvol = 2;
       if (postStatus == fWorldBoundary ) inwvol = 3; 
       chname_nlevel = step->GetPreStepPoint()->GetTouchable()->GetVolume()->GetName();
       noCopy = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1);
	   
	   G4Track* track = step->GetTrack();
       const G4ParticleDefinition* particle = track->GetParticleDefinition();
	   if (particle->GetParticleType() == "opticalphoton") itrtyp = 7;
	   Materials* mates = Materials::Instance();
       const G4Material* numed_ = step->GetPreStepPoint()->GetMaterial();
       G4Material* n_detmate_ = mates->GetMatByIndex(n_detmate); 
	   auto touch = step->GetPreStepPoint()->GetTouchable();
	   nlevel = touch->GetHistoryDepth();
       number = step->GetPreStepPoint()->GetTouchable()->GetVolume(nlevel)->GetCopyNo(); 
	   tofg = track->GetGlobalTime()/s;
	   edep = step->GetTotalEnergyDeposit();
	   if(this->GetName() == "HSNG")idtype = 1;
       if(this->GetName() == "PMT")idtype = 2;
	   if(this->GetName() == "ENDV")idtype = 3;
	   istop = (track->GetTrackStatus() == fAlive) ? 0 : 1;
	   G4ThreeVector pos = step->GetPostStepPoint()->GetPosition();
       vect[0] = pos.x()/cm;
       vect[1] = pos.y()/cm;
       vect[2] = pos.z()/cm;
   
       if(idtype == 2)ghipmt(step);
	   if(idtype == 1)ghidet(step);
		  
       if(savehit)
	     {
	      DRAGONHit* hit = new DRAGONHit();
	      hit->Setx_fngr(hits[0]);
	      hit->Sety_fngr(hits[1]);
	      hit->Setz_fngr(hits[2]);
	      hit->Settofg(hits[3]);
	      hit->Setedep(hits[4]);
	      hit->SetPMTNumber(noCopy);
	    
              fHitsCollection->insert(hit);
	   
	      in_new_vol = 0;
		  savehit = false;
		  }

       return true;
}

void DRAGONSensitiveDetector::ghipmt(const G4Step*)
{
	   //From ghipmt.f
	   ///////////////////////////////////////////////

          if(in_new_vol = 1 && itrtyp == 7)
	        {
	         if(chname_nlevel == "PMT")
  		       {
          	    hits[0] = DRAGON_det->x_fngr[number];
				hits[1] = DRAGON_det->y_fngr[number];
				hits[2] = DRAGON_det->z_fngr[number];
				hits[3] = tofg;
				hits[4] = 1.0;
				savehit = true;
			    }
		    }	        
		///////////////////////////////////////////////		
}
	

void DRAGONSensitiveDetector::ghidet(const G4Step*)
{
 	   //From ghidet.f
	   ///////////////////////////////////////////////
          if(numed == n_detmate && itrtyp != 7)
	        {
			 if(idtype <= 0)return;
	         if(inwvol == 1 && istop != 0)return;
	   
//C.
//C *** INWVOL = 1     store entrance quantities
//C.
	   
	         if (inwvol == 1)
	            {
                 edep = 0.;
                 tofin = tofg;
		         ucopy(vect,x1,3);
		         gmtod(x1,xd,1,chname_nlevel);
		         }

             if(inwvol == 0)
	           {
		        if(istop != 0)
	              {
	               if(inwvol == 1)
			         {				 
			   	      edep = 0.;
			  	      tofin = tofg;
				      }
                   hits[0] = (x1[0]+x2[0])/2.;
                   hits[1] = (x1[1]+x2[1])/2.;
                   hits[2] = (x1[2]+x2[2])/2.;
                   hits[3] = (tofin + tofout)/2. * 1.E9;
                   hits[4] = 1000. * edep;
				   
				   if(edep > 1.E-6)savehit = true;
				   }
	            }
	        } 
	   ///////////////////////////////////////////////
	
}


void DRAGONSensitiveDetector::EndOfEvent(G4HCofThisEvent* /*hce*/) 
{
}

}


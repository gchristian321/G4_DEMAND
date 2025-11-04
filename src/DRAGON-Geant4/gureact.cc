
#include "G4VProcess.hh"                //Geant4
#include "G4ParticleChange.hh"
#include "G4StackManager.hh"
#include "G4ParticleTable.hh"
#include "G4EventManager.hh"
#include "G4RunManager.hh"

#include "DRAGONSteppingAction.hh"      //local
#include "DRAGONIon.hh"
#include "DRAGONRunAction.hh"

#include <cstdlib> 

//C.
//C======================================================================C
//C                                                                      C
//C                KINEMATICS OF CAPTURE GAMMA REACTION                  C
//C                                                                      C
//C======================================================================C
//C.

namespace DRAGON {



class G4DRAGONRESProcess : public G4VProcess {
public:
    G4DRAGONRESProcess() : G4VProcess("RESR", fUserDefined) {}
    ~G4DRAGONRESProcess() override {}

    G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&) override {return &fParticleChange;}

    G4VParticleChange* AlongStepDoIt(const G4Track&, const G4Step&) override {return &fParticleChange;}

    G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&) override {return &fParticleChange;}

    G4double PostStepGetPhysicalInteractionLength(const G4Track&, G4double, G4ForceCondition* condition) override 
            {
             *condition = NotForced;
             return DBL_MAX;
             }

    G4double AlongStepGetPhysicalInteractionLength(const G4Track&, G4double, G4double, G4double&, G4GPILSelection*) override {return DBL_MAX;}

    G4double AtRestGetPhysicalInteractionLength(const G4Track&, G4ForceCondition* condition) override 
             {
              *condition = NotForced;     
              return DBL_MAX;
              }

private:
    G4ParticleChange fParticleChange;
};


void DRAGONSteppingAction::gureact(const G4Step* step)
{
 ngkine = 1;

 kcase = "RESR";

 ipart = 81;
  
 auto recoilDef = G4ParticleTable::GetParticleTable()->FindParticle(ipart); 
 G4String napart = recoilDef->GetParticleName();
 G4String itrtyp = recoilDef->GetParticleType();
 G4double charge = recoilDef->GetPDGCharge();
 G4double tlife = recoilDef->GetPDGLifeTime();

 std::cout << "Creating particle " << ipart << std::endl;
 std::cout << ipart << ", " << napart << ", " << itrtyp << ", " << newm << ", " << charge << ", " << tlife << ", " << std::endl;
 
 G4ThreeVector pos = G4ThreeVector(vect[0]*cm,vect[1]*cm,vect[2]*cm);
 G4ThreeVector p = G4ThreeVector((vect[3]*vect[6])*GeV,(vect[4]*vect[6])*GeV,(vect[5]*vect[6])*GeV);
 G4double time = step->GetTrack()->GetGlobalTime();
 
 G4Track* track = step->GetTrack();

 // Posición actual (PostStepPoint)
 G4ThreeVector pos1 = track->GetPosition();
 G4ThreeVector mom = track->GetMomentum();
 G4double E_tot = track->GetTotalEnergy(); 
 G4double p_mag = mom.mag(); 

 G4cout << "Real track Position: x=" << pos1.x()/cm 
       << " cm, y=" << pos1.y()/cm 
       << " cm, z=" << pos1.z()/cm << " cm" << G4endl;

 G4cout << "Real track Momentum: px=" << mom.x()/GeV 
       << " GeV/c, py=" << mom.y()/GeV 
       << " GeV/c, pz=" << mom.z()/GeV << " GeV/c" 
       << G4endl;
       
 G4cout << "Real momentum magnitude: " << p_mag/GeV << " GeV/c" << G4endl;

 G4cout << "Real kinetic Energy: " << track->GetKineticEnergy()/GeV << " GeV" << G4endl;

 G4cout << "Real total Energy: " << E_tot/GeV << " GeV" << G4endl;

 std::cout << "x: " << vect[0] << " cm" << std::endl;
 std::cout << "y: " << vect[1] << " cm" << std::endl;
 std::cout << "z: " << vect[2] << " cm" << std::endl;
 std::cout << "I_px: " << vect[3] << std::endl;
 std::cout << "I_py: " << vect[4] << std::endl;
 std::cout << "I_pz: " << vect[5] << std::endl;
 std::cout << "px: " << vect[3]*vect[6] << " GeV/c" << std::endl;
 std::cout << "py: " << vect[4]*vect[6] << " GeV/c" << std::endl;
 std::cout << "pz: " << vect[5]*vect[6] << " GeV/c" << std::endl;

 G4DynamicParticle* recoilDyn = new G4DynamicParticle(recoilDef, p);
 recoilDyn->SetMass(newm*GeV);
 G4cout << "New dynamic particle mass: " << recoilDyn->GetMass()/GeV << " GeV/c2" << G4endl;

 G4Track* newTrack = new G4Track(recoilDyn, time, pos);
 newTrack->SetParentID(step->GetTrack()->GetTrackID());
 static G4DRAGONRESProcess* DRAGONRESProcess = new G4DRAGONRESProcess();
 newTrack->SetCreatorProcess(DRAGONRESProcess);
                         
 G4StackManager* stackManager = G4EventManager::GetEventManager()->GetStackManager();
 stackManager->PushOneTrack(newTrack);
 }

}



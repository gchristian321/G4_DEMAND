#ifndef DemandPrimaryGeneratorMessenger_h
#define DemandPrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcommand;
class DemandPrimaryGeneratorAction;

class DemandPrimaryGeneratorMessenger: public G4UImessenger
{
public:
	DemandPrimaryGeneratorMessenger(DemandPrimaryGeneratorAction*);
	virtual ~DemandPrimaryGeneratorMessenger();
	virtual void SetNewValue(G4UIcommand*, G4String);

private:
	G4UIdirectory*        fPrimaryDir;
	G4UIcmdWithAnInteger* fBeamZ;
	G4UIcmdWithAnInteger* fBeamA;
	G4UIcmdWithAnInteger* fSourceZ;
	G4UIcmdWithAnInteger* fSourceA;
	G4UIcmdWithADoubleAndUnit *fBeamEnergy;
	G4UIcmdWithAString *fReactionCmd;
	G4UIcmdWithAString *fAngularDistributionCmd;
	G4UIcmdWithADoubleAndUnit *fExRecoil;
	G4UIcmdWithADoubleAndUnit *fBeamSigmaX;
	G4UIcmdWithADoubleAndUnit *fBeamSigmaY;
	G4UIcmdWithADoubleAndUnit *fBeamSigmaThetaX;
	G4UIcmdWithADoubleAndUnit *fBeamSigmaThetaY;
	G4UIcmdWith3VectorAndUnit *fSourceEnergies;

	DemandPrimaryGeneratorAction* fPrimary;
};

#endif

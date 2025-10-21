#ifndef DemandRunMessenger_h
#define DemandRunMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcommand;
class DemandRunAction;

class DemandRunMessenger: public G4UImessenger
{
public:
	DemandRunMessenger(DemandRunAction*);
	virtual ~DemandRunMessenger();
	virtual void SetNewValue(G4UIcommand*, G4String);

private:
	G4UIdirectory* fRunDir;
	G4UIcmdWithAString* fGeant3SetupCmd;

	DemandRunAction* fRun;
};

#endif

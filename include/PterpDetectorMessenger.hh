#ifndef PterpDetectorMessenger_h
#define PterpDetectorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class PterpDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcommand;
class G4UIcmdWithABool;
class G4UIcmdWithADouble;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;

class PterpDetectorMessenger: public G4UImessenger
{
  public:

    PterpDetectorMessenger(PterpDetectorConstruction*);
    virtual ~PterpDetectorMessenger();
 
    virtual void SetNewValue(G4UIcommand*, G4String);
 
  private:

    PterpDetectorConstruction*   fPterpDetector;
    G4UIdirectory*               fDetectorDir;
    G4UIcmdWithADoubleAndUnit*   fDimensionsCmdX;
    G4UIcmdWithADoubleAndUnit*   fDimensionsCmdY;
    G4UIcmdWithADoubleAndUnit*   fDimensionsCmdZ;
    G4UIcmdWithADoubleAndUnit*   fPosCmdX;
    G4UIcmdWithADoubleAndUnit*   fPosCmdY;
    G4UIcmdWithADoubleAndUnit*   fPosCmdZ;
	  G4UIcmdWithADoubleAndUnit*   fDistanceCmd;
	  G4UIcmdWithADoubleAndUnit*   fThetaXCmd;
	  G4UIcmdWithADoubleAndUnit*   fPhiCmd;
	  G4UIcmdWithADoubleAndUnit*   fThetaCmd;
  	G4UIcmdWithADoubleAndUnit*   fThresholdCmd;
    G4UIcmdWithAnInteger*        fNxCmd;
    G4UIcmdWithAnInteger*        fNyCmd;
    G4UIcmdWithAnInteger*        fNzCmd;
	  G4UIcmdWithAString *         fReadoutTypeCmd;
	  G4UIcmdWithABool *           fRotateCmd;
	  G4UIcmdWithADoubleAndUnit *fTargetThicknessCmd;
	  G4UIcmdWithAString *fTargetMaterialCmd;
	  G4UIcmdWithoutParameter* fAddModuleCmd;

	  size_t fCurrentModuleNumber;
};

#endif

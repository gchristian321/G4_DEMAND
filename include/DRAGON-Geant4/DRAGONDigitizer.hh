#ifndef DRAGONDigitizer_h
#define DRAGONDigitizer_h 1

#include "G4VDigitizerModule.hh"
#include "DRAGONDigi.hh"

namespace DRAGON
{

class DRAGONDigitizerMessenger;   
class DRAGONEventAction;
class DRAGONDetectorConstruction;

class DRAGONDigitizer : public G4VDigitizerModule {
public:
    DRAGONDigitizer(const G4String& name);
    virtual ~DRAGONDigitizer();

    virtual void Digitize() override;
    void SetANAL(G4double val){E_threshold = val;}
    void SetTHLD(const G4String& input);
    void uvinit();
    void SetEventAction(DRAGONEventAction* EventAction) { fEventAction = EventAction;}
    void SetDetConst(DRAGONDetectorConstruction* DetectorConst) {fDetector = DetectorConst;}

    void gudigi(); 
    void anadet(G4int itrig);
    void genclu(G4int nelem, G4int jelem[][500], G4int& nclu, G4int nhit[], G4int jhit[][10]);
    void dhpmt(); 
    void interac();

      static const G4int mclu = 10;
      static const G4int mhits = 10;

//C.    nclu              Number of energy clusters found

      G4int nclu;

//C.    nhit[iclu]:       Number of neighbours on hit list
//C.    jhit[ihit][iclu]:  Volume # of entry on hit list

      G4int nhit[mclu], jhit[mhits][mclu];

      G4double eclu[mclu], xclu[mclu], yclu[mclu], zclu[mclu];

      G4double dir_clu[3][mclu];

	  G4double raddeg = 180.0/M_PI;
      G4int nhits;
      G4double e_detect, z_react;  //OJO no se de donde viene z_react
      G4int neff[11][29], ifngr_max;
      G4double gammatof;
      G4double melem_gbox;

      const static G4int mhit_gbox = 20;
      G4int jelem_gbox[1][mhit_gbox];
      G4double energy_gbox[mhit_gbox];
      const static G4int max_hexagon = 30;
      const static G4int Nn = 10;
      G4int n_fngr[Nn][max_hexagon];   //OJO quitar de aqui y tomar de DetConst
      G4int iswit[10];                 //OJO quitar de aqui y tomar de RunAction
      
    G4double GetE_threshold(){return E_threshold;}
     

	
private:
    DRAGONDigitizerMessenger* fDRAGONDigitizerMessenger = nullptr;
    DRAGONEventAction* fEventAction = nullptr;
    DRAGONDetectorConstruction* fDetector = nullptr;
    G4double E_threshold, pmt_thrshld, tot_thrshld;


};

}

#endif

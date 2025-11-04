#ifndef GLOBAL_VARIABLES_HH
#define GLOBAL_VARIABLES_HH

#include "G4Types.hh"

#include "mitray_setup.hh"

namespace DRAGON
{


class GlobalVariables{
      public:
             static GlobalVariables* GetInstance() {
        if (!Instance) {
            Instance = new GlobalVariables();
        }
        return Instance;
    }
             
    G4bool GetDialog() const { return ldiag; }
    void SetDialog(G4bool value) { ldiag = value; }  
    
    G4int GetIswit(G4int index) const { return iswit[index]; }        
    void SetIswit(G4int index,G4int value) { iswit[index] = value; }   
    
    G4int GetNo() const { return no; }
    void SetNo(G4int value) {no = value;} 
    
    G4double Getbscale() const { return bscale; }
    void Setbscale(G4double value) { bscale = value; }
    
    G4double Getescale() const { return escale; }
    void Setescale(G4double value) { escale = value; }
    
    void FillItitle(const G4String values[nmax]){
    for (G4int i = 0; i < nmax; ++i) {
        ititle[i] = values[i];
    }
    }

    void FillIdata(const G4int values[nmax]){
    for (G4int i = 0; i < nmax; ++i) {
        idata[i] = values[i];
    }
    }
    
    void FillData(const G4double values[mmax][nmax]) {
    for (G4int i = 0; i < mmax; ++i) {     
        for (G4int j = 0; j < nmax; ++j) { 
            data[i][j] = values[i][j];      
        }
    }
}

  G4String GetItitle(G4int index){return ititle[index];}
  G4int GetIdata(G4int index){return idata[index];}
  G4double GetData(G4int row,G4int col){return data[row][col];}


      public:
            GlobalVariables() : ldiag(false) {}
    static GlobalVariables* Instance;
    
		G4bool ldiag;
		G4int iswit[6]; 
		G4int no;
		G4String ititle[nmax];
		G4int idata[nmax];
		G4double data[mmax][nmax];
		G4double escale = 3.27479124;
		G4double bscale = 1.7610544;
	
		std::map<G4int,int> IDMap;
   
      };

}

#endif

#ifndef MATERIALS_HH
#define MATERIALS_HH 1

#include "G4NistManager.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"

#include <map>

class Materials{
      public:
             static Materials* Instance();
	         void DefineMaterials();
			 				  
             void ugmate();
             void ugmate_trgt();
             void ugstmed();
             void ugstmed_trgt();

             G4int Getatarg(){return atarg;}
             G4int Getmtarg(){return mtarg;}
             G4int* Getment(){return ment;}
             G4double Gettargetl(){return targetl;}
             G4double Getexitdens(){return exitdens;}
             G4double Getentdens(){return entdens;}
			
             void SetTARG(G4int targ){atarg = targ;}
             
             G4Material* Hydrogen;
             G4Material* Deuterium;
             G4Material* Helium;
             G4Material* Lithium;
             G4Material* Beryllium;
             G4Material* Carbon;
             G4Material* Nitrogen;
             G4Material* Neon;
             G4Material* Aluminium;
             G4Material* Iron;
             G4Material* Copper;
             G4Material* Tungsten;
             G4Material* Lead;
             G4Material* Uranium;
             G4Material* Air;
             G4Material* Vacuum;
             G4Material* Scintillator;
             G4Material* BariumFluoride;
             G4Material* CesiumFluoride;
             G4Material* SodiumIodide;
             G4Material* CesiumIodide;
             G4Material* BGO;
             G4Material* LSO;
             G4Material* MgO;
             G4Material* Glass;
			 G4Material* Silicon;
			 G4Material* Target_Carbon;
			 G4Material* ModifiedSilicon;
			 G4Material* MCP_Carbon;
             G4Material* StainlessSteel;
             G4Material* Target;
             G4Material* CentralVacuum;
             G4Material* Entrance1;
             G4Material* Entrance2;
             G4Material* Entrance3;
             G4Material* Ext1;
             G4Material* Ext2;
             G4Material* Ext3;
             G4Material* Baseline; 

             G4int atarg, mcent, mtarg, ment[6], mex[6], mbox;
             G4double targetl, exitdens, entdens;
             G4double bulk_absorption, paint_absorption;
			 
	     std::map<int, G4Material*> materialMap;
	     G4Material* GetMatByIndex(int index) const;


      private:
              Materials(){DefineMaterials();}
              
              static Materials* instance;
      };

#endif

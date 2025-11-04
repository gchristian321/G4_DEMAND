#include "G4SystemOfUnits.hh"       

#include "Materials.hh"                          //local

void Materials::ugmate()
     {
      //From ugmate.f  
      //************************************************************************
      //*                                                                      *
      //*                 Routine to define tracking material                  *
      //*                                                                      * 
      //************************************************************************    
      
      G4NistManager* nistManager = G4NistManager::Instance();
      // Scintillator                                                                        //17
      G4Element* C = nistManager->FindOrBuildElement("C");                                    
      G4Element* H = nistManager->FindOrBuildElement("H"); 
      Scintillator = new G4Material("SCINTILLATOR", 1.032*g/cm3, 2);
      Scintillator->AddElement(C, 0.915);
      Scintillator->AddElement(H, 0.085);
	  materialMap[9] = Scintillator;

      // Barium Fluoride (BaF2)                                                              //18
      G4Element* Ba = nistManager->FindOrBuildElement("Ba");
      G4Element* F = nistManager->FindOrBuildElement("F");
      BariumFluoride = new G4Material("BARIUM FLORIDE BAF2", 4.89*g/cm3, 2);
      BariumFluoride->AddElement(Ba, 0.783);
      BariumFluoride->AddElement(F, 0.217);
	  materialMap[10] = BariumFluoride;

      // Cesium Fluoride (CsF)                                                               //19
      G4Element* Cs = nistManager->FindOrBuildElement("Cs");
      CesiumFluoride = new G4Material("CESIUM FLORIDE CSF", 4.64*g/cm3, 2);
      CesiumFluoride->AddElement(Cs, 0.875);
      CesiumFluoride->AddElement(F, 0.125);
	  materialMap[11] = CesiumFluoride;

      // Sodium Iodide (NaI:Tl)                                                              //20
      G4Element* Na = nistManager->FindOrBuildElement("Na");
      G4Element* I = nistManager->FindOrBuildElement("I");
      SodiumIodide = new G4Material("SODIUM IODIDE NAI:TL", 3.67*g/cm3, 2);
      SodiumIodide->AddElement(Na, 0.153);
      SodiumIodide->AddElement(I, 0.847);
	  materialMap[12] = SodiumIodide;

      // Cesium Iodide (CsI:Tl)                                                              //21
      CesiumIodide = new G4Material("CESIUM IODIDE CSI:TL", 4.51*g/cm3, 2);
      CesiumIodide->AddElement(Cs, 0.512);
      CesiumIodide->AddElement(I, 0.488);
	  materialMap[13] = CesiumIodide;

      // BGO (Bi4Ge3O12)                                                                     //22
      G4Element* Bi = nistManager->FindOrBuildElement("Bi");
      G4Element* Ge = nistManager->FindOrBuildElement("Ge");
      G4Element* O = nistManager->FindOrBuildElement("O");
      BGO = new G4Material("BGO BI4GE3O12", 7.13*g/cm3, 3);
      BGO->AddElement(Bi, 0.671);
      BGO->AddElement(Ge, 0.175);
      BGO->AddElement(O, 0.154);
	  materialMap[14] = BGO;

      // LSO (Lu2(SiO4)O:Ce)                                                                 //23
      G4Element* Lu = nistManager->FindOrBuildElement("Lu");
      G4Element* Si = nistManager->FindOrBuildElement("Si");
      LSO = new G4Material("LSO LU2(SI04)O:CE", 7.4*g/cm3, 3);
      LSO->AddElement(Lu, 0.764);
      LSO->AddElement(Si, 0.061);
      LSO->AddElement(O, 0.175);
	  materialMap[15] = LSO;
      
       // MgO (Powder)                                                                       //24
      G4Element* Mg = nistManager->FindOrBuildElement("Mg");
      MgO = new G4Material("MGO (POWDER)", 1.87*g/cm3, 2);
      MgO->AddElement(Mg, 0.603);
      MgO->AddElement(O, 0.397);
	  materialMap[16] = MgO;

      // Glass                                                                               //25
      Glass = new G4Material("GLASS", 1.032*g/cm3, 2);
      Glass->AddElement(C, 0.915);
      Glass->AddElement(H, 0.085);
	  materialMap[17] = Glass;

      //SILICON
      Silicon = new G4Material("SILICON", 14.0, 28.08*g/mole, 2.33*g/cm3);                   //50
	  materialMap[19] = Silicon;
      
      //TARGET CARBON      
      Target_Carbon = new G4Material("TARGET CARBON", 6.0, 12.0*g/mole, 0.0225*g/cm3);       //62    ! TARGET CARBON density 100 x less dense than C.
      materialMap[21] = Target_Carbon;
	  
      //Modified SILICON
      ModifiedSilicon = new G4Material("Modified SILICON", 14.0, 28.08*g/mole, 0.0233*g/cm3);//51    ! Modified SILICON 100 x less dense than Si.
      materialMap[23] = ModifiedSilicon;

      //MCP Carbon
      MCP_Carbon = new G4Material("MCP Carbon", 6.0, 12.0*g/mole, 0.00225*g/cm3);            //52    ! MCP Carbon density 1000 x less dense than C. 
      materialMap[24] = MCP_Carbon;
	  
      //Stainless Steel                                                                      //26
      G4Element* Cr = nistManager->FindOrBuildElement("Cr");
      G4Element* Ni = nistManager->FindOrBuildElement("Ni");
      G4Element* Mn = nistManager->FindOrBuildElement("Mn");
      G4Element* Fe = nistManager->FindOrBuildElement("Fe");
      
      StainlessSteel = new G4Material("STAINLESS STEEL", 7.705*g/cm3, 6);
      StainlessSteel->AddElement(Fe, 0.709);
      StainlessSteel->AddElement(Cr, 0.180);
      StainlessSteel->AddElement(Ni, 0.080);
      StainlessSteel->AddElement(Mn, 0.020); 
      StainlessSteel->AddElement(Si, 0.010);
      StainlessSteel->AddElement(C, 0.002);
	  materialMap[20] = StainlessSteel;
	  materialMap[26] = StainlessSteel;

      ugmate_trgt();
}










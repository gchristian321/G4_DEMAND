#include "G4SystemOfUnits.hh"       

#include "Materials.hh"                          //local

void Materials::ugstmed() 
     {
     
//C.
//************************************************************************
//*                                                                      *
//*                   Routine to define tracking media                   *
//*                                                                      *
//************************************************************************
//C.
      const G4int npckov = 32;
      G4double absco[npckov], effic[npckov], rindex[npckov]; 
      

      G4double ppckov[npckov] = {2.038*eV, 2.072*eV, 2.107*eV, 2.143*eV, 2.181*eV,
                                 2.220*eV, 2.260*eV, 2.302*eV, 2.346*eV, 2.391*eV,
                                 2.438*eV, 2.486*eV, 2.537*eV, 2.590*eV, 2.645*eV,
                                 2.702*eV, 2.763*eV, 2.825*eV, 2.891*eV, 2.960*eV,
                                 3.032*eV, 3.108*eV, 3.188*eV, 3.271*eV, 3.360*eV,
                                 3.453*eV, 3.552*eV, 3.656*eV, 3.767*eV, 3.884*eV,
                                 4.010*eV, 4.144*eV};
                                 
      G4double absco_scnt[npckov] = { 344.8*cm,  408.2*cm,  632.9*cm,  917.4*cm, 1234.6*cm, 1388.9*cm,
                                     1515.2*cm, 1724.1*cm, 1886.8*cm, 2000.0*cm, 2631.6*cm, 3571.4*cm,
                                     4545.5*cm, 4761.9*cm, 5263.2*cm, 5263.2*cm, 5555.6*cm, 5263.2*cm,
                                     5263.2*cm, 4761.9*cm, 4545.5*cm, 4166.7*cm, 3703.7*cm, 3333.3*cm,
                                     3000.0*cm, 2850.0*cm, 2700.0*cm, 2450.0*cm, 2200.0*cm, 1950.0*cm,
                                     1750.0*cm, 1450.0*cm};

      G4double rindex_scnt[npckov] = {1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82,
                                      1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82,
                                      1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82,
                                      1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82,
                                      1.82, 1.82, 1.82, 1.82};
                        
      //G4double effic_pmt[npckov] = {0.005, 0.01, 0.02,  0.03, 0.04,  0.05,  0.06,  0.07,
      //                               0.08, 0.09, 0.10, 0.115, 0.13,  0.15,  0.16,  0.18,
      //                              0.195, 0.22, 0.23,  0.24, 0.25, 0.255,  0.26, 0.265,
      //                               0.26, 0.25, 0.24, 0.215,0.175,  0.14, 0.085,  0.0};

      G4double effic_pmt[npckov] = { 0.02, 0.025, 0.03, 0.035,  0.04, 0.05, 0.075, 0.09,
                                     0.12,  0.14, 0.15, 0.175, 0.185, 0.20,  0.21, 0.22,
                                     0.25,  0.26, 0.27,  0.28,  0.30, 0.30, 0.295, 0.29,
                                    0.285,  0.28, 0.26,  0.20, 0.175, 0.10,  0.05, 0.0};
                                    
      auto createMaterialPropertiesTable = [&](G4double* rindex, G4double* absco, G4double* effic) 
           {
            G4MaterialPropertiesTable* MPT = new G4MaterialPropertiesTable();
            MPT->AddProperty("RINDEX", ppckov, rindex, npckov);
            MPT->AddProperty("ABSLENGTH", ppckov, absco, npckov);
            MPT->AddProperty("EFFICIENCY", ppckov, effic, npckov);
            return MPT;
            };

      std::copy(std::begin(absco_scnt), std::end(absco_scnt), absco);   //! dielectric - scintillator                              
      std::fill_n(absco, npckov, bulk_absorption*cm); 
      std::copy(std::begin(effic_pmt), std::end(effic_pmt), effic);
      std::copy(std::begin(rindex_scnt), std::end(rindex_scnt), rindex);
    
      G4MaterialPropertiesTable* scintillatorMPT = createMaterialPropertiesTable(rindex, absco, effic);
      Scintillator->SetMaterialPropertiesTable(scintillatorMPT);
      BariumFluoride->SetMaterialPropertiesTable(scintillatorMPT);
      CesiumFluoride->SetMaterialPropertiesTable(scintillatorMPT);
      SodiumIodide->SetMaterialPropertiesTable(scintillatorMPT);
      CesiumIodide->SetMaterialPropertiesTable(scintillatorMPT);
      BGO->SetMaterialPropertiesTable(scintillatorMPT);
      LSO->SetMaterialPropertiesTable(scintillatorMPT);

      std::fill_n(effic, npckov, 0);                                    //! dielectric - air and vacuum
      std::fill_n(absco, npckov, 1.E10*cm);                            
      std::fill_n(rindex, npckov, 1.00);
   
      G4MaterialPropertiesTable* vacuumMPT = createMaterialPropertiesTable(rindex, absco, effic);
      Vacuum->SetMaterialPropertiesTable(vacuumMPT);
      Air->SetMaterialPropertiesTable(vacuumMPT);

      std::fill_n(effic, npckov, 0);                                    //! dielectric - glass
      std::fill_n(absco, npckov, 10.0*cm);                             
      std::fill_n(rindex, npckov, 1.50);                               
   
      G4MaterialPropertiesTable* glassMPT = createMaterialPropertiesTable(rindex, absco, effic);
      Glass->SetMaterialPropertiesTable(glassMPT);
                            
      std::fill_n(effic, npckov, 0);                                    //! metal - Al and MgO
      std::fill_n(rindex, npckov, 0);
      std::fill_n(absco, npckov, paint_absorption);      

      G4MaterialPropertiesTable* paintMPT = createMaterialPropertiesTable(rindex, absco, effic);    
      Aluminium->SetMaterialPropertiesTable(paintMPT);
      MgO->SetMaterialPropertiesTable(paintMPT);
      }














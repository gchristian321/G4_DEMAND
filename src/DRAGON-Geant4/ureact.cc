//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//C                                                                    C
//C               Define defaults of run-time variables                C
//C                                                                    C
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

#include "DRAGONPhysicsList.hh"

#include "G4DecayTable.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4ProcessManager.hh"
#include "G4Decay.hh"

#include "rescom.hh"
#include "geant3functions.hh"

#include "DRAGONIon.hh"
#include "Materials.hh"


namespace DRAGON 
{
void DRAGONPhysicsList::ureact()
     {
  
      G4int ireaction; 
      
      G4int ubuf[2];
      
      const G4double hbar = 6.582122E-22;
      const G4double amugev = 0.93149432;
      const G4double hmass  = 1.007825032*amugev;
      const G4double hemass = 4.002603250*amugev;
      const G4double deutmass = 2.0141018*amugev;
      const G4double c12mass = 12.0*amugev;
      
      G4double devmass,aamass,aamass1, elevel, aprod,  tlif;
      
      G4int mode[10] = {0};
      G4double brat[10] = {0.};
      G4String targtyp;
      G4String num[] = {"","1",  "2",  "3",  "4",  "5",  "6",  "7",  "8",  "9", "10",
                        "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
                        "21", "22", "23", "24", "25", "26", "27", "28", "29", "30"};

//C
//C======================================================================C
//C                                                                      C
//C      Initial beam and reaction information passed via FKIN card      C
//C                                                                      C
//C      LKINE reaction number                                           C
//C                                                                      C
//C      ( 1) 13N(p,g)14O                                                C
//C      ( 2) 15O(alpha,g)19Ne                                           C
//C      ( 3) 25Al(p,g)26Si                                              C
//C      ( 4) 17F(p,g)18Ne                                               C
//C      ( 5) 18F(p,g)19Ne                                               C
//C      ( 6) 19Ne(p,g)20Na                                              C
//C      ( 7) 20Na(p,g)21Mg 532 3/2-                                     C
//C      ( 8) 21Na(p,g)22Mg(220)                                         C
//C      ( 9) 23Mg(p,g)24Al                                              C
//C      (10) 26mAl(p,g)27Si                                             C
//C      (11) 7Be(p,g)8B                                                 C
//C      (12) 21Na(d,n)22Mg                                              C
//C      (13) 23Na(d,n)24Mg 2+                                           C
//C      (14) 23Na(p,g)24Mg 2+                                           C
//C      (15) 20Ne(p,g)21Na nonres                                       C
//C      (16) 20Na(p,g)21Mg131                                           C
//C      (17) 22Ne(alpha,g)26Mg                                          C
//C      (18) 21Na(p,g)22Mg(825)                                         C
//C      (19) 12C(a,g)16O                                                C
//C======================================================================C
//C
//C
//C-- MOD 10/06/03 C.Ruiz
//C-- Gamma cascades (isotropic emission) are included for
//C-- product nuclei through up to 32 states (including final 
//C-- ground state).
//C--
//C-- Each product energy level is defined as a GEANT particle.
//C-- Mass in GeV, lifetime and branching ratios to the lower states
//C-- are specified for each level.
//C--
//C-- Within this routine the array  BRAT and MODE   specify
//C-- the decay branches. MODE = part1 + 100*part2
//C-- where part1 is the GEANT particle number correspnding to the
//C-- energy level or gamma (=1).
//C--
//C--                   LABEL     IDPART
//C--  resonance                   81
//C--                              82
//C--                              83
//C--                              84
//C--                              85
//C--                              86
//C--                              .
//C--                              .
//C--  ground state              irecoil < 100
//C--
//C=======================================================================
//C       
//C     .                                        //From uginit.f
//C     .-->     Define Ne19++++ particle
//C 
      ipart = 61;
      aamass = 19.*0.93149432;    
      ubuf[0] = 4;
      tlif = 1000;
      DRAGONIon* ion = new DRAGONIon(
                    "Ne19++++",                  
                    aamass*GeV,      
                    0.0*keV,                    
                    ubuf[0]*eplus,          
                    5, +1, 0, 0, 0, 0,       
                    "nucleus",              
                    0, 0, ipart,               
                    false,                  
                    tlif*s,                 
                    nullptr,                
                    false,
                    10,
                    19            
                    );


/*
/////////////////////////////////////////////////////////////////////////////////////////////////////////  //OJO eliminar
// Define Co60copia (Cobalt-60)
DRAGONIon* Co60copia = new DRAGONIon(
    "Co60copia", 55.934936 * GeV, 0.0, 27. * eplus,
    0, 0, 0, 0, 0, 0,
    "nucleus", 0, 0, 16838,
    false, 100. * ps, nullptr, false,
    27, 60
);

// Define Ni60copia ground state (Nickel-60)
DRAGONIon* Ni60copia = new DRAGONIon(
    "Ni60copia", 55.942132 * GeV, 0.0, 28. * eplus,
    0, 0, 0, 0, 0, 0,
    "nucleus", 0, 0, 65460,
    true, -1.0, NULL, false,
    28, 60
);

// Define Ni60 excited states
DRAGONIon* Ni60_1173copia = new DRAGONIon(
    "Ni60_1173copia", 55.942132 * GeV + 1.173 * MeV, 0.0, 28. * eplus,
    0, 0, 0, 0, 0, 0,
    "nucleus", 0, 0, 65461,
    false, 1.0e-9 * second, nullptr, true, // Short-lived
    28, 60
);

DRAGONIon* Ni60_1332copia = new DRAGONIon(
    "Ni60_1332copia", 55.942132 * GeV + 1.332 * MeV, 0.0, 28. * eplus,
    0, 0, 0, 0, 0, 0,
    "nucleus", 0, 0, 65462,
    false, 1.0e-9 * second, nullptr, true,
    28, 60
);

DRAGONIon* Ni60_2505copia = new DRAGONIon(
    "Ni60_2505copia", 55.942132 * GeV + 2.505 * MeV, 0.0, 28. * eplus,
    0, 0, 0, 0, 0, 0,
    "nucleus", 0, 0, 65463,
    false, 1.0e-9 * second, nullptr, true,
    28, 60
);

// Decay table for Co60copia
G4DecayTable* co60DecayTable = new G4DecayTable();
co60DecayTable->Insert(new G4PhaseSpaceDecayChannel("Co60copia", 0.3333, 1, "Ni60_1173copia"));
co60DecayTable->Insert(new G4PhaseSpaceDecayChannel("Co60copia", 0.3333, 1, "Ni60_1332copia"));
co60DecayTable->Insert(new G4PhaseSpaceDecayChannel("Co60copia", 0.3334, 1, "Ni60_2505copia"));
Co60copia->SetDecayTable(co60DecayTable);

// Decay tables for Ni60 excited states
G4DecayTable* ni60_1173DecayTable = new G4DecayTable();
ni60_1173DecayTable->Insert(new G4PhaseSpaceDecayChannel("Ni60_1173copia", 1.0, 1, "Ni60copia"));
Ni60_1173copia->SetDecayTable(ni60_1173DecayTable);

G4DecayTable* ni60_1332DecayTable = new G4DecayTable();
ni60_1332DecayTable->Insert(new G4PhaseSpaceDecayChannel("Ni60_1332copia", 1.0, 1, "Ni60copia"));
Ni60_1332copia->SetDecayTable(ni60_1332DecayTable);

G4DecayTable* ni60_2505DecayTable = new G4DecayTable();
ni60_2505DecayTable->Insert(new G4PhaseSpaceDecayChannel("Ni60_2505copia", 1.0, 1, "Ni60copia"));
Ni60_2505copia->SetDecayTable(ni60_2505DecayTable);

PrintIonProperties(Co60copia);
PrintIonProperties(Ni60_1332copia);
PrintIonProperties(Ni60_2505copia);
PrintIonProperties(Ni60copia);

G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
G4ParticleDefinition* ionCo60copia = pTable->FindParticle("Co60copia");

if (ionCo60copia) {
    std::cout << "IOIOIOIOIO" << std::endl;
    G4ProcessManager* pManager = ionCo60copia->GetProcessManager();
    if (pManager) {
        // Agregar el proceso de decaimiento simple
        G4Decay* decay = new G4Decay();
        pManager->AddProcess(decay);
        pManager->SetProcessOrdering(decay, idxPostStep);
        pManager->SetProcessOrdering(decay, idxAtRest);
    }
}
else std::cout << "TTRTRTRTR" << std::endl;
/////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

//C     .
//C     .-->   Define radioactive ion reactions
//C     .
 
      ireaction = std::abs(lkine); 

      switch (ireaction) 
             {
              case 0:  
                     break;
              case 1:
                    {
//C
//C       ' (1) 13N(p,g)14O '
//C                     
                     zbeam =  7.;
                     abeam = 13.;
                     atarg = 1.;
                     aprod = atarg + abeam;
//C
                     resenerg = 0.526;
                     reswidth = 0.000037;
//C
                     elevel = 0.0;
//C                     
                     std::cout << "|**** 13N(p,g)14O reaction NOT implemented yet ****|" << std::endl;  
                     std::cout << "Resonance energy " << resenerg << " MeV " << ", " << "width " << reswidth << " MeV " << ", " << "level " << elevel << " MeV" << std::endl;   
                     break;
                     }
              case 2:
                    {
//C
//C       ' (2) 15O(alpha,g)19Ne '
//C
//C
                     ipart = 80;
                     zbeam =  8.;
                     abeam = 15.;                     
                     atarg = 4.;
                     zprod = 10.;
                     aprod = atarg + abeam;
//C
//C--     create beam particle --> idpart = 80
//C
                     devmass = 2855.4E-6;
                     aamass  = abeam*amugev + devmass;
                     tlif    = 122.2;
                     ubuf[0] = fkine[0];

                     DRAGONIon* ion = new DRAGONIon(
                                                 "O15",                  
                                                 aamass*GeV,      
                                                 0.0,                    
                                                 ubuf[0]*eplus,          
                                                 1, -1, 0, 0, 0, 0,       
                                                 "nucleus",              
                                                 0, +15, ipart,               
                                                 false,                  
                                                 tlif*s,                 
                                                 nullptr,                
                                                 false,
                                                 zbeam,
                                                 abeam               
                                                 );

//C
//C--     resonant level populated --> idpart = 81
//C                     
                     ipart = 81;
                     resenerg = 0.5036;
                     reswidth = 1.E-11;
                     aamass   = aamass + resenerg/1000. + hemass;
                     tlif     = hbar/reswidth;
                     ubuf[0]     = fkine[1];         
                     
                     ion = new DRAGONIon(
                           "res_Ne19_3/2+",                  
                           aamass*GeV,      
                           0.0*keV,                    
                           ubuf[0]*eplus,          
                           3, +1, 0, 0, 0, 0,       
                           "nucleus",              
                           0, 0, ipart,               
                           false,                  
                           tlif*s,                 
                           nullptr,                
                           false,
                           zbeam+2,
                           abeam+4                
                           );
                              
//C
//C--     Define states in Ne19 gamma decays from resonance
//C
                     ipart = 82;
                     devmass = 1751.1E-6;
                     prodm   = aprod*amugev + devmass;
//C
                     elevel  = 1.536;
                     aamass  = prodm + elevel/1000.;
                     tlif    = 2.8E-11;
//C                  
                     ion = new DRAGONIon(
                           "3_Ne19_3/2+",                  
                           aamass*GeV,      
                           0.0*keV,                    
                           ubuf[0]*eplus,          
                           3, +1, 0, 0, 0, 0,       
                           "nucleus",              
                           0, 0, ipart,               
                           false,                  
                           tlif*s,                 
                           nullptr,                
                           false,
                           zbeam+2,
                           abeam+4                
                           );
//C
                     ipart = 83;
                     elevel = .275;
                     aamass = prodm + elevel/1000.;
                     tlif   = 6.3E-11;
//C
                     ion = new DRAGONIon(
                           "2_Ne19_1/2",                  
                           aamass*GeV,      
                           0.0*keV,                    
                           ubuf[0]*eplus,          
                           1, +1, 0, 0, 0, 0,       
                           "nucleus",              
                           0, 0, ipart,               
                           false,                  
                           tlif*s,                 
                           nullptr,                
                           false,
                           zbeam+2,
                           abeam+4                                          
                           );
                   
//C
                     ipart = 84;
                     elevel =.238;
                     aamass = prodm + elevel/1000.;
                     tlif   = 2.6E-8;
//C
                     ion = new DRAGONIon(
                           "1_Ne19_5/2+",                  
                           aamass*GeV,      
                           0.0*keV,                    
                           ubuf[0]*eplus,          
                           5, +1, 0, 0, 0, 0,       
                           "nucleus",              
                           0, 0, ipart,               
                           false,                  
                           tlif*s,                 
                           nullptr,                
                           false,
                           zbeam+2,
                           abeam+4                                       
                           );
                     
//C
//C--     ground state --> idpart = 85
                     ipart = 85;
                     irecoil = 85;
//C
                     tlif = 1000.;                    

                     ion = new DRAGONIon(
                           "Ne19_1/2+",                  
                           prodm*GeV,      
                           0.0*keV,                    
                           ubuf[0]*eplus,          
                           1, +1, 0, 0, 0, 0,       
                           "nucleus",              
                           0, 0, ipart,               
                           false,                  
                           tlif*s,                 
                           nullptr,                
                           false,
                           zbeam+2,
                           abeam+4                                         
                           );                                   
//C
//C--     branch info -- resonance decays
//C
                     brat[0] = 80.;
                     mode[0] = 85;
                     brat[1] = 15.;
                     mode[1] = 83;
                     brat[2] = 5.;
                     mode[2] = 82;
        
                     G4DecayTable* decayTable1 = new G4DecayTable();
                     
                     for (int i = 0; i< 10; i++)
                         { 
                          if(mode[i])
                            {
                             G4VDecayChannel* decayChannel = new G4PhaseSpaceDecayChannel(
                                                     G4ParticleTable::GetParticleTable()->FindParticle(81)->GetParticleName(),     
                                                     brat[i]/100.,                        
                                                     1,                          
                                                     G4ParticleTable::GetParticleTable()->FindParticle(mode[i])->GetParticleName()                
                                                     );
                             decayTable1->Insert(decayChannel);
                             }
                          } 

                     DRAGONIon* part = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(81));
                     part->SetDecayTable(decayTable1);
                     decayTable1->DumpInfo();

                     vzero(brat,10);
                     vzero(mode,10);
//C
//C--     3/2+ state
//C
                     brat[0] = 95.;
                     mode[0] = 84;
                     brat[1] = 5.;
                     mode[1] = 83;
                     
                     G4DecayTable* decayTable2 = new G4DecayTable();

                     for (int i = 0; i< 10; i++)
                         { 
                         if(mode[i])
                           {
                            G4VDecayChannel* decayChannel = new G4PhaseSpaceDecayChannel(
                                                     G4ParticleTable::GetParticleTable()->FindParticle(82)->GetParticleName(),     
                                                     brat[i]/100.,                        
                                                     1,                          
                                                     G4ParticleTable::GetParticleTable()->FindParticle(mode[i])->GetParticleName()                
                                                     );
                            decayTable2->Insert(decayChannel);
                            }
                          }
                     DRAGONIon* part2 = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(82));
 
                     part2->SetDecayTable(decayTable2);
                     decayTable2->DumpInfo();                    
//C
                     vzero(brat,10);
                     vzero(mode,10);
//C
                     brat[0] = 100.;
                     mode[0] = 95;                   
                     
                     G4DecayTable* decayTable3 = new G4DecayTable();
        
                     for (int i = 0; i< 10; i++)
                         { 
                         if(mode[i])
                         {
                          G4VDecayChannel* decayChannel = new G4PhaseSpaceDecayChannel(
                                                     G4ParticleTable::GetParticleTable()->FindParticle(83)->GetParticleName(),     
                                                     brat[i]/100.,                        
                                                     1,                          
                                                     G4ParticleTable::GetParticleTable()->FindParticle(mode[i])->GetParticleName()                
                                                     );
                          decayTable3->Insert(decayChannel);
                          }
                          }
                     DRAGONIon* part3 = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(83));
                     part3->SetDecayTable(decayTable3);
                     decayTable3->DumpInfo(); 

                   G4DecayTable* decayTable4 = new G4DecayTable();

                   for (int i = 0; i < decayTable3->entries(); ++i) {
                       G4VDecayChannel* channel = decayTable3->GetDecayChannel(i);
                       G4VDecayChannel* clonedChannel = new G4PhaseSpaceDecayChannel(*dynamic_cast<G4PhaseSpaceDecayChannel*>(channel));
                       decayTable4->Insert(clonedChannel);
                       }

                   DRAGONIon* part4 = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(84));
                  part4->SetDecayTable(decayTable4);
                  decayTable4->DumpInfo();
//C                                    
                     std::cout << "|**** 15O(alpha,gamma)19Ne reaction ****|" << std::endl; 
                     std::cout << "Resonance energy " << resenerg << " MeV " << ", " << " width " << reswidth << " MeV " << ", " << "level " << level << " MeV" << std::endl; 
                     break; 
                     }
              case 3:
                    {
//C
//C       ' (3) 25Al(p,g)26Si '
//C
//C
                     ipart = 80;
                     zbeam = 13;
                     abeam = 25;                     
                     atarg = 1.;
                     zprod = 14.;
                     aprod = atarg + abeam;
//C
//C--     create beam particle --> idpart = 80
//C
                     devmass = -8915.7E-6;
                     aamass  = abeam*amugev + devmass;
                     tlif    = 7.18;
                     ubuf[0]    = fkine[0];
                     
                     DRAGONIon* ion = new DRAGONIon(
                                                 "Al25",                  
                                                 aamass*GeV,      
                                                 0.0*keV,                    
                                                 ubuf[0]*eplus,          
                                                 5, +1, 0, 0, 0, 0,       
                                                 "nucleus",              
                                                 0, 0, ipart,               
                                                 false,                  
                                                 tlif*s,                 
                                                 nullptr,                
                                                 false,
                                                 zbeam,
                                                 abeam                                                              
                                                 );
                               
//C
//C--     resonant level populated --> idpart = 81
//C
                     ipart = 80;
                     resenerg = 0.452;   //!calculated state...see C. Illiadis Phys
//C                                        Rev C53(1)(1995)475
                     reswidth = 0.00006;
                     aamass   = aamass + resenerg/1000. + hmass;
                     tlif     = hbar/reswidth;
                     ubuf[0]  = fkine[0]; 
//C 
                     ion = new DRAGONIon(
                           "res_Si26_3+",                  
                           aamass*GeV,      
                           0.0*keV,                    
                           ubuf[0]*eplus,          
                           6, +1, 0, 0, 0, 0,       
                           "nucleus",              
                           0, 0, ipart,               
                           false,                  
                           tlif*s,                 
                           nullptr,                
                           false,
                           zbeam+1,
                           abeam+1                                        
                           );
                                       
//C
//C--     Define states in Si26 gamma decays from resonance
//C 
                     ipart = 82;
                     devmass = -7145.E-6;
                     prodm   = aprod*amugev + devmass;
//C
                     elevel = 4.183;
                     aamass = prodm + elevel/1000.;
                     tlif    = 150.E-15;
//C 
                     ion = new DRAGONIon(
                           "4_Si26_3+",                  
                           aamass*GeV,      
                           0.0*keV,                    
                           ubuf[0]*eplus,          
                           6, +1, 0, 0, 0, 0,       
                           "nucleus",              
                           0, 0, ipart,               
                           false,                  
                           tlif*s,                 
                           nullptr,                
                           false,
                           zbeam+1,
                           abeam+1                                         
                           );
            
//C 
                     ipart = 83;
                     elevel = 3.756;
                     aamass = prodm + elevel/1000.;
                     tlif   = 700.E-15;
//C
                     ion = new DRAGONIon(
                           "3_Si26_3+",                  
                           aamass*GeV,      
                           0.0*keV,                    
                           ubuf[0]*eplus,          
                           6, +1, 0, 0, 0, 0,       
                           "nucleus",              
                           0, 0, ipart,               
                           false,                  
                           tlif*s,                 
                           nullptr,                
                           false,
                           zbeam+1,
                           abeam+1                                           
                           );
                     
//C
                     ipart = 84;
                     elevel = 2.783;
                     aamass = prodm + elevel/1000.;
                     tlif   = 210.E-15;
//C
                     ion = new DRAGONIon(
                                                 "2_Si26_2+",                  
                                                 aamass*GeV,      
                                                 0.0*keV,                    
                                                 ubuf[0]*eplus,          
                                                 4, +1, 0, 0, 0, 0,       
                                                 "nucleus",              
                                                 0, 0, ipart,               
                                                 false,                  
                                                 tlif*s,                 
                                                 nullptr,                
                                                 false,
                                                 zbeam+1,
                                                 abeam+1               
                                                 );

//C
                     ipart = 85;
                     elevel = 1.796;
                     aamass = prodm + elevel/1000.;
                     tlif   = 620.E-15;
//C 
                     ion = new DRAGONIon(
                           "1_Si26_2+",                  
                           aamass*GeV,      
                           0.0*keV,                    
                           ubuf[0]*eplus,          
                           4, +1, 0, 0, 0, 0,       
                           "nucleus",              
                           0, 0, ipart,               
                           false,                  
                           tlif*s,                 
                           nullptr,                
                           false,
                           zbeam+1,
                           abeam+1                                          
                           );
                     
//C
//C--     ground state --> idpart = 86
                     ipart = 86;
                     irecoil = 86;
//C
                     tlif = 1000.;
//C
                     ion = new DRAGONIon(
                           "Si26_0+",                  
                           prodm*GeV,      
                           0.0*keV,                    
                           ubuf[0]*eplus,          
                           0, +1, 0, 0, 0, 0,       
                           "nucleus",              
                           0, 0, 96,                    //OJO
                           false,                  
                           tlif*s,                 
                           nullptr,                
                           false,
                           zbeam+1,
                           abeam+1                                          
                           );
                    
//C
//C--     branch info -- resonance decays
//C
                     brat[0] = 87.;
                     mode[0] = 82;
                     brat[1] = 8.;
                     mode[1] = 83;
                     brat[1] = 5.;
                     mode[1] = 85;
        
                     G4DecayTable* decayTable1 = new G4DecayTable();
                     
                     for (int i = 0; i< 10; i++)
                         { 
                         if(mode[i])
                         {
                          G4VDecayChannel* decayChannel = new G4PhaseSpaceDecayChannel(
                                                     G4ParticleTable::GetParticleTable()->FindParticle(81)->GetParticleName(),     
                                                     brat[i]/100.,                        
                                                     1,                          
                                                     G4ParticleTable::GetParticleTable()->FindParticle(mode[i])->GetParticleName()                
                                                     );
                     
                          decayTable1->Insert(decayChannel);
                          }
                          }
                     DRAGONIon* part1 = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(81));
                     part1->SetDecayTable(decayTable1);
                     decayTable1->DumpInfo();
//C
                     vzero(brat,10);
                     vzero(mode,10);
//C
                     brat[0] = 47.;
                     mode[0] = 84;
                     brat[1] = 53.;
                     mode[1] = 85;
                     
                     G4DecayTable* decayTable2 = new G4DecayTable();
                            
                     for (int i = 0; i< 10; i++)
                         { 
                         if(mode[i])
                         {
                          G4VDecayChannel* decayChannel = new G4PhaseSpaceDecayChannel(
                                                     G4ParticleTable::GetParticleTable()->FindParticle(82)->GetParticleName(),     
                                                     brat[i]/100.,                        
                                                     1,                          
                                                     G4ParticleTable::GetParticleTable()->FindParticle(mode[i])->GetParticleName()                
                                                     );
                     
                          decayTable2->Insert(decayChannel);
                          }
                          }
                     DRAGONIon* part2 = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(82));
                     part2->SetDecayTable(decayTable2);
                     decayTable2->DumpInfo();
//C
                     vzero(brat,10);
                     vzero(mode,10);
//C
                     brat[0] = 70.;
                     mode[0] = 85;
                     brat[1] = 30.;
                     mode[1] = 84;
                     
                     G4DecayTable* decayTable3 = new G4DecayTable();
                            
                     for (int i = 0; i< 10; i++)
                         { 
                         if(mode[i])
                         {
                          G4VDecayChannel* decayChannel = new G4PhaseSpaceDecayChannel(
                                                     G4ParticleTable::GetParticleTable()->FindParticle(83)->GetParticleName(),     
                                                     brat[i]/100.,                        
                                                     1,                          
                                                     G4ParticleTable::GetParticleTable()->FindParticle(mode[i])->GetParticleName()               
                                                     );
                     
                          decayTable3->Insert(decayChannel);
                          }
                          }
                     DRAGONIon* part3 = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(83));
                     part3->SetDecayTable(decayTable3);
                     decayTable3->DumpInfo();
//C
                     vzero(brat,10);
                     vzero(mode,10);
//C
                     brat[0] = 69.;
                     mode[0] = 85;
                     brat[1] = 31.;
                     mode[1] = 86;
                     
                     G4DecayTable* decayTable4 = new G4DecayTable();
                            
                     for (int i = 0; i< 10; i++)
                         { 
                         if(mode[i])
                         {
                          G4VDecayChannel* decayChannel = new G4PhaseSpaceDecayChannel(
                                                     G4ParticleTable::GetParticleTable()->FindParticle(84)->GetParticleName(),     
                                                     brat[i]/100.,                        
                                                     1,                          
                                                     G4ParticleTable::GetParticleTable()->FindParticle(mode[i])->GetParticleName()                
                                                     );
                     
                          decayTable4->Insert(decayChannel);
                          }
                          }
                     DRAGONIon* part4 = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(84));
                     part4->SetDecayTable(decayTable4);
                     decayTable4->DumpInfo();
//C
                     vzero(brat,10);
                     vzero(mode,10);
//C
                     brat[0] = 100.;
                     mode[0] = 86;
                     
                     G4DecayTable* decayTable5 = new G4DecayTable();
   
                     for (int i = 0; i< 10; i++)
                         { 
                         if(mode[i])
                         {
                          G4VDecayChannel* decayChannel = new G4PhaseSpaceDecayChannel(
                                                     G4ParticleTable::GetParticleTable()->FindParticle(85)->GetParticleName(),     
                                                     brat[i]/100.,                        
                                                     1,                          
                                                     G4ParticleTable::GetParticleTable()->FindParticle(mode[i])->GetParticleName()                
                                                     );
                     
                          decayTable5->Insert(decayChannel);
                          }
                          }
                     DRAGONIon* part5 = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(85));
                     part5->SetDecayTable(decayTable5);
                     decayTable5->DumpInfo();
//C
                     std::cout << "|**** 25Al(p,g)26Si reaction ****|" << std::endl; 
                     std::cout << "Resonance energy " << resenerg << " MeV " << ", " << " width " << reswidth << " MeV " << ", " << "level " << level << " MeV" << std::endl; 
                     break;
                     }
              case 4:
                    {
//C
//C       ' (4) 17F(p,g)18Ne '
//C
//C
                     zbeam =  9;
                     abeam = 17;                     
                     atarg = 1.;
                     aprod = abeam +atarg;
//C
                     resenerg = 0.64;
                     reswidth = 0.;
//C
                     elevel = 2.67;                    
//C--     create beam particle --> idpart = 80
//C
                     std::cout << "|**** 17F(p,g)18Ne reaction NOT implemented yet ****|" << std::endl;  
                     std::cout << "Resonance energy " << resenerg << " MeV " << ", " << "width " << reswidth << " MeV " << ", " << "level " << elevel << " MeV" << std::endl; 
                     break;
                     }
              case 5:
                    {
//C
//C       ' (5) 18F(p,g)19Ne '
//C
//C
                     zbeam =  9;
                     abeam = 18;                       
                     atarg = 1.;
//C
                     resenerg = 0.3308;
                     reswidth = 0.0;
//C
                     elevel = 6.742;
//C--     create beam particle --> idpart = 80
//C
                     std::cout << "|**** 18F(p,g)19Ne reaction NOT implemented yet ****|" << std::endl;  
                     std::cout << "Resonance energy " << resenerg << " MeV " << ", " << "width " << reswidth << " MeV " << ", " << "level " << elevel << " MeV" << std::endl; 
                     break;
                     }
              case 6:
                    {
//C
//C       ' (5) 19Ne(p,g)20Na '
//C
//C
                     ipart = 80;
                     zbeam =  10;
                     abeam = 19;                       
                     atarg = 1.;
                     zprod = 11.;
                     aprod = atarg + abeam;
//C
//C--     create beam particle --> idpart = 80
//C
                     devmass = 1751.E-6;
                     aamass  = abeam*amugev + devmass;
                     tlif    = 1.;
                     ubuf[0] = fkine[0];
//C                     
                     DRAGONIon* ion = new DRAGONIon(
                                                 "Ne19",                  
                                                 aamass*GeV,      
                                                 0.0*keV,                    
                                                 ubuf[0]*eplus,          
                                                 1, +1, 0, 0, 0, 0,       
                                                 "nucleus",              
                                                 0, 0, ipart,               
                                                 false,                  
                                                 tlif*s,                 
                                                 nullptr,                
                                                 false,
                                                 zbeam,
                                                 abeam                
                                                 );
                             
//C
//C--     resonant level populated --> idpart = 81
//C
                     ipart = 81;
                     resenerg = 0.451;
                     reswidth = 7.E-9;       //!w=1 don't know spins
                     aamass   = aamass + resenerg/1000. + hmass;
                     tlif     = hbar/reswidth;
                     ubuf[0]  = fkine[2];
//C                     
                     ion = new DRAGONIon(
                           "res_Na20",                  
                           aamass*GeV,      
                           0.0*keV,                    
                           ubuf[0]*eplus,          
                           1, +1, 0, 0, 0, 0,       
                           "nucleus",              
                           0, 0, ipart,               
                           false,                  
                           tlif*s,                 
                           nullptr,                
                           false,
                           zbeam+1,
                           abeam+1                                        
                           );                        
                     
//C
//C--     Define states in Na20 gamma decays from resonance
//C
                     elevel = 2.660;
	             aamass = aamass - elevel/1000.;
//C
//C--     ground state --> idpart = 82
//C
                     ipart = 82;
                     irecoil = 82;
                     tlif = 1.;         //! made up
//C
                     ion = new DRAGONIon(
                           "gs_Na20",                  
                           prodm*GeV,      
                           0.0*keV,                    
                           ubuf[0]*eplus,          
                           4, +1, 0, 0, 0, 0,       
                           "nucleus",              
                           0, 0, ipart,               
                           false,                  
                           tlif*s,                 
                           nullptr,                
                           false,
                           zbeam+1,
                           abeam+1                                          
                           );                                    
//C
                     vzero(brat,10);
                     vzero(mode,10);
//C
//C--     branch info -- resonance decays
//C
                     brat[0] = 100.;
                     mode[0] = 82;
                                        
                     G4DecayTable* decayTable = new G4DecayTable();
                     
                     for (int i = 0; i< 10; i++)
                         { 
                         if(mode[i])
                         {
                          G4VDecayChannel* decayChannel = new G4PhaseSpaceDecayChannel(
                                                     G4ParticleTable::GetParticleTable()->FindParticle(81)->GetParticleName(),     
                                                     brat[i]/100.,                        
                                                     1,                          
                                                     G4ParticleTable::GetParticleTable()->FindParticle(mode[i])->GetParticleName()                
                                                     );
                     
                          decayTable->Insert(decayChannel);
                          }
                          }
                     DRAGONIon* part = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(81));
                     part->SetDecayTable(decayTable);
                     decayTable->DumpInfo();
                     
                     std::cout << "|**** 19Ne(p,gamma)20Na reaction ****|" << std::endl;  
                     std::cout << "Resonance energy " << resenerg << " MeV " << ", " << "width " << reswidth << " MeV " << ", " << "level " << elevel << " MeV" << std::endl; 
                     break;
                     }
              case 7:
                    {
//C
//C       ' (7) 20Na(p,g)21Mg536'
//C
//C
                     ipart = 80;
                     zbeam = 11;
                     abeam = 20;                       
                     atarg = 1.;
                     zprod = 12;
                     aprod = atarg + abeam;
//C--     create beam particle --> idpart = 80
//C
       	             devmass = 6845E-6;
       	             aamass  = abeam*amugev + devmass;
                     tlif    = 0.03;
                     ubuf[0] = fkine[0];
//C
                     DRAGONIon* ion = new DRAGONIon(
                                                 "Na20",                  
                                                 aamass*GeV,      
                                                 0.0*keV,                    
                                                 ubuf[0]*eplus,          
                                                 4, +1, 0, 0, 0, 0,       
                                                 "nucleus",              
                                                 0, 0, ipart,               
                                                 false,                  
                                                 tlif*s,                 
                                                 nullptr,                
                                                 false,
                                                 zbeam,
                                                 abeam                                                                
                                                 );
                                      
//C
//C--     resonant level populated --> idpart = 81
//C
                     ipart = 81;
                     resenerg = 0.536;
                     reswidth = hbar/tlif;
                     aamass   = aamass + resenerg/1000. + hmass;
                     ubuf[0]  = fkine[1];
                     tlif     = 4.0E-14;
//C                     
                     ion = new DRAGONIon(
                           "res_Mg21_536",                  
                           aamass*GeV,      
                           0.0*keV,                    
                           ubuf[0]*eplus,          
                           5, +1, 0, 0, 0, 0,       
                           "nucleus",              
                           0, 0, ipart,               
                           false,                  
                           tlif*s,                 
                           nullptr,                
                           false,
                           zbeam+1,
                           abeam+1                                         
                           );
                           
//C    Define states for Mg21 gamma decays from resonance state
                     ipart = 82;
                     devmass = 10912E-6;
                     prodm   = aprod*amugev + devmass;
//C
                     elevel = 0.0;
                     amass = prodm +elevel/1000.;           //OJO no tengo claro el origen
                     tlif    = 1000.;
                     irecoil = 82;
//C
                     ion = new DRAGONIon(
                           "Mg21_gs",                  
                           aamass*GeV,      
                           0.0*keV,                    
                           ubuf[0]*eplus,          
                           5, +1, 0, 0, 0, 0,       
                           "nucleus",              
                           0, 0, ipart,               
                           false,                  
                           tlif*s,                 
                           nullptr,                
                           false,
                           zbeam+1,
                           abeam+1               
                           );
                               
//C
                     vzero(brat,10);
                     vzero(mode,10);
//C                   
//C--     branch info -- resonance decays
//C
                     brat[0] = 100.;
                     mode[0] = 82;
                                        
                     G4DecayTable* decayTable = new G4DecayTable();
                     
                     for (int i = 0; i< 10; i++)
                         { 
                         if(mode[i])
                         {
                          G4VDecayChannel* decayChannel = new G4PhaseSpaceDecayChannel(
                                                     G4ParticleTable::GetParticleTable()->FindParticle(81)->GetParticleName(),     
                                                     brat[i]/100.,                        
                                                     1,                          
                                                     G4ParticleTable::GetParticleTable()->FindParticle(mode[i])->GetParticleName()                
                                                     );
                     
                          decayTable->Insert(decayChannel);
                          }
                          }
                     DRAGONIon* part = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(81));
                     part->SetDecayTable(decayTable);
                     decayTable->DumpInfo();

                     std::cout << "|**** 20Na(p,g)21Mg reaction ****|" << std::endl;  
                     std::cout << "Resonance energy " << resenerg << " MeV " << ", " << "width " << reswidth << " MeV " << ", " << "level " << elevel << " MeV" << std::endl; 
                     break;
                     }
              case 8:
                    {
//C
//C       ' (8) 21Na(p,g)22Mg(220) '
//C
//C
                     ipart = 80;
                     zbeam = 11;
                     abeam = 21;                       
                     atarg = 1.;
                     zprod = 12;
                     aprod = atarg + abeam;
//C
//C--     create beam particle --> idpart = 80
//C
       	             devmass = -2184.3E-6;
       	             aamass  = abeam*amugev + devmass;
                     tlif    = 0.03;
                     ubuf[0] = fkine[0];
//C
                     DRAGONIon* ion = new DRAGONIon(
                                                 "Na21",                  
                                                 aamass*GeV,      
                                                 0.0*keV,                    
                                                 ubuf[0]*eplus,          
                                                 3, +1, 0, 0, 0, 0,       
                                                 "nucleus",              
                                                 0, 0, ipart,               
                                                 false,                  
                                                 tlif*s,                 
                                                 nullptr,                
                                                 false,
                                                 zbeam,
                                                 abeam                                                 
                                                 );
                                          
//C
//C--     resonant level populated --> idpart = 81
//C
                     ipart    = 81;
                     tlif     = 4.E-14;
                     resenerg = 0.2124;
                     reswidth = hbar/tlif;
                     aamass   = aamass + resenerg/1000. + hmass;
                     ubuf[0]  = fkine[1];
//C
                     ion = new DRAGONIon(
                           "res_Mg22_2+",                  
                           aamass*GeV,      
                           0.0*keV,                    
                           ubuf[0]*eplus,          
                           4, +1, 0, 0, 0, 0,       
                           "nucleus",              
                           0, 0, ipart,               
                           false,                  
                           tlif*s,                 
                           nullptr,                
                           false,
                           zbeam+1,
                           abeam+1                                          
                           );  

//C
//C--     Define states in Mg22 gamma decays from resonance
//C
                     ipart   = 82; 
                     devmass = -396.8E-6;
                     prodm   = aprod*amugev + devmass;
//C
                     elevel  = 1.246;
                     aamass  = prodm + elevel/1000.;
                     tlif    = 3.E-11;
//C
                     ion = new DRAGONIon(
                           "Mg22_2+",                  
                           aamass*GeV,      
                           0.0*keV,                    
                           ubuf[0]*eplus,          
                           0, +1, 0, 0, 0, 0,       
                           "nucleus",              
                           0, 0, ipart,               
                           false,                  
                           tlif*s,                 
                           nullptr,                
                           false,
                           zbeam+1,
                           abeam+1                                         
                           );

//C
//C--     ground state --> idpart = 83
//C
                     ipart = 83;
                     irecoil = 83;
                     tlif   = 1000.;
//C
                     ion = new DRAGONIon(
                           "Mg22_0++",                  
                           prodm*GeV,      
                           0.0*keV,                    
                           ubuf[0]*eplus,          
                           0, +1, 0, 0, 0, 0,       
                           "nucleus",              
                           0, 0, ipart,               
                           false,                  
                           tlif*s,                 
                           nullptr,                
                           false,
                           zbeam+1,
                           abeam+1                                        
                           );
//C
                     vzero(brat,10);
                     vzero(mode,10);
//C
//C--     branch info -- resonance decays
//C
                     brat[0] = 87.;
                     mode[0] = 82;
                     brat[0] = 13.;
                     mode[0] = 83;                     
                                        
                     G4DecayTable* decayTable1 = new G4DecayTable();
                     
                     for (int i = 0; i< 10; i++)
                         { 
                         if(mode[i])
                         {
                          G4VDecayChannel* decayChannel = new G4PhaseSpaceDecayChannel(
                                                     G4ParticleTable::GetParticleTable()->FindParticle(81)->GetParticleName(),     
                                                     brat[i]/100.,                        
                                                     1,                          
                                                     G4ParticleTable::GetParticleTable()->FindParticle(mode[i])->GetParticleName()                
                                                     );
                     
                          decayTable1->Insert(decayChannel);
                          }
                          }
                     DRAGONIon* part1 = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(81));
                     part1->SetDecayTable(decayTable1);
                     decayTable1->DumpInfo();
//C
                     vzero(brat,10);
                     vzero(mode,10);
//C
                     brat[0] = 100.;
                     mode[0] = 83;
                     
                     G4DecayTable* decayTable2 = new G4DecayTable();
                     
                     for (int i = 0; i< 10; i++)
                         { 
                         if(mode[i])
                         {
                          G4VDecayChannel* decayChannel = new G4PhaseSpaceDecayChannel(
                                                     G4ParticleTable::GetParticleTable()->FindParticle(82)->GetParticleName(),     
                                                     brat[i]/100.,                        
                                                     1,                          
                                                     G4ParticleTable::GetParticleTable()->FindParticle(mode[i])->GetParticleName()                
                                                     );
                     
                          decayTable2->Insert(decayChannel);
                          }
                          }
                     DRAGONIon* part2 = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(82));
                     part2->SetDecayTable(decayTable2);
                     decayTable2->DumpInfo();
//C
                     std::cout << "|**** 21Na(p,g)22Mg reaction ****|" << std::endl;  
                     std::cout << "Resonance energy " << resenerg << " MeV " << ", " << "width " << reswidth << " MeV " << ", " << "level " << elevel << " MeV" << std::endl; 
                     break;
                     }
              case 9:
                    {
//C
//C       ' (9) 23Mg(p,g)24Al '
//C
//C
                     zbeam = 12;
                     abeam = 23;                       
                     atarg = 1.;
//C
                     resenerg = 0.51;
                     reswidth = 0.0;
//C
                     elevel = 2.38;
//C
                     std::cout << "|**** 23Mg(p,g)24Al reaction NOT implemented yet ****|" << std::endl;  
                     std::cout << "Resonance energy " << resenerg << " MeV " << ", " << "width " << reswidth << " MeV " << ", " << "level " << elevel << " MeV" << std::endl; 
                     break;
                     }
              case 10:
                     {
//C
//C       ' (10) 26mAl(p,g)27Si '
//C
//C
                     zbeam = 13;
                     abeam = 26;                       
                     atarg = 1.;
//C
                     resenerg = 0.201;
                     reswidth = 0.0;
//C
                     elevel = 7.893;           
//C
                     std::cout << "|**** 26mAl(p,g)27Si reaction NOT implemented yet ****|" << std::endl;  
                     std::cout << "Resonance energy " << resenerg << " MeV " << ", " << "width " << reswidth << " MeV " << ", " << "level " << elevel << " MeV" << std::endl; 
                     break;
                     }
              case 11:
                     {
//C
//C       ' (11) 7Be(p,g)8B '
//C
//C
                     zbeam = 4;
                     abeam = 7;                       
                     atarg = 1.;
//C
                     resenerg = 0.2;
                     reswidth = 0.;
//C
                     elevel = 0.338;
                     std::cout << "|**** 7Be(p,g)8B reaction NOT implemented yet ****|" << std::endl;  
                     std::cout << "Resonance energy " << resenerg << " MeV " << ", " << "width " << reswidth << " MeV " << ", " << "level " << elevel << " MeV" << std::endl; 
                     break;
                     }
              case 12:
                     {
//C
//C       ' (12) 21Na(d,n)22Mg '
//C
                     ipart = 80;
                     zbeam = 11;
                     abeam = 21;                       
                     atarg = 2.;
                     zprod = 12;
                     aprod = atarg + abeam;
//C
//C--     create beam particle --> idpart = 80
//C
       	             devmass = -2184.3E-6;
       	             aamass  = abeam*amugev + devmass;
                     tlif    = 0.03;
                     ubuf[0] = fkine[0];
//C
                     DRAGONIon* ion = new DRAGONIon(
                                                 "Na21",                  
                                                 aamass*GeV,      
                                                 0.0*keV,                    
                                                 ubuf[0]*eplus,          
                                                 3, +1, 0, 0, 0, 0,       
                                                 "nucleus",              
                                                 0, 0, ipart,               
                                                 false,                  
                                                 tlif*s,                 
                                                 nullptr,                
                                                 false,
                                                 zbeam,
                                                 abeam                                             
                                                 );                     

//C
//C--     nonresonant level populated --> idpart = 81
//C
                     ipart    = 81;
                     tlif     = 1e-20;  //!10 KeV width
                     resenerg = 0.00;
                     reswidth = hbar/tlif;
                     aamass = sqrt(std::pow(aamass+deutmass,2.0) + 2.*deutmass*beamenerg*.001);
                     ubuf[0]= fkine[1];
//C
                     ion = new DRAGONIon(
                           "nonres_Mg23_",                  
                           aamass*GeV,      
                           0.0*keV,                    
                           ubuf[0]*eplus,          
                           3, +1, 0, 0, 0, 0,       
                           "nucleus",              
                           0, 0, ipart,               
                           false,                  
                           tlif*s,                 
                           nullptr,                
                           false,
                           zbeam+1,
                           abeam                                      
                           );
 
//C
//C--     Define states in Mg22  from neutron decays from resonance
//C
                     ipart   = 82;
                     devmass = -396.8E-6;
                     prodm   = (aprod-1)*amugev + devmass; //!gs of 22Mg
//C

                     elevel  = 4.401;

                     aamass  = prodm + elevel/1000.;
                     tlif    = 3.E-11;
//C
                     ion = new DRAGONIon(
                           "Mg22_2+_2",                  
                           aamass*GeV,      
                           0.0*keV,                    
                           ubuf[0]*eplus,          
                           4, +1, 0, 0, 0, 0,       
                           "nucleus",              
                           0, 0, ipart,               
                           false,                  
                           tlif*s,                 
                           nullptr,                
                           false,
                           zbeam+1,
                           abeam+2                                        
                           );                     

//C
                     ipart   = 83;
                     elevel  = 1.246;

                     aamass  = prodm + elevel/1000.;
                     tlif    = 3.E-11;
//C
                     ion = new DRAGONIon(
                           "Mg22_2+_1",                  
                           aamass*GeV,      
                           0.0*keV,                    
                           ubuf[0]*eplus,          
                           4, +1, 0, 0, 0, 0,       
                           "nucleus",              
                           0, 0, ipart,               
                           false,                  
                           tlif*s,                 
                           nullptr,                
                           false,
                           zbeam+1,
                           abeam                                           
                           );                     

//C
//C--     ground state --> idpart = 84;
//C
                     ipart  = 84;
                     irecoil= 84;
                     tlif   = 1000.;
//C
                     ion = new DRAGONIon(
                           "Mg22_0+",                  
                           prodm*GeV,      
                           0.0*keV,                    
                           ubuf[0]*eplus,          
                           0, +1, 0, 0, 0, 0,       
                           "nucleus",              
                           0, 0, ipart,               
                           false,                  
                           tlif*s,                 
                           nullptr,                
                           false,
                           zbeam+1,
                           abeam                                           
                           );                     
                    
//C
                     vzero(brat,10);
                     vzero(mode,10);
//C
//C--     branch info -- resonance decays
//C
                     brat[0] = 100.;
                     mode[0] = 82;
                                       
                     G4DecayTable* decayTable1 = new G4DecayTable();
                     
                     for (int i = 0; i< 10; i++)
                         { 
                         if(mode[i])
                         {
                          G4VDecayChannel* decayChannel = new G4PhaseSpaceDecayChannel(
                                                     G4ParticleTable::GetParticleTable()->FindParticle(81)->GetParticleName(),     
                                                     brat[i]/100.,                        
                                                     1,                          
                                                     G4ParticleTable::GetParticleTable()->FindParticle(mode[i])->GetParticleName()                
                                                     );
                     
                          decayTable1->Insert(decayChannel);
                          }
                          }
                     DRAGONIon* part1 = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(81));
                     part1->SetDecayTable(decayTable1);
                     decayTable1->DumpInfo();                  
                   
//C
                     vzero(brat,10);
                     vzero(mode,10);
//C
                     brat[0] = 87.;
                     mode[0] = 83;
                     brat[1] = 13.;
                     mode[1] = 84;

                     G4DecayTable* decayTable2 = new G4DecayTable();
                     
                     for (int i = 0; i< 10; i++)
                         { 
                         if(mode[i])
                         {
                          G4VDecayChannel* decayChannel = new G4PhaseSpaceDecayChannel(
                                                     G4ParticleTable::GetParticleTable()->FindParticle(82)->GetParticleName(),     
                                                     brat[i]/100.,                        
                                                     1,                          
                                                     G4ParticleTable::GetParticleTable()->FindParticle(mode[i])->GetParticleName()                
                                                     );
                     
                          decayTable2->Insert(decayChannel);
                          }
                          }
                     DRAGONIon* part2 = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(82));
                     part2->SetDecayTable(decayTable2); 
                     decayTable2->DumpInfo();                  
//C
                     vzero(brat,10);
                     vzero(mode,10);
//C
                     brat[0] = 100.;
                     mode[0] = 84;
                     
                     G4DecayTable* decayTable3 = new G4DecayTable();
                     
                     for (int i = 0; i< 10; i++)
                         { 
                         if(mode[i])
                         {
                          G4VDecayChannel* decayChannel = new G4PhaseSpaceDecayChannel(
                                                     G4ParticleTable::GetParticleTable()->FindParticle(83)->GetParticleName(),     
                                                     brat[i]/100.,                        
                                                     1,                          
                                                     G4ParticleTable::GetParticleTable()->FindParticle(mode[i])->GetParticleName()                
                                                     );
                     
                          decayTable3->Insert(decayChannel);
                          }
                          }
                     DRAGONIon* part3 = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(83));
                     part3->SetDecayTable(decayTable3);  
                     decayTable3->DumpInfo();                   
//C
                     std::cout << "|**** 21Na(d,ng)22Mg reaction ****|" << std::endl;
                     std::cout << " 100% to gs  + neutron" << std::endl;
                     break;
                     }
              case 13:
                     {
//C
//C       ' (13) 23Na(d,n)24Mg 2+ '
//C
                     ipart = 80;
                     zbeam = 11;
                     abeam = 23;                       
                     atarg = 2.;
                     zprod = 12;
                     aprod = atarg + abeam;
//C
//C--     create beam particle --> idpart = 80
//C
       	             devmass = -9530.0E-6;
       	             aamass  = abeam*amugev + devmass;
                     tlif    = 1000.;
                     ubuf[0] = fkine[0];
//C
                     DRAGONIon* ion = new DRAGONIon(
                                                 "Na23",                  
                                                 aamass*GeV,      
                                                 0.0*keV,                    
                                                 ubuf[0]*eplus,          
                                                 3, +1, 0, 0, 0, 0,       
                                                 "nucleus",              
                                                 0, 0, ipart,               
                                                 false,                  
                                                 tlif*s,                 
                                                 nullptr,                
                                                 false,
                                                 zbeam,
                                                 abeam                                          
                                                 );                     
 
//C
//C--     nonresonant level populated --> idpart = 81
//C
                     ipart    = 81;
                     tlif     = 1e-20;  //!10 KeV width
                     resenerg = 0.00;
                     reswidth = hbar/tlif;
                     aamass = std::sqrt(std::pow(aamass+deutmass,2.0) + 2.*deutmass*beamenerg*.001);
                     ubuf[0]= fkine[1];
//C
                     ion = new DRAGONIon(
                           "nonres_Mg25_",                  
                           aamass*GeV,      
                           0.0*keV,                    
                           ubuf[0]*eplus,          
                           5, +1, 0, 0, 0, 0,       
                           "nucleus",              
                           0, 0, ipart,               
                           false,                  
                           tlif*s,                 
                           nullptr,                
                           false,
                           zbeam+1,
                           abeam                                           
                           );                     

//C
//C--     Define states in Mg24  from neutron decays from resonance
//C
                     ipart   = 82;
                     devmass = -13930.7E-6;
                     prodm   = (aprod-1)*amugev + devmass; //!gs of 24Mg
//C

                     elevel  = 1.368;

                     aamass1  = prodm + elevel/1000.;
                     if (aamass < aamass1)
                     std::cout << "Beam energy below 2+ threshold" << std::endl;
                     tlif    = 3.E-11;
//C
                     ion = new DRAGONIon(
                           "Mg24_2+_2",                  
                           aamass1*GeV,      
                           0.0*keV,                    
                           ubuf[0]*eplus,          
                           4, +1, 0, 0, 0, 0,       
                           "nucleus",              
                           0, 0, ipart,               
                           false,                  
                           tlif*s,                 
                           nullptr,                
                           false,
                           zbeam+1,
                           abeam                                         
                           );                     
 
//C
//C--     ground state --> idpart = 93
//C
                     ipart  = 83;
                     irecoil= 83;
                     tlif   = 1000.;
//C
                     ion = new DRAGONIon(
                           "Mg24_0+",                  
                           prodm*GeV,      
                           0.0*keV,                    
                           ubuf[0]*eplus,          
                           0, +1, 0, 0, 0, 0,       
                           "nucleus",              
                           0, 0, ipart,               
                           false,                  
                           tlif*s,                 
                           nullptr,                
                           false,
                           zbeam+1,
                           abeam                                          
                           );  

//C
                     vzero(brat,10);
                     vzero(mode,10);
//C
//C--     branch info -- resonance decays
//C
                     brat[0] = 100.;
                     mode[0] = 82;
    
                     G4DecayTable* decayTable1 = new G4DecayTable();
                     
                     for (int i = 0; i< 10; i++)
                         { 
                         if(mode[i])
                         {
                          G4VDecayChannel* decayChannel = new G4PhaseSpaceDecayChannel(
                                                     G4ParticleTable::GetParticleTable()->FindParticle(81)->GetParticleName(),     
                                                     brat[i]/100.,                        
                                                     1,                          
                                                     G4ParticleTable::GetParticleTable()->FindParticle(mode[i])->GetParticleName()                
                                                     );
                     
                          decayTable1->Insert(decayChannel);
                          }
                          }
                     DRAGONIon* part1 = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(81));
                     part1->SetDecayTable(decayTable1); 
                     decayTable1->DumpInfo();
//C
                     vzero(brat,10);
                     vzero(mode,10);
//C
                     brat[0] = 100.;
                     mode[0] = 83;
                     
                     G4DecayTable* decayTable2 = new G4DecayTable();
                     
                     for (int i = 0; i< 10; i++)
                         { 
                         if(mode[i])
                         {
                          G4VDecayChannel* decayChannel = new G4PhaseSpaceDecayChannel(
                                                     G4ParticleTable::GetParticleTable()->FindParticle(82)->GetParticleName(),     
                                                     brat[i]/100.,                        
                                                     1,                          
                                                     G4ParticleTable::GetParticleTable()->FindParticle(mode[i])->GetParticleName()                
                                                     );
                     
                          decayTable2->Insert(decayChannel);
                          }
                          }
                     DRAGONIon* part2 = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(82));
                     part2->SetDecayTable(decayTable2); 
                     decayTable2->DumpInfo();
//C
                     std::cout << "|**** 23Na(d,ng)24Mg ****|" << std::endl;  
                     std::cout << "Resonance energy " << resenerg << " MeV " << ", " << "width " << reswidth << " MeV " << ", " << "level " << elevel << " MeV" << std::endl; 
                     break; 
                     }
              case 14:
                     {
//C
//C       ' (14) 23Na(p,g)24Mg '
//C
//C
                      ipart = 80;
                      zbeam = 11;
                      abeam = 23;                       
                      atarg = 1.;
                      zprod = 12;
                      aprod = atarg + abeam;
//C
//C--     create beam particle --> idpart = 90
//C
       	              devmass = -9530.0E-6;
         	      aamass  = abeam*amugev + devmass;
                      tlif    = 0.03;
                      ubuf[0] = fkine[0];
//C
                      DRAGONIon* ion = new DRAGONIon(
                                                 "Na23",                  
                                                 aamass*GeV,      
                                                 0.0*keV,                    
                                                 ubuf[0]*eplus,          
                                                 3, +1, 0, 0, 0, 0,       
                                                 "nucleus",              
                                                 0, 0, ipart,               
                                                 false,                  
                                                 tlif*s,                 
                                                 nullptr,                
                                                 false,
                                                 zbeam,
                                                 abeam                                            
                                                 );  

//C
//C--     resonant level populated --> idpart = 91
//C
                      ipart = 81;
                      tlif     = 4.E-14;
                      elevel = 1.368;
                      devmass = -13930.7E-6;
                      prodm   = aprod*amugev + devmass;
                      resenerg = (prodm + elevel/1000. -aamass -hmass)*1000.;
                      if (resenerg <= 0.)
                         {
                          resenerg = 0.0;
                          elevel = (aamass + hmass -prodm)*1000.;
                          std::cout << " resonance energy negative - set to 0!" << std::endl;
                          }
                      aamass = aamass+ hmass + resenerg/1000.;
                      reswidth = hbar/tlif;
                      ubuf[0]  = fkine[1];
// C
                      ion = new DRAGONIon(
                            "res_Mg24_2+",                  
                            aamass*GeV,      
                            0.0*keV,                    
                            ubuf[0]*eplus,          
                            4, +1, 0, 0, 0, 0,       
                            "nucleus",              
                            0, 0, ipart,               
                            false,                  
                            tlif*s,                 
                            nullptr,                
                            false,
                            zbeam+1,
                            abeam+1                                          
                            );  

//C
//C--     Define states in Mg24 gamma decays from resonance
//C
//C
//C
//C--     ground state --> idpart = 82
//C
                      ipart = 82;
                      irecoil=82;
                      tlif   = 1000.;
//C
                      ion = new DRAGONIon(
                            "Mg24_0+",                  
                            aamass*GeV,      
                            0.0*keV,                    
                            ubuf[0]*eplus,          
                            0, +1, 0, 0, 0, 0,       
                            "nucleus",              
                            0, 0, ipart,               
                            false,                  
                            tlif*s,                 
                            nullptr,                
                            false,
                            zbeam+1,
                            abeam+1                                           
                            );  

//C
                     vzero(brat,10);
                     vzero(mode,10);
//C
//C--     branch info -- resonance decays
//C                      
                     brat[0] = 100.;
                     mode[0] = 82;
    
                     G4DecayTable* decayTable = new G4DecayTable();
                     
                     for (int i = 0; i< 10; i++)
                         { 
                         if(mode[i])
                         {
                          G4VDecayChannel* decayChannel = new G4PhaseSpaceDecayChannel(
                                                     G4ParticleTable::GetParticleTable()->FindParticle(81)->GetParticleName(),     
                                                     brat[i]/100.,                        
                                                     1,                          
                                                     G4ParticleTable::GetParticleTable()->FindParticle(mode[i])->GetParticleName()                
                                                     );
                     
                          decayTable->Insert(decayChannel);
                          }
                          }
                     DRAGONIon* part = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(81));
                     part->SetDecayTable(decayTable); 
                     decayTable->DumpInfo();                     
//C
                     std::cout << "|**** 23Na(p,g)24Mg reaction ****|" << std::endl;  
                     std::cout << "Resonance energy " << resenerg << " MeV " << ", " << "width " << reswidth << " MeV " << ", " << "level " << elevel << " MeV" << std::endl;                      
                     break;
                     }
              case 18:
                     {
//C
//C       ' (18) 21Na(p,g)22Mg(822) '
//C
                     ipart = 80;
                     zbeam = 11;
                     abeam = 21;                       
                     atarg = 1.;
                     zprod = 12;
                     aprod = atarg + abeam;
//C
//C--     create beam particle --> idpart = 80
//C
       	             devmass = -2184.3E-6;
                     aamass  = abeam*amugev + devmass;
                     tlif    = 0.03;
                     ubuf[0] = fkine[0];
//C
                     DRAGONIon* ion = new DRAGONIon(
                                                 "Na21",                  
                                                 aamass*GeV,      
                                                 0.0*keV,                    
                                                 ubuf[0]*eplus,          
                                                 3, +1, 0, 0, 0, 0,       
                                                 "nucleus",              
                                                 0, 0, ipart,               
                                                 false,                  
                                                 tlif*s,                 
                                                 nullptr,                
                                                 false,
                                                 zbeam,
                                                 abeam                                                
                                                 );  
                  
//C
//C--     resonant level populated --> idpart = 81
//C
                     ipart    = 81;
                     tlif     = 84.0e-19;
                     resenerg = 0.8225;
                     reswidth = hbar/tlif/2; //!This is gamma/2 =HWHM for BW
                     aamass   = aamass + resenerg/1000. + hmass;
                     ubuf[0]  = fkine[1];
//C
                     ion = new DRAGONIon(
                           "res_Mg22_822",                  
                           aamass*GeV,      
                           0.0*keV,                    
                           ubuf[0]*eplus,          
                           0, +1, 0, 0, 0, 0,       
                           "nucleus",              
                           0, 0, ipart,               
                           false,                  
                           tlif*s,                 
                           nullptr,                
                           false,
                           zbeam+1,
                           abeam+1                                            
                           );  

//C
//C--     Define states in Mg22 gamma decays from resonance
//C
                     ipart   = 82;
                     devmass = -396.8E-6;
                     prodm   = aprod*amugev + devmass;
//C
                     elevel = (aamass -prodm)*1000.;
                     std::cout << "|**** 21Na(p,g)22Mg(822) reaction ****|" << std::endl;  
                     std::cout << "Resonance energy " << resenerg << " MeV " << ", " << "width " << reswidth << " MeV " << ", " << "level " << elevel << " MeV" << std::endl; 
//C
                     elevel  = 4.401;
                     aamass  = prodm + elevel/1000.;
                     tlif    = 3.E-14;
//C
                     ion = new DRAGONIon(
                           "Mg22_4401+",                  
                           aamass*GeV,      
                           0.0*keV,                    
                           ubuf[0]*eplus,          
                           0, +1, 0, 0, 0, 0,       
                           "nucleus",              
                           0, 0, ipart,               
                           false,                  
                           tlif*s,                 
                           nullptr,                
                           false,
                           zbeam+1,
                           abeam+1                
                           );  

//C
//C
                     ipart   = 83;
                     elevel  = 3.308;
                     aamass  = prodm + elevel/1000.;
                     tlif    = 3.E-14;
//C
                     ion = new DRAGONIon(
                           "Mg22_3308",                  
                           aamass*GeV,      
                           0.0*keV,                    
                           ubuf[0]*eplus,          
                           0, +1, 0, 0, 0, 0,       
                           "nucleus",              
                           0, 0, ipart,               
                           false,                  
                           tlif*s,                 
                           nullptr,                
                           false,
                           zbeam+1,
                           abeam+1                                          
                           );  

//C
//C
                     ipart   = 84;
                     elevel  = 1.246;
                     aamass  = prodm + elevel/1000.;
                     tlif    = 3.E-11;
//C
                     ion = new DRAGONIon(
                           "Mg22_2+",                  
                           aamass*GeV,      
                           0.0*keV,                    
                           ubuf[0]*eplus,          
                           4, +1, 0, 0, 0, 0,       
                           "nucleus",              
                           0, 0, ipart,               
                           false,                  
                           tlif*s,                 
                           nullptr,                
                           false,
                           zbeam+1,
                           abeam+1                                          
                           );  

//C
//C--     ground state --> idpart = 85
//C
                     ipart  = 85;
                     irecoil= 85;
                     tlif   = 1000.;
//C
                     ion = new DRAGONIon(
                           "Mg22_0+",                  
                           aamass*GeV,      
                           0.0*keV,                    
                           ubuf[0]*eplus,          
                           0, +1, 0, 0, 0, 0,       
                           "nucleus",              
                           0, 0, ipart,               
                           false,                  
                           tlif*s,                 
                           nullptr,                
                           false,
                           zbeam+1,
                           abeam+1                                           
                           );  

//C
                     vzero(brat,10);
                     vzero(mode,10);
//C
//C--     branch info -- resonance decays
//C
                     brat[0] = 50.;
                     mode[0] = 84;
                     brat[0] = 50.;
                     mode[0] = 85;
    
                     G4DecayTable* decayTable1 = new G4DecayTable();
                     
                     for (int i = 0; i< 10; i++)
                         { 
                         if(mode[i])
                         {
                          G4VDecayChannel* decayChannel = new G4PhaseSpaceDecayChannel(
                                                     G4ParticleTable::GetParticleTable()->FindParticle(81)->GetParticleName(),     
                                                     brat[i]/100.,                        
                                                     1,                          
                                                     G4ParticleTable::GetParticleTable()->FindParticle(mode[i])->GetParticleName()                
                                                     );
                     
                          decayTable1->Insert(decayChannel);
                          }
                          }
                     DRAGONIon* part1 = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(81));
                     part1->SetDecayTable(decayTable1); 
                     decayTable1->DumpInfo();
//C
                     vzero(brat,10);
                     vzero(mode,10);
//C
                     brat[0] = 87.;
                     mode[0] = 84;
                     brat[0] = 13.;
                     mode[0] = 85;
                     
                     G4DecayTable* decayTable2 = new G4DecayTable();
                     
                     for (int i = 0; i< 10; i++)
                         { 
                         if(mode[i])
                         {
                          G4VDecayChannel* decayChannel = new G4PhaseSpaceDecayChannel(
                                                     G4ParticleTable::GetParticleTable()->FindParticle(82)->GetParticleName(),     
                                                     brat[i]/100.,                        
                                                     1,                          
                                                     G4ParticleTable::GetParticleTable()->FindParticle(mode[i])->GetParticleName()                
                                                     );
                     
                          decayTable2->Insert(decayChannel);
                          }
                          }
                     DRAGONIon* part2 = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(82));
                     part2->SetDecayTable(decayTable2); 
                     decayTable2->DumpInfo();
//C
                     vzero(brat,10);
                     vzero(mode,10);
//C
                     brat[0] = 100.;
                     mode[0] = 84;
                     
                     G4DecayTable* decayTable3 = new G4DecayTable();
                     
                     for (int i = 0; i< 10; i++)
                         { 
                         if(mode[i])
                         {
                          G4VDecayChannel* decayChannel = new G4PhaseSpaceDecayChannel(
                                                     G4ParticleTable::GetParticleTable()->FindParticle(83)->GetParticleName(),     
                                                     brat[i]/100.,                        
                                                     1,                          
                                                     G4ParticleTable::GetParticleTable()->FindParticle(mode[i])->GetParticleName()                
                                                     );
                     
                          decayTable3->Insert(decayChannel);
                          }
                          }
                     DRAGONIon* part3 = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(83));
                     part3->SetDecayTable(decayTable3); 
                     decayTable3->DumpInfo();
//C
                     vzero(brat,10);
                     vzero(mode,10);
//C
                     brat[0] = 100.;
                     mode[0] = 85;
                     
                     G4DecayTable* decayTable4 = new G4DecayTable();
                     
                     for (int i = 0; i< 10; i++)
                         { 
                         if(mode[i])
                         {
                          G4VDecayChannel* decayChannel = new G4PhaseSpaceDecayChannel(
                                                     G4ParticleTable::GetParticleTable()->FindParticle(84)->GetParticleName(),     
                                                     brat[i]/100.,                        
                                                     1,                          
                                                     G4ParticleTable::GetParticleTable()->FindParticle(mode[i])->GetParticleName()                
                                                     );
                     
                          decayTable4->Insert(decayChannel);
                          }
                          }
                     DRAGONIon* part4 = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(84));
                     part4->SetDecayTable(decayTable4); 
                     decayTable4->DumpInfo();
                     break;
                     }
              case 19:
                     {
//C      
//C       ' (19) 12C(alpha,g)16O '
//C
//C
                      ipart = 80;
                      zbeam = 6;
                      abeam = 12;                        
                      atarg = 4.;
                      zprod = 8.;
                      aprod = atarg + abeam;
//C
//C--     create beam particle --> idpart = 80
//C
                      devmass = 0.0E-6;
                      aamass  = abeam*amugev + devmass;
                      tlif    = 1000.;
                      ubuf[0] = fkine[0];
//C
                      DRAGONIon* ion = new DRAGONIon(
                                                  "C12gs",                  
                                                  aamass*GeV,      
                                                  0.0*keV,                    
                                                  ubuf[0]*eplus,          
                                                  0, +1, 0, 0, 0, 0,       
                                                  "nucleus",              
                                                  0, 0, ipart,               
                                                  false,                  
                                                  tlif*s,                 
                                                  nullptr,                
                                                  false,
                                                  zbeam,
                                                  abeam                                                                  
                                                  );  

//C
//C--     resonant level populated --> idpart = 81
//C
                      ipart = 81;
                      resenerg = 4.358;
                      reswidth = 7.E-05;
                      aamass   = aamass + resenerg/1000. + hemass;
                      std::cout << aamass << std::endl;
                      tlif     = hbar/reswidth;
                      ubuf[0]  = fkine[1];
//C
                      ion = new DRAGONIon(
                      "res_O16_2+",                  
                      aamass*GeV,      
                      0.0*keV,                    
                      ubuf[0]*eplus,          
                      4, +1, 0, 0, 0, 0,       
                      "nucleus",              
                      0, 0, ipart,               
                      false,                  
                      tlif*s,                 
                      nullptr,                
                      false,
                      zbeam+2,
                      abeam+4                                     
                      );  

//C
                      elevel = 11.520;
                      std::cout << "Resonant mass: " << aamass << std::endl;
                      rmass = aamass;
                      
                      std::cout << "|**** 12C(alpha,gamma)16O reaction ****|" << std::endl;  
                      std::cout << "Resonance energy " << resenerg << " MeV " << ", " << "width " << reswidth << " MeV " << ", " << "level " << elevel << " MeV" << std::endl; 
//C
//C--     Set common block variables for cross-section
//C
                      m1 = abeam*amugev + devmass;
                      m2 = hemass;
                      z1 = zbeam;
                      z2 = 2.;
                      er = resenerg;
                      gp = 70./1000.;
                      gg = 0.007/1000.;
                      omg = 5.;
                      ell = 1.;
                      ires = 81;
//C
//C--     Define states in O16 gamma decays from resonance
//C
                      ipart = 82;
                      devmass = -4736.998E-6;
                      prodm   = aprod*amugev + devmass;
//C
                      elevel  = 11.260;
                      aamass  = prodm + elevel/1000.;
                      tlif    = hbar/2500.E-06;
//C
                      ion = new DRAGONIon(
                            "12_O16_0+",                  
                            aamass*GeV,      
                            0.0*keV,                    
                            ubuf[0]*eplus,          
                            0, +1, 0, 0, 0, 0,       
                            "nucleus",              
                            0, 0, ipart,               
                            false,                  
                            tlif*s,                 
                            nullptr,                
                            false,
                            zbeam+2,
                            abeam+4                                             
                            );  

//C
                      ipart  = 83; 
                      elevel = 11.0967;
                      aamass = prodm + elevel/1000.;
                      tlif   = hbar/0.28E-06;
//C
                      ion = new DRAGONIon(
                            "11_O16_4+",                  
                            aamass*GeV,      
                            0.0*keV,                    
                            ubuf[0]*eplus,          
                            8, +1, 0, 0, 0, 0,       
                            "nucleus",              
                            0, 0, ipart,               
                            false,                  
                            tlif*s,                 
                            nullptr,                
                            false,
                            zbeam+2,
                            abeam+4                                            
                            );  

//C
                      ipart = 84;
                      elevel = 11.080;
                      aamass = prodm + elevel/1000.;
                      tlif   = hbar/12.E-06;
//C
                      ion = new DRAGONIon(
                            "10_O16_3+",                  
                            aamass*GeV,      
                            0.0*keV,                    
                            ubuf[0]*eplus,          
                            6, +1, 0, 0, 0, 0,       
                            "nucleus",              
                            0, 0, ipart,               
                            false,                  
                            tlif*s,                 
                            nullptr,                
                            false,
                            zbeam+2,
                            abeam+4                                            
                            );  

//C
                      ipart  = 85;
                      elevel = 10.957;
                      aamass = prodm + elevel/1000.;
                      tlif = 5.5E-15;
//C
                      ion = new DRAGONIon(
                            "9_O16_0-",                  
                            aamass*GeV,      
                            0.0*keV,                    
                            ubuf[0]*eplus,          
                            0, -1, 0, 0, 0, 0,       
                            "nucleus",              
                            0, 0, ipart,               
                            false,                  
                            tlif*s,                 
                            nullptr,                
                            false,
                            zbeam+2,
                            abeam+4                                           
                            );  

//C
                      ipart  = 86;
                      elevel = 10.356;
                      aamass = prodm + elevel/1000.;
                      tlif = hbar/26.E-06;
//C
                      ion = new DRAGONIon(
                            "8_O16_4+",                  
                            aamass*GeV,      
                            0.0*keV,                    
                            ubuf[0]*eplus,          
                            8, +1, 0, 0, 0, 0,       
                            "nucleus",              
                            0, 0, ipart,               
                            false,                  
                            tlif*s,                 
                            nullptr,                
                            false,
                            zbeam+2,
                            abeam+4                 
                            );  

//C
                      ipart  = 87;
                      elevel = 9.8445;
                      aamass = prodm + elevel/1000.;
                      tlif = hbar/0.62E-06;
//C
                      ion = new DRAGONIon(
                            "7_O16_2+",                  
                            aamass*GeV,      
                            0.0*keV,                    
                            ubuf[0]*eplus,          
                            4, +1, 0, 0, 0, 0,       
                            "nucleus",              
                            0, 0, ipart,               
                            false,                  
                            tlif*s,                 
                            nullptr,                
                            false,
                            zbeam+2,
                            abeam+4                                           
                            );  

//C
                      ipart  = 88;
                      elevel = 9.585;
                      aamass = prodm + elevel/1000.;
                      tlif = hbar/420.E-06;
//C
                      ion = new DRAGONIon(
                            "6_O16_1-",                  
                            aamass*GeV,      
                            0.0*keV,                    
                            ubuf[0]*eplus,          
                            2, -1, 0, 0, 0, 0,       
                            "nucleus",              
                            0, 0, ipart,               
                            false,                  
                            tlif*s,                 
                            nullptr,                
                            false,
                            zbeam+2,
                            abeam+4                                             
                            );  

//C
                      ipart  = 89; 
                      elevel = 8.8719;
                      aamass = prodm + elevel/1000.;
                      tlif = 125.E-15;
//C
                      ion = new DRAGONIon(
                            "5_O16_2-",                  
                            aamass*GeV,      
                            0.0*keV,                    
                            ubuf[0]*eplus,          
                            4, -1, 0, 0, 0, 0,       
                            "nucleus",              
                            0, 0, ipart,               
                            false,                  
                            tlif*s,                 
                            nullptr,                
                            false,
                            zbeam+2,
                            abeam+4                                        
                            );  

//C
                      ipart  = 90;
                      elevel = 7.11685;
                      aamass = prodm + elevel/1000.;
                      tlif = 8.3E-15;
//C
                      ion = new DRAGONIon(
                            "4_O16_1-",                  
                            aamass*GeV,      
                            0.0*keV,                    
                            ubuf[0]*eplus,          
                            0, +1, 0, 0, 0, 0,       
                            "nucleus",              
                            0, 0, ipart,               
                            false,                  
                            tlif*s,                 
                            nullptr,                
                            false,
                            zbeam+2,
                            abeam+4                                           
                            );  

//C
                      ipart  = 91;
                      elevel = 6.9171;
                      aamass = prodm + elevel/1000.;
                      tlif = 4.7E-15;
//C
                      ion = new DRAGONIon(
                            "3_O16_2+",                  
                            aamass*GeV,      
                            0.0*keV,                    
                            ubuf[0]*eplus,          
                            4, +1, 0, 0, 0, 0,       
                            "nucleus",              
                            0, 0, ipart,               
                            false,                  
                            tlif*s,                 
                            nullptr,                
                            false,
                            zbeam+2,
                            abeam+4                                            
                            );  

//C
                      ipart  = 92;
                      elevel = 6.12989;
                      aamass = prodm + elevel/1000.;
                      tlif = 18.4E-12;
//C
                      ion = new DRAGONIon(
                            "2_O16_3-",                  
                            aamass*GeV,      
                            0.0*keV,                    
                            ubuf[0]*eplus,          
                            6, -1, 0, 0, 0, 0,       
                            "nucleus",              
                            0, 0, ipart,               
                            false,                  
                            tlif*s,                 
                            nullptr,                
                            false,
                            zbeam+2,
                            abeam+4                                           
                            );  

//C
                      ipart = 93;
                      elevel = 6.0494;
                      aamass = prodm + elevel/1000.;
                      tlif = 67.E-12;
//C
                      ion = new DRAGONIon(
                            "1_O16_0+",                  
                            aamass*GeV,      
                            0.0*keV,                    
                            ubuf[0]*eplus,          
                            0, +1, 0, 0, 0, 0,       
                            "nucleus",              
                            0, 0, ipart,               
                            false,                  
                            tlif*s,                 
                            nullptr,                
                            false,
                            zbeam+2,
                            abeam+4              
                            );  

//C
//C--     ground state --> idpart = 104
                      ipart   = 94;
                      irecoil = 94;
//C
                      tlif = 1000.;
//C
                      ion = new DRAGONIon(
                            "gs_O16_0+",                  
                            aamass*GeV,      
                            0.0*keV,                    
                            ubuf[0]*eplus,          
                            0, +1, 0, 0, 0, 0,       
                            "nucleus",              
                            0, 0, ipart,               
                            false,                  
                            tlif*s,                 
                            nullptr,                
                            false,
                            zbeam+2,
                            abeam+4               
                            );  

                      std::cout << "Ground state mass: " << prodm << std::endl;
//C
                     vzero(brat,10);
                     vzero(mode,10);
//C
//C--     branch info -- resonance decays
//C
                     brat[0] = 90.99;
                     mode[0] = 94;
                     brat[1] = 4.19;
                     mode[1] = 93;
                     brat[2] = 4.;
                     mode[2] = 91;
                     brat[3] = 0.82;
                     mode[3] = 90;
    
                     G4DecayTable* decayTable1 = new G4DecayTable();
                     
                     for (int i = 0; i< 10; i++)
                         { 
                         if(mode[i])
                         {
                          G4VDecayChannel* decayChannel = new G4PhaseSpaceDecayChannel(
                                                     G4ParticleTable::GetParticleTable()->FindParticle(81)->GetParticleName(),     
                                                     brat[i]/100.,                        
                                                     1,                          
                                                     G4ParticleTable::GetParticleTable()->FindParticle(mode[i])->GetParticleName()                
                                                     );
                     
                          decayTable1->Insert(decayChannel);
                          }
                          }
                     DRAGONIon* part1 = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(81));
                     part1->SetDecayTable(decayTable1); 
                     decayTable1->DumpInfo();
//C
//C--     12th ex. state
//C
//C
                     vzero(brat,10);
                     vzero(mode,10);                     
//C
                     G4DecayTable* decayTable2 = new G4DecayTable();

                     for (int i = 0; i < decayTable1->entries(); ++i) {
                       G4VDecayChannel* channel = decayTable1->GetDecayChannel(i);
                       G4VDecayChannel* clonedChannel = new G4PhaseSpaceDecayChannel(*dynamic_cast<G4PhaseSpaceDecayChannel*>(channel));
                       decayTable2->Insert(clonedChannel);
                       }
          
                     DRAGONIon* part2 = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(82));
                     part2->SetDecayTable(decayTable2); 
                     decayTable2->DumpInfo();
//C
//C--     11th ex. state
//C
                     brat[0] = 55.25;
                     mode[0] = 92;
                     brat[1] = 44.75;
                     mode[1] = 91;
                     
                     G4DecayTable* decayTable3 = new G4DecayTable();
                     
                     for (int i = 0; i< 10; i++)
                         { 
                         if(mode[i])
                         {
                          G4VDecayChannel* decayChannel = new G4PhaseSpaceDecayChannel(
                                                     G4ParticleTable::GetParticleTable()->FindParticle(83)->GetParticleName(),     
                                                     brat[i]/100.,                        
                                                     1,                          
                                                     G4ParticleTable::GetParticleTable()->FindParticle(mode[i])->GetParticleName()                
                                                     );
                     
                          decayTable3->Insert(decayChannel);
                          }
                          }
                     DRAGONIon* part3 = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(83));
                     part3->SetDecayTable(decayTable3); 
                     decayTable3->DumpInfo();
//C
//C--     10th ex. state
//C
                     vzero(brat,10);
                     vzero(mode,10);
//C
                     G4DecayTable* decayTable4 = new G4DecayTable();

                     for (int i = 0; i < decayTable3->entries(); ++i) {
                       G4VDecayChannel* channel = decayTable3->GetDecayChannel(i);
                       G4VDecayChannel* clonedChannel = new G4PhaseSpaceDecayChannel(*dynamic_cast<G4PhaseSpaceDecayChannel*>(channel));
                       decayTable4->Insert(clonedChannel);
                       }
          
                     DRAGONIon* part4 = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(82));
                     part4->SetDecayTable(decayTable4);
                     decayTable4->DumpInfo();
//C
//C--     9th ex. state
//C
                     brat[0] = 100.;
                     mode[0] = 90;
                     
                     G4DecayTable* decayTable5 = new G4DecayTable();
                     
                     for (int i = 0; i< 10; i++)
                         { 
                         if(mode[i])
                         {
                          G4VDecayChannel* decayChannel = new G4PhaseSpaceDecayChannel(
                                                     G4ParticleTable::GetParticleTable()->FindParticle(85)->GetParticleName(),     
                                                     brat[i]/100.,                        
                                                     1,                          
                                                     G4ParticleTable::GetParticleTable()->FindParticle(mode[i])->GetParticleName()                
                                                     );
                     
                          decayTable5->Insert(decayChannel);
                          }
                          }
                     DRAGONIon* part5 = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(85));
                     part5->SetDecayTable(decayTable5); 
                     decayTable5->DumpInfo();
//C
                     vzero(brat,10);
                     vzero(mode,10);
//C
//C--     8th ex. state
//C
                     brat[0] = 1.57;
                     mode[0] = 92;
                     brat[1] = 98.43;
                     mode[1] = 91;
                     
                     G4DecayTable* decayTable6 = new G4DecayTable();
                     
                     for (int i = 0; i< 10; i++)
                         { 
                         if(mode[i])
                         {
                          G4VDecayChannel* decayChannel = new G4PhaseSpaceDecayChannel(
                                                     G4ParticleTable::GetParticleTable()->FindParticle(86)->GetParticleName(),     
                                                     brat[i]/100.,                        
                                                     1,                          
                                                     G4ParticleTable::GetParticleTable()->FindParticle(mode[i])->GetParticleName()                
                                                     );
                     
                          decayTable6->Insert(decayChannel);
                          }
                          }
                     DRAGONIon* part6 = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(86));
                     part6->SetDecayTable(decayTable6); 
                     decayTable6->DumpInfo();
//C
                     vzero(brat,10);
                     vzero(mode,10);
//C
//C--     7th ex. state
//C
                     brat[0] = 60.98;
                     mode[0] = 94;
                     brat[1] = 18.29;
                     mode[1] = 93;
                     brat[2] = 20.73;
                     mode[2] = 91;
                     
                     G4DecayTable* decayTable7 = new G4DecayTable();
                     
                     for (int i = 0; i< 10; i++)
                         { 
                         if(mode[i])
                         {
                          G4VDecayChannel* decayChannel = new G4PhaseSpaceDecayChannel(
                                                     G4ParticleTable::GetParticleTable()->FindParticle(87)->GetParticleName(),     
                                                     brat[i]/100.,                        
                                                     1,                          
                                                     G4ParticleTable::GetParticleTable()->FindParticle(mode[i])->GetParticleName()                
                                                     );
                     
                          decayTable7->Insert(decayChannel);
                          }
                          }
                     DRAGONIon* part7 = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(87));
                     part7->SetDecayTable(decayTable7); 
                     decayTable7->DumpInfo();
//C
                     vzero(brat,10);
                     vzero(mode,10);
//C
//C--     6th ex. state
//C
                     brat[0] = 89.29;
                     mode[0] = 94;
                     brat[1] = 10.71;
                     mode[1] = 91;
                     
                     G4DecayTable* decayTable8 = new G4DecayTable();
                     
                     for (int i = 0; i< 10; i++)
                         { 
                         if(mode[i])
                         {
                          G4VDecayChannel* decayChannel = new G4PhaseSpaceDecayChannel(
                                                     G4ParticleTable::GetParticleTable()->FindParticle(88)->GetParticleName(),     
                                                     brat[i]/100.,                        
                                                     1,                          
                                                     G4ParticleTable::GetParticleTable()->FindParticle(mode[i])->GetParticleName()                
                                                     );
                     
                          decayTable8->Insert(decayChannel);
                          }
                          }
                     DRAGONIon* part8 = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(88));
                     part8->SetDecayTable(decayTable8); 
                     decayTable8->DumpInfo();
//C
                     vzero(brat,10);
                     vzero(mode,10);
//C
//C--     5th ex. state
//C
                     brat[0] = 7.22;
                     mode[0] = 94;
                     brat[1] = 0.12;
                     mode[1] = 93;
                     brat[0] = 3.57;
                     mode[0] = 91;
                     brat[1] = 11.42;
                     mode[1] = 90;
                     
                     G4DecayTable* decayTable9 = new G4DecayTable();
                     
                     for (int i = 0; i< 10; i++)
                         { 
                         if(mode[i])
                         {
                          G4VDecayChannel* decayChannel = new G4PhaseSpaceDecayChannel(
                                                     G4ParticleTable::GetParticleTable()->FindParticle(89)->GetParticleName(),     
                                                     brat[i]/100.,                        
                                                     1,                          
                                                     G4ParticleTable::GetParticleTable()->FindParticle(mode[i])->GetParticleName()                
                                                     );
                     
                          decayTable9->Insert(decayChannel);
                          }
                          }
                     DRAGONIon* part9 = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(89));
                     part9->SetDecayTable(decayTable9); 
                     decayTable9->DumpInfo();
//C
                     vzero(brat,10);
                     vzero(mode,10);
//C
//C--     4th ex. state
//C
                     brat[0] = 99.93;
                     mode[0] = 94;
                     brat[1] = 0.07;
                     mode[1] = 92;
                     
                     G4DecayTable* decayTable10 = new G4DecayTable();
                     
                     for (int i = 0; i< 10; i++)
                         { 
                         if(mode[i])
                         {
                          G4VDecayChannel* decayChannel = new G4PhaseSpaceDecayChannel(
                                                     G4ParticleTable::GetParticleTable()->FindParticle(90)->GetParticleName(),     
                                                     brat[i]/100.,                        
                                                     1,                          
                                                     G4ParticleTable::GetParticleTable()->FindParticle(mode[i])->GetParticleName()                
                                                     );
                     
                          decayTable10->Insert(decayChannel);
                          }
                          }
                     DRAGONIon* part10 = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(90));
                     part10->SetDecayTable(decayTable10); 
                     decayTable10->DumpInfo();
//C
                     vzero(brat,10);
                     vzero(mode,10);
//C
//C--     3rd ex. state
//C
                     brat[0] = 99.97;
                     mode[0] = 94;
                     brat[1] = 0.027;
                     mode[1] = 93;
                     brat[2] = 0.008;
                     mode[2] = 92;
                     
                     G4DecayTable* decayTable11 = new G4DecayTable();
                     
                     for (int i = 0; i< 10; i++)
                         { 
                         if(mode[i])
                         {
                          G4VDecayChannel* decayChannel = new G4PhaseSpaceDecayChannel(
                                                     G4ParticleTable::GetParticleTable()->FindParticle(91)->GetParticleName(),     
                                                     brat[i]/100.,                        
                                                     1,                          
                                                     G4ParticleTable::GetParticleTable()->FindParticle(mode[i])->GetParticleName()                
                                                     );
                     
                          decayTable11->Insert(decayChannel);
                          }
                          }
                     DRAGONIon* part11 = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(91));
                     part11->SetDecayTable(decayTable11); 
                     decayTable11->DumpInfo();
//C
                     vzero(brat,10);
                     vzero(mode,10);
//C
//C--     2nd ex. state
//C
                     brat[0] = 100.;
                     mode[0] = 94;
                     
                     G4DecayTable* decayTable12 = new G4DecayTable();
                     
                     for (int i = 0; i< 10; i++)
                         { 
                         if(mode[i])
                         {
                          G4VDecayChannel* decayChannel = new G4PhaseSpaceDecayChannel(
                                                     G4ParticleTable::GetParticleTable()->FindParticle(92)->GetParticleName(),     
                                                     brat[i]/100.,                        
                                                     1,                          
                                                     G4ParticleTable::GetParticleTable()->FindParticle(mode[i])->GetParticleName()                
                                                     );
                     
                          decayTable12->Insert(decayChannel);
                          }
                          }
                     DRAGONIon* part12 = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(92));
                     part12->SetDecayTable(decayTable12); 
                     decayTable12->DumpInfo();
//C
//C--     1st ex. state
//C  
                     brat[0] = 100.;
                     mode[0] = 94;
//C
                     vzero(brat,10);
                     vzero(mode,10);
                     
                     G4DecayTable* decayTable13 = new G4DecayTable();

                   for (int i = 0; i < decayTable12->entries(); ++i) {
                       G4VDecayChannel* channel = decayTable3->GetDecayChannel(i);
                       G4VDecayChannel* clonedChannel = new G4PhaseSpaceDecayChannel(*dynamic_cast<G4PhaseSpaceDecayChannel*>(channel));
                       decayTable12->Insert(clonedChannel);
                       }
//C                     
                     DRAGONIon* part13 = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(93));
                     part13->SetDecayTable(decayTable13); 
                     decayTable13->DumpInfo();
                     break;
                     }

//C
//C--> Case 20: reaction drawn from user input card
//C    only implemented for (a,g) or (p,g) reactions at present
//C    CR: 22.07.2003
//C    Now implemented for (c12,g) reactions.
//C    JS: 15.02.2004      
              case 20:
                     {
                      ReadUserReactionNameList();   
//C
//C    Setup beam particle
//C
                      ipart = 80;
                      aprod = atarg + abeam;
                      devmass = beam_mass_excess;
                      aamass  = abeam*amugev + devmass;
                      tlif    = beamlifetime;
                      ubuf[0] = fkine[0];

     G4cout << "RRRRRRRRRRRRRRRRR" << G4endl;
     G4cout << "atarg: "   << atarg << G4endl;
     G4cout << "abeam: "   << abeam << G4endl;
     G4cout << "aprod: "   << aprod << G4endl;
     G4cout << "devmass: " << devmass << G4endl;
     G4cout << "amugev: "  << amugev << G4endl;
     G4cout << "aamass: "  << aamass << G4endl;
     G4cout << "tlif: "    << tlif << G4endl;
     G4cout << "ubuf[0]: " << ubuf[0] << G4endl;
     G4cout << "RRRRRRRRRRRRRRRRR" << G4endl;
     
//C
                      DRAGONIon* ion = new DRAGONIon(
                                                 beamtyp+"gs",                  
                                                 aamass*GeV,      
                                                 0.0*keV,                    
                                                 ubuf[0]*eplus,          
                                                 0, +1, 0, 0, 0, 0,       
                                                 "nucleus",              
                                                 0, 0, ipart,               
                                                 false,                  
                                                 tlif*s,                 
                                                 nullptr,                
                                                 false,
                                                 zbeam,
                                                 abeam                
                                                 ); 

                      PrintIonProperties(ion);

//C
//C    Create resonant particle
//C
//C    Note: Resonant particle energy based on beam mass, target mass and
//C          variable resenerg.  Not based on energy levels specified in INPUT
//C          file.  (However, energy of states below resonance are based on
//C          eneryg levels specified in INPUT)
//C

                      if(atarg == 4.)
                        { 
                         targtyp = "a";
                         targmass = hemass;
                         aamass   = aamass + resenerg/1000. + hemass;
                         std::cout << aamass << ", " << hemass << ", " << resenerg << std::endl;
                         }
                      else if(atarg == 1.)
                             {
                              targtyp = "p";
                              targmass = hmass;
                              aamass   = aamass + resenerg/1000. + hmass;
                              }
                      else if(atarg == 12.)
                             {
                              targtyp = "12C";
                              targmass = c12mass;
                              aamass   = aamass + resenerg/1000. + c12mass; 
                              }
              
                      tlif    = hbar/((part_width+gam_width)/1000.);
                      ubuf[0] = fkine[1];
//C
                      ipart = 81;
                      ion = new DRAGONIon(
                            "res_"+rectyp+"_"+num[rstate],                  
                            aamass*GeV,      
                            0.0*keV,                    
                            ubuf[0]*eplus,          
                            0, +1, 0, 0, 0, 0,       
                            "nucleus",              
                            0, 0, ipart,               
                            false,                  
                            tlif*s,                 
                            nullptr,                
                            false,
                            ztarg+zbeam,
                            atarg+abeam               
                            );  
                    
                      PrintIonProperties(ion);

                      ires = 81;
                      elevel = level[rstate];
                      std::cout << "Resonant mass: " << aamass << std::endl;
//C                   rmass = aamass;
                      resmass = aamass;
      
                      std::cout << "|**** " << beamtyp << "C" << targtyp << ",gamma)" << rectyp << "reaction ****|" << std::endl;  
                      std::cout << "Resonance energy " << resenerg << " MeV " << ", " << "level " << elevel << " MeV" << std::endl; 

//C
//C    Setup excited states below resonance
//C
                      prodm   = aprod*amugev + recoil_mass_excess;
                      for (int i = 0; i<rstate; i++)
                          {
                           elevel = level[i];
                           aamass = prodm + elevel/1000.;
                           tlif = life[i];
//C           
                           ipart = 81+rstate-i;
                           DRAGONIon* ion = new DRAGONIon(
                                                       rectyp+"_"+num[i],                  
                                                       aamass*GeV,      
                                                       0.0*keV,                    
                                                       ubuf[0]*eplus,          
                                                       0, +1, 0, 0, 0, 0,       
                                                       "nucleus",              
                                                       0, 0, ipart,               
                                                       false,                  
                                                       tlif*s,                 
                                                       nullptr,                
                                                       false,
                                                       ztarg+zbeam,
                                                       atarg+abeam                
                                                       );  
                           PrintIonProperties(ion);
                           } 
                      irecoil = 81+rstate;
//C
//C    Setup decay branching ratios and modes
//C
                      for (int i = 0; i < rstate; i++)
                          {
                           G4DecayTable* decayTable = new G4DecayTable();
                           vzero(brat,10);
                           vzero(mode,10);
                           for (int j = 0; j < 10; j++)
                               {
                                if(br[i][j] != 0)
                                  {
                                   brat[j] = br[i][j];
                                   mode[j] = irecoil-md[i][j];
                                   
                                   //std::cout << "i: " << i << " j: " << j << " "<< br[i][j] << " " << md[i][j] << std::endl;
                                   //std::cout << "i: " << i << " j: " << j << " "<< brat[j] << " " << mode[j] << std::endl;
                             
                                   int a = 81+rstate-i;
                                        
                                   std::cout << "Padre: " << 81+rstate-i << " " << G4ParticleTable::GetParticleTable()->FindParticle(81+rstate-i)->GetParticleName() << std::endl;
                                   //std::cout << "Mode " << mode[j] << std::endl;
   
                                   G4VDecayChannel* decayChannel = new G4PhaseSpaceDecayChannel(
                                                                 G4ParticleTable::GetParticleTable()->FindParticle(81+rstate-i)->GetParticleName(),     
                                                                 brat[j]/100.,                        
                                                                 1,                          
                                                                 G4ParticleTable::GetParticleTable()->FindParticle(mode[j])->GetParticleName()                
                                                                 ); 
                                   //std::cout << "Padre " << G4ParticleTable::GetParticleTable()->FindParticle(81+rstate-i)->GetParticleName() << std::endl;  
                                   //std::cout << "Hijo " << G4ParticleTable::GetParticleTable()->FindParticle(mode[j])->GetParticleName() << std::endl;
                                   //std::cout << "BR " << brat[j]/100. << std::endl;                            
                                   
                                   decayTable->Insert(decayChannel);   
                                   }
                                }
                            
                                if (decayTable->entries() > 0) {
        decayTable->DumpInfo();  
    }
    
     DRAGONIon* part = dynamic_cast<DRAGONIon*>(G4ParticleTable::GetParticleTable()->FindParticle(81 + rstate - i));
    if (part) {
        part->SetDecayTable(decayTable);  
    } else {
        std::cout << "Error: La partícula no es un DRAGONIon." << std::endl;  
    }
                          } 
//C
//C    Set cross-section variables
//C
                      gp = part_width;
                      gg = gam_width;
                      er = resenerg;
                      omg = spin_stat_fac;
                      ell = ell;
                      m1 = abeam*amugev + beam_mass_excess;
                      std::cout << "m1" << std::endl;
    
                      if(atarg == 1.)
                        {
                         m2 = hmass;
                         std::cout << "++++++" << m2 << std::endl;
                         }
                      else if(atarg == 4.)
                             {m2 = hemass;}
                      else if(atarg == 12.)
                             {m2 = c12mass;}
            
                      z1 = zbeam;
                      z2 = ztarg;  
                      
                      std::cout << "YUYUYUY" << ipart << std::endl;
                      break;
                      }
//C     alpha source
              case 21:
                     {
                      alpha = true;
                      atarg = 1;
                      aamass = hemass;
                      ubuf[0] = fkine[0];
                      ubuf[1] = fkine[1];
                      tlif = 1000.;
                      ipart = 80;
                      
                      DRAGONIon* ion = new DRAGONIon(
                                                  "Alpha",                  
                                                  aamass*GeV,      
                                                  0.0*keV,                    
                                                  ubuf[0]*eplus,          
                                                  0, +1, 0, 0, 0, 0,       
                                                  "nucleus",              
                                                  0, 0, ipart,               
                                                  false,                  
                                                  tlif*s,                 
                                                  nullptr,                
                                                  false,
                                                  2,
                                                  4              
                                                  );   
                      break;  
                      }   
                               
              Materials* materials = Materials::Instance();  
              materials->atarg = atarg;
              }
     } 

}













#include "G4SystemOfUnits.hh"       

#include "Materials.hh"                          //local

void Materials::ugmate_trgt()
{
      //From ugmate_trgt.f
      G4NistManager* nistManager = G4NistManager::Instance();

      G4double a, z, dens, radl, pdens, ptarg;
      G4double centralvac_factor;

      //C     Target materials (1 atm 20 deg C)
      a = atarg;
      if (atarg < 1.2) 
         {
          //C     Hydrogen
          z = 1;
          a *= g/mole; 
          dens = 8.38e-5 * g/cm3 ;   
          radl = 752300 * cm;    
          } 
      else 
         {
          //C     Helium
          z = 2;
          a = 4.0026 * g/mole;
          dens = 1.6586e-4 * g/cm3 ;
          radl = 568686 * cm;
          }

     //C
     //C     Target pressure
     ptarg = 5.0/760.;
     targetl = 5.5 * cm;
     exitdens = targetl;
     entdens = targetl;
	 
     std::cout << "AAAAAAAAAAAAAAAAAAAAAAAA" << std::endl;
	 std::cout << " z: " << z << std::endl;
	 std::cout << " a: " << a/(g/mole) << std::endl;
	 std::cout << " atarg: " << atarg/cm << std::endl;
     std::cout << " targetl: " << targetl/cm << std::endl;
	 std::cout << " entdens: " << entdens/cm << std::endl;
	 std::cout << " exitdens: " << exitdens/cm << std::endl;


     //C
     //C     Solid target, still want to set mtarg
      if (atarg == 12.)
         {
          mtarg = 62;                                                                          //62
          Target = Target_Carbon; 
          }
      else
         {
          mtarg = 27;   
          G4Material* Gas_Target = new G4Material("target", z, a, ptarg*dens);             //27     		  
          Target = Gas_Target;                                
		  materialMap[27] = Gas_Target;
          }
    
     // Central vacuum material (1/500th of central pressure)
     mcent = 28;                                                                               //28
     centralvac_factor = 0.002;
     pdens = ptarg * centralvac_factor;
     CentralVacuum = new G4Material("centralvacuum", z, a, pdens * dens);
	 materialMap[22] = CentralVacuum;
	 materialMap[28] = CentralVacuum;

      //C     gas pressures for inside stepped collimator sections and eff.length
      ment[0] = 29;                                                                            //29
      pdens = ptarg*0.05; //!! ratio from Knudsen equation. 
	  std::cout << " entdens0: " << entdens/cm << std::endl;
      entdens = entdens + pdens*5.08/ptarg *cm; //!! (5.08cm section length)
      Entrance1 = new G4Material("ent1", z, a,  pdens * dens);
	  materialMap[29] = Entrance1;

   
	 std::cout << " pdens: " << pdens << std::endl;
	 std::cout << " ptarg: " << ptarg << std::endl;
	 std::cout << " entdens1: " << entdens/cm << std::endl;
	 std::cout << " exitdens: " << exitdens/cm << std::endl;
	  
      //C
      ment[1] = 30;                                                                            //30
      pdens = ptarg*0.036; //!! ratio from Knudsen equation.
      entdens = entdens + pdens*5.08/ptarg *cm;
      Entrance2 = new G4Material("ent2", z, a, pdens*dens);
	  materialMap[30] = Entrance2;

	 std::cout << " pdens: " << pdens << std::endl;	  
	 std::cout << " entdens: " << entdens/cm << std::endl;
	 std::cout << " exitdens: " << exitdens/cm << std::endl;
      //C
      ment[2] = 31;                                                                            //31
      pdens = ptarg*0.016; //!! ratio from Knudsen equation.
      entdens = entdens + pdens*5.08/ptarg *cm;
      Entrance3 = new G4Material("ent3", z, a, pdens*dens);
	  materialMap[31] = Entrance3;

	 std::cout << " pdens: " << pdens << std::endl;	  
	 std::cout << " entdens: " << entdens/cm << std::endl;
	 std::cout << " exitdens: " << exitdens/cm << std::endl;
	 
      //C
      mex[0] = 32;                                                                             //32
      pdens = ptarg*0.05; //!! ratio from Knudsen equation. 
      exitdens = exitdens + pdens*5.08/ptarg *cm; //!! (5.08cm section length)
      Ext1 = new G4Material("ext1", z, a, pdens*dens);
	  materialMap[29] = Ext1;

	 std::cout << " pdens: " << pdens << std::endl;	  
	 std::cout << " entdens: " << entdens/cm << std::endl;
	 std::cout << " exitdens: " << exitdens/cm << std::endl;
	 
      //C
      mex[1] = 33;                                                                             //33
      pdens = ptarg*0.036; //!! ratio from Knudsen equation. 
      exitdens = exitdens + pdens*5.08/ptarg *cm; //!! (5.08cm section length)
      Ext2 = new G4Material("ext2", z, a, pdens*dens);
	  materialMap[29] = Ext2;

	 std::cout << " pdens: " << pdens << std::endl;	  
	 std::cout << " entdens: " << entdens/cm << std::endl;
	 std::cout << " exitdens: " << exitdens/cm << std::endl;
	 
      //C
      mex[2] = 34;                                                                             //34
      pdens = ptarg*0.016; //!! ratio from Knudsen equation. 
      exitdens = exitdens + pdens*5.08/ptarg *cm; //!! (5.08cm section length)
      Ext3 = new G4Material("ext3", z, a, pdens*dens);
	  materialMap[29] = Ext3;

	 std::cout << " pdens: " << pdens << std::endl;	  
	 std::cout << " entdens: " << entdens/cm << std::endl;
	 std::cout << " exitdens: " << exitdens/cm << std::endl;
	 
      //C
      //C     Box gas baseline pressure
      mbox = 35;                                                                               //35
      pdens = ptarg*0.056; //!! DAH Gas target profile report June 2002
      entdens = entdens + pdens*2.2/ptarg *cm;
      exitdens = exitdens + pdens*2.2/ptarg *cm;
      Baseline = new G4Material("base", z, a, pdens*dens); 
	  materialMap[35] = Baseline;

	 std::cout << " pdens: " << pdens << std::endl;	  
	 std::cout << " entdens: " << entdens/cm << std::endl;
	 std::cout << " exitdens: " << exitdens/cm << std::endl;
	 
	 std::cout << "AAAAAAAAAAAAAAAAAAAAAAAA" << std::endl;
}






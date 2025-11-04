#include "DRAGONEMField.hh"
#include "G4RunManager.hh"


namespace DRAGON {

void DRAGONEMField::mitray_field(G4String devname, G4double xpos[3], G4double bfld[3], G4double efld[3]) const
     { 
      /*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      C                                                                      C
      C     This subroutine takes as input the name of the element and the   C
      C     position of the particle, and returns the magnetic and electric  C
      C     field components due to the specified element.                   C
      C                                                                      C
      C     INPUT:                                                           C
      C     -------                                                          C
      C     DEVNAME         CHAR*12   Name of element whose B/E field        C
      C                               is to be evaluated;  must match one    C
      C                               of the names in ITITLE(i),i=1,NO       C
      C                                                                      C
      C     XPOS(i)         REAL*8(3) Vector containing the coordinates      C
      C                               of the point where the B/E field       C
      C                               components must be evaluated, in       C
      C                               the RAYTRACE A-coordinate system of    C
      C                               the element. [cm]                      C
      C                                                                      C
      C     OUTPUT:                                                          C
      C     -------                                                          C
      C     BFLD(i)         REAL*8(3) B-field components Bx, By, Bz in       C
      C                               BFLD(1), BFLD(2), BFLD(3). [Tesla]     C
      C                                                                      C
      C     EFLD(i)         REAL*8(3) E-field components Ex, Ey, Ez in       C
      C                               EFLD(1), EFLD(2), EFLD(3).             C
      C                                                                      C
      C     Original Geant3 Code: S. Yen (TRIUMF)  e-mail  STAN@TRIUMF.CA    C                     
      C     Adapted into Geant4: Ben Marlow Nov-19-2024 bmarlow@triumf.ca    C
      C                                                                      C
      CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/

     G4int i, ino, itype;
     G4double devdata[mmax];
     
     for(i=0;i<3;i++) 
	{
	 bfld[i] = 0.;
	 efld[i] = 0.;
	 }
	
//C
//C *** Set flag for diagnostic mode
//C
     ldiag = false;                 
     GlobalVariables::GetInstance()->SetIswit(5,0);     //OJO sacar de aqui e implementarlo en el BeginRunAction
     if(GlobalVariables::GetInstance()->GetIswit(5)==1) ldiag = true;
     if(ldiag) std::cout << "Entering subroutine MITRAY_FIELD";
  
//C
//C *** Now find which element number in the beamline this corresponds to.
//C
     G4bool found = false;
 
     no = GlobalVariables::GetInstance()->GetNo();
 
     for (i = 0; i < no; i++) 
         {
	  if (devname == GlobalVariables::GetInstance()->GetItitle(i)) 
	     {
	      ino = i;
	      found = true;
	      std::cout << "KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK " << std::endl;
	      std::cout << devname << " " << GlobalVariables::GetInstance()->GetItitle(i) << std::endl;
	      std::cout << "KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK " << std::endl;
	      break;
	      }
	  }
 
     if(!found) 
       {
	std::cout << "MITRAY_FIELD called with unknown device name " << devname << std::endl;
	std::cout << "All B- and E-field components set to zero." << std::endl;
	std::cout << "!!! Abort current event !!!" << std::endl;
	//istop = 1;     //OJO
	//ieotri = 1;    //OJO
	G4RunManager::GetRunManager()->AbortEvent();
	
	return;
	}
	
//C *** INO now contains the number code for the device, where INO=1
//C     corresponds to the first device of the beam transport system, etc.
//C
//C *** IDATA(INO) is the type code for this device, e.g. 2=dipole, etc.

      itype = GlobalVariables::GetInstance()->GetIdata(ino); 
     		
      for (i = 0; i < mmax; i++) 
          {	
	   devdata[i] = GlobalVariables::GetInstance()->GetData(i,ino);  //devdata[i] = data[i][ino];
//	   std::cout << "i " << i << " devdata[i]  " << devdata[i] << "\n";
	   }
	 
	//dipole magnet
      if(itype==2)mitray_dipole(devdata,xpos,bfld);

	//electrostatic deflector
      else if(itype==7)mitray_edipol(devdata,xpos,efld);

	//solenoid magnet
           else if(itype==14)mitray_solnd(devdata,xpos,bfld);

	//multipole magnet
	        else if(itype==9)mitray_poles(devdata,xpos,bfld);

	//SASP dipole magnet
	             else if(itype==20)mitray_sasp(devdata,xpos,bfld);

	                  else 
	                     {
		              std::cout << "**ERROR** IN MITRAY_FIELD:" << std::endl;
		              std::cout << "NO FIELD CALC DEFINED FOR ITYPE= " << itype << std::endl;
	                      }
}

}

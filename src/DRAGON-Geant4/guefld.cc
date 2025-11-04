#include "G4SystemOfUnits.hh"        //!geant
#include "DRAGONEMField.hh"
#include "G4LogicalVolumeStore.hh"

#include "geom_dipole.hh"
#include "geom_edipol.hh"
#include "geom_mpole.hh"
#include "geom_sole.hh"

#include "geant3functions.hh"

namespace DRAGON {


void DRAGONEMField::guefld(G4double vector[7], G4double time, G4double* bfield, G4double* efield) const
     {
//c     GEANT mother volume:
//c     --------------------
//c
//c     XM(3) is in the beam direction (down stream)
//c     XM(2) is UP
//c     XM(1) completes a right-handed coordinate system

//c     Local variables    

      G4int k, irot;
      G4int iyes, ifld; 
      G4int iglvolu;

      G4double xm[3], xd[3], bfm[3], bfd[3], efm[3], efd[3];
      G4double B_field, E_field, vmod;

      G4double xdd[3], bfdd[3], efdd[3];

      G4String chname, new_chname;
      G4String kdname;

      G4int nlevel_old,number_old[2],lvolum_old[2];
      G4double gonly_old[15];
      G4int nlevel = 1,number[2],lvolum[2],gonly[2];
      G4String names[2],names_old[2];

 G4ThreeVector center(vector[0],vector[1],vector[2]);   
 G4VPhysicalVolume* physicalVolume = fNavigator->LocateGlobalPointAndSetup(center*cm);
 
 std::cout << "YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY" << std::endl;
 std::cout << physicalVolume->GetName() << std::endl;
 std::cout << "YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY" << std::endl;

 G4String devname = physicalVolume->GetName();

if(devname != "WRLD")
  {
   names[nlevel] = devname ;
   
      for (G4int i=0;i<DRAGONEMelements.size();i++)
          {
//std::cout << "DRAGONEMelements[i].name " << DRAGONEMelements[i].name << std::endl;
//std::cout << "names[nlevel] " << names[nlevel] << std::endl;
           if(names[nlevel] == DRAGONEMelements[i].name)
             {
              number[nlevel] = DRAGONEMelements[i].multiplicity;
              ifld = DRAGONEMelements[i].ifld;
              gonly[nlevel] = DRAGONEMelements[i].manyORonly;
              kdname = DRAGONEMelements[i].type;
              irot = DRAGONEMelements[i].irot;
              theta[0] = DRAGONEMelements[i].theta1;
              theta[1] = DRAGONEMelements[i].theta2;
              theta[2] = DRAGONEMelements[i].theta3;
              phi[0] = DRAGONEMelements[i].phi1;
              phi[1] = DRAGONEMelements[i].phi2;
              phi[2] = DRAGONEMelements[i].phi3;  
                        
//std::cout << "number[nlevel] " << number[nlevel] << std::endl;
//std::cout << "ifld " << ifld << std::endl;
//std::cout << "gonly[nlevel] " << gonly[nlevel] << std::endl;
//std::cout << "kdname " << kdname << std::endl;
              break;
              }
           }
 
      static G4int istep = 0;
std::cout << "VVVVVVVVVVVVVVVVVVGEGINVVVVVVVVVVVVVVVVVVVV"<< std::endl;
      iglvolu = 0;
std::cout << "iglvolu " << iglvolu << std::endl;    
      istep = istep + 1;
std::cout << "AFUERA istep " << istep << std::endl;
 
      if(istep == 1)
        {
std::cout << "-------------------------------------"<< std::endl;
std::cout << "DENTRO istep " << istep << std::endl;        
         nlevel_old = nlevel;
std::cout << "nlevel " << nlevel << std::endl;
std::cout << "chname " << names[nlevel] << std::endl;
         names_old[nlevel] = names[nlevel];
std::cout << "chname " << names_old[nlevel] << std::endl;
std::cout << "number " << number[nlevel] << std::endl;
         number_old[nlevel] = number[nlevel];
std::cout << "number " << number_old[nlevel] << std::endl;
         lvolum_old[nlevel] = lvolum[nlevel];
std::cout << "gonly " << gonly[nlevel] << std::endl;
         gonly_old[nlevel] = gonly[nlevel];
std::cout << "gonly " << gonly_old[nlevel] << std::endl;
std::cout << "-------------------------------------"<< std::endl;
         }
   
      if (ifield == 1) 
         {
          if (istep % 3 == 0) 
             {istep = 0;}
          } 
      else 
         {istep = 0;} 
//c
//c     Initialize the field values
//c
      vzero(bfield,3);
      vzero(efield,3);

//c
//c *** Obtain the field for the present volume
//c
      k = number[nlevel];

      chname = names[nlevel];

std::cout << "ifield " << ifld << std::endl;
std::cout << "k " << k << std::endl;
std::cout << "chname " << chname << std::endl;
std::cout << "kdname " << kdname << std::endl;
std::cout << "vector " << vector[0] << " " << vector[1] << " " << vector[2] << std::endl;
std::cout << "xm " << xm[0] << " " << xm[1] << " " << xm[2] << std::endl;
      ucopy(vector,xm,3); 
std::cout << "vector " << vector[0] << " " << vector[1] << " " << vector[2] << std::endl;
std::cout << "xm " << xm[0] << " " << xm[1] << " " << xm[2] << std::endl;
std::cout << "xd " << xd[0] << " " << xd[1] << " " << xd[2] << std::endl;
      gmtod(xm,xd,1,chname);
std::cout << "xm " << xm[0] << " " << xm[1] << " " << xm[2] << std::endl;
std::cout << "xd " << xd[0] << " " << xd[1] << " " << xd[2] << std::endl;
std::cout << "vector " << vector[0] << " " << vector[1] << " " << vector[2] << std::endl;


//std::cout << "PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP\n";
//std::cout << "Q";
//std::cout << " dx_mpole " << dx_mpole[0][1] << " " << dx_mpole[1][1] << " " << dx_mpole[2][1] << std::endl;
//std::cout << " dx_mpole " << dx_mpole[0][2] << " " << dx_mpole[1][2] << " " << dx_mpole[2][2] << std::endl;
//std::cout << " dx_mpole " << dx_mpole[0][3] << " " << dx_mpole[1][3] << " " << dx_mpole[2][3] << std::endl;
//std::cout << " dx_mpole " << dx_mpole[0][4] << " " << dx_mpole[1][4] << " " << dx_mpole[2][4] << std::endl;
//std::cout << "E";
//std::cout << " dx_edipol " << dx_edipol[0][1] << " " << dx_edipol[1][1] << " " << dx_edipol[2][1] << std::endl;
//std::cout << " dx_edipol " << dx_edipol[0][2] << " " << dx_edipol[1][2] << " " << dx_edipol[2][2] << std::endl;
//std::cout << "D";
//std::cout << " dx_dipole " << dx_dipole[0][1] << " " << dx_dipole[1][1] << " " << dx_dipole[2][1] << std::endl;
//std::cout << " dx_dipole " << dx_dipole[0][2] << " " << dx_dipole[1][2] << " " << dx_dipole[2][2] << std::endl;
//std::cout << "PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP\n";

     if(kdname == "D")
        {  //Problema aqui
         G4double dx_pole[3] = {dx_dipole[0][k-1],dx_dipole[1][k-1],dx_dipole[2][k-1]};
         std::cout << "VBVBVBVBVBVB" << std::endl;   
         std::cout << "dx_pole " << dx_pole[0] << " " << dx_pole[1] << " " << dx_pole[2] << std::endl;
         std::cout << "theta " << theta[0] << " " << theta[1] << " " << theta[2] << std::endl;
         std::cout << "phi " << phi[0] << " " << phi[1] << " " << phi[2] << std::endl;
         std::cout << "irot " << irot << std::endl;
         std::cout << "xd " << xd[0] << " " << xd[1] << " " << xd[2] << std::endl;                 
         gitran(xd,dx_pole,irot,xd,theta,phi); 
         std::cout << "xd " << xd[0] << " " << xd[1] << " " << xd[2] << std::endl;  
         std::cout << "VBVBVBVBVBVB" << std::endl;    
         }
      else if(kdname == "E")
            {
std::cout << "EEEEEEEEEEEEEEE 1" << std::endl;            
             G4double dx_pole[3] = {dx_edipol[0][k-1],dx_edipol[1][k-1],dx_edipol[2][k-1]};
std::cout << "irot " << irot << std::endl;
std::cout << "xd " << xd[0] << " " << xd[1] << " " << xd[2] << std::endl;
std::cout << "dx_edipol " << dx_pole[0] << " " << dx_pole[1] << " " << dx_pole[2] << std::endl;
             gitran(xd,dx_pole,irot,xd,theta,phi);
std::cout << "xd " << xd[0] << " " << xd[1] << " " << xd[2] << std::endl;
std::cout << "EEEEEEEEEEEEEEE 2" << std::endl; 
             }
           else if (kdname == "Q")
                   {
                    irot = 0;
                    G4double dx_pole[3] = {dx_mpole[0][k-1],dx_mpole[1][k-1],dx_mpole[2][k-1]};
                    gitran(xd,dx_pole,0,xd);
                    }
                else if (kdname == "S")
                        {
                         irot = 0;
                         gitran(xd,dx_sole[k-1],0,xd);
                         }
  
      xdd[0] = xd[0];
      xdd[1] = xd[1];
      xdd[2] = xd[2];

std::cout << "xd " << xd[0] << " " << xd[1] << " " << xd[2] << std::endl;
std::cout << "xdd " << xdd[0] << " " << xdd[1] << " " << xdd[2] << std::endl;

      devname = chname;

std::cout << "devname " << devname << std::endl;
std::cout << "chname " << chname << std::endl;

      mitray_field(devname, xdd, bfdd, efdd);   

//std::cout << "==================================" << std::endl;
//std::cout << "bfdd " << bfdd[0] << " " << bfdd[1] << " " << bfdd[2] << std::endl;
//std::cout << "efdd " << efdd[0] << " " << efdd[1] << " " << efdd[2] << std::endl;
//std::cout << "==================================" << std::endl;

      bfd[0] = bfdd[0];
      bfd[1] = bfdd[1];
      bfd[2] = bfdd[2];

std::cout << "==================================" << std::endl;
std::cout << "xdd " << xdd[0] << " " << xdd[1] << " " << xdd[2] << std::endl;
std::cout << "bfdd " << bfdd[0] << " " << bfdd[1] << " " << bfdd[2] << std::endl;
std::cout << "efdd " << efdd[0] << " " << efdd[1] << " " << efdd[2] << std::endl;
std::cout << "==================================" << std::endl;

      efd[0] = efdd[0];
      efd[1] = efdd[1];
      efd[2] = efdd[2];

      B_field = 1.E+01*Vmod(bfd,3);
      E_field = 1.E-06*Vmod(efd,3);

std::cout << "B_field " << B_field << std::endl;
std::cout << "E_field " << E_field << std::endl;

      if(B_field > 0.0)
        {
std::cout << "bfd " << bfd[0] << " " << bfd[1] << " " << bfd[2] << std::endl; 
         vunit(bfd,3);   //ver como funciona esta
std::cout << "bfd " << bfd[0] << " " << bfd[1] << " " << bfd[2] << std::endl;         
         if(irot != 0)
           {
std::cout << "irot " << irot << std::endl;
            bfm[0] = bfd[0]*std::sin(theta[0]*deg)*std::cos(phi[0]*deg)+
                     bfd[1]*std::sin(theta[1]*deg)*std::cos(phi[1]*deg)+ 
                     bfd[2]*std::sin(theta[2]*deg)*std::cos(phi[2]*deg);
            bfm[1] = bfd[0]*std::sin(theta[0]*deg)*std::sin(phi[0]*deg)+
                     bfd[1]*std::sin(theta[1]*deg)*std::sin(phi[1]*deg)+ 
                     bfd[2]*std::sin(theta[2]*deg)*std::sin(phi[2]*deg);
            bfm[2] = bfd[0]*std::cos(theta[0]*deg)+
                     bfd[1]*std::cos(theta[1]*deg)+ 
                     bfd[2]*std::cos(theta[2]*deg);

std::cout << "bfm " << bfm[0] << " " << bfm[1] << " " << bfm[2] << std::endl; 
std::cout << "bfd " << bfd[0] << " " << bfd[1] << " " << bfd[2] << std::endl; 
            ucopy(bfm,bfd,3);
std::cout << "bfm " << bfm[0] << " " << bfm[1] << " " << bfm[2] << std::endl; 
std::cout << "bfd " << bfd[0] << " " << bfd[1] << " " << bfd[2] << std::endl; 
            }
std::cout << "gdtom 4444444444444444444444444 " << std::endl;  
std::cout << "bfd " << bfd[0] << " " << bfd[1] << " " << bfd[2] << std::endl;          
std::cout << "bfm " << bfm[0] << " " << bfm[1] << " " << bfm[2] << std::endl; 
        gdtom(bfd,bfm,2,devname);
std::cout << "bfd " << bfd[0] << " " << bfd[1] << " " << bfd[2] << std::endl;          
std::cout << "bfm " << bfm[0] << " " << bfm[1] << " " << bfm[2] << std::endl; 
std::cout << "gdtom 4444444444444444444444444 " << std::endl; 

std::cout << "vscale 555555555555555555555555 " << std::endl;  
std::cout << "bfm " << bfm[0] << " " << bfm[1] << " " << bfm[2] << std::endl; 
std::cout << "B_field " << B_field << std::endl; 
std::cout << "bfield " << bfield[0] << " " << bfield[1] << " " << bfield[2] << std::endl; 
         vscale(bfm,B_field,bfield,3);
std::cout << "bfm " << bfm[0] << " " << bfm[1] << " " << bfm[2] << std::endl; 
std::cout << "B_field " << B_field << std::endl; 
std::cout << "bfield " << bfield[0] << " " << bfield[1] << " " << bfield[2] << std::endl; 
std::cout << "vscale 555555555555555555555555 " << std::endl; 
         }

      if(E_field > 0.0)
        {
std::cout << "efd " << efd[0] << " " << efd[1] << " " << efd[2] << std::endl;  
         vunit(efd,3);
std::cout << "efd " << efd[0] << " " << efd[1] << " " << efd[2] << std::endl; 
         if(irot != 0)
           {
std::cout << "irot " << irot << std::endl;
std::cout << "theta " << theta[0] << " " << theta[1] << " " << theta[2] << std::endl; 
std::cout << "phi " << phi[0] << " " << phi[1] << " " << phi[2] << std::endl; 

std::cout << "efm " << efm[0] << " " << efm[1] << " " << efm[2] << std::endl; 
            efm[0] = efd[0]*std::sin(theta[0]*deg)*std::cos(phi[0]*deg)+
                     efd[1]*std::sin(theta[1]*deg)*std::cos(phi[1]*deg)+
                     efd[2]*std::sin(theta[2]*deg)*std::cos(phi[2]*deg);
            efm[1] = efd[0]*std::sin(theta[0]*deg)*std::sin(phi[0]*deg)+
                     efd[1]*std::sin(theta[1]*deg)*std::sin(phi[1]*deg)+ 
                     efd[2]*std::sin(theta[2]*deg)*std::sin(phi[2]*deg);
            efm[2] = efd[0]*std::cos(theta[0]*deg)+
                     efd[1]*std::cos(theta[1]*deg)+ 
                     efd[2]*std::cos(theta[2]*deg);    ///OJO ERROR EN EL VALOR COMPONENTE Z

std::cout << "efm " << efm[0] << " " << efm[1] << " " << efm[2] << std::endl; 
std::cout << "efd " << efd[0] << " " << efd[1] << " " << efd[2] << std::endl; 
            ucopy(efm,efd,3); 
std::cout << "efm " << efm[0] << " " << efm[1] << " " << efm[2] << std::endl; 
std::cout << "efd " << efd[0] << " " << efd[1] << " " << efd[2] << std::endl;
            } 
std::cout << "efm " << efm[0] << " " << efm[1] << " " << efm[2] << std::endl; 
std::cout << "efd " << efd[0] << " " << efd[1] << " " << efd[2] << std::endl;
         gdtom(efd,efm,2,devname);
         vscale(efm,E_field,efield,3);
std::cout << "efm " << efm[0] << " " << efm[1] << " " << efm[2] << std::endl; 
std::cout << "efd " << efd[0] << " " << efd[1] << " " << efd[2] << std::endl;
std::cout << "E_field " << E_field << std::endl; 
         }
//c
//c *** Add all overlapping fields
//c
std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;

      if(gonly[nlevel] == 0)	
        {
std::cout << "gonly " << gonly[nlevel] << std::endl;
         for (G4int i=0;i<18;i++)
             { 
             // for (G4int i=0;i<18;i++){std::cout << i << " " << DRAGONEMelements[i].name << std::endl;}
             
              new_chname = DRAGONEMelements[i].name;
std::cout << "chname " << chname << std::endl;
std::cout << "new_chname " << new_chname << std::endl;

              if(new_chname == chname)
                {
//c
//c       Restore original COMMON/GCVOLU
//c             
                 if(iglvolu > 0)
                   {
std::cout << "Dentro de iglvolu" << std::endl;
                    nlevel = nlevel_old; 
                    names[nlevel] = names_old[nlevel];
                    number[nlevel] = number_old[nlevel];
                    lvolum[nlevel] = lvolum_old[nlevel];
                    gonly[nlevel] = gonly_old[nlevel];
                    }
                 continue;
                 }
               
std::cout << "I AM INSIDE " << std::endl;  
              ifld = DRAGONEMelements[i].ifld;
              number[nlevel] = DRAGONEMelements[i].multiplicity;
              kdname = DRAGONEMelements[i].type;
              irot = DRAGONEMelements[i].irot;
              theta[0] = DRAGONEMelements[i].theta1;
              theta[1] = DRAGONEMelements[i].theta2;
              theta[2] = DRAGONEMelements[i].theta3;
              phi[0] = DRAGONEMelements[i].phi1;
              phi[1] = DRAGONEMelements[i].phi2;
              phi[2] = DRAGONEMelements[i].phi3; 
 
              if(ifld > 0)			
                {
std::cout << "ifld " << ifld << std::endl;
                 iglvolu = iglvolu + 1; 
std::cout << "iglvolu1 " << iglvolu << std::endl;
//std::cout << "xm " << xm[0] << " " << xm[1] << " " << xm[2] << std::endl; 
//std::cout << "xd " << xd[0] << " " << xd[1] << " " << xd[2] << std::endl; 
                 gmtod(xm,xd,1,new_chname);
std::cout << "xm " << xm[0] << " " << xm[1] << " " << xm[2] << std::endl; 
std::cout << "xd " << xd[0] << " " << xd[1] << " " << xd[2] << std::endl;  
             
                  G4VSolid* solid = G4LogicalVolumeStore::GetInstance()->GetVolume(new_chname)->GetSolid();

                  if (solid->Inside(G4ThreeVector(xd[0]*cm,xd[1]*cm,xd[2]*cm)) == kInside)
                      {iyes = 1;}
                  else
                      {iyes = 0;} 

                  std::cout << "iyes " << iyes << std::endl; 
                  
                  if (iyes > 0)
                     {
                      std::cout << "A"<< std::endl; 
                      std::cout << "A"<< std::endl; 
                      std::cout << "A"<< std::endl; 
                      std::cout << "A"<< std::endl;  
                      devname = new_chname;

                      iglvolu = iglvolu + 1;
std::cout << "devname " << devname << std::endl;
std::cout << "iglvolu2 " << iglvolu << std::endl;
std::cout << "kdname " << kdname << std::endl;

std::cout << "xm " << xm[0] << " " << xm[1] << " " << xm[2] << std::endl; 
std::cout << "xd " << xd[0] << " " << xd[1] << " " << xd[2] << std::endl; 
                      gmtod(xm,xd,1,devname);
std::cout << "xm " << xm[0] << " " << xm[1] << " " << xm[2] << std::endl; 
std::cout << "xd " << xd[0] << " " << xd[1] << " " << xd[2] << std::endl; 
             
                      k = number[nlevel];
             
                      if(kdname == "D")
                        {
std::cout << "xd " << xd[0] << " " << xd[1] << " " << xd[2] << std::endl; 
                         G4double dx_pole[3] = {dx_dipole[0][k-1],dx_dipole[1][k-1],dx_dipole[2][k-1]};
                         gitran(xd,dx_pole,irot,xd,theta,phi);
std::cout << "xd " << xd[0] << " " << xd[1] << " " << xd[2] << std::endl;                       
std::cout << "theta " << theta[0] << " " << theta[1] << " " << theta[2] << std::endl; 
std::cout << "phi " << phi[0] << " " << phi[1] << " " << phi[2] << std::endl; 
std::cout << "dx_pole " << dx_pole[0] << " " << dx_pole[1] << " " << dx_pole[2] << std::endl; 
std::cout << "irot " << irot << std::endl;
                         }
                      else if(kdname == "E")
                             {
                              G4double dx_pole[3] = {dx_edipol[0][k-1],dx_edipol[1][k-1],dx_edipol[2][k-1]};
                              gitran(xd,dx_pole,irot,xd,theta,phi);
                              }
                           else if (kdname == "Q")
                                   {
                                    irot = 0;
                                    G4double dx_pole[3] = {dx_mpole[0][k-1],dx_mpole[1][k-1],dx_mpole[2][k-1]};
                                    gitran(xd,dx_pole,0,xd);
                                    }
                                else if (kdname == "S")
                                        {
                                         irot = 0;
                                         gitran(xd,dx_sole[k-1],0,xd);
                                         }

                      xdd[0] = xd[0];
                      xdd[1] = xd[1];
                      xdd[2] = xd[2];
 
std::cout << "xdd " << xdd[0] << " " << xdd[1] << " " << xdd[2] << std::endl;
             
                      mitray_field(devname, xdd, bfdd, efdd);
              
                      bfd[0] = bfdd[0];
                      bfd[1] = bfdd[1];
                      bfd[2] = bfdd[2]; 
std::cout << "==================" << std::endl;
std::cout << "bfd " << bfd[0] << " " << bfd[1] << " " << bfd[2] << std::endl;

                      efd[0] = efdd[0];
                      efd[1] = efdd[1];
                      efd[2] = efdd[2]; 

std::cout << "efd " << efd[0] << " " << efd[1] << " " << efd[2] << std::endl;
std::cout << "==================" << std::endl;

                      B_field = 1.E+01*Vmod(bfd,3);
                      E_field = 1.E-06*Vmod(efd,3);

std::cout << "B_field " << B_field << std::endl;
std::cout << "E_field " << E_field << std::endl;
      
                      if(B_field > 0.0)
                        {
std::cout << "Dentro de B_field" << std::endl;
                         vunit(bfd,3);
std::cout << "bfd " << bfd[0] << " " << bfd[1] << " " << bfd[2] << std::endl;        
                         if(irot != 0)
                           {
                            bfm[0] = bfd[0]*std::sin(theta[0]*deg)*std::cos(phi[0]*deg)+
                                     bfd[1]*std::sin(theta[1]*deg)*std::cos(phi[1]*deg)+ 
                                     bfd[2]*std::sin(theta[2]*deg)*std::cos(phi[2]*deg);
                            bfm[1] = bfd[0]*std::sin(theta[0]*deg)*std::sin(phi[0]*deg)+
                                     bfd[1]*std::sin(theta[1]*deg)*std::sin(phi[1]*deg)+ 
                                     bfd[2]*std::sin(theta[2]*deg)*std::sin(phi[2]*deg);
                            bfm[2] = bfd[0]*std::cos(theta[0]*deg)+
                                     bfd[1]*std::cos(theta[1]*deg)+ 
                                     bfd[2]*std::cos(theta[2]*deg);
std::cout << "bfm " << bfm[0] << " " << bfm[1] << " " << bfm[2] << std::endl;  
std::cout << "bfd " << bfd[0] << " " << bfd[1] << " " << bfd[2] << std::endl;                            
                            ucopy(bfm,bfd,3); 
std::cout << "bfm " << bfm[0] << " " << bfm[1] << " " << bfm[2] << std::endl;  
std::cout << "bfd " << bfd[0] << " " << bfd[1] << " " << bfd[2] << std::endl; 
                            }
std::cout << "bfm " << bfm[0] << " " << bfm[1] << " " << bfm[2] << std::endl;  
std::cout << "bfd " << bfd[0] << " " << bfd[1] << " " << bfd[2] << std::endl;  
                        gdtom(bfd,bfm,2,devname);
std::cout << "bfm " << bfm[0] << " " << bfm[1] << " " << bfm[2] << std::endl;  
std::cout << "bfd " << bfd[0] << " " << bfd[1] << " " << bfd[2] << std::endl;  
std::cout << "B_field " << B_field << std::endl;
                        vscale(bfm,B_field,bfm,3);
std::cout << "bfm " << bfm[0] << " " << bfm[1] << " " << bfm[2] << std::endl;  
std::cout << "B_field" << B_field << std::endl;
std::cout << "bfield " << bfield[0] << " " << bfield[1] << " " << bfield[2] << std::endl; 
                        vadd(bfield,bfm,bfield,3);
std::cout << "bfield " << bfield[0] << " " << bfield[1] << " " << bfield[2] << std::endl; 
std::cout << "bfm " << bfm[0] << " " << bfm[1] << " " << bfm[2] << std::endl;  
                        } 

                      if(E_field > 0.0)
                        {
                         vunit(efd,3);
                         if(irot != 0)
                           {
                            efm[0] = efd[0]*std::sin(theta[0]*deg)*std::cos(phi[0]*deg)+
                                     efd[1]*std::sin(theta[1]*deg)*std::cos(phi[1]*deg)+
                                     efd[2]*std::sin(theta[2]*deg)*std::cos(phi[2]*deg);
                            efm[1] = efd[0]*std::sin(theta[0]*deg)*std::sin(phi[0]*deg)+
                                     efd[1]*std::sin(theta[1]*deg)*std::sin(phi[1]*deg)+ 
                                     efd[2]*std::sin(theta[2]*deg)*std::sin(phi[2]*deg);
                            efm[2] = efd[0]*std::cos(theta[0]*deg)+
                                     efd[1]*std::cos(theta[1]*deg)+ 
                                     efd[2]*std::cos(theta[2]*deg); 
                            ucopy(efm,efd,3);
                            }
                         gdtom(efd,efm,2,devname);
                         vscale(efm,E_field,efm,3);
                         vadd(efield,efm,efield,3); 
                         }
                     } 
                 }
             }
//c
//c       Restore original COMMON/GCVOLU
//c  
              if(iglvolu > 0)
                {
                 std::cout << "Dentro de iglvolu" << std::endl;
                 nlevel = nlevel_old; 
                 names[nlevel] = names_old[nlevel];
                 number[nlevel] = number_old[nlevel];
                 lvolum[nlevel] = lvolum_old[nlevel];
                 gonly[nlevel] = gonly_old[nlevel];
                 }
//OJO estas variables se definen en el primary, ver si esto se puede llevar al constructor mejor de esta clase
              G4double bscale = GlobalVariables::GetInstance()->Getbscale();
              G4double escale = GlobalVariables::GetInstance()->Getescale();
std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
//    Rescale the fields from the reference tune according to the reaction
std::cout << "Dentro de vscale" << std::endl;
std::cout << "bfield " << bfield[0] << " " << bfield[1] << " " << bfield[2] << std::endl;
std::cout << "efield " << efield[0] << " " << efield[1] << " " << efield[2] << std::endl; 
std::cout << "bscale " << bscale << std::endl; 
              vscale(bfield,bscale,bfield,3);
              vscale(efield,escale,efield,3);
std::cout << "bfield " << bfield[0] << " " << bfield[1] << " " << bfield[2] << std::endl;
std::cout << "efield " << efield[0] << " " << efield[1] << " " << efield[2] << std::endl; 
std::cout << "escale " << escale << std::endl;              
        }
std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;    
std::cout << "VVVVVVVVVVVVVVVVVVENDVVVVVVVVVVVVVVVVVVVV"<< std::endl;

  }
}

}



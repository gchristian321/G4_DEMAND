#include "geant3functions.hh"


namespace DRAGON {

//Geant3 helper functions transferred into Geant4
//Made by Ben Marlow, DRAGON Co-op Student, on 10 December, 2024

//===========================================================

//VECTOR FUNCTIONS

//===========================================================


void ucopy(G4double* array1, G4double* array2, G4int dimension) 
     {
      if(dimension <= 0) return;
      for(G4int i = 0; i < dimension; i++) 
         {array2[i] = array1[i];}
      }

void ucopy(G4String* names1, G4String* names2, G4int dimension)
     {  
      if(dimension <= 0) return;
      for(G4int i = 0; i < dimension; i++) 
         {names2[i] = names1[i];}
      }

void ucopy(G4int* array1, G4int* array2, G4int dimension) 
     {
      if(dimension <= 0) return;
      for(G4int i = 0; i < dimension; i++) 
         {array2[i] = array1[i];}
      }

void ucopy(G4double* array1, G4double array2[][max_photon], G4int rows, G4int col_index)
     {
      for (G4int i = 0; i < rows; ++i)
           array2[i][col_index] = array1[i];
      }

void vzero(G4double* array, G4int dimension) 
          {
	   if(dimension <= 0) return;
	   for(G4int i = 0; i < dimension; i++) 
	      {array[i] = 0.0;}
           }

void vzero(G4int* array, G4int dimension) 
          {
	   if(dimension <= 0) return;
	   for(G4int i = 0; i < dimension; i++) 
	      {array[i] = 0;}
           }

void vscale(G4double* input_array, G4double scale, G4double* output_array, G4int dimension) 
     {
      if(dimension <= 0) return;
      for(G4int i = 0; i < dimension; i++) 
         {output_array[i] = scale * input_array[i];}
      }

void vunit(G4double* input_array, G4int dimension) 
     {
      if(dimension <= 0) return;
      G4double norm = 0.0;
      for(G4int i = 0; i < dimension; i++) 
         {norm += input_array[i]*input_array[i];}
      if(norm <= 0) return;
      norm = std::sqrt(norm);
      for(G4int i = 0; i < dimension; i++) 
         {input_array[i] = input_array[i]/norm;}
      }

void vadd(G4double* array1, G4double* array2, G4double* array3, G4int dimension) {
	if(dimension <= 0) return;
	for(G4int i = 0; i< dimension; i++) {
		array3[i] = array1[i] + array2[i];
	}
}

G4double Vmod(G4double* vector, G4int n) 
         {
          G4double sum_squares = 0.0;
          for (G4int i = 0; i < n; ++i) 
              {sum_squares += vector[i] * vector[i];}

          return std::sqrt(sum_squares);  
          }

void vfill(G4double* array, G4int dimension, G4double value)
{
    for(G4int i=0;i<dimension;i++) 
		array[i]=value;
}

void vfill(G4int* array, G4int dimension, G4int value)
{
    for (G4int i = 0; i < dimension; ++i)
        array[i] = value;
}

void ublank(G4String* array, G4int dimension)
{
    if (dimension <= 0) return;
    for (G4int i = 0; i < dimension; ++i)
        array[i] = " ";
}

G4double vdotn(const G4double vec1[3][10], const G4double vec2[3][10], G4int col) {
    return vec1[0][col]*vec2[0][col] +
           vec1[1][col]*vec2[1][col] +
           vec1[2][col]*vec2[2][col];
}


void sortzv(G4double* a, G4int* index, G4int n, G4int mode, G4int nway, G4int nsort)
{
    for (int i = 0; i < n; i++)
        index[i] = i;  

    bool descending = (mode < 0);


    auto compare = [&](int i, int j) {
        if (descending)
            return a[i] > a[j]; 
        else
            return a[i] < a[j];  
    };

    std::sort(index, index + n, compare);
}






void gitran(G4double x[3], G4double dx[3], G4int irot, G4double xnew[3], G4double theta[3], G4double phi[3]) 
     {
/*
std::cout << "HOLA HOLA " << std::endl;
std::cout << "x " << x[0] << " " << x[1] << " " << x[2] << std::endl;
std::cout << "dx " << dx[0] << " " << dx[1] << " " << dx[2] << std::endl;
std::cout << "irot " << irot << std::endl;

std::cout << "theta[0] "  << theta[0] << std::endl;
//std::cout << "phi[0] "  << phi[0] << std::endl;
std::cout << "theta[1] "  << theta[1] << std::endl;
//std::cout << "phi[1] "  << phi[1] << std::endl;
std::cout << "theta[2] "  << theta[2] << std::endl;
std::cout << "phi[2] " << phi[2] << std::endl;
*/
G4double pi = 3.141592653589793;

      if(irot == 0)    
	{
	 xnew[0] = x[0] - dx[0];
	 xnew[1] = x[1] - dx[1];
	 xnew[2] = x[2] - dx[2];
	 }
      else
         {
          G4double xl0 = x[0] - dx[0];
	  G4double xl1 = x[1] - dx[1];
	  G4double xl2 = x[2] - dx[2];

/*
std::cout << "theta[0] " << std::setprecision(7)  << theta[0] << std::endl;
std::cout << "phi[0] " << std::setprecision(7)  << phi[0] << std::endl;
std::cout << "theta[1] " << std::setprecision(7)  << theta[1] << std::endl;
std::cout << "phi[1] " << std::setprecision(7)  << phi[1] << std::endl;
std::cout << "theta[2] " << std::setprecision(7)  << theta[2] << std::endl;
std::cout << "phi[2] " << std::setprecision(7)  << phi[2] << std::endl;


std::cout << "x[0] " << std::setprecision(7)  << x[0] << std::endl;
std::cout << "x[1] " << std::setprecision(7)  << x[1] << std::endl;
std::cout << "x[2] " <<  std::setprecision(7)  <<  x[2] << std::endl;

std::cout << "dx[0] " << std::setprecision(7)  << dx[0] << std::endl;
std::cout << "dx[1] " << std::setprecision(7)  << dx[1] << std::endl;
std::cout << "dx[2] " << std::setprecision(7)  << dx[2] << std::endl;

std::cout << "xl0 " << std::setprecision(7) << xl0 << std::endl;
std::cout << "xl1 " <<  std::setprecision(7)  << xl1 << std::endl;
std::cout << "xl2 " <<  std::setprecision(7)  << xl2 << std::endl;

std::cout << "theta[0] " << theta[0] << std::endl;
std::cout << "theta[1] " << theta[1] << std::endl;
std::cout << "theta[2] " << theta[2] << std::endl;
std::cout << "phi[0] " << phi[0] << std::endl;
std::cout << "phi[1] " << phi[1] << std::endl;
std::cout << "phi[2] " << phi[2] << std::endl;
*/	   
	   xnew[0] = xl0*std::sin(theta[0]*pi/180)*std::cos(phi[0]*pi/180) + xl1*std::sin(theta[0]*pi/180)*std::sin(phi[0]*pi/180)+xl2*std::cos(theta[0]*pi/180);
	   xnew[1] = xl0*std::sin(theta[1]*pi/180)*std::cos(phi[1]*pi/180) + xl1*std::sin(theta[1]*pi/180)*std::sin(phi[1]*pi/180)+xl2*std::cos(theta[1]*pi/180);
	   xnew[2] = xl0*std::sin(theta[2]*pi/180)*std::cos(phi[2]*pi/180) + xl1*std::sin(theta[2]*pi/180)*std::sin(phi[2]*pi/180)+xl2*std::cos(theta[2]*pi/180);

/*
std::cout << "Q(JR+1) " << std::setprecision(7)  << std::sin(theta[0]*pi/180.)*std::cos(phi[0]*pi/180.)  << std::endl;
std::cout << "Q(JR+2) " << std::setprecision(7)  << std::sin(theta[0]*pi/180.)*std::sin(phi[0]*pi/180.) << std::endl;
std::cout << "Q(JR+3) " << std::setprecision(7)  << std::cos(theta[0]*pi/180.) << std::endl;

std::cout << "Q(JR+4) " << std::setprecision(7)  << std::sin(theta[1]*pi/180.)*std::cos(phi[1]*pi/180.)  << std::endl;
std::cout << "Q(JR+5) " << std::setprecision(7)  << std::sin(theta[1]*pi/180.)*std::sin(phi[1]*pi/180.) << std::endl;
std::cout << "Q(JR+6) " << std::setprecision(7)  << std::cos(theta[1]*pi/180.) << std::endl;

std::cout << "Q(JR+7) " << std::setprecision(7)  << std::sin(theta[2]*pi/180.)*std::cos(phi[2]*pi/180.)  << std::endl;
std::cout << "Q(JR+8) " << std::setprecision(7)  << std::sin(theta[2]*pi/180.)*std::sin(phi[2]*pi/180.) << std::endl;
std::cout << "Q(JR+9) " << std::setprecision(7)  << std::cos(theta[2]*pi/180.) << std::endl;

std::cout << "xnew[0] " << std::setprecision(7)  << xnew[0] << std::endl;
std::cout << "xnew[1] " << std::setprecision(7)  << xnew[1] << std::endl;
std::cout << "xnew[2] " << std::setprecision(7)  << xnew[2] << std::endl;
*/


           }

//std::cout << "xnew " << xnew[0] << " " << xnew[1] << " " << xnew[2] << std::endl;
//std::cout << "dx " << dx[0] << " " << dx[1] << " " << dx[2] << std::endl;
//std::cout << "HOLA HOLA " << std::endl;
 
     }

void gdtom(G4double xd[3], G4double xm[3], G4int iflag, G4String devname) {
	
	G4double grmatN[10], gtranN[3];
	
    if(devname == "WRLD") {
		grmatN[0] = 1.; grmatN[1] = 0.; grmatN[2] = 0.;
		grmatN[3] = 0.; grmatN[4] = 1.; grmatN[5] = 0.;
		grmatN[6] = 0.; grmatN[7] = 0.; grmatN[8] = 1.;
		grmatN[9] = 0.;
		gtranN[0] = 0.; gtranN[1] = 0.; gtranN[2] = 0.;
	}
	else if(devname == "Q1") {
		grmatN[0] = 1.; grmatN[1] = 0.; grmatN[2] = 0.;
		grmatN[3] = 0.; grmatN[4] = 1.; grmatN[5] = 0.;
		grmatN[6] = 0.; grmatN[7] = 0.; grmatN[8] = 1.;
		grmatN[9] = 1.;
		gtranN[0] = 0.; gtranN[1] = 0.; gtranN[2] = 119.5;
	}
	else if(devname == "Q2") {
		grmatN[0] = 1.; grmatN[1] = 0.; grmatN[2] = 0.;
		grmatN[3] = 0.; grmatN[4] = 1.; grmatN[5] = 0.;
		grmatN[6] = 0.; grmatN[7] = 0.; grmatN[8] = 1.;
		grmatN[9] = 1.;
		gtranN[0] = 0.; gtranN[1] = 0.; gtranN[2] = 174.497513;
	}
	else if(devname == "D1") {
		grmatN[0] = -0.42261824; grmatN[1] = 0.; grmatN[2] = 0.906307757;
		grmatN[3] = 0.906307757; grmatN[4] = 0.; grmatN[5] = 0.42261824;
		grmatN[6] = 0.; grmatN[7] = 1.; grmatN[8] = 0.;
		grmatN[9] = 1.;
		gtranN[0] = -9.56029606; gtranN[1] = 0.; gtranN[2] = 297.259338;
	}
	else if(devname == "TST1") {
		grmatN[0] = 0.642787576; grmatN[1] = 0.; grmatN[2] = 0.766044438;
		grmatN[3] = -4.37113883E-08; grmatN[4] = 1.; grmatN[5] = -4.37113883E-08;
		grmatN[6] = -0.766044438; grmatN[7] = 0.; grmatN[8] = 0.642787635;
		grmatN[9] = 1.;
		gtranN[0] = -62.363575; gtranN[1] = 0.; gtranN[2] = 354.0466;
	}
	else if(devname == "Q3") {
		grmatN[0] = 0.642787576; grmatN[1] = 0.; grmatN[2] = 0.766044438;
		grmatN[3] = -4.37113883E-08; grmatN[4] = 1.; grmatN[5] = -4.37113883E-08;
		grmatN[6] = -0.766044438; grmatN[7] = 0.; grmatN[8] = 0.642787635;
		grmatN[9] = 1.;
		gtranN[0] = -121.268555; gtranN[1] = 0.; gtranN[2] = 403.473755;
	}
	else if(devname == "Q4") {
		grmatN[0] = 0.642787576; grmatN[1] = 0.; grmatN[2] = 0.766044438;
		grmatN[3] = -4.37113883E-08; grmatN[4] = 1.; grmatN[5] = -4.37113883E-08;
		grmatN[6] = -0.766044438; grmatN[7] = 0.; grmatN[8] = 0.642787635;
		grmatN[9] = 1.;
		gtranN[0] = -153.599457; gtranN[1] = 0.; gtranN[2] = 430.6026;
	}
	else if(devname == "Q5") {
		grmatN[0] = 0.642787576; grmatN[1] = 0.; grmatN[2] = 0.766044438;
		grmatN[3] = -4.37113883E-08; grmatN[4] = 1.; grmatN[5] = -4.37113883E-08;
		grmatN[6] = -0.766044438; grmatN[7] = 0.; grmatN[8] = 0.642787635;
		grmatN[9] = 1.;
		gtranN[0] = -195.731888; gtranN[1] = 0.; gtranN[2] = 465.955902;
	}
	else if(devname == "Q6") {
		grmatN[0] = 0.642787576; grmatN[1] = 0.; grmatN[2] = 0.766044438;
		grmatN[3] = -4.37113883E-08; grmatN[4] = 1.; grmatN[5] = -4.37113883E-08;
		grmatN[6] = -0.766044438; grmatN[7] = 0.; grmatN[8] = 0.642787635;
		grmatN[9] = 1.;
		gtranN[0] = -237.864319; gtranN[1] = 0.; gtranN[2] = 501.309204;
	}
	else if(devname == "Q7") {
		grmatN[0] = 0.642787576; grmatN[1] = 0.; grmatN[2] = 0.766044438;
		grmatN[3] = -4.37113883E-08; grmatN[4] = 1.; grmatN[5] = -4.37113883E-08;
		grmatN[6] = -0.766044438; grmatN[7] = 0.; grmatN[8] = 0.642787635;
		grmatN[9] = 1.;
		gtranN[0] = -270.195221; gtranN[1] = 0.; gtranN[2] = 528.438049;
	}
	else if(devname == "E1") {
		grmatN[0] = -0.866025388; grmatN[1] = 0.; grmatN[2] = 0.5;
		grmatN[3] = 0.5; grmatN[4] = 0.; grmatN[5] = 0.866025388;
		grmatN[6] = 0.; grmatN[7] = 1.; grmatN[8] = 0.;
		grmatN[9] = 1.;
		gtranN[0] = -392.964996; gtranN[1] = 0.; gtranN[2] = 563.123169;
	}
	else if(devname == "TST2") {
		grmatN[0] = 0.342020094; grmatN[1] = 0.; grmatN[2] = 0.939692616;
		grmatN[3] = -3.83236625E-08; grmatN[4] = 1.; grmatN[5] = -1.05293395E-07;
		grmatN[6] = -0.939692616; grmatN[7] = 0.; grmatN[8] = 0.342020184;
		grmatN[9] = 1.;
		gtranN[0] = -501.589172; gtranN[1] = 0.; gtranN[2] = 658.362976;
	}
	else if(devname == "Q8") {
		grmatN[0] = 0.342020094; grmatN[1] = 0.; grmatN[2] = 0.939692616;
		grmatN[3] = -3.83236625E-08; grmatN[4] = 1.; grmatN[5] = -1.05293395E-07;
		grmatN[6] = -0.939692616; grmatN[7] = 0.; grmatN[8] = 0.342020184;
		grmatN[9] = 1.;
		gtranN[0] = -585.597717; gtranN[1] = 0.; gtranN[2] = 688.939636;
	}
	else if(devname == "Q9") {
		grmatN[0] = 0.342020094; grmatN[1] = 0.; grmatN[2] = 0.939692616;
		grmatN[3] = -3.83236625E-08; grmatN[4] = 1.; grmatN[5] = -1.05293395E-07;
		grmatN[6] = -0.939692616; grmatN[7] = 0.; grmatN[8] = 0.342020184;
		grmatN[9] = 1.;
		gtranN[0] = -637.280823; gtranN[1] = 0.; gtranN[2] = 707.750732;
	}
	else if(devname == "Q10") {
		grmatN[0] = 0.342020094; grmatN[1] = 0.; grmatN[2] = 0.939692616;
		grmatN[3] = -3.83236625E-08; grmatN[4] = 1.; grmatN[5] = -1.05293395E-07;
		grmatN[6] = -0.939692616; grmatN[7] = 0.; grmatN[8] = 0.342020184;
		grmatN[9] = 1.;
		gtranN[0] = -677.170776; gtranN[1] = 0.; gtranN[2] = 722.26947;
	}
	else if(devname == "D2") {
		grmatN[0] = -0.953716993; grmatN[1] = 0.; grmatN[2] = -0.300705791;
		grmatN[3] = -0.300705791; grmatN[4] = 0.; grmatN[5] = 0.953716993;
		grmatN[6] = 0.; grmatN[7] = 1.; grmatN[8] = 0.;
		grmatN[9] = 1.;
		gtranN[0] = -781.07251; gtranN[1] = 0.; gtranN[2] = 741.305786;
	}
	else if(devname == "Q11") {
		grmatN[0] = -0.819152057; grmatN[1] = 0.; grmatN[2] = 0.57357645;
		grmatN[3] = -1.21985684E-08; grmatN[4] = 1.; grmatN[5] = -1.61318837E-07;
		grmatN[6] = -0.573576391; grmatN[7] = 0.; grmatN[8] = -0.819151998;
		grmatN[9] = 1.;
		gtranN[0] = -901.007996; gtranN[1] = 0.; gtranN[2] = 600.789185;
	}
	else if(devname == "Q12") {
		grmatN[0] = -0.819152057; grmatN[1] = 0.; grmatN[2] = 0.57357645;
		grmatN[3] = -1.21985684E-08; grmatN[4] = 1.; grmatN[5] = -1.61318837E-07;
		grmatN[6] = -0.573576391; grmatN[7] = 0.; grmatN[8] = -0.819151998;
		grmatN[9] = 1.;
		gtranN[0] = -925.356323; gtranN[1] = 0.; gtranN[2] = 566.016174;
	}
	else if(devname == "E2") {
		grmatN[0] = -0.30070582; grmatN[1] = 0.; grmatN[2] = -0.953716993;
		grmatN[3] = -0.953716993; grmatN[4] = 0.; grmatN[5] = 0.30070582;
		grmatN[6] = 0.; grmatN[7] = 1.; grmatN[8] = 0.;
		grmatN[9] = 1.;
		gtranN[0] = -941.914368; gtranN[1] = 0.; gtranN[2] = 397.334045;
	}
	else if(devname == "Q13") {
		grmatN[0] = -1.; grmatN[1] = 0.; grmatN[2] = 4.29742499E-08;
		grmatN[3] = 4.86795244E-08; grmatN[4] = 1.; grmatN[5] = -1.50584398E-07;
		grmatN[6] = 5.76116967E-08; grmatN[7] = 0.; grmatN[8] = -1.;
		grmatN[9] = 1.;
		gtranN[0] = -1025.10291; gtranN[1] = 0.; gtranN[2] = 201.410583;
	}
	else if(devname == "Q14") {
		grmatN[0] = -1.; grmatN[1] = 0.; grmatN[2] = 4.29742499E-08;
		grmatN[3] = 4.86795244E-08; grmatN[4] = 1.; grmatN[5] = -1.50584398E-07;
		grmatN[6] = 5.76116967E-08; grmatN[7] = 0.; grmatN[8] = -1.;
		grmatN[9] = 1.;
		gtranN[0] = -1025.10291; gtranN[1] = 0.; gtranN[2] = 134.810577;
	}
	else if(devname == "HOL0") {
		grmatN[0] = 1.; grmatN[1] = 0.; grmatN[2] = -8.59484999E-08;
		grmatN[3] = -8.67542241E-15; grmatN[4] = 1.; grmatN[5] = 0.;
		grmatN[6] = -1.15223393E-07; grmatN[7] = 0.; grmatN[8] = 1.;
		grmatN[9] = 1.;
		gtranN[0] = -1025.10291; gtranN[1] = 0.; gtranN[2] = 8.79128265;
	}
	else if(devname == "MCP0") {
		grmatN[0] = -1.; grmatN[1] = 0.; grmatN[2] = 1.2892275E-07;
		grmatN[3] = 4.86795351E-08; grmatN[4] = 1.; grmatN[5] = -1.50584398E-07;
		grmatN[6] = 1.72835087E-07; grmatN[7] = 0.; grmatN[8] = -1.;
		grmatN[9] = 1.;
		gtranN[0] = -1025.10291; gtranN[1] = 0.; gtranN[2] = 8.79128265;
	}
    else if(devname == "HOL1") {
		grmatN[0] = 1.; grmatN[1] = 0.; grmatN[2] = -8.59484999E-08;
		grmatN[3] = -8.67542241E-15; grmatN[4] = 1.; grmatN[5] = 0.;
		grmatN[6] = -1.15223393E-07; grmatN[7] = 0.; grmatN[8] = 1.;
		grmatN[9] = 1.;
		gtranN[0] = -1025.10291; gtranN[1] = 0.; gtranN[2] = -50.2087173;
	}
	else if(devname == "MCP1") {
		grmatN[0] = -1.; grmatN[1] = 0.; grmatN[2] = 1.2892275E-07;
		grmatN[3] = 4.86795351E-08; grmatN[4] = 1.; grmatN[5] = -1.50584398E-07;
		grmatN[6] = 1.72835087E-07; grmatN[7] = 0.; grmatN[8] = -1.;
		grmatN[9] = 1.;
		gtranN[0] = -1025.10291; gtranN[1] = 0.; gtranN[2] = -50.2087173;
	}
	if(iflag == 1) {
		if(grmatN[9] != 0) {
			xm[0] = gtranN[0] + grmatN[0] * xd[0] + grmatN[3] * xd[1] + grmatN[6] * xd[2];
			xm[1] = gtranN[1] + grmatN[1] * xd[0] + grmatN[4] * xd[1] + grmatN[7] * xd[2];
			xm[2] = gtranN[2] + grmatN[2] * xd[0] + grmatN[5] * xd[1] + grmatN[8] * xd[2];
		}
		
		else {
			xm[0] = gtranN[0] + xd[0];
			xm[1] = gtranN[1] + xd[1];
			xm[2] = gtranN[2] + xd[2];
		}
	}
	
	else {
		if(grmatN[9] != 0) {
			xm[0] = grmatN[0] * xd[0] + grmatN[3] * xd[1] + grmatN[6] * xd[2];
			xm[1] = grmatN[1] * xd[0] + grmatN[4] * xd[1] + grmatN[7] * xd[2];
			xm[2] = grmatN[2] * xd[0] + grmatN[5] * xd[1] + grmatN[8] * xd[2];
		}
		
		else {
			xm[0] = xd[0];
			xm[1] = xd[1];
			xm[2] = xd[2];
		}
	}
}

void gmtod(G4double xm[3], G4double xd[3], G4int iflag, G4String devname) 
     {
      G4double grmatN[10], gtranN[3];
	
      if(devname == "Q1") 
        {
	 grmatN[0] = 1.; grmatN[1] = 0.; grmatN[2] = 0.;
	 grmatN[3] = 0.; grmatN[4] = 1.; grmatN[5] = 0.;
	 grmatN[6] = 0.; grmatN[7] = 0.; grmatN[8] = 1.;
	 grmatN[9] = 1.;
	 gtranN[0] = 0.; gtranN[1] = 0.; gtranN[2] = 119.5;
	 std::cout << "Ejecutando gmtod para Q1" << std::endl;
	 }
      else if(devname == "Q2") 
             {
	      grmatN[0] = 1.; grmatN[1] = 0.; grmatN[2] = 0.;
	      grmatN[3] = 0.; grmatN[4] = 1.; grmatN[5] = 0.;
	      grmatN[6] = 0.; grmatN[7] = 0.; grmatN[8] = 1.;
	      grmatN[9] = 1.;
	      gtranN[0] = 0.; gtranN[1] = 0.; gtranN[2] = 174.497513;
	      }
	   else if(devname == "D1") 
	          {
		   grmatN[0] = -0.42261824; grmatN[1] = 0.; grmatN[2] = 0.906307757;
		   grmatN[3] = 0.906307757; grmatN[4] = 0.; grmatN[5] = 0.42261824;
		   grmatN[6] = 0.; grmatN[7] = 1.; grmatN[8] = 0.;
		   grmatN[9] = 1.;
		   gtranN[0] = -9.56029606; gtranN[1] = 0.; gtranN[2] = 297.259338;
	           }
	        else if(devname == "Q3") 
	               {
		        grmatN[0] = 0.642787576; grmatN[1] = 0.; grmatN[2] = 0.766044438;
		        grmatN[3] = -4.37113883E-08; grmatN[4] = 1.; grmatN[5] = -4.37113883E-08;
		        grmatN[6] = -0.766044438; grmatN[7] = 0.; grmatN[8] = 0.642787635;
		        grmatN[9] = 1.;
		        gtranN[0] = -121.268555; gtranN[1] = 0.; gtranN[2] = 403.473755;
	                }
 	             else if(devname == "Q4") 
 	                    {
		             grmatN[0] = 0.642787576; grmatN[1] = 0.; grmatN[2] = 0.766044438;
		             grmatN[3] = -4.37113883E-08; grmatN[4] = 1.; grmatN[5] = -4.37113883E-08;
	 	             grmatN[6] = -0.766044438; grmatN[7] = 0.; grmatN[8] = 0.642787635;
		             grmatN[9] = 1.;
		             gtranN[0] = -153.599457; gtranN[1] = 0.; gtranN[2] = 430.6026;
	                     }
	                  else if(devname == "Q5")
	                         {
		                  grmatN[0] = 0.642787576; grmatN[1] = 0.; grmatN[2] = 0.766044438;
		                  grmatN[3] = -4.37113883E-08; grmatN[4] = 1.; grmatN[5] = -4.37113883E-08;
		                  grmatN[6] = -0.766044438; grmatN[7] = 0.; grmatN[8] = 0.642787635;
		                  grmatN[9] = 1.;
		                  gtranN[0] = -195.731888; gtranN[1] = 0.; gtranN[2] = 465.955902;
	                          }
	                       else if(devname == "Q6") 
	                              {
		                       grmatN[0] = 0.642787576; grmatN[1] = 0.; grmatN[2] = 0.766044438;
		                       grmatN[3] = -4.37113883E-08; grmatN[4] = 1.; grmatN[5] = -4.37113883E-08;
	 	                       grmatN[6] = -0.766044438; grmatN[7] = 0.; grmatN[8] = 0.642787635;
		                       grmatN[9] = 1.;
		                       gtranN[0] = -237.864319; gtranN[1] = 0.; gtranN[2] = 501.309204;
	                               }
	                            else if(devname == "Q7") 
	                                   {
		                            grmatN[0] = 0.642787576; grmatN[1] = 0.; grmatN[2] = 0.766044438;
		                            grmatN[3] = -4.37113883E-08; grmatN[4] = 1.; grmatN[5] = -4.37113883E-08;
		                            grmatN[6] = -0.766044438; grmatN[7] = 0.; grmatN[8] = 0.642787635;
		                            grmatN[9] = 1.;
		                            gtranN[0] = -270.195221; gtranN[1] = 0.; gtranN[2] = 528.438049;
	                                    }
	                                  else if(devname == "E1") 
	                                         {
		                                  grmatN[0] = -0.866025388; grmatN[1] = 0.; grmatN[2] = 0.5;
		                                  grmatN[3] = 0.5; grmatN[4] = 0.; grmatN[5] = 0.866025388;
	 	                                  grmatN[6] = 0.; grmatN[7] = 1.; grmatN[8] = 0.;
		                                  grmatN[9] = 1.;
		                                  gtranN[0] = -392.964996; gtranN[1] = 0.; gtranN[2] = 563.123169;
	                                          }
	                                       else if(devname == "Q8") 
	                                              {
		                                       grmatN[0] = 0.342020094; grmatN[1] = 0.; grmatN[2] = 0.939692616;
		                                       grmatN[3] = -3.83236625E-08; grmatN[4] = 1.; grmatN[5] = -1.05293395E-07;
		                                       grmatN[6] = -0.939692616; grmatN[7] = 0.; grmatN[8] = 0.342020184;
		                                       grmatN[9] = 1.;
		                                       gtranN[0] = -585.597717; gtranN[1] = 0.; gtranN[2] = 688.939636;
	                                               }
	                                            else if(devname == "Q9") 
	                                                   {
		                                            grmatN[0] = 0.342020094; grmatN[1] = 0.; grmatN[2] = 0.939692616;
		                                            grmatN[3] = -3.83236625E-08; grmatN[4] = 1.; grmatN[5] = -1.05293395E-07;
		                                            grmatN[6] = -0.939692616; grmatN[7] = 0.; grmatN[8] = 0.342020184;
		                                            grmatN[9] = 1.;
		                                            gtranN[0] = -637.280823; gtranN[1] = 0.; gtranN[2] = 707.750732;
	                                                    }
	                                                 else if(devname == "Q10")
	                                                        {
		                                                 grmatN[0] = 0.342020094; grmatN[1] = 0.; grmatN[2] = 0.939692616;
		                                                 grmatN[3] = -3.83236625E-08; grmatN[4] = 1.; grmatN[5] = -1.05293395E-07;
		                                                 grmatN[6] = -0.939692616; grmatN[7] = 0.; grmatN[8] = 0.342020184;
		                                                 grmatN[9] = 1.;
		                                                 gtranN[0] = -677.170776; gtranN[1] = 0.; gtranN[2] = 722.26947;
	                                                         }
	                                                      else if(devname == "D2") 
	                                                             {
		                                                      grmatN[0] = -0.953716993; grmatN[1] = 0.; grmatN[2] = -0.300705791;
		                                                      grmatN[3] = -0.300705791; grmatN[4] = 0.; grmatN[5] = 0.953716993;
		                                                      grmatN[6] = 0.; grmatN[7] = 1.; grmatN[8] = 0.;
		                                                      grmatN[9] = 1.;
		                                                      gtranN[0] = -781.07251; gtranN[1] = 0.; gtranN[2] = 741.305786;
	                                                              }
	                                                           else if(devname == "Q11") 
	                                                                  {
		                                                           grmatN[0] = -0.819152057; grmatN[1] = 0.; grmatN[2] = 0.57357645;
		                                                           grmatN[3] = -1.21985684E-08; grmatN[4] = 1.; grmatN[5] = -1.61318837E-07;
		                                                           grmatN[6] = -0.573576391; grmatN[7] = 0.; grmatN[8] = -0.819151998;
		                                                           grmatN[9] = 1.;
		                                                           gtranN[0] = -901.007996; gtranN[1] = 0.; gtranN[2] = 600.789185;
	                                                                   }
	                                                                else if(devname == "Q12")
	                                                                       {
		                                                                grmatN[0] = -0.819152057; grmatN[1] = 0.; grmatN[2] = 0.57357645;
		                                                                grmatN[3] = -1.21985684E-08; grmatN[4] = 1.; grmatN[5] = -1.61318837E-07;
		                                                                grmatN[6] = -0.573576391; grmatN[7] = 0.; grmatN[8] = -0.819151998;
		                                                                grmatN[9] = 1.;
		                                                                gtranN[0] = -925.356323; gtranN[1] = 0.; gtranN[2] = 566.016174;
	                                                                        }
	                                                                     else if(devname == "E2") 
	                                                                            {
		                                                                     grmatN[0] = -0.30070582; grmatN[1] = 0.; grmatN[2] = -0.953716993;
		                                                                     grmatN[3] = -0.953716993; grmatN[4] = 0.; grmatN[5] = 0.30070582;
		                                                                     grmatN[6] = 0.; grmatN[7] = 1.; grmatN[8] = 0.;
		                                                                     grmatN[9] = 1.;
		                                                                     gtranN[0] = -941.914368; gtranN[1] = 0.; gtranN[2] = 397.334045;
	                                                                             }
	                                                                          else if(devname == "Q13")
	                                                                                 {
		                                                                          grmatN[0] = -1.; grmatN[1] = 0.; grmatN[2] = 4.29742499E-08;
		                                                                          grmatN[3] = 4.86795244E-08; grmatN[4] = 1.; grmatN[5] = -1.50584398E-07;
		                                                                          grmatN[6] = 5.76116967E-08; grmatN[7] = 0.; grmatN[8] = -1.;
		                                                                          grmatN[9] = 1.;
		                                                                          gtranN[0] = -1025.10291; gtranN[1] = 0.; gtranN[2] = 201.410583;
	                                                                                  }
	                                                                               else if(devname == "Q14") 
	                                                                                      {
		                                                                               grmatN[0] = -1.; grmatN[1] = 0.; grmatN[2] = 4.29742499E-08;
		                                                                               grmatN[3] = 4.86795244E-08; grmatN[4] = 1.; grmatN[5] = -1.50584398E-07;
		                                                                               grmatN[6] = 5.76116967E-08; grmatN[7] = 0.; grmatN[8] = -1.;
		                                                                               grmatN[9] = 1.;
		                                                                               gtranN[0] = -1025.10291; gtranN[1] = 0.; gtranN[2] = 134.810577;
	                                                                                       }
	
      if(iflag == 1) 
        {
	 if(grmatN[9] != 0) 
	   {
	    G4double t1 = xm[0] - gtranN[0];
            G4double t2 = xm[1] - gtranN[1];
            G4double t3 = xm[2] - gtranN[2];
            xd[0] = grmatN[0]*t1 + grmatN[1]*t2 + grmatN[2]*t3;
            xd[1] = grmatN[3]*t1 + grmatN[4]*t2 + grmatN[5]*t3;
            xd[2] = grmatN[6]*t1 + grmatN[7]*t2 + grmatN[8]*t3;
	    }
	 else 
	    {
	     xd[0] = xm[0] - gtranN[0];
             xd[1] = xm[1] - gtranN[1];
             xd[2] = xm[2] - gtranN[2];
	     }
	 }
      else 
         {
	  if(grmatN[9] != 0) 
	    {
	     xd[0] = grmatN[0] * xm[0] + grmatN[1] * xm[1] + grmatN[2] * xm[2];
	     xd[1] = grmatN[3] * xm[0] + grmatN[4] * xm[1] + grmatN[5] * xm[2];
	     xd[2] = grmatN[6] * xm[0] + grmatN[7] * xm[1] + grmatN[8] * xm[2];
	     }
	  else 
	     {
	      xd[0] = xm[0];
	      xd[1] = xm[1];
	      xd[2] = xm[2];
	      }
	}
     }

}

#include "DRAGONDetectorConstruction.hh"     //local
#include "mitray_setup.hh"         
#include "geom_dipole.hh"
#include "geom_mpole.hh"
#include "geom_edipol.hh"
#include "geom_sole.hh"

#include <fstream>                           //std        


G4int ndipole;                        
G4int irot_dipole[max_dipole];        
G4double dx_dipole[3][max_dipole];     
G4double gap_dipole[max_dipole];      
G4double phi_dipole[max_dipole];      
G4double r_dipole[max_dipole];        
G4double dr_dipole[max_dipole];       
G4double alpha_dipole[max_dipole];    
G4double beta_dipole[max_dipole];    
G4double z11_dipole[max_dipole];      
G4double z22_dipole[max_dipole];      
G4int jcol_dipole[2][max_dipole];      
G4double xcol_dipole[2][max_dipole];  
G4double ycol_dipole[2][max_dipole];  
G4double dxcol_dipole[2][max_dipole];  
G4double dycol_dipole[2][max_dipole];  

G4int nedipol;    
G4double gap_edipol[max_edipol];
G4double phi_edipol[max_edipol];
G4double r_edipol[max_edipol];
G4double dr_edipol[max_edipol];
G4double z11_edipol[max_edipol];
G4double z22_edipol[max_edipol];
G4int jcol_edipol[2][max_edipol];      
G4double xcol_edipol[2][max_edipol];   
G4double ycol_edipol[2][max_edipol];   
G4double dxcol_edipol[2][max_edipol];  
G4double dycol_edipol[2][max_edipol];
G4double dx_edipol[3][max_edipol];  
G4int irot_edipol[max_edipol];   

G4int nmpole; 
G4double dx_mpole[3][max_mpole];       
G4double efblength_mpole[max_mpole];   
G4double r_mpole[max_mpole];           
G4double z11_mpole[max_mpole];         
G4double z22_mpole[max_mpole];         
G4int jcol_mpole[2][max_mpole];        
G4double xcol_mpole[2][max_mpole];     
G4double ycol_mpole[2][max_mpole];     
G4double dxcol_mpole[2][max_mpole];    
G4double dycol_mpole[2][max_mpole];    

G4int nsole;                          
G4double pos_sole[3][max_sole];       
G4double dx_sole[3][max_sole];        
G4double efblength_sole[max_sole];    
G4double r_sole[max_sole];            
G4double z11_sole[max_sole];          
G4double z22_sole[max_sole];          
G4int jcol_sole[2][max_sole];         
G4double xcol_sole[2][max_sole];      
G4double ycol_sole[2][max_sole];      
G4double dxcol_sole[2][max_sole];     
G4double dycol_sole[2][max_sole];     

G4double pos_dipole[3][max_dipole] = {
        {-9.56029606, 0.0, 297.259338},
        {-781.07251, 0.0, 741.305786}  
    };

G4double pos_edipol[3][max_dipole] = {
        {-392.964996, 0.0, 563.123169},
        {-941.914368, 0.0, 397.334045}  
    };

G4double pos_mpole[max_mpole][3] = {
        {0., 0., 119.5},
        {0., 0., 174.497513},
        {-121.268555, 0., 403.473755},   
        {-153.599457, 0., 430.6026},
        {-195.731888, 0., 465.955902},
        {-237.864319, 0., 501.309204},
        {-270.195221, 0., 528.438049},
        {-585.597717, 0., 688.939636},
        {-637.280823, 0., 707.750732},
        {-677.170776, 0., 722.26947},
        {-901.007996, 0., 600.789185},
        {-925.356323, 0., 566.016174},
        {-1025.10291, 0., 201.410583},
        {-1025.10291, 0., 134.810577}      
    };

G4double col_pos[] = {         0.,     0.,    217.590012,  //RC9    C1     
				      -42.5028305,     0.,    371.648407,  //QSLT   C2
				       -76.249176,     0.,    331.431061,          
				      -59.3760033,  26.25,    351.539734,  //       C3
				      -59.3760033, -26.25,    351.539734,   												                   
				      -79.7144775,     0.,    368.605743,  //RC12   C4 
				       -299.00235,     0.,    552.610107,  //RC23   C5
				      -320.428619,     0.,    570.588928,  //RC25   C6
				      -326.180939,     0.,    588.384033,  //FC1    C7
				      -339.036713,     0.,     573.06311,  
				      -332.608826,    28.,    580.723572,  //       C8
					  -332.608826,   -28.,    580.723572, 
					  -404.458099,     0.,    633.651978,  //FC2    C9
					  -411.298492,     0.,    614.858154, 
					  -407.878296,    28.,    624.255066,  //       C10
					  -407.878296,   -28.,    624.255066, 
	                    -453.4534,     0.,    640.843018,  //RC27   C11
					  -495.384613,     0.,    664.884155,  //MSLT   C12
					  -501.027924,     0.,    649.379272,  
					  -498.206268,  13.75,    657.131714,  //       C13
					  -498.206268, -13.75,    657.131714, 
					  -523.108154,     0.,    666.195251,  //RC28   C14
						-561.0177,     0.,    679.993225,  //RC30   C15
					  -859.309509,     0.,    660.340637,  //RC40   C16
					  -877.100159,     0.,    634.933044,  //RC42   C17
					  -947.123535,     0.,    534.929321,  //RC47   C18
					  -964.431213,     0.,    510.211395,  //RC49   C19
					  -982.902954,     0.,    501.110413,  //FC3    C20
					  -966.519897,     0.,    489.638855, 
					  -974.711426,    30.,    495.374634,  //       C21
				  	  -974.711426,   -30.,    495.374634,
	                  -1035.10291,     0.,    335.835571,  //FC4    C22
					  -1015.10291,     0.,    335.835571,
					  -1025.10291,    30.,    335.835571,  //       C23
					  -1025.10291,   -30.,    335.835571,									                   
					  -1025.10291,     0.,    318.585571,  //RC51   C24
					  -1025.10291,     0.,    282.435577,  //RC53   C25
					  -1025.10291,     0.,    249.710587,  //RC55   C26
					  -1025.10291,     0.,    52.6259308,  //RC58   C27
					  -1049.85291,     0.,   -6.20871639,  //FSLT   C28
					  -1000.35291,     0.,    -6.2087183, 
                      -1025.10291,  24.75,   -6.20872116,  //       C29
                      -1025.10291, -24.75,   -6.20872116, 
                      -1025.10291,     0.,   -59.3437157,  //FSB2   C30
                      -1025.10291,     0.,   -68.7087173}; //ICS1   C31
                                             
G4double col_rot_angles[] = {0.,       0.,    0.,  //RC9    C1     
				         0.,       0.,    0.,  //QSLT   C2
				         0.,       0.,    0.,         
				         0.,       0.,    0.,  //       C3
				         0.,      50.,    0.,   												                   
				         0.,      50.,    0.,  //RC12   C4 
				         0.,      50.,    0.,  //RC23   C5
				         0.,      50.,    0.,  //RC25   C6
				         0.,       0.,    0.,  //FC1    C7
				         0.,       0.,    0.,  
				         0.,       0.,    0.,  //       C8
					     0.,      50.,    0., 
					     0.,       0.,    0.,  //FC2    C9
					     0.,       0.,    0., 
					     0.,       0.,    0.,  //       C10
					     0.,      70.,    0., 
	                     0.,      70.,    0.,  //RC27   C11
					     0.,       0.,    0.,  //MSLT   C12
					     0.,       0.,    0.,  
					     0.,       0.,    0.,  //       C13
					     0.,      70.,    0., 
					     0.,      70.,    0.,  //RC28   C14
						 0.,      70.,    0.,  //RC30   C15
					     0.,     145.,    0.,  //RC40   C16
					     0.,     145.,    0.,  //RC42   C17
					     0.,     145.,    0.,  //RC47   C18
					     0.,     145.,    0.,  //RC49   C19
					     0.,       0.,    0.,  //FC3    C20
					     0.,       0.,    0., 
					     0.,       0.,    0.,  //       C21
				  	     0.,     145.,    0.,
	                     0.,       0.,    0.,  //FC4    C22
					     0.,       0.,    0.,
					     0.,       0.,    0.,  //       C23
					     0.,     180.,    0.,									                   
					     0.,     180.,    0.,  //RC51   C24
					     0.,     180.,    0.,  //RC53   C25
					     0.,     180.,    0.,  //RC55   C26
					     0.,     180.,    0.,  //RC58   C27
					     0.,     180.,    0.,  //FSLT   C28
					     0.,     180.,    0., 
                         0.,     180.,    0.,  //       C29
                         0.,     180.,    0., 
                         0.,     180.,    0.,  //FSB2   C30
                         0.,     180.,    0.}; //ICS1   C31

G4double dipole_pos[] = {0.,       0.,    0.,  //VV1
				         0.,       0.,    0.,  //VV2
				         0.,       0.,    0.,  //VV3
				         0.,       0.,    0.,  //VV4
				         0.,       0.,    0.,  //VV5
				         0.,       0.,    0.,  //VV6
				         0.,       0.,    0.,  //VV7
				         0.,       0.,    0.,  //VV8
				         0.,       0.,    0.,  //VV9
				         0.,       0.,    0.,  //VV10
				         0.,       0.,    0.,  //VV11
				         0.,       0.,    0.,  //VV12
				         0.,       0.,    0.,  //VV13
				         0.,       0.,    0.,  //VV14
				         0.,       0.,    0.,  //VV15
				         0.,       0.,    0.,  //TRAP
				         };  
                 
G4double dipole_rot_angles[] = {0.,       0.,    0.,  //VV1
				                0.,       0.,    0.,  //VV2
				                0.,       0.,    0.,  //VV3
				                0.,       0.,    0.,  //VV4
				                0.,       0.,    0.,  //VV5
				                0.,       0.,    0.,  //VV6
				                0.,       0.,    0.,  //VV7
				                0.,       0.,    0.,  //VV8
				                0.,       0.,    0.,  //VV9
				                0.,       0.,    0.,  //VV10
				                0.,       0.,    0.,  //VV11
				                0.,       0.,    0.,  //VV12
				                0.,       0.,    0.,  //VV13
				                0.,       0.,    0.,  //VV14
				                0.,       0.,    0.,  //VV15
				                0.,       0.,    0.,  //TRAP
				                };  

G4double mpole_rot_angles[] = {0.,      0.,    0.,  //Q1
			                   0.,      0.,    0.,  //Q2
				               0.,     50.,    0.,  //Q3
				               0.,     50.,    0.,  //Q4
				               0.,     50.,    0.,  //Q5
				               0.,     50.,    0.,  //Q6
				               0.,     50.,    0.,  //Q7
				               0.,     70.,    0.,  //Q8
				               0.,     70.,    0.,  //Q9
				               0.,     70.,    0.,  //Q10
				               0.,     145,    0.,  //Q11   
				               0.,     145,    0.,  //Q12   
				               0.,      0.,    0.,  //Q13
				               0.,      0.,    0.,  //Q14
				               };

namespace DRAGON
{
void DRAGONDetectorConstruction::mitray_setup()
     {
	  //CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      //C                                                                      C
      //C	This subroutine reads the input file for MIT-RAYTRACE              C
      //C	and sets up the array of parameters specifying the beam            C
      //C       transport system.                                              C
      //C                                                                      C
      //C       It was adapted from the Aug 1989 version of MIT-RAYTRACE by    C
      //C       S. Yen (TRIUMF)  e-mail  STAN@TRIUMF.CA       and              C
      //C       P. Gumplinger (TRIUMF) e-mail GUM@TRIUMF.CA                    C
      //C                                                                      C
      //CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	  
	  G4int i, J, ielement;

      G4String nwd, mitray, jtitle;
         
      for (i=0;i<nmax;i++)
           idata[i]=0.0;
      for (i=0;i<mmax;i++)
          {
	 	  for (J=0;J<nmax;J++)
	 	       data[i][J] = 0.0;
	       }
      //C *** Now open the file containing the specification
      //C *** of all the beam transport elements   
      
      mitray = "dragon_2014_IC.dat";
      std::ifstream infile(mitray.c_str());
   
      if (!infile.is_open()) 
         {
          mitray = "mitray.dat";
          infile.clear(); 
          infile.open(mitray.c_str());
          }  
    
      //C Output data file for ascii data, EOC   
      std::ofstream file20("TP7.dat");
      std::ofstream file21("chargeslit.dat");
      std::ofstream file23("massslit.dat");
      std::ofstream file24("finalslit.dat");
      
      std::ofstream lout("Errors and warnings.dat");
    
      if (!infile.is_open())
         {
          G4cout << " Error opening MITRAY file: " << mitray << G4endl;
          std::exit(EXIT_FAILURE);
          }
	 
	  no = -1;
	  
	  G4String line;	
	  	  
	  for(ielement=0;ielement<nmax;ielement++)
	        {
	         std::getline(infile, line);
   	         std::istringstream iss(line);
	         iss >> nwd >> jtitle;
	                             
	         if(nwd == "COMM")
			   {;}   //C         This is a comment, so do nothing!
	         else
	            {
	             if(nwd == "DIPO")
			       { 
		   		    //C          DIPOLE  LENS           TYPE = 2	  
			        no = no + 1;
                    ititle[no] = jtitle;
                    idata[no] = 2;
					   
			        for (G4int J=0;J<6;J++)       //! RECORD 2
			             infile >> data[J][no]; 
			   	    for (G4int J=10;J<15;J++)     //! RECORD 3
			             infile >> data[J][no]; 
			        for (G4int J=15;J<18;J++)     //! RECORD 4
			             infile >> data[J][no]; 
				    for (G4int J=18;J<22;J++)     //! RECORD 5
			             infile >> data[J][no]; 
			        for (G4int J=24;J<28;J++)     //! RECORD 6
			             infile >> data[J][no]; 
				    for (G4int J=28;J<34;J++)     //! RECORD 7
			             infile >> data[J][no]; 
			        for (G4int J=34;J<40;J++)     //! RECORD 8
			             infile >> data[J][no]; 
				    for (G4int J=40;J<46;J++)     //! RECORD 9
			             infile >> data[J][no]; 
				    for (G4int J=46;J<50;J++)     //! RECORD 10
			             infile >> data[J][no]; 
				    for (G4int J=50;J<57;J++)     //! RECORD 11
			             infile >> data[J][no]; 
				    for (G4int J=57;J<64;J++)     //! RECORD 12
			             infile >> data[J][no]; 
			        std::getline(infile, line); 
			        }  
			     else
			        {
			   	     if(nwd == "EINZ")	  
				       {
			    	    //C          EINZEL  LENS           TYPE = 3	  
			            no = no + 1;
                        ititle[no] = jtitle;
                        idata[no] = 3;  
						  
					    for (G4int J=0;J<2;J++)       //! RECORD 2
			                 infile >> data[J][no]; 
			     	    for (G4int J=9;J<14;J++)      //! RECORD 3
			                 infile >> data[J][no]; 
			            for (G4int J=14;J<16;J++)     //! RECORD 4
			                 infile >> data[J][no]; 
			            for (G4int J=16;J<22;J++)     //! RECORD 5
			                 infile >> data[J][no]; 		           					 
					    std::getline(infile, line); 
					    }
			         else 
			            {
					     if(nwd == "EDIP")	 
					       {
			    	        //C          ELECTROSTATIC DEFLECTOR  TYPE = 7
			                no = no + 1;
                            ititle[no] = jtitle;
                            idata[no] = 7;  
						  
					        for (G4int J=0;J<4;J++)       //! RECORD 2
			                     infile >> data[J][no]; 
			     	        for (G4int J=10;J<15;J++)     //! RECORD 3
			                     infile >> data[J][no]; 
			                infile >> data[15][no];       //! RECORD 4
			                for (G4int J=16;J<20;J++)     //! RECORD 5
			                     infile >> data[J][no]; 
			                for (G4int J=24;J<28;J++)     //! RECORD 6
			                     infile >> data[J][no];
			                for (G4int J=28;J<34;J++)     //! RECORD 7
			                     infile >> data[J][no]; 
			                for (G4int J=34;J<40;J++)     //! RECORD 8
			                     infile >> data[J][no]; 		          		 
				            std::getline(infile, line); 
					        }
					     else 
					        {
						     if(nwd == "VELS")	 
					           {
			    	            //C          VELOCITY SELECTOR      TYPE = 8
			                    no = no + 1;
                                ititle[no] = jtitle;
                                idata[no] = 8;  
						  
					            for (G4int J=0;J<4;J++)       //! RECORD 2
			                         infile >> data[J][no]; 
			     	            for (G4int J=6;J<11;J++)      //! RECORD 3
			                         infile >> data[J][no]; 
			                    for (G4int J=11;J<13;J++)     //! RECORD 4
			                         infile >> data[J][no]; 
			                    for (G4int J=15;J<19;J++)     //! RECORD 5
			                         infile >> data[J][no];
			                    for (G4int J=19;J<23;J++)     //! RECORD 6
			                         infile >> data[J][no]; 
			                    for (G4int J=23;J<27;J++)     //! RECORD 7
			                         infile >> data[J][no]; 	
			                    for (G4int J=27;J<33;J++)     //! RECORD 8
			                         infile >> data[J][no]; 
			                    for (G4int J=33;J<39;J++)     //! RECORD 9
			                         infile >> data[J][no]; 
			                    for (G4int J=39;J<45;J++)     //! RECORD 10
			                         infile >> data[J][no]; 
			                    for (G4int J=45;J<51;J++)     //! RECORD 11
			                         infile >> data[J][no]; 
		                        std::getline(infile, line); 
		                        }
						      else  
    			                 {
    			                  if(nwd == "POLE")	 
					                {
			    	                 //C          MULTIPOLE (POLES)      TYPE =  9
			                         no = no + 1;
                                     ititle[no] = jtitle;
                                     idata[no] = 9;  
						    
					                 for (G4int J=0;J<3;J++)       //! RECORD 2
			                              infile >> data[J][no]; 
			                         for (G4int J=9;J<13;J++)      //! RECORD 3
			                              infile >> data[J][no]; 
			                         for (G4int J=13;J<18;J++)     //! RECORD 4
			                              infile >> data[J][no]; 
			                         for (G4int J=18;J<22;J++)     //! RECORD 5
			                              infile >> data[J][no]; 
			                         for (G4int J=22;J<28;J++)     //! RECORD 6
			                              infile >> data[J][no]; 
			                         for (G4int J=28;J<34;J++)     //! RECORD 7
			                              infile >> data[J][no]; 
			                         for (G4int J=34;J<42;J++)     //! RECORD 8
			                             infile >> data[J][no]; 
			                         std::getline(infile, line); 
			                         }
			                      else
			                         {
								      if(nwd == "MULT")	 
					                    {
			    	                     //C          MULTIPOLE (POLES)      TYPE =  10
			                             no = no + 1;
                                         ititle[no] = jtitle;
                                         idata[no] = 10;  
						    
					                     for (G4int J=0;J<2;J++)       //! RECORD 2
			                                  infile >> data[J][no]; 
			                             for (G4int J=9;J<15;J++)      //! RECORD 3
			                                  infile >> data[J][no]; 
			                             for (G4int J=15;J<17;J++)     //! RECORD 4
			                                  infile >> data[J][no]; 
			                             for (G4int J=19;J<25;J++)     //! RECORD 5
			                                  infile >> data[J][no]; 
			                             for (G4int J=25;J<28;J++)     //! RECORD 6
			                                  infile >> data[J][no]; 
			                             std::getline(infile, line); 
			                             }
			                          else
			                             {
									      if(nwd == "SHRT")	 
					                        {
			    	                         //C          SHIFT AND ROTATE       TYPE = 11
			                                 no = no + 1;
                                             ititle[no] = jtitle;
                                             idata[no] = 11;  
						    
					                         for (G4int J=0;J<6;J++)          //! RECORD 2
			                                      infile >> data[J][no]; 
			                                 std::getline(infile, line); 
			                                 }
			                               else
			                                  {
											   if(nwd == "DRIF")	 
					                             {
			    	                              //C          DRIFT                  TYPE = 12
			                                      no = no + 1;
                                                  ititle[no] = jtitle;
                                                  idata[no] = 12;  
						    
					                              infile >> data[0][no];      //! RECORD 2
			                                      std::getline(infile, line); 
			                                      }
			                                    else
			                                       {
											   	    if(nwd == "FCUP")	 
					                                  {
			    	                                   //C          FARRADAY CUP           TYPE = 21
			                                           no = no + 1;
                                                       ititle[no] = jtitle;
                                                       idata[no] = 21;  
						    
						                               for (G4int J=0;J<3;J++) 
					                                        infile >> data[J][no];   //! RECORD 2
			                                           std::getline(infile, line); 
			                                           }
			                                         else
			                                            {
														 if(nwd == "COLL")
														   {
															//C          COLLIMATOR             TYPE = 13
															no = no + 1;
                                                            ititle[no] = jtitle;
                                                            idata[no] = 13;  
						    
						                                    for (G4int J=0;J<5;J++) 
					                                             infile >> data[J][no];   //! RECORD 2
			                                                std::getline(infile, line);  
														    }
													     else
													        { 
														     if(nwd == "RCOL")	 
					                                           {
			    	                                            //C          COLLIMATOR             TYPE = 17
			                                                    no = no + 1;
                                                                ititle[no] = jtitle;
                                                                idata[no] = 17;  
						    
						                                        for (G4int J=0;J<6;J++) 
					                                                 infile >> data[J][no];   //! RECORD 2
			                                                    std::getline(infile, line); 
			                                                    }
			                                                 else
			                                                    {
														         if(nwd == "SOLE")	 
					                                               {
			    	                                                //C          SOLENOID               TYPE = 14
			                                                        no = no + 1;
                                                                    ititle[no] = jtitle;
                                                                    idata[no] = 14;  
						    
						                                            infile >> data[0][no];
						                                            for (G4int J=9;J<14;J++) 
					                                                     infile >> data[J][no];   //! RECORD 2
			                                                        for (G4int J=14;J<16;J++) 
					                                                     infile >> data[J][no];   //! RECORD 2
			                                                        std::getline(infile, line); 
			                                                        }   
			                                                     else
			                                                        {
																     if(nwd == "LENS")	 
					                                                   {
			    	                                                    //C          LENS                   TYPE = 15
			                                                            no = no + 1;
                                                                        ititle[no] = jtitle;
                                                                        idata[no] = 15;  
						    
						                                                for (G4int J=0;J<8;J++) 
					                                                         infile >> data[J][no];   //! RECORD 2
			                                                            for (G4int J=8;J<11;J++) 
					                                                         infile >> data[J][no];   //! RECORD 2
			                                                            std::getline(infile, line); 
			                                                            }    
			                                                         else
			                                                            { 
																	     if(nwd == "ACCE")	 
					                                                       {
			    	                                                        //C           ACCELERATOR            TYPE = 16
			                                                                no = no + 1;
                                                                            ititle[no] = jtitle;
                                                                            idata[no] = 16;  
						    
						                                                    for (G4int J=0;J<3;J++) 
					                                                             infile >> data[J][no];   //! RECORD 2
			                                                                for (G4int J=9;J<13;J++) 
					                                                             infile >> data[J][no];   //! RECORD 3
			                                                                infile >> data[13][no];
			                                                                for (G4int J=14;J<18;J++) 
					                                                             infile >> data[J][no];   //! RECORD 4
			                                                                for (G4int J=18;J<24;J++) 
					                                                             infile >> data[J][no];   //! RECORD 5
			                                                                for (G4int J=24;J<30;J++) 
					                                                             infile >> data[J][no];   //! RECORD 6
			                                                                std::getline(infile, line); 
			                                                                } 
			                                                             else
			                                                                {
																		     if(nwd == "SASP")	 
					                                                           {
			    	                                                            //C          Added by S.Yen - SASP clamshell dipole at TRIUMF
                                                                                //C          (same as normal DIPOLE except that sagging of the
                                                                                //C           field at the high field end is taken into account)

                                                                                //C          TRIUMF SASP DIPOLE     TYPE = 20
			                                                                    no = no + 1;
                                                                                ititle[no] = jtitle;
                                                                                idata[no] = 20;  
						    
						                                                        for (G4int J=0;J<6;J++) 
					                                                                 infile >> data[J][no];   //! RECORD 2
			                                                                    for (G4int J=10;J<15;J++) 
					                                                                 infile >> data[J][no];   //! RECORD 3
			                                                                    for (G4int J=15;J<18;J++) 
					                                                                 infile >> data[J][no];   //! RECORD 4
			                                                                    for (G4int J=18;J<22;J++) 
					                                                                 infile >> data[J][no];   //! RECORD 5
			                                                                    for (G4int J=24;J<28;J++) 
					                                                                 infile >> data[J][no];   //! RECORD 6
					                                                            for (G4int J=28;J<34;J++) 
					                                                                 infile >> data[J][no];   //! RECORD 7
					                                                            for (G4int J=34;J<40;J++) 
					                                                                 infile >> data[J][no];   //! RECORD 8
					                                                            for (G4int J=40;J<46;J++) 
					                                                                 infile >> data[J][no];   //! RECORD 9
					                                                            for (G4int J=46;J<50;J++) 
					                                                                 infile >> data[J][no];   //! RECORD 10
					                                                            for (G4int J=50;J<57;J++) 
					                                                                 infile >> data[J][no];   //! RECORD 11
					                                                            for (G4int J=57;J<64;J++) 
					                                                                 infile >> data[J][no];   //! RECORD 12
			                                                                    std::getline(infile, line); 
			                                                                    } 
			                                                                 else
			                                                                    {
																			     if(nwd == "ENDV")	 
					                                                               {
			    	                                                                //C          END CONTROL VOLUME     TYPE = 18
			                                                                        no = no + 1;
                                                                                    ititle[no] = jtitle;
                                                                                    idata[no] = 18;  
			                                                                        }    
																			     else
																			        {
																				     if(nwd == "TEST")	 
					                                                                   {
			    	                                                                    //C          TEST VOLUME            TYPE = 22
			                                                                            no = no + 1;
                                                                                        ititle[no] = jtitle;
                                                                                        idata[no] = 22;  
			                                                                            }
			                                                                         else
			                                                                            {
																				   	     if(nwd == "MCPF")	 
					                                                                       {
			    	                                                                        //C          MCP                   TYPE = 24
			                                                                                no = no + 1;
                                                                                            ititle[no] = jtitle;
                                                                                            idata[no] = 24;  
                                                                                          
                                                                                            for (G4int J=0;J<2;J++) 
					                                                                             infile >> data[J][no];   //! RECORD 12
					                                                                        std::getline(infile, line);
			                                                                                }	  
			                                                                              else
			                                                                                 {
																							  if(nwd == "TUBS")	 
					                                                                            {
			    	                                                                             //C          ANGLED TUBE            TYPE = 23
			                                                                                     no = no + 1;
                                                                                                 ititle[no] = jtitle;
                                                                                                 idata[no] = 23;  
                                                                                          
                                                                                                 for (G4int J=0;J<5;J++) 
					                                                                                  infile >> data[J][no];   //! RECORD 12
					                                                                             std::getline(infile, line);
			                                                                                     }
			                                                                                   else
			                                                                                      {
																							  	   if(nwd == "STRV")	 
					                                                                                 {
			    	                                                                                  //C          START CONTROL VOLUME   TYPE = 19
			                                                                                          no = no + 1;
                                                                                                      ititle[no] = jtitle;
                                                                                                      idata[no] = 19;  
                                                                                          			  }
                                                                                          		   else
                                                                                          		      {
																								  	   if(nwd == "SENT")	 
					                                                                                     {
			    	                                                                                      //SYSTEM END             TYPE = 1
			                                                                                              no = no + 1;
                                                                                                          ititle[no] = jtitle;
                                                                                                          idata[no] = 1;  
                                                                                                          //C	   Close the input file that has been opened
                                                                                                          
                                                                                                          infile.close();
                                                                                                          
                                                                                                          //C      INITIALIZE Geometry
                                                                                                          ugeom_setup();                                                                                                          
                                                                                                          break;
                                                                                          			      } 
                                                                                          			    else 
                                                                                          			        G4cout << " UNKNOWN ELEMENT TYPE " << nwd << " IGNORED" << G4endl;
																									   }
																							      }	
																							}
																					   }
																			       }
																		       }
																	      } 
																        } 
															       }  
														        }
													      }
													  }
											     }
											}
			                           }
							       }
						       }
    			          }
			    	 }
			      }
    		   }
    	 }
        
	  if (lout.is_open())
	      lout.close();  

     GlobalVariables::GetInstance()->SetNo(no); 
     GlobalVariables::GetInstance()->FillItitle(ititle);
     GlobalVariables::GetInstance()->FillIdata(idata);   
     GlobalVariables::GetInstance()->FillData(data);
	  }
	

	
void DRAGONDetectorConstruction::ugeom_setup()
     {
	  //CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      //C                                                                      C
      //C     Subroutine which initializes the necessary ugeom COMMON          C
      //C     blocks from the RAYTRACE input file                              C
      //C                                                                      C
      //CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	
      G4int ielement;

      G4int ncol, mcol, jcol[2], last_type, ntest, nmcp;

      G4double d[3], p[3], rmat[10];
      G4double xcol[2], ycol[2], dxcol[2], dycol[2];

      G4double par[5];
      
      rmat[0] = 1.;
      rmat[4] = 1.;
      rmat[8] = 1.;
      
      ndipole = 0;
      nedipol = 0;
      nmpole  = 0;
      ntest = 0;
      nmcp = 0;
      nsole = 0;
      ncol = 1;
     
      for(ielement = 0;ielement < no;ielement++)
         { 
		  if(idata[ielement] == 21)
		    {
		     //C          FCUP                    TYPE = 21
     	     p[0] = -495.041718;
		     p[1] = 0.;
		     p[2] = 668.750061; 
	         
		     ugeo_fcup(p);
	         }
	      if(idata[ielement] == 23)
	        {
	         //C          TUBS                  TYPE = 23
             par[0] = data[0][22];
             par[1] = data[1][22];
             par[2] = data[2][22];
             par[3] = data[3][22];
             par[4] = data[4][22];
	         }
	      if(idata[ielement] == 18)
	        { 
			 //C          ENDV                   TYPE =  18
		 
			 p[0] = -1025.10291;
		     p[1] = 0.;
		     p[2] = -71.9687195;
		     
			 ugeo_end(p);
		     }
		  else
		     {
		      if(idata[ielement] == 24)
		        {
		         //C          MCP                   TYPE =  24
				 
	             if(nmcp == 0)
	               {
	                p[0] = -1025.10291;
		            p[1] = 0.;
		            p[2] = 8.79128265;
			        }
			     if(nmcp == 1)
	               {
	                p[0] = -1025.10291;
		            p[1] = 0.;
		            p[2] = -50.2087173;
			        }
	             
	             G4double mcp_data[] = {data[0][ielement],data[1][ielement]};       
	             ugeo_mcp(p,mcp_data,ititle[ielement],nmcp);
	             nmcp = nmcp + 1;
	             }
	           else
	              { 
				   if(idata[ielement] == 22)
		             { 
		              //C          TEST                   TYPE =  22
		              G4double rot_angles[3];
		              ntest = ntest + 1;
            
		              if(ntest == 1)
		                {
		                 p[0] = -62.363575;
		                 p[1] = 0.;
		                 p[2] = 354.0466;
		                 
		                 rot_angles[0] = 0.;
		                 rot_angles[1] = 50.;
		                 rot_angles[2] = 0.;
		                 
		                 ugeo_test(ntest, p, rot_angles);
				         }
				      if(ntest == 2)
		                {
		                 p[0] = -501.589172;
		                 p[1] = 0.;
		                 p[2] = 658.362976;
		                 
		                 rot_angles[0] = 0.;
		                 rot_angles[1] = 70.;
		                 rot_angles[2] = 0.;
		                 
		                 ugeo_test(ntest, p, rot_angles);		                 
				         }
				      if(ntest == 3)
		                {
		                 p[0] = -1025.10291;          
    	                 p[1] = 0.;
		                 p[2] = -71.9187164;
		                 
		                 ugeo_test(ntest, p, rot_angles);		                 
		                 }
		               } 
		           else
		              {
					   if(idata[ielement] == 19)
						 {
					      //C          STRV                   TYPE =  19
					   
					      p[0] = 0.0;          
    	                  p[1] = 0.0;
		                  p[2] = 0.0;
					     
					      ugeo_start(p);
					      }
					    else
					       {
							if(idata[ielement] == 12)
							  {
							   //C          DRIFT                  TYPE = 12
					  		   
					  		   d[0] = 0.0;
		                       d[1] = 0.0;
		                       d[2] = data[0][ielement];
		               	       }
		               	     else
		               	        {
								 if(idata[ielement] == 11)
								   {
									//C          SHIFT AND ROTATE       TYPE = 11
									
									d[0] = data[0][ielement]; 
		                            d[1] = data[1][ielement]; 
		                            d[2] = data[2][ielement]; 

								    if(data[3][ielement] == 0)   //! Rotation around x-axis
								      {
									   G4RotationMatrix *irotmat = new G4RotationMatrix();
									   irotmat->rotateX(-data[3][ielement]*deg);	
									   }									
									if(data[4][ielement] == 0)   //! Rotation around y-axis
								      {
									   G4RotationMatrix *irotmat = new G4RotationMatrix();
									   irotmat->rotateY(data[4][ielement]*deg);	
									   }   
									if(data[5][ielement] == 0)   //! Rotation around z-axis
								      {
									   G4RotationMatrix *irotmat = new G4RotationMatrix();
									   irotmat->rotateZ(data[5][ielement]*deg);	
									   }     
									}
							      else
							         {
								      if(idata[ielement] == 13)
								        {
										 //C          COLLIMATOR             TYPE = 13
									
										 ncol = ncol + 1;        
                                         mcol = ncol % 2;
										 if(mcol == 1)mcol = 0;     
										 if(mcol == 0)mcol = 1;
										 
										 jcol[mcol] = data[0][ielement];
                                         xcol[mcol] = data[1][ielement];
                                         ycol[mcol] = data[2][ielement];
                                         dxcol[mcol] = data[3][ielement];
                                         dycol[mcol] = data[4][ielement];
										 	  
									     if(mcol == 1)
									       {
										    if(last_type == 2)
										      {
                                               jcol_dipole[1][ndipole] =  jcol[mcol];
                                               xcol_dipole[1][ndipole] =  xcol[mcol];
                                               ycol_dipole[1][ndipole] =  ycol[mcol];
                                               dxcol_dipole[1][ndipole] = dxcol[mcol];
                                               dycol_dipole[1][ndipole] = dycol[mcol];
									           }
									        else
									           {
											    if(last_type == 7)
											  	  {
											       jcol_edipol[1][nedipol] =  jcol[mcol];
                                                   xcol_edipol[1][nedipol] =  xcol[mcol];
                                                   ycol_edipol[1][nedipol] =  ycol[mcol];
                                                   dxcol_edipol[1][nedipol] = dxcol[mcol];
                                                   dycol_edipol[1][nedipol] = dycol[mcol];
										           }
										        else
										           {
											  	    if(last_type == 9)
												      {
											           jcol_mpole[1][nmpole] =  jcol[mcol];
                                                       xcol_mpole[1][nmpole] =  xcol[mcol];
                                                       ycol_mpole[1][nmpole] =  ycol[mcol];
                                                       dxcol_mpole[1][nmpole] = dxcol[mcol];
                                                       dycol_mpole[1][nmpole] = dycol[mcol];
										               } 
										            else
										               {
													    if(last_type == 14)
												          {
											               jcol_sole[1][nsole] =  jcol[mcol];
                                                           xcol_sole[1][nsole] =  xcol[mcol];
                                                           ycol_sole[1][nsole] =  ycol[mcol];
                                                           dxcol_sole[1][nsole] = dxcol[mcol];
                                                           dycol_sole[1][nsole] = dycol[mcol];
										                   }  	 
													    }
												    }
										        }
									        }
									     }
									  else
									     {
										  if(idata[ielement] == 17)
										    {
											 //C          REAL COLLIMATOR        TYPE = 17

											 G4double col_data[] = {data[0][ielement],data[1][ielement],data[2][ielement],data[3][ielement],data[4][ielement],data[5][ielement]};  
											                             
											 ugeo_col(col_pos,col_rot_angles,col_data,ititle[ielement]);       
											 }	  
										  else
											 { 
											  if(idata[ielement] == 2 || idata[ielement] == 20)
											    {
												 //C          DIPOLE  LENS           TYPE =  2
                                                 //C          TRIUMF SASP DIPOLE     TYPE = 20
												  
												 ndipole = ndipole + 1;
                                                 ndipole = std::min(ndipole,max_dipole);
                                                 
                                                 if(mcol == 1)
                                                   {
                                                    jcol_dipole[0][ndipole] =  jcol[mcol];
                                                    xcol_dipole[0][ndipole] =  xcol[mcol];
                                                    ycol_dipole[0][ndipole] =  ycol[mcol];
                                                    dxcol_dipole[0][ndipole] = dxcol[mcol];
                                                    dycol_dipole[0][ndipole] = dycol[mcol];
                                                    }
           
                                                 last_type = idata[ielement];

                                                 gap_dipole[ndipole] = data[12][ielement];
                                                 phi_dipole[ndipole] = data[15][ielement];
                                                 r_dipole[ndipole] = data[13][ielement];
                                                 dr_dipole[ndipole]  = (data[48][ielement]+data[49][ielement])/2.;
                                                 alpha_dipole[ndipole] = data[16][ielement];
                                                 beta_dipole[ndipole] = data[17][ielement];
                                                 z11_dipole[ndipole] = data[24][ielement];
                                                 z22_dipole[ndipole] = data[27][ielement];

                                                 ugeo_dipole(ndipole);
											     }
											  else
											     {
												  if(idata[ielement] == 7)  
												    {
												     //C          ELECTROSTATIC DEFLECTOR  TYPE = 7
												     nedipol = nedipol + 1;
                                                     nedipol = std::min(nedipol,max_edipol);
                                                     
                                                     if(mcol == 1)
                                                       {
                                                        jcol_edipol[0][nedipol] =  jcol[mcol];
                                                        xcol_edipol[0][nedipol] =  xcol[mcol];
                                                        ycol_edipol[0][nedipol] =  ycol[mcol];
                                                        dxcol_edipol[0][nedipol] = dxcol[mcol];
                                                        dycol_edipol[0][nedipol] = dycol[mcol];
                                                        }
												     last_type = idata[ielement];
												  
											         gap_edipol[nedipol] = data[12][ielement];
                                                     phi_edipol[nedipol] = data[15][ielement];
                                                     r_edipol[nedipol] = data[13][ielement];
                                                     dr_edipol[nedipol]  = (data[1][ielement]+data[2][ielement])/2.;
                                                     z11_edipol[nedipol] = data[24][ielement];
                                                     z22_edipol[nedipol] = data[27][ielement];
                                                     
                                                     ugeo_edipol(nedipol);
											         }
											      else
											         {
													  if(idata[ielement] == 9)
													    {
													     //C          MULTIPOLE (POLES)      TYPE =  9
                                                         nmpole = nmpole + 1;
                                                         nmpole = std::min(nmpole,max_mpole);

														 if(mcol == 1)
														   {
														    jcol_mpole[0][nmpole] =  jcol[mcol];
                                                            xcol_mpole[0][nmpole] =  xcol[mcol];
                                                            ycol_mpole[0][nmpole] =  ycol[mcol];
                                                            dxcol_mpole[0][nmpole] = dxcol[mcol];
                                                            dycol_mpole[0][nmpole] = dycol[mcol];
														    }
                                                         last_type = idata[ielement];
                                                         
                                                         efblength_mpole[nmpole] = data[11][ielement];
                                                         r_mpole[nmpole] = data[12][ielement];
                                                         z11_mpole[nmpole] = data[18][ielement];
                                                         z22_mpole[nmpole] = data[21][ielement];
														  
														 ugeo_mpole(nmpole,mpole_rot_angles); 
														 }
												       else
														  {
														   if(idata[ielement] == 14)
														     {
													          //C          SOLENOID               TYPE = 14
													          nsole = nsole + 1;
                                                              nsole = std::min(nsole,max_sole);
													          
													          if(mcol == 1)
			  										            {
                                                                 jcol_sole[0][nsole] =  jcol[mcol];
                                                                 xcol_sole[0][nsole] =  xcol[mcol];
                                                                 ycol_sole[0][nsole] =  ycol[mcol];
                                                                 dxcol_sole[0][nsole] = dxcol[mcol];
                                                                 dycol_sole[0][nsole] = dycol[mcol];
                                                                 }
                                                              last_type = idata[ielement];
                                                              
                                                              efblength_sole[nsole] = data[11][ielement];
                                                              r_sole[nsole] = data[12][ielement]/2.;
                                                              z11_sole[nsole] = data[14][ielement];
                                                              z22_sole[nsole] = data[15][ielement];
													          
													          ugeo_sole(nsole);          		   
													 	      }
														   }
												      }
												  }	
											  }
										  }
									     
									  }
							     }    
							}
		              }
		           }
	          }
          }
      }	 
		 
	}	 
		 
		 
		 
	

//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software 
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
// 
// If you use this example, please cite the following publication:
// Rad. Prot. Dos. 133 (2009) 2-11

#ifndef DRAGONEMField_h
#define DRAGONEMField_h 1

#include "G4ElectroMagneticField.hh"
#include "G4TransportationManager.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4SystemOfUnits.hh"

#include "global_variables.hh"
#include <vector>

class G4MagIntegratorStepper;
class G4EqMagElectricField;

namespace DRAGON
{

class DRAGONEMFieldMessenger;

class DRAGONEMField : public G4ElectroMagneticField
{  
  public:
    DRAGONEMField();
    ~DRAGONEMField();
  
    void GetFieldValue( const  double Point[4], double *Bfield ) const;
    
	G4double GetBScale() const {return bscale;}
	G4double GetEscale() const {return escale;}
	
	G4double GetfMinStep() const {return fMinStep;}
	G4double GetfDeltaChord() const {return fDeltaChord;}
	G4double GetfDeltaOneStep() const {return fDeltaOneStep;}
	G4double GetfDeltaIntersection() const {return fDeltaIntersection;}
	G4double GetfEpsMin() const {return fEpsMin;}
	G4double GetfEpsMax() const {return fEpsMax;}
        G4MagIntegratorStepper* GetStepper() const {return fStepper;}
        G4EqMagElectricField* GetEquation() const {return fEquation;}

        void SetStepperType(G4int StepperType) { fStepperType = StepperType; }
        void SetStepper();
        void SetMinStep(G4double MinStep) { fMinStep = MinStep; }
        void SetDeltaChord(G4double DeltaChord) { fDeltaChord = DeltaChord; }
        void SetDeltaOneStep(G4double DeltaOneStep) { fDeltaOneStep = DeltaOneStep; }
        void SetDeltaIntersection(G4double DeltaIntersection) { fDeltaIntersection = DeltaIntersection; }
        void SetEpsMin(G4double EpsMin) { fEpsMin = EpsMin; }
        void SetEpsMax(G4double EpsMax) { fEpsMax = EpsMax; }

    G4bool DoesFieldChangeEnergy() const {return false;}
    void SetSCAL(const G4String& input);
    
    //guefld.cc
  	void guefld(G4double vector[7], G4double time, G4double* bfield, G4double* efield) const;
  	
  	//mitray_field.cc
  	void mitray_field(G4String devname, G4double xpos[3], G4double bfld[3], G4double efld[3]) const;
 	
 	//mitray_solnd.cc
 	void mitray_solnd(G4double DATA[75], G4double XPOS[3], G4double BFLD[3]) ;  
 	
 	//mitray_sasp.cc
 	void mitray_sasp(G4double data[75], G4double xpos[3], G4double bfld[3]);	

  	//mitray_poles.cc
	void mitray_poles(G4double* DATA, G4double* XPOS, G4double* BFLD) const;
  	void mitray_bpoles() const;
  	void mitray_bpls(G4int IGP, G4double D, G4double S, G4double& RE, G4double& G1, G4double& G2, G4double& G3, G4double& G4, G4double& G5, G4double& G6) const;
   
        //mitray_zone.cc
        void mitray_zone(G4double ZB, G4double ZC, G4double Z11, G4double Z12, G4double Z21, G4double Z22, G4int& IZONE) const;
           
   	//mitray_dipole.cc
  	void mitray_dipo_AtoB(G4double XA, G4double YA, G4double ZA, G4double& XB, G4double& YB, G4double& ZB) const;
  	void mitray_dipo_BtoC(G4double XB, G4double YB, G4double ZB, G4double& XC, G4double& YC, G4double& ZC) const;
  	void mitray_dipo_BBtoBA(G4double BXB, G4double BYB, G4double BZB, G4double& BXA, G4double& BYA, G4double& BZA) const;
  	void mitray_dipo_BCtoBB(G4double BXC, G4double BYC, G4double BZC, G4double& BXB,  G4double& BYB, G4double& BZB) const;
  	void mitray_dipole(G4double* DATA, G4double* XPOS, G4double* BFLD) const;
  	void mitray_bdip() const;
  	void mitray_bdpp(G4double& BFLD, G4double Z, G4double X, G4bool BDPPXflag = false) const;
  	void mitray_bdppx() const;
  	void mitray_ndip() const;
  	void mitray_ndpp(G4double& BFLD, G4double Z, G4double X, G4double DR, G4bool NDPPXflag = false) const;
  	void mitray_bpretz() const;
  	void mitray_sdip(G4double X, G4double Z, G4bool SIJflag = false) const;
  	void mitray_bdmp(G4double& BZZ, G4double Z, G4double X) const;
	void mitray_fmap();
	void mitray_dmap(G4int II);
	G4double xmitray_zefb(G4double XP) const;
	G4double xmitray_dzdx(G4double XP) const;
	G4double xmitray_dzdx2(G4double XP) const;
	
	//mitray_edipol.cc
	void mitray_edipol(G4double* DATA, G4double* XPOS, G4double* EFLD) const;
	void mitray_edip() const;
	void mitray_edpp(G4double D, G4double S, G4double& RE, G4double& G1, G4double& G2, G4double& G3, G4double& G4, G4double& G5, G4double& G6) const;
	void mitray_edip_ECtoEA(G4double EXC, G4double EYC, G4double EZC, G4double& EXA, G4double& EYA, G4double& EZA) const;
  
        void mitray_sasp(G4double* DATA, G4double* XPOS, G4double* BFLD) const;
        void mitray_saspratio(G4double BMEAS0, G4double XA, G4double& RATIO) const; 
        
        void mitray_solnd(G4double* DATA, G4double* XPOS, G4double* BFLD) const;
        void mitray_bsol() const;
        void mitray_FB02AD(G4double CAYSQ, G4double SINP, G4double COSP, G4double& E, G4double& F) const; 
        void mitray_FB03AD(G4double GN, G4double CACA, G4double& P) const;
        void mitray_FB01AD(G4double C, G4double& VK, G4double& VE) const;

  private:
        //mutable G4int irot_dipole[2],irot_edipol[2];
  
        mutable G4double GRAD1,GRAD2,GRAD3,GRAD4,GRAD5;
      
        mutable G4double DH, DO, DD, DDD, DSH, DSO, DSD, DSDD;
  
  	const G4double deg_rad = 180. / CLHEP::pi;
  	
        G4double bscale, escale;
  	
  	//from COMMON MITRAY10
  	mutable G4double BX, BY, BZ, K, TC[6], DTC[6];
  	
  	//from COMMON MITRAY11
  	mutable G4double EX, EY, EZ, QMC, IVEC;
  	
  	//from COMMON MITRAY20
  	mutable G4double EC2, EC4, WE, WC, NDX, BET1, GAMA, DELT;
  	
  	//from COMMON MITRAY21
  	mutable G4double RCA, DELS, BR, S2, S3, S4, S5, S6, S7, S8;
  	
  	//from COMMON MITRAY22
  	mutable G4double D, DG, S, BF, BT, EF, ET, DUM, WDIP;
  	
  	//from COMMON MITRAY23
        mutable G4double C0, C1, C2, C3, C4, C5;
     
        //from COMMON MITRAY24 
        mutable G4double RB, XC_OFFSET, ZC_OFFSET;
    
        //from COMMON MITRAY25
        mutable G4int IN, MTYP, NSRF, IR, IDUM1, IDUM2, IMAP;
    
        //from COMMON MITRAY26  
        G4double JMAP[5], IX, IZ, IDUM, BZMAP[101][101][2][5];
    
        //from COMMON MITRAY_AXES
        mutable G4double XA, YA, ZA, XB, YB, ZB, XC, YC, ZC;
    
        //from COMMON MITRAY_BOUNDS
        G4double XBMIN, XBMAX, XCMIN, XCMAX;
     
        //from COMMON MITRAY_EDIPO
        mutable G4double A, B, PHI;
  	
  	mutable G4bool ldiag;
  	  	
  	mutable G4double RE, G1, G2, G3, G4, G5, G6;
  	
  	//from COMMON MITRAYBLSDIP 
  	mutable G4double XO, ZO, SS, DCS, DSN;
  	
  	//from COMMON MITRAY_DIPO
  	mutable G4double ALPHA, BETA, XCR1, XCR2;
  	
  	mutable G4double B0, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12;
  	
  	mutable G4int IZONE;
  	
  	mutable G4double AL, RAD;
  	
  	mutable G4int ifield;
  	
        mutable G4double theta[3],phi[3];
  	
        struct DRAGONEMFilter
               {
                G4String type;
                G4String name;  
                G4int multiplicity;
                G4int manyORonly; 
                G4int ifld;
                G4int irot;  
                G4double theta1,theta2,theta3;
                G4double phi1,phi2,phi3;         
                
                DRAGONEMFilter() : type(" "), name(""), multiplicity(0), manyORonly(1), ifld(1), irot(0), 
                                   theta1(90.0), phi1(0.0), theta2(90.0), phi2(90.0), theta3(90.0), phi3(0.0) {}
                DRAGONEMFilter(G4String j, G4String n, G4int i, G4int f, G4int o, G4int iirot, G4double theta_1, G4double phi_1, 
                               G4double theta_2, G4double phi_2, G4double theta_3, G4double phi_3)
            : type(j), name(n), multiplicity(i), manyORonly(f), ifld(o), irot(iirot),theta1(theta_1),phi1(phi_1),
              theta2(theta_2),phi2(phi_2),theta3(theta_3),phi3(phi_3){}


    };  
    
           G4Navigator* fNavigator;
           G4EqMagElectricField* fEquation = nullptr;
           G4MagIntegratorStepper* fStepper = nullptr;
           G4int fStepperType = 8; 
           G4double fMinStep = 0.001 * cm;
           G4double fDeltaChord = 0.3 * cm;
           G4double fDeltaOneStep = 0.001 * cm;
           G4double fDeltaIntersection = 0.01 * cm;
           G4double fEpsMin = 2.5e-7;
           G4double fEpsMax = 0.001; 

      public:
        //DRAGONEMFilter DRAGONEMelements[18]; 
        //DRAGONEMFilter* GetAllDRAGONEMFilters() {return DRAGONEMelements;}    
        std::vector<DRAGONEMFilter> DRAGONEMelements;
        const std::vector<DRAGONEMFilter>& GetAllDRAGONEMFilters() const { return DRAGONEMelements;}
                
        DRAGONEMFieldMessenger* fEMFieldMessenger = nullptr;
   
  	
};

}
#endif

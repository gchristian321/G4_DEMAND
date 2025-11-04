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
//
// Based on purging magnet advanced example.
//

#include "DRAGONEMField.hh"
#include "DRAGONEMFieldMessenger.hh"
#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4CashKarpRKF45.hh"
#include "G4EqMagElectricField.hh"

#include <fstream>



namespace DRAGON
{

DRAGONEMField::DRAGONEMField()
{
 fEMFieldMessenger = new DRAGONEMFieldMessenger(this);

 std::ifstream inputFile("DRAGONEMFilters.data");
 std::string line;
 G4String type, name;
 G4int multiplicity, manyORonly, ifld, irot;
 G4double theta1, theta2, theta3, phi1, phi2, phi3;
       
 if (inputFile.is_open()) 
    {
     while (std::getline(inputFile, line)) 
           {
            std::istringstream ss(line);
            ss >> type >> name >> multiplicity >> manyORonly >> ifld >> irot >> theta1 >> phi1 >> theta2 >> phi2 >> theta3 >> phi3;
            DRAGONEMelements.push_back(DRAGONEMFilter(type, name, multiplicity, manyORonly, ifld, irot, theta1, phi1, theta2, phi2, theta3, phi3));         
            //DRAGONEMelements[index] = DRAGONEMFilter(type, name, multiplicity, manyORonly, ifld, irot, theta1, phi1, theta2, phi2, theta3, phi3);
            }
     inputFile.close();
     }
 else
    {
     std::cout << "UNABLE TO READ DRAGONEMFilters.data file" << std::endl;
     exit(EXIT_FAILURE); 
     }

 G4Navigator* theNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
 if (theNavigator->GetWorldVolume())
    {
     fNavigator = new G4Navigator();
     fNavigator->SetWorldVolume(theNavigator->GetWorldVolume());
     }   
}

DRAGONEMField::~DRAGONEMField() 
{
 delete fEMFieldMessenger;
 delete fNavigator;
 delete fStepper;
 delete fEquation;
 }

void DRAGONEMField::GetFieldValue(const double point[4], double *Bfield ) const 
     { 
//c     Initialize the field values

// Magnetic field
      Bfield[0] = Bfield[1] = Bfield[2] = 0.;
// Electric field
      Bfield[3] = Bfield[4] = Bfield[5] = 0.;
 
      double bfield[3] = {Bfield[0], Bfield[1], Bfield[2]};
      double efield[3] = {Bfield[3], Bfield[4], Bfield[5]};
// Position from mm (Geant4) to cm (Geant3)
      double vector[7] = {point[0]/10., point[1]/10., point[2]/10., 0.0, 0.0, 0.0, 0.0};
 
      std::cout << "DENTROGETVALUECOOOOOOOOOOOGEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE" << std::endl;
      guefld(vector, 0., bfield, efield);
      std::cout << "FUERAGETVALUECOOOOOOOOOOOGEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE" << std::endl;
      
      Bfield[0] = bfield[0]*0.00017;               //OJO 
      Bfield[1] = bfield[1]*0.00017; 
      Bfield[2] = bfield[2]*0.00017; 
// Electric field from V/cm (Geant3) to V/m (Geant4)
      Bfield[3] = efield[0]*100.;
      Bfield[4] = efield[1]*100.;
      Bfield[5] = efield[2]*100.;

      std::cout << "FINAL Bfield " << Bfield[0] << " " << Bfield[1] << " " << Bfield[2] << std::endl;  
      std::cout << "FINAL Efield " << Bfield[3] << " " << Bfield[4] << " " << Bfield[5] << std::endl;       
      std::cout << "=================" << std::endl; 
      }

void DRAGONEMField::SetSCAL(const G4String& input) 
     {
      std::istringstream iss(input);
      iss >> bscale;
      iss >> escale;
      }

void DRAGONEMField::SetStepper()
{
  delete fStepper;
  delete fEquation;
  
  fEquation = new G4EqMagElectricField(this);

  switch (fStepperType) {
    case 0:
      //      fStepper = new G4ExplicitEuler( fEquation, 8 ); // no spin tracking
      fStepper = new G4ExplicitEuler(fEquation, 12);  // with spin tracking
      G4cout << "G4ExplicitEuler is called" << G4endl;
      break;
    case 1:
      //      fStepper = new G4ImplicitEuler( fEquation, 8 ); // no spin tracking
      fStepper = new G4ImplicitEuler(fEquation, 12);  // with spin tracking
      G4cout << "G4ImplicitEuler is called" << G4endl;
      break;
    case 2:
      //      fStepper = new G4SimpleRunge( fEquation, 8 ); // no spin tracking
      fStepper = new G4SimpleRunge(fEquation, 12);  // with spin tracking
      G4cout << "G4SimpleRunge is called" << G4endl;
      break;
    case 3:
      //      fStepper = new G4SimpleHeum( fEquation, 8 ); // no spin tracking
      fStepper = new G4SimpleHeum(fEquation, 12);  // with spin tracking
      G4cout << "G4SimpleHeum is called" << G4endl;
      break;
    case 4:
      //      fStepper = new G4ClassicalRK4( fEquation, 8 ); // no spin tracking
      fStepper = new G4ClassicalRK4(fEquation, 12);  // with spin tracking
      G4cout << "G4ClassicalRK4 (default) is called" << G4endl;
      break;
    case 5:
      //      fStepper = new G4CashKarpRKF45( fEquation, 8 ); // no spin tracking
      fStepper = new G4CashKarpRKF45(fEquation, 12);  // with spin tracking
      G4cout << "G4CashKarpRKF45 is called" << G4endl;
      break;
    default:
      fStepper = nullptr;
  }
}

}

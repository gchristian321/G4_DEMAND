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
/// \file hadronic/Hadr01/src/PhysicsList.cc
/// \brief Implementation of the PhysicsList class
//
//
//
/////////////////////////////////////////////////////////////////////////
//
// DRAGONPhysicsList
//
// Created: 31.04.2006 V.Ivanchenko
//
// Modified:
// 04.06.2006 Adoptation of Hadr01 (V.Ivanchenko)
// 26.04.2007 Physics according to 8.3 Physics List (V.Ivanchenko)
// 16.10.2012 Renamed used classes (A.Ribon)
//
////////////////////////////////////////////////////////////////////////
// 
#include "G4EmPenelopePhysics.hh"                //Geant4
#include "G4EmLivermorePhysics.hh"
#include "G4EmStandardPhysicsGS.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"

#include "G4DecayPhysics.hh"
#include "G4IonElasticPhysics.hh"
#include "G4HadronPhysicsFTFP_BERT.hh"
#include "G4HadronPhysicsFTFP_BERT_HP.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4ProcessManager.hh"

#include "G4StepLimiterPhysics.hh"
#include "G4OpticalPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4Ions.hh"

#include "G4EmParameters.hh"
#include "G4EmConfigurator.hh"
#include "G4LossTableManager.hh"
#include "G4BraggIonGasModel.hh"
#include "G4BetheBlochIonGasModel.hh"
#include "G4IonFluctuations.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4UniversalFluctuation.hh"

#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "G4DecayTable.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4Threading.hh"

#include "G4PhysicalConstants.hh"
#include "G4PhysicsListHelper.hh"
#include "G4hMultipleScattering.hh"
#include "G4ionIonisation.hh"
#include "G4NuclearStopping.hh"
#include "G4SystemOfUnits.hh"

#include "DRAGONDetectorConstruction.hh"         //Local  
#include "DRAGONPhysicsList.hh"
#include "DRAGONPhysicsListMessenger.hh"
#include "DRAGONIon.hh"

#include <regex>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
namespace DRAGON
{

DRAGONPhysicsList::DRAGONPhysicsList(DRAGONDetectorConstruction* det) : G4VModularPhysicsList(),
DRAGON_det(det),CARDNAME("c12ag.dat"),alpha(false),beamenerg(0)
{
  SetDefaultCutValue(0.7*CLHEP::mm);
  uvinit();
 
  verboseLevel = 1;

  // Messengers
  fMessenger = new DRAGONPhysicsListMessenger(this);

  // EM physics
  fEmPhysicsList = new G4EmStandardPhysics(verboseLevel);
  
  // Step limiter physics (limits max step size per volume or particle)
  fStepLimiterPhysics = new G4StepLimiterPhysics();
 
  // Particles
  fDecPhysicsList = new G4DecayPhysics(verboseLevel);
  
  // Radioactive decay physics
  //fRadDecPhysicsList = new G4RadioactiveDecayPhysics(verboseLevel);
 
  // Heavy ion physics
  fIonPhysicsList = new G4IonPhysics(verboseLevel);
  
  // Optical physics (Cherenkov, scintillation, absorption, etc.)
//  fOptPhysicsList = new G4OpticalPhysics(verboseLevel);   
  // (Opcional) configura efectos ópticos si quieres
  // fOptPhysicsList->SetScintillationYieldFactor(1.0);
  // fOptPhysicsList->SetTrackSecondariesFirst(kCerenkov, true);
  // fOptPhysicsList->SetTrackSecondariesFirst(kScintillation, true);


 
  // Hadronic physics
  fHadElastPhysics =  new G4HadronElasticPhysics(verboseLevel);        // Elastic scattering
  //fHadPhysicsFTFP_BERT = new G4HadronPhysicsFTFP_BERT(verboseLevel);   // Inelastic interactions
  fHadPhysicsFTFP_BERT = new G4HadronPhysicsFTFP_BERT_HP(verboseLevel);   // Inelastic interactions
	
  // Specialized ion gas models for accurate dE/dx of heavy ions in gas
  AddIonGasModels(); 

  lkine     =   2;           //! full 15O(alpha,g)19Ne simulation
  fkine[0] =   8.0;          //! effective charge of the BEAM in vacuum
  fkine[1] =   4.0;          //! effective charge of the PRODUCT in vacuum
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

DRAGONPhysicsList::~DRAGONPhysicsList()
{
  delete fMessenger;
  delete fDecPhysicsList;
  //delete fRadDecPhysicsList;
  delete fEmPhysicsList;
  delete fIonPhysicsList;
//  delete fOptPhysicsList;
  delete fStepLimiterPhysics;
  delete fHadElastPhysics;
  delete fHadPhysicsFTFP_BERT;
  }


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DRAGONPhysicsList::ConstructParticle()
{
   std::cout << "-FISICA-Particulas" << std::endl;

    // Ahora construimos las partículas estándar
    fEmPhysicsList->ConstructParticle();
    fDecPhysicsList->ConstructParticle();
    fIonPhysicsList->ConstructParticle();
//   fOptPhysicsList->ConstructParticle();
    fHadElastPhysics->ConstructParticle();
    fHadPhysicsFTFP_BERT->ConstructParticle();
#if 0
    if (G4Threading::IsMasterThread()) {
        static G4bool isDefined = false;
        if (!isDefined) {
            G4ParticleDefinition* ion = DefineCo60cascade();
			PrintIonProperties(ion);
            isDefined = true;
        }
    }
#endif

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DRAGONPhysicsList::ConstructProcess()
{
    std::cout << "-FISICA-Procesos" << std::endl;

    AddTransportation();

    fEmPhysicsList->ConstructProcess();
    fDecPhysicsList->ConstructProcess();
    fIonPhysicsList->ConstructProcess();
//    fOptPhysicsList->ConstructProcess();
    fStepLimiterPhysics->ConstructProcess();
    fHadElastPhysics->ConstructProcess();
    fHadPhysicsFTFP_BERT->ConstructProcess();
	
    AddIonGasModels();   

    G4ParticleDefinition* ionCo60copia = G4ParticleTable::GetParticleTable()->FindParticle("Co60copia");
    if (ionCo60copia) {
        std::cout << "IOIOIOIOIO" << std::endl;
        G4ProcessManager* pManager = ionCo60copia->GetProcessManager();
        if (pManager) {
            std::cout << "1212121212" << std::endl;
            G4Decay* decay = new G4Decay();
            pManager->AddProcess(decay);
            pManager->SetProcessOrdering(decay, idxPostStep);
            pManager->SetProcessOrdering(decay, idxAtRest);
        }
    } else {
        std::cout << "TTRTRTRTR" << std::endl;
    }
	
}


G4ParticleDefinition* DRAGONPhysicsList::DefineCo60cascade()
{
    G4ParticleDefinition* existing = G4ParticleTable::GetParticleTable()->FindParticle("Co60copia");
    if (existing) return existing;	
	
    // Masa aproximada
    G4double mass = 59.933822 * g / mole * amu_c2;

    // --- Estado base (Co60 estable) ---
    G4ParticleDefinition* Co60 = new G4ParticleDefinition(
        "Co60", mass, 0.0, 0.0, 0, +1, 0, 0, 0, 0,
        "nucleus", 0, +1, 1000270600,
        true, 0.0, nullptr, false, "nucleus", 0, 0.0);

    // --- Estado excitado intermedio a 1173 keV ---
    G4ParticleDefinition* Co60_1173 = new G4ParticleDefinition(
        "Co60_1173", mass + 1173.*keV, 0.0, 0.0, 0, +1, 0, 0, 0, 0,
        "nucleus", 0, +1, 1000270601,  // diferente encoding
        false, 0.0, nullptr, false, "nucleus", 0, 1173.*keV);

    // --- Estado inicial "copia" ---
    G4ParticleDefinition* Co60copia = new G4ParticleDefinition(
        "Co60copia", mass + (1173.+1332.)*keV, 0.0, 0.0, 0, +1, 0, 0, 0, 0,
        "nucleus", 0, +1, 1000270602,
        false, 0.0, nullptr, false, "nucleus", 0, (1173.+1332.)*keV);

    // --- Tabla de decaimiento de Co60copia ---
    G4DecayTable* table_copia = new G4DecayTable();
    G4VDecayChannel* decay1 = new G4PhaseSpaceDecayChannel("Co60copia", 1.0, 2, "Co60_1173", "gamma");
    table_copia->Insert(decay1);
    Co60copia->SetDecayTable(table_copia);

    // --- Tabla de decaimiento de Co60_1173 ---
    G4DecayTable* table_1173 = new G4DecayTable();
    G4VDecayChannel* decay2 = new G4PhaseSpaceDecayChannel("Co60_1173", 1.0, 2, "Co60", "gamma");
    table_1173->Insert(decay2);
    Co60_1173->SetDecayTable(table_1173);
    
   G4IonTable* ionTable = G4IonTable::GetIonTable();
   ionTable->Insert(Co60copia);

    return Co60copia;
}



void DRAGONPhysicsList::SetCuts()
{
  std::cout << "-FISICA-Cortes" << std::endl;
  SetCutValue(cutgam, "gamma");
  SetCutValue(cutele, "e-");
  SetCutValue(cutele, "e+");

  DumpCutValuesTable();   
}

void DRAGONPhysicsList::uvinit()
{
  //From uvinit.f
  //C.
  //C.    *************** RUN DEFAULTS ***************
  //C.
  cutgam = defaultCutValue;     //! Photon energy tracking cut
  cutele = defaultCutValue;     //! Electron kin. energy tracking cut
}

void DRAGONPhysicsList::SetFKIN(const G4String& input) {
    std::istringstream iss(input);
    iss >> lkine;
    iss >> fkine[0];
    iss >> fkine[1];
}

void DRAGONPhysicsList::SetTUNE(const G4String& input) {
    std::istringstream iss(input);
    iss >> refenerg;
    iss >> refatno;
    iss >> refq; 
}


void DRAGONPhysicsList::SetCUTs(const G4String& input) {
    std::istringstream iss(input);
    iss >> cutgam >> cutele;
}


void DRAGONPhysicsList::SetMTUN(const G4String& input) {
    std::istringstream iss(input);
    for (G4int i = 0; i < 5; ++i) {
        if (!(iss >> offset[i])) {
            std::cout << "Invalid input string in method SetMTUN. Expected 5 numbers.";
        }
    }
    iss >> energscale;
}

void DRAGONPhysicsList::SetInputCardname(const G4String& input) {
    std::istringstream iss(input);
    iss >> CARDNAME;
}

void DRAGONPhysicsList::ReadUserReactionNameList()
{
    if (CARDNAME == " ") {
        CARDNAME = "c12ag.dat";
    }

    std::ifstream inFile(CARDNAME);
    if (!inFile.is_open()) {
        std::cerr << "No input reaction card! (LKINE set to 20)" << std::endl;
        return;
    }

    std::string line;
    std::regex string_or_number_pattern(R"(^\s*([a-zA-Z_][a-zA-Z0-9_]*)\s*=\s*([^\s#]+))");
    std::regex array1d_pattern(R"(^\s*([a-zA-Z_][a-zA-Z0-9_]*)\s*\(\s*(\d+)\s*\)\s*=\s*(-?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][-+]?[0-9]+)?))");
    std::regex array2d_pattern(R"(^\s*([a-zA-Z_][a-zA-Z0-9_]*)\((\d+),(\d+)\)\s*=\s*(-?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][-+]?[0-9]+)?)\s*$)");
    std::regex is_number_pattern(R"(^-?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][-+]?[0-9]+)?$)");

    std::cout << "=== Abriendo archivo de entrada: " << CARDNAME << " ===" << std::endl;

    while (std::getline(inFile, line)) {
        if (line.empty() || line[0] == '#' || line[0] == '$') continue;

        std::smatch match;

        // Matriz 2D
        if (std::regex_search(line, match, array2d_pattern)) {
            std::string var = match[1];
            int i = std::stoi(match[2]);
            int j = std::stoi(match[3]);

            if (var == "br" && i < 15 && j < 10) {
                double val = std::stod(match[4]);
                br[i-1][j-1] = val;
                std::cout << "br[" << i-1 << "][" << j-1 << "] = " << br[i-1][j-1] << std::endl;
            }
            else if (var == "md" && i < 15 && j < 10) {
                int val = std::stoi(match[4]);
                md[i-1][j-1] = val;
                std::cout << "md[" << i-1 << "][" << j-1 << "] = " << md[i-1][j-1] << std::endl;
            }
        }

        // Arreglo 1D
        else if (std::regex_search(line, match, array1d_pattern)) {
            std::string var = match[1];
            int idx = std::stoi(match[2]);
            double val = std::stod(match[3]);

            if (var == "level" && idx < 16) level[idx] = val;
            else if (var == "life" && idx < 16) life[idx] = val;
      
        
        }

        // Escalares y strings
        else if (std::regex_search(line, match, string_or_number_pattern)) {
            std::string var = match[1];
            std::string val = match[2];

            if (std::regex_match(val, is_number_pattern)) {
                double dval = std::stod(val);

                if (var == "beam_mass_excess") beam_mass_excess = dval;
                else if (var == "recoil_mass_excess") recoil_mass_excess = dval;
                else if (var == "part_width") part_width = dval;
                else if (var == "gam_width") gam_width = dval;
                else if (var == "spin_stat_fac") spin_stat_fac = dval;
                else if (var == "zbeam") zbeam = dval;
                else if (var == "ztarg") ztarg = dval;
                else if (var == "atarg") atarg = dval;
                else if (var == "zprod") zprod = dval;
                else if (var == "abeam") abeam = dval;
                else if (var == "beamlifetime") beamlifetime = dval;
                else if (var == "resenerg") resenerg = dval;
                else if (var == "ell") ell = dval;
                else if (var == "rstate") rstate = dval;
            }
            else {
                if (var == "beamtyp") beamtyp = val;
                else if (var == "rectyp") rectyp = val;
            }
        }
    }

    inFile.close();

std::cout << "DESDE "<< std::endl;
for(int i =0 ; i<15;i++)
{
for(int j =0 ; j<10;j++)
{
std::cout << "i: " << i << " j: " << j << " br[" << i << "][" << j << "] = " << br[i][j]
          << " md[" << i << "][" << j << "] = " << md[i][j] << std::endl;
}
}
std::cout << "DESDE "<< std::endl;
/*
    std::cout << "Imprimir valores para depuración "<< std::endl;
    std::cout << "beamtyp: " << beamtyp << std::endl;
    std::cout << "rectyp: " << rectyp << std::endl;
    std::cout << "beam_mass_excess: " << beam_mass_excess << std::endl;
    std::cout << "recoil_mass_excess: " << recoil_mass_excess << std::endl;
    std::cout << "part_width: " << part_width << std::endl;
    std::cout << "gam_width: " << gam_width << std::endl;
    std::cout << "spin_stat_fac: " << spin_stat_fac << std::endl;
    std::cout << "zbeam: " << zbeam << std::endl;
    std::cout << "ztarg: " << ztarg << std::endl;
    std::cout << "atarg: " << atarg << std::endl;
    std::cout << "zprod: " << zprod << std::endl;
    std::cout << "abeam: " << abeam << std::endl;
    std::cout << "beamlifetime: " << beamlifetime << std::endl;
    std::cout << "resenerg: " << resenerg << std::endl;
    std::cout << "ell: " << ell << std::endl;
    std::cout << "rstate: " << rstate << std::endl;

for(int i =0 ; i<16;i++)
std::cout << "i " << i << " level "<< level[i] << std::endl;
for(int i =0 ; i<16;i++)
std::cout << "i " << i << " life "<< life[i] << std::endl;
*/

}


void DRAGONPhysicsList::PrintIonProperties(G4ParticleDefinition* ion)
{
   if (!ion) {
        G4cout << "[ERROR] Null ion pointer." << G4endl;
        return;
    }

    G4cout << "===== Ion Properties =====" << G4endl;
    G4cout << "Name: " << ion->GetParticleName() << G4endl;
    G4cout << "Type: " << ion->GetParticleType() << G4endl;
    G4cout << "SubType: " << ion->GetParticleSubType() << G4endl;
    G4cout << "PDG Encoding: " << ion->GetPDGEncoding() << G4endl;
    G4cout << "Mass: " << G4BestUnit(ion->GetPDGMass(), "Energy") << G4endl;
    G4cout << "Charge: " << ion->GetPDGCharge() / eplus << " e" << G4endl;
    G4cout << "Spin: " << ion->GetPDGSpin() << G4endl;
    G4cout << "Magnetic Moment: " << ion->GetPDGMagneticMoment() << G4endl;
    G4cout << "Parity: " << ion->GetPDGiParity() << G4endl;
    G4cout << "Isospin: " << ion->GetPDGIsospin() << G4endl;
    G4cout << "Isospin3: " << ion->GetPDGIsospin3() << G4endl;
    G4cout << "Stable: " << (ion->GetPDGStable() ? "Yes" : "No") << G4endl;
    G4cout << "Lifetime: ";
    if (ion->GetPDGStable())
        G4cout << "Stable" << G4endl;
    else
        G4cout << G4BestUnit(ion->GetPDGLifeTime(), "Time") << G4endl;

    // Solo para núcleos
    if (ion->GetParticleType() == "nucleus") {
        G4cout << "Atomic Number (Z): " << ion->GetAtomicNumber() << G4endl;
        G4cout << "Atomic Mass (A): " << ion->GetAtomicMass() << G4endl;

        if (auto ionPtr = dynamic_cast<G4Ions*>(ion)) {
            G4cout << "Excitation Energy: " << G4BestUnit(ionPtr->GetExcitationEnergy(), "Energy") << G4endl;
        } else {
            G4cout << "Excitation Energy: [N/A - not a G4Ions]" << G4endl;
        }
    }

    // Mostrar tabla de decaimiento si existe
    if (!ion->GetPDGStable()) {
        G4DecayTable* decayTable = ion->GetDecayTable();
        if (decayTable) {
            G4cout << "--- Decay Table ---" << G4endl;
            G4int nDecays = decayTable->entries();
            for (G4int i = 0; i < nDecays; ++i) {
                G4VDecayChannel* channel = decayTable->GetDecayChannel(i);
                G4cout << "Mode " << i + 1 << ": "
                       << channel->GetVerboseLevel() << " "
                       << channel->GetBR() * 100.0 << "% -> ";
                for (G4int j = 0; j < channel->GetNumberOfDaughters(); ++j) {
                    G4cout << channel->GetDaughterName(j);
                    if (j < channel->GetNumberOfDaughters() - 1) G4cout << " + ";
                }
                G4cout << G4endl;
            }
        } else {
            G4cout << "[No decay table available]" << G4endl;
        }
    }

    G4cout << "==========================" << G4endl;
}

void DRAGONPhysicsList::AddIonGasModels()
{
  G4EmConfigurator* em_config =
    G4LossTableManager::Instance()->EmConfigurator();
  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while ((*particleIterator)())
  {
    G4ParticleDefinition* particle = particleIterator->value();
    G4String partname = particle->GetParticleName();
    if(partname == "alpha" || partname == "He3" || partname == "GenericIon") {
      G4BraggIonGasModel* mod1 = new G4BraggIonGasModel();
      G4BetheBlochIonGasModel* mod2 = new G4BetheBlochIonGasModel();
      G4double eth = 2.*MeV*particle->GetPDGMass()/proton_mass_c2;
      em_config->SetExtraEmModel(partname,"ionIoni",mod1,"",0.0,eth,
                                 new G4IonFluctuations());
      em_config->SetExtraEmModel(partname,"ionIoni",mod2,"",eth,100*TeV,
                                 new G4UniversalFluctuation());

    }
  }
}




}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


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
/// \file DRAGON/src/DRAGONHistoManager.cc
/// \brief Implementation of the DRAGONHistoManager class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



#include "DRAGONHistoManager.hh"
#include "DRAGONPrimaryGeneratorAction.hh"
#include "DRAGONDetectorConstruction.hh"

#include "G4AnalysisManager.hh"
#include "global_variables.hh"
#include <iostream>


namespace DRAGON{
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



DRAGONHistoManager::DRAGONHistoManager(DRAGONDetectorConstruction* detector, DRAGONPrimaryGeneratorAction* primary, DRAGONPhysicsList* physics)
:fDet(detector),fPrim(primary),fPhys(physics)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DRAGONHistoManager::~DRAGONHistoManager()    
{

}

void DRAGONHistoManager::ScanHitos()
{
    auto analysisManager = G4AnalysisManager::Instance();

    std::cout << "*************************************************************" << std::endl;
    std::cout << "               Registered histograms                        " << std::endl;
    std::cout << "*************************************************************" << std::endl;

    std::cout << "---------------" << std::endl;
    std::cout << " 1D histograms " << std::endl;
    std::cout << "---------------" << std::endl;

    G4int nH1 = analysisManager->GetNofH1s();
    for (G4int i = 1; i < nH1 + 1; ++i) {
        if (analysisManager->GetH1Name(i).size() > 0) {
            G4int nbins = analysisManager->GetH1Nbins(i);
            G4double xmin = analysisManager->GetH1Xmin(i);
            G4double xmax = analysisManager->GetH1Xmax(i);

            std::cout << "H1 ID: " << i
                      << " Name: " << analysisManager->GetH1Name(i)
                      << " Title: " << analysisManager->GetH1Title(i)
                      << " bins: " << nbins
                      << " xmin: " << xmin
                      << " xmax: " << xmax;

            if (xmin >= xmax) {
                std::cout << "   <<< WARNING: min >= max!";
            }

            std::cout << std::endl;
        }
    }

    std::cout << "---------------" << std::endl;
    std::cout << " 2D histograms " << std::endl;
    std::cout << "---------------" << std::endl;

    G4int nH2 = analysisManager->GetNofH2s();
    for (G4int i = 1; i < nH2 + 1; ++i) {
        if (analysisManager->GetH2Name(i).size() > 0) {
            G4int nx = analysisManager->GetH2Nxbins(i);
            G4double xmin = analysisManager->GetH2Xmin(i);
            G4double xmax = analysisManager->GetH2Xmax(i);

            G4int ny = analysisManager->GetH2Nybins(i);
            G4double ymin = analysisManager->GetH2Ymin(i);
            G4double ymax = analysisManager->GetH2Ymax(i);

            std::cout << "H2 ID: " << i
                      << " Name: " << analysisManager->GetH2Name(i)
                      << " Title: " << analysisManager->GetH2Title(i)
                      << " nbinsX: " << nx << " xmin: " << xmin << " xmax: " << xmax
                      << " nbinsY: " << ny << " ymin: " << ymin << " ymax: " << ymax;

            if (xmin >= xmax || ymin >= ymax) {
                std::cout << "   <<< WARNING: min >= max!";
            }

            std::cout << std::endl;
        }
    }

    std::cout << "*************************************************************" << std::endl;
}

void DRAGONHistoManager::PrintIDMap()
{
    std::cout << "================= IDMap Contents =================" << std::endl;
    for (const auto& pair : GlobalVariables::Instance->IDMap) {
        std::cout << "UserKey: " << pair.first
                  << "  -->  HistID: " << pair.second << std::endl;
    }
    std::cout << "==================================================" << std::endl;
}


}


















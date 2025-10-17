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
// $Id$
//
/// \file DemandAnalysis.hh
/// \brief Selection of the analysis technology

#ifndef DemandAnalysis_h
#define DemandAnalysis_h 1
#include <string>

namespace CLHEP { class HepLorentzVector; }
namespace g4gen { class ReactionKinematics; }

class DemandAnalysis {
private:
	DemandAnalysis();
public:
	static DemandAnalysis* Instance();
	virtual ~DemandAnalysis();
	void OpenFile(const std::string& filename);
	void CloseFile();
	void Write();
	void Clear();
	void AddHit(double edep, double time, double xpos, double ypos, double zpos, int pA, int pZ);
	void SetFirstInteraction(double,double,double,double);
	void Analyze();
	long GetEventsAboveThreshold() const;
	long GetEventsCrossingDetector() const;
	void AddEventCrossingDetector();
	void SetGeneratedNeutron(const CLHEP::HepLorentzVector& p);
	void SetGeneratedRecoil(const CLHEP::HepLorentzVector& p);
	void FillGenTree();
private:
	void CalculateReaction(g4gen::ReactionKinematics*);
};


#endif

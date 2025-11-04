#include "G4AnalysisManager.hh"     //Geant4
#include "G4Run.hh" 
#include "Randomize.hh"

#include "DRAGONRunAction.hh"     //Local


namespace DRAGON {


void DRAGONRunAction::uglast(const G4Run* run)
{
//C.
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//C                                                                      C
//C             UGLAST terminates the GEANT/USER  program                C
//C                                                                      C
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//C.
    G4int nbeam;

    auto TIMEND = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(TIMEND - TIMINT);

    G4cout << "Time elapsed after initialization = " << duration.count() << " seconds" << G4endl;

    ievent = run->GetNumberOfEvent();
    nbeam = ievent - nreact;

    idevt = idevt - ievent; // idevt = sum of all events + sum of good events

    G4cout << "Total number of events generated " << ievent << G4endl;
    G4cout << "Number of events that went to output " << idevt << G4endl;

     CLHEP::HepRandom::saveEngineStatus("engine_after_run.rndm");

     const long* seeds = CLHEP::HepRandom::getTheSeeds();
    G4cout << "Random Number at the beginning of last event" << G4endl;
    G4cout << "**** " << seeds[0] << ", " << seeds[1] << " ****" << G4endl;

    G4cout << "Beam particles entering gas cell " << ngascell << G4endl;

    if (ngascell > 0) {
        G4cout << "Reactions occurring " << nreact << " "
               << 100. * nreact / ngascell << "%  +- "
               << 100. * std::sqrt(nreact * (1. - nreact / static_cast<float>(ngascell))) / ngascell << "%" 
               << G4endl;
    }

    if (nreact > 0) {
        G4cout << "\n** Recoil efficiencies are with respect to the number of Reactions **" << G4endl;

        G4cout << "Recoils exiting target " << ntargexit << " "
               << 100. * ntargexit / nreact << "%  +- "
               << 100. * std::sqrt(ntargexit * (1. - ntargexit / static_cast<float>(nreact))) / nreact << "%"
               << G4endl;
//C
//C MT adds counters for the number of recoils that make it
//C   to Q3 (Sext 1) and to Q8 (Quad 6).
//C
        G4cout << "Recoils reaching Q3 (Sext 1) " << Num_Recoils_Q3 << " "
               << 100. * Num_Recoils_Q3 / nreact << "%  +- "
               << 100. * std::sqrt(Num_Recoils_Q3 * (1. - Num_Recoils_Q3 / static_cast<float>(nreact))) / nreact << "%"
               << G4endl;

        G4cout << "Recoils reaching Q8 (Quad 6) " << Num_Recoils_Q8 << " "
               << 100. * Num_Recoils_Q8 / nreact << "%  +- "
               << 100. * std::sqrt(Num_Recoils_Q8 * (1. - Num_Recoils_Q8 / static_cast<float>(nreact))) / nreact << "%"
               << G4endl;

        G4cout << "Recoils reaching end detector " << nend << " "
               << 100. * nend / nreact << "%  +- "
               << 100. * std::sqrt(nend * (1. - nend / static_cast<float>(nreact))) / nreact << "%"
               << G4endl;
//C
//C MT adds counter for the number of beam particles that reach
//C   the end detector.
//C
        G4cout << "Beam particles exiting target " << nbeamout << " "
               << 100. * nbeamout / nbeam << "%  +- "
               << 100. * std::sqrt(nbeamout * (1. - nbeamout / static_cast<float>(nbeam))) / nbeam << "%"
               << G4endl;

        G4cout << "Beam particles reaching FCM2 " << nfcm2 << " "
               << 100. * nfcm2 / nbeam << "%  +- "
               << 100. * std::sqrt(nfcm2 * (1. - nfcm2 / static_cast<float>(nbeam))) / nbeam << "%"
               << G4endl;

        G4cout << "Beam particles reaching end detector " << Num_BeamPart_ENDV << " "
               << 100. * Num_BeamPart_ENDV / nbeam << "%  +- "
               << 100. * std::sqrt(Num_BeamPart_ENDV * (1. - Num_BeamPart_ENDV / static_cast<float>(nbeam))) / nbeam << "%"
               << G4endl;
    }

//C.
//C.-->   Save histograms
//C.
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->IsActive() ) {
    analysisManager->Write();
    analysisManager->CloseFile();
  } 
}

}

#include <vector>

#include "globals.hh"
#include "G4Neutron.hh"
#include "G4RunManager.hh"
#include "Randomize.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4EmCalculator.hh"
#include "G4LorentzVector.hh"

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "DemandAnalysis.hh"
#include "DemandPrimaryGeneratorAction.hh"
#include "DemandDetectorConstruction.hh"
#include "ReactionKinematics.hh"
#include "DemandSD.hh"

using std::string;
using std::vector;
using namespace CLHEP;

namespace {

const double FWHM = 1/(2*sqrt(2*log(2)));
const double TIME_RES = 0.6 * ns * FWHM;
// double ResEnergy(double energy) {
// 	// Use resolution from NIMA 792, p. 74 (2015)
// 	// (Eq. 3, NOTE it's given in % in the paper)
// 	const double a=0.93, b=7.68, c=0.2;
// 	const double res_fwhm = energy*sqrt(pow(a,2) + pow(b,2)/energy + pow(c/energy,2))/100;
// 	return res_fwhm * FWHM;
// }

TFile* fFile = 0;
TTree* fTree = 0;
TTree* fGenTree = 0;

vector<double> *fEdep = 0;
vector<double> *fEdep_noquench = 0;
vector<double> *fTime = 0;
vector<double> *fXpos = 0;
vector<double> *fYpos = 0;
vector<double> *fZpos = 0;
vector<int> *fDetno = 0;
vector<int> *fpA = 0;
vector<int> *fpZ = 0;
vector<int> *fpID = 0;

int fNumHits = 0;
double fEkin = 0;
double fTheta = 0;
double fPhi = 0;
double fEx = 0;

double fRx = 0;
double fRy = 0;
double fRz = 0;

long fEventsAboveThreshold = 0;
long fEventsCrossingDetector = 0;
bool fDetected = false;

TLorentzVector *fFirstInteraction = 0;

TVector3       *fReacPos = 0;
TLorentzVector *fNeutronMomentum = 0;
TLorentzVector *fRecoilMomentum = 0;
Int_t fCrossedDetector = 0;

std::vector<double> *fPrimaryScatterX = 0;
std::vector<double> *fPrimaryScatterY = 0;
std::vector<double> *fPrimaryScatterZ = 0;
std::vector<double> *fPrimaryScatterT = 0;
std::vector<double> *fPrimaryScatterE = 0;
std::vector<std::string> *fPrimaryScatterVolumeName = 0;

double fTimeHit0;


// step tree
bool fSaveStepTree = false;
TTree* fStepTree = 0;

G4int fIsRecdet = 0;
G4int fNumSteps = 0;
std::vector<G4int>    *fStepDetno = 0;
std::vector<G4double> *fStepEdep = 0;
std::vector<G4double> *fStepEdepQuenched = 0;
std::vector<G4double> *fStepLen = 0;
std::vector<G4int>    *fStepParticleID = 0;
std::vector<std::string>  *fStepParticleName = 0;
}

DemandAnalysis::DemandAnalysis() { }

DemandAnalysis::~DemandAnalysis(){ }

DemandAnalysis* DemandAnalysis::Instance()
{
	static DemandAnalysis* instance = 0;
	if(!instance) {
		instance = new DemandAnalysis();
	}
	return instance;
}

void DemandAnalysis::OpenFile(const string& filename)
{
	// static int ntimes = 0;
	// if(++ntimes > 1) {
	// 	throw std::invalid_argument(
	// 		"DemandAnalysis::OpenFile -- Called more than once");
	// }
	fFile = new TFile(filename.c_str(), "recreate");
	fTree = new TTree("DemandTree", "Para-Terphenyl detector tree");

	fTree->Branch("edep",&fEdep);
	fTree->Branch("edep_noquench",&fEdep_noquench);
	fTree->Branch("time",&fTime);
	fTree->Branch("xpos",&fXpos);
	fTree->Branch("ypos",&fYpos);
	fTree->Branch("zpos",&fZpos);
	fTree->Branch("detno",&fDetno);

	fTree->Branch("particleA",&fpA);
	fTree->Branch("particleZ",&fpZ);
	fTree->Branch("particleID",&fpID);

	fTree->Branch("numHits",&fNumHits,"numHits/I");
	fTree->Branch("ekin",&fEkin,"ekin/D");
	fTree->Branch("theta",&fTheta,"theta/D");
	fTree->Branch("phi",&fPhi,"phi/D");
	fTree->Branch("ex",&fEx,"ex/D");

	fTree->Branch("firstInteraction", "TLorentzVector", &fFirstInteraction);

	fTree->Branch("reacpos","TVector3",&fReacPos);
	fTree->Branch("pneut","TLorentzVector",&fNeutronMomentum);
	fTree->Branch("precoil","TLorentzVector",&fRecoilMomentum);

	fTree->Branch("xreact",&fRx,"xreact/D");
	fTree->Branch("yreact",&fRy,"yreact/D");
	fTree->Branch("zreact",&fRz,"zreact/D");

	fTree->Branch("scatterX",&fPrimaryScatterX);
	fTree->Branch("scatterY",&fPrimaryScatterY);
	fTree->Branch("scatterZ",&fPrimaryScatterZ);
	fTree->Branch("scatterT",&fPrimaryScatterT);
	fTree->Branch("scatterE",&fPrimaryScatterE);
	fTree->Branch("scatterVolumeName",&fPrimaryScatterVolumeName);

	fTree->Branch("timeHit0", &fTimeHit0);
	

	fEventsAboveThreshold = 0;
	fEventsCrossingDetector = 0;

	fGenTree = new TTree("GenTree", "Tree of generated neutrons");

	fGenTree->Branch("reacpos","TVector3",&fReacPos);
	fGenTree->Branch("pneut","TLorentzVector",&fNeutronMomentum);
	fGenTree->Branch("precoil","TLorentzVector",&fRecoilMomentum);
	fGenTree->Branch("crossed_detector",&fCrossedDetector);
	fGenTree->Branch("detected",&fDetected);
	fGenTree->Branch("recdet",&fIsRecdet);

	// fGenTree->Branch("scatterX",&fPrimaryScatterX);
	// fGenTree->Branch("scatterY",&fPrimaryScatterY);
	// fGenTree->Branch("scatterZ",&fPrimaryScatterZ);
	// fGenTree->Branch("scatterT",&fPrimaryScatterT);
	// fGenTree->Branch("scatterE",&fPrimaryScatterE);
	// fGenTree->Branch("scatterVolumeName",&fPrimaryScatterVolumeName);

	// fGenTree->Branch("timeHit0",&fTimeHit0);


	if(fSaveStepTree) {
		fStepTree = new TTree("StepTree", "Tree of all steps");

		fStepTree->Branch("fIsRecdet", &fIsRecdet);
		fStepTree->Branch("fNumSteps", &fNumSteps);
		fStepTree->Branch("fStepDetno", &fStepDetno);
		fStepTree->Branch("fStepEdep", &fStepEdep);
		fStepTree->Branch("fStepEdepQuenched", &fStepEdepQuenched);
		fStepTree->Branch("fStepLen", &fStepLen);
		fStepTree->Branch("fStepParticleName", &fStepParticleName);
		fStepTree->Branch("fStepParticleID", &fStepParticleID);
	}
}

void DemandAnalysis::CloseFile()
{
	fFile->Close();
	fEventsAboveThreshold = 0;
	fEventsCrossingDetector = 0;
}

void DemandAnalysis::Write()
{
	fFile->cd();
	fTree->Write();
	fGenTree->Write();
	if(fSaveStepTree && fStepTree){
		fStepTree->Write();
	}
}

void DemandAnalysis::ClearPrimaryScatters()
{
	fPrimaryScatterX->clear();
	fPrimaryScatterY->clear();
	fPrimaryScatterZ->clear();
	fPrimaryScatterT->clear();
	fPrimaryScatterE->clear();
	fPrimaryScatterVolumeName->clear();
}

void DemandAnalysis::Clear()
{
	fTimeHit0 = -1;
	
	fEdep->clear();
	fEdep_noquench->clear();
	fTime->clear();
	fXpos->clear();
	fYpos->clear();
	fZpos->clear();
	fDetno->clear();

	fpA->clear();
	fpZ->clear();
	fpID->clear();

	fNumHits = 0;
	fEkin = 0;
	fTheta = 0;
	fPhi = 0;
	fEx = 0;

	fFirstInteraction->SetXYZT(0,0,0,0);
//	fIsRecdet = 0;  //handled.mamnually.in EndOfEventAction
	
	// step tree
	if(fSaveStepTree){
		fNumSteps = 0;
		fStepDetno->clear();
		fStepEdep->clear();
		fStepEdepQuenched->clear();
		fStepLen->clear();
		fStepParticleName->clear();
		fStepParticleID->clear();
	};
}

void DemandAnalysis::SetFirstInteraction(double time, double x, double y, double z)
{
	fFirstInteraction->SetXYZT(
		x,y,z,time);
}

void DemandAnalysis::AddPrimaryScatter(
	const CLHEP::Hep3Vector& pos, double time, double energy, const G4String& vname)
{
	fPrimaryScatterX->push_back(pos.x());
	fPrimaryScatterY->push_back(pos.y());
	fPrimaryScatterZ->push_back(pos.z());
	fPrimaryScatterT->push_back(time);
	fPrimaryScatterE->push_back(energy);
	fPrimaryScatterVolumeName->push_back(vname);
}

void DemandAnalysis::SortPrimaryScatters()
{
	// time sort
	vector<Int_t> indx(fPrimaryScatterT->size());
	TMath::Sort(int(indx.size()), fPrimaryScatterT->data(), &indx[0], false);
	auto X = *fPrimaryScatterX;
	auto Y = *fPrimaryScatterY;
	auto Z = *fPrimaryScatterZ;
	auto T = *fPrimaryScatterT;
	auto E = *fPrimaryScatterE;
	auto N = *fPrimaryScatterVolumeName;
	for(size_t i=0; i< indx.size(); ++i){
		fPrimaryScatterX->at(i) = X.at(indx.at(i));
		fPrimaryScatterY->at(i) = Y.at(indx.at(i));
		fPrimaryScatterZ->at(i) = Z.at(indx.at(i));
		fPrimaryScatterE->at(i) = E.at(indx.at(i));
		fPrimaryScatterT->at(i) = T.at(indx.at(i));
		fPrimaryScatterVolumeName->at(i) = N.at(indx.at(i));
	}
}

void DemandAnalysis::AddHit(
	double edep, double edep_noquench, double time,
	double xpos, double ypos, double zpos,
	int pA, int pZ, int pID, int detno)
{
	// record TRUE time of earliest hit
	if(fTimeHit0 < 0 || time < fTimeHit0){
		fTimeHit0 = time;
	}
	
	// add resolutions
	time += G4RandGauss::shoot(0, TIME_RES/FWHM);
	// energy resolution already done in DemandEventAction
//	edep += G4RandGauss::shoot(0, ResEnergy(edep));
// 	edep = DemandSD::CalculateEnergyResolution(edep);

	if(edep > 0){
		// time sort
		auto it = lower_bound(fTime->begin(), fTime->end(), time);
		auto dit = it - fTime->begin();
		fTime->emplace(it, time);
		fEdep->emplace(fEdep->begin() + dit, edep);
		fEdep_noquench->emplace(fEdep_noquench->begin() + dit, edep_noquench);
		fXpos->emplace(fXpos->begin() + dit, xpos);
		fYpos->emplace(fYpos->begin() + dit, ypos);
		fZpos->emplace(fZpos->begin() + dit, zpos);
		fDetno->emplace(fDetno->begin() + dit, detno);
		fpA->emplace(fpA->begin() + dit, pA);
		fpZ->emplace(fpZ->begin() + dit, pZ);
		fpID->emplace(fpID->begin() + dit, pID);
		fNumHits++;
	}
}

void DemandAnalysis::Analyze()
{
	if(fNumHits == 0) { return; }

	++fEventsAboveThreshold;
	fDetected = true;

	// first hit
	const double M0 = G4Neutron::Definition()->GetPDGMass();

	TVector3 hitPos(fXpos->at(0), fYpos->at(0), fZpos->at(0));
	const double hitTime = fTime->at(0);
	const double hitVel = hitPos.Mag() / hitTime;
	const double hitBeta = hitVel / (TMath::C()/1e6);
	fEkin = (1/sqrt(1 - hitBeta*hitBeta) - 1) * M0;
	fTheta = hitPos.Theta()/deg;
	fPhi = hitPos.Phi()/deg;

	g4gen::ReactionKinematics* reaction =
		dynamic_cast<const DemandPrimaryGeneratorAction&>(
			*(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction())
			).GetReactionKinematics();
	if(reaction) {
		CalculateReaction(reaction);
	}

	const DemandPrimaryGeneratorAction* genAction =
    static_cast<const DemandPrimaryGeneratorAction*>(
        G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());

	G4ThreeVector pos = genAction->GetReactionPosition();
	fRx = pos.x();
	fRy = pos.y();
	fRz = pos.z();

	SortPrimaryScatters();
	
	fFile->cd();
	fTree->Fill();
}

long DemandAnalysis::GetEventsAboveThreshold() const
{
	return fEventsAboveThreshold;
}

long DemandAnalysis::GetEventsCrossingDetector() const
{
	return fEventsCrossingDetector;
}

void DemandAnalysis::AddEventCrossingDetector()
{
	if(fCrossedDetector == 0) ++fEventsCrossingDetector;
	++fCrossedDetector;
}

void DemandAnalysis::CalculateReaction(g4gen::ReactionKinematics* reaction)
{
	G4LorentzVector beam;
	G4LorentzVector target = reaction->GetTarget();

	auto beamDefinition = dynamic_cast<const DemandPrimaryGeneratorAction&>(
		*(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction())
		).GetBeamDefinition();
	const G4double incident_beam_energy =
		dynamic_cast<const DemandPrimaryGeneratorAction&>(
			*(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction())
			).GetBeamEnergy();

	G4double ebeam = incident_beam_energy;

	// energy loss through half target
	if(dynamic_cast<const DemandDetectorConstruction&>(
			 *(G4RunManager::GetRunManager()->GetUserDetectorConstruction())).
		 GetHaveTarget())
	{
		G4LogicalVolume* trgtLV =
			G4LogicalVolumeStore::GetInstance()->GetVolume("Target");
		if ( trgtLV ) {
			double thickness = dynamic_cast<G4Box&>(*(trgtLV->GetSolid())).
				GetZHalfLength()*2.;

			static G4EmCalculator* emCalc = 0;
			if(!emCalc) { emCalc = new G4EmCalculator(); }
			G4double dedx = emCalc->ComputeTotalDEDX(
				incident_beam_energy, beamDefinition, trgtLV->GetMaterial());
			G4double dx = thickness/2;
			ebeam = incident_beam_energy - dedx*dx;
//			G4cout << dedx*dx*2 << G4endl;
		}
		else { throw trgtLV; }
	}
	const G4double mbeam = beamDefinition->GetPDGMass();
	const G4double pbeam = sqrt(pow(ebeam+mbeam,2) - pow(mbeam,2));
	const G4double thbeam = 0;
	const G4double phbeam = 0;
	beam = G4LorentzVector(
		pbeam*sin(thbeam)*cos(phbeam),
		pbeam*sin(thbeam)*sin(phbeam),
		pbeam*cos(thbeam),
		mbeam + ebeam);



	double m3 = reaction->GetProduct(0).m();
	double p3 = sqrt(pow(m3+fEkin,2) - m3*m3);
	G4LorentzVector ejectile(
		p3*sin(fTheta*deg)*cos(fPhi*deg),
		p3*sin(fTheta*deg)*cos(fPhi*deg),
		p3*cos(fTheta*deg),
		m3+fEkin);

	auto recoilDefinition = dynamic_cast<const DemandPrimaryGeneratorAction&>(
		*(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction())
		).GetRecoilDefinition();
	if(!recoilDefinition) {
		throw std::logic_error(
			"DemandAnalysis::CalculateReaction :: NULL Recoil Definition");
	}

	double m4 = recoilDefinition->GetPDGMass();
	G4LorentzVector recoil = beam + target - ejectile;
	fEx = recoil.m() - m4;
}

void DemandAnalysis::SetReacPos(const G4ThreeVector& p)
{
	fReacPos->SetXYZ(p.x(),p.y(),p.z());
}

void DemandAnalysis::SetGeneratedNeutron(const G4LorentzVector& p)
{
	fNeutronMomentum->SetPxPyPzE(
		p.px(),p.py(),p.pz(),p.e());
}

void DemandAnalysis::SetGeneratedRecoil(const G4LorentzVector& p)
{
	fRecoilMomentum->SetPxPyPzE(
		p.px(),p.py(),p.pz(),p.e());
}

void DemandAnalysis::FillGenTree()
{
	fGenTree->Fill();
	fCrossedDetector = 0;
	fDetected = false;
	fReacPos->SetXYZ(0,0,0);\
	fNeutronMomentum->SetPxPyPzE(0,0,0,0);
	fRecoilMomentum->SetPxPyPzE(0,0,0,0);
}


std::vector<long> DemandAnalysis::GetEventsAboveSoftwareCut(
	std::vector<G4double> cuts) const
{
	std::vector<long> output;
	for(const auto& c : cuts){
		long N = fTree->GetPlayer()->GetEntries(
			Form("edep[0] > %.6E && particleID[0] == 2212", c) );
		output.push_back(N);
	}
	return output;
}


// Step Tree //
void DemandAnalysis::SetSaveStepTree(bool save)
{ fSaveStepTree = save; }

bool DemandAnalysis::GetSaveStepTree() const
{ return fSaveStepTree; }

void DemandAnalysis::FillStepTree()
{
	if(fSaveStepTree && fStepTree && fNumSteps > 0) {
		fStepTree->Fill();
	}
}

void DemandAnalysis::SetRecdet(G4int recdet)
{ fIsRecdet = recdet; }

void DemandAnalysis::AddStep(
	G4int detno, G4double edep, G4double edep_quenched, G4double stepLen,
	G4String particleName, G4int particleID)
{
	if(fSaveStepTree){
		++fNumSteps;
		fStepDetno->push_back(detno);
		fStepEdep->push_back(edep);
		fStepEdepQuenched->push_back(edep_quenched);
		fStepLen->push_back(stepLen);
		fStepParticleName->push_back(particleName);
		fStepParticleID->push_back(particleID);
	}
}

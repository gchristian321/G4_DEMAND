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

#include "PterpAnalysis.hh"
#include "PterpPrimaryGeneratorAction.hh"
#include "PterpDetectorConstruction.hh"
#include "ReactionKinematics.hh"

using std::string;
using std::vector;
using namespace CLHEP;

namespace {

const double FWHM = 1/(2*sqrt(2*log(2)));
const double TIME_RES = 0.6 * ns * FWHM;
double ResEnergy(double energy) {
	// Use resolution from NIMA 792, p. 74 (2015)
	// (Eq. 3, NOTE it's given in % in the paper)
	const double a=0.93, b=7.68, c=0.2;
	const double res_fwhm = energy*sqrt(pow(a,2) + pow(b,2)/energy + pow(c/energy,2))/100;
	return res_fwhm * FWHM;
}

TFile* fFile = 0;
TTree* fTree = 0;
TTree* fGenTree = 0;

vector<double> *fEdep = 0;
vector<double> *fTime = 0;
vector<double> *fXpos = 0;
vector<double> *fYpos = 0;
vector<double> *fZpos = 0;
vector<int> *fpA = 0;
vector<int> *fpZ = 0;

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

TLorentzVector *fNeutronMomentum = 0;
TLorentzVector *fRecoilMomentum = 0;
Int_t fCrossedDetector = 0;
}

PterpAnalysis::PterpAnalysis() { }

PterpAnalysis::~PterpAnalysis(){ }

PterpAnalysis* PterpAnalysis::Instance()
{
	static PterpAnalysis* instance = 0;
	if(!instance) {
		instance = new PterpAnalysis();
	}
	return instance;
}

void PterpAnalysis::OpenFile(const string& filename)
{
	// static int ntimes = 0;
	// if(++ntimes > 1) {
	// 	throw std::invalid_argument(
	// 		"PterpAnalysis::OpenFile -- Called more than once");
	// }
	fFile = new TFile(filename.c_str(), "recreate");
	fTree = new TTree("PterpTree", "Para-Terphenyl detector tree");

	fTree->Branch("edep",&fEdep);
	fTree->Branch("time",&fTime);
	fTree->Branch("xpos",&fXpos);
	fTree->Branch("ypos",&fYpos);
	fTree->Branch("zpos",&fZpos);

	fTree->Branch("particleA",&fpA);
	fTree->Branch("particleZ",&fpZ);

	fTree->Branch("numHits",&fNumHits,"numHits/I");
	fTree->Branch("ekin",&fEkin,"ekin/D");
	fTree->Branch("theta",&fTheta,"theta/D");
	fTree->Branch("phi",&fPhi,"phi/D");
	fTree->Branch("ex",&fEx,"ex/D");
	
	fTree->Branch("firstInteraction", "TLorentzVector", &fFirstInteraction);

	fTree->Branch("pneut","TLorentzVector",&fNeutronMomentum);
	fTree->Branch("precoil","TLorentzVector",&fRecoilMomentum);

	fTree->Branch("xreact",&fRx,"xreact/D");
	fTree->Branch("yreact",&fRy,"yreact/D");
	fTree->Branch("zreact",&fRz,"zreact/D");


	fEventsAboveThreshold = 0;
	fEventsCrossingDetector = 0;

	fGenTree = new TTree("GenTree", "Tree of generated neutrons");

	fGenTree->Branch("pneut","TLorentzVector",&fNeutronMomentum);
	fGenTree->Branch("precoil","TLorentzVector",&fRecoilMomentum);
	fGenTree->Branch("crossed_detector",&fCrossedDetector);
	fGenTree->Branch("detected",&fDetected);
}

void PterpAnalysis::CloseFile()
{
	fFile->Close();
	fEventsAboveThreshold = 0;
	fEventsCrossingDetector = 0;
}

void PterpAnalysis::Write()
{
	fFile->cd();
	fTree->Write();
	fGenTree->Write();
}

void PterpAnalysis::Clear()
{
	fEdep->clear();
	fTime->clear();
	fXpos->clear();
	fYpos->clear();
	fZpos->clear();

	fpA->clear();
	fpZ->clear();

	fNumHits = 0;
	fEkin = 0;
	fTheta = 0;
	fPhi = 0;
	fEx = 0;

	fFirstInteraction->SetXYZT(0,0,0,0);
}

void PterpAnalysis::SetFirstInteraction(double time, double x, double y, double z)
{
	fFirstInteraction->SetXYZT(
		x,y,z,time);
}

void PterpAnalysis::AddHit(	
	double edep, double time, double xpos, double ypos, double zpos, int pA, int pZ)
{
	// add resolutions
	time += G4RandGauss::shoot(0, TIME_RES/FWHM);
	edep += G4RandGauss::shoot(0, ResEnergy(edep));
	
	// time sort
	auto it = lower_bound(fTime->begin(), fTime->end(), time);
	auto dit = it - fTime->begin();
	fTime->emplace(it, time);
	fEdep->emplace(fEdep->begin() + dit, edep);
	fXpos->emplace(fXpos->begin() + dit, xpos);
	fYpos->emplace(fYpos->begin() + dit, ypos);
	fZpos->emplace(fZpos->begin() + dit, zpos);
	fpA->emplace(fpA->begin() + dit, pA);
	fpZ->emplace(fpZ->begin() + dit, pZ);
	fNumHits++;
}

void PterpAnalysis::Analyze()
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
		dynamic_cast<const PterpPrimaryGeneratorAction&>(
			*(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction())
			).GetReactionKinematics();
	if(reaction) {
		CalculateReaction(reaction);
	}

	const PterpPrimaryGeneratorAction* genAction =
    static_cast<const PterpPrimaryGeneratorAction*>(
        G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());

	G4ThreeVector pos = genAction->GetReactionPosition();
	fRx = pos.x();
	fRy = pos.y();
	fRz = pos.z();
	
	fFile->cd();
	fTree->Fill();
}

long PterpAnalysis::GetEventsAboveThreshold() const
{
	return fEventsAboveThreshold;
}

long PterpAnalysis::GetEventsCrossingDetector() const
{
	return fEventsCrossingDetector;
}

void PterpAnalysis::AddEventCrossingDetector()
{
	if(fCrossedDetector == 0) ++fEventsCrossingDetector;
	++fCrossedDetector;
}

void PterpAnalysis::CalculateReaction(g4gen::ReactionKinematics* reaction)
{
	G4LorentzVector beam;
	G4LorentzVector target = reaction->GetTarget();

	auto beamDefinition = dynamic_cast<const PterpPrimaryGeneratorAction&>(
		*(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction())
		).GetBeamDefinition();
	const G4double incident_beam_energy =
		dynamic_cast<const PterpPrimaryGeneratorAction&>(
			*(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction())
			).GetBeamEnergy();

	G4double ebeam = incident_beam_energy;
	
	// energy loss through half target
	if(dynamic_cast<const PterpDetectorConstruction&>(
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

	auto recoilDefinition = dynamic_cast<const PterpPrimaryGeneratorAction&>(
		*(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction())
		).GetRecoilDefinition();
	if(!recoilDefinition) {
		throw std::logic_error(
			"PterpAnalysis::CalculateReaction :: NULL Recoil Definition");
	}
	
	double m4 = recoilDefinition->GetPDGMass();
	G4LorentzVector recoil = beam + target - ejectile;
	fEx = recoil.m() - m4;
}

void PterpAnalysis::SetGeneratedNeutron(const G4LorentzVector& p)
{
	fNeutronMomentum->SetPxPyPzE(
		p.px(),p.py(),p.pz(),p.e());
}

void PterpAnalysis::SetGeneratedRecoil(const G4LorentzVector& p)
{
	fRecoilMomentum->SetPxPyPzE(
		p.px(),p.py(),p.pz(),p.e());
}

void PterpAnalysis::FillGenTree()
{
	fGenTree->Fill();
	fCrossedDetector = 0;
	fDetected = false;
	fNeutronMomentum->SetPxPyPzE(0,0,0,0);
	fRecoilMomentum->SetPxPyPzE(0,0,0,0);
}

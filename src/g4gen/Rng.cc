#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>

#include <cassert>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdexcept>

#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "g4gen/Error.hh"
#include "g4gen/Rng.hh"


namespace { 

void handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
} }

void g4gen::CheckMaxTries::operator() (G4int& n, const char* fct, G4double* low, G4double* high)
{
	if(n++ == kMaxTries) {
		G4GWAR << "g4gen::Rng::" << fct << ":: Failed to generate RNG fulfilling requirements";
		if(low)  { 
			G4cerr << "(above " << *low; 
		}
		if(high) { 
			if(low) { G4cerr << " and"; } else { G4cerr << "("; }
			G4cerr << " below " << *high << ")"; 
		}
		G4cerr << " after " << kMaxTries << " attempts. Enter 1 to abort program "
					 << "(and print stack trace), or 0 to continue..." << G4endl;
		G4int val;
		G4cin >> val;
		if(val != 0) { handler(0); exit(1); }
	}
}



void g4gen::SetRngSeed(G4int seed)
{
	CLHEP::HepRandom::setTheSeed( seed );
}

G4int g4gen::GetRngSeed()
{
	return CLHEP::HepRandom::getTheSeed( );
}


// RNG (BASE CLASS) //

G4double g4gen::Rng::GenerateAbove(G4double low)
{
	G4int n = 0;
	G4double value;
	do { 
		value = Generate();
		g4gen::CheckMaxTries() (n, "GenerateAbove", &low, 0);
	} while (value < low);
	return value;
}

G4double g4gen::Rng::GenerateBelow(G4double high)
{
	G4int n = 0;
	G4double value;
	do {
		value = Generate();
		g4gen::CheckMaxTries() (n, "GenerateBelow", 0, &high);
	} while (value >= high);
	return value;
}

G4double g4gen::Rng::GenerateBetween(G4double low, G4double high)
{
	G4int n = 0;
	G4double value;
	do { 
		value = Generate();
		g4gen::CheckMaxTries() (n, "GenerateBetween", &low, &high);
	} while (value < low || value >= high);
	return value;
}


// RNG UNIFORM //

g4gen::RngUniform::RngUniform(G4double low, G4double high):
	mLow(low), mHigh(high)
{ }

g4gen::RngUniform::~RngUniform()
{ }

G4double g4gen::RngUniform::DoGenerate()
{
	G4double rndm = G4UniformRand();
	return mLow + (mHigh-mLow)*rndm;
}


// RNG GAUS //

g4gen::RngGaus::RngGaus(G4double mean, G4double sigma):
	mMean(mean), mSigma(sigma)
{ }

g4gen::RngGaus::~RngGaus()
{ }

G4double g4gen::RngGaus::DoGenerate()
{
	return mSigma > 1e-9 ? CLHEP::RandGauss::shoot(mMean, mSigma) : mMean;
}


// RNG BREIT WIGNER //

g4gen::RngBreitWigner::RngBreitWigner(G4double mean, G4double sigma):
	mMean(mean), mSigma(sigma)
{ }

g4gen::RngBreitWigner::~RngBreitWigner()
{ }

G4double g4gen::RngBreitWigner::DoGenerate()
{
	return mSigma > 1e-9 ? CLHEP::RandGauss::shoot(mMean, mSigma) : mMean;
}


namespace { void do_file_init(const G4String& filename,
															std::vector<G4double>& mXlow,
															std::vector<G4double>& mCdf,
															std::vector<G4double>* mPdf)
{
	if(filename == "NULL") return;
	
	std::ifstream ifs(filename.c_str());
  if(!ifs.good()) {
		G4cerr<< "ERROR:: Invalid filename passed to g4gen::RngCustom ctor:: "
					<< filename << G4endl;
		throw std::invalid_argument(filename);
	}

  std::string line;
  std::vector<double> pdf;
  double integral = 0;

  double x, w;
  while(std::getline(ifs, line)) {
    std::stringstream iss(line);
    if (!(iss >> x >> w)) { continue; } // skip comments
		mXlow.push_back(x);
		pdf.push_back(w);
    integral += w;
  }

  double cdfi = 0;
  mCdf.resize(pdf.size());
  for(size_t i=0; i< pdf.size(); ++i) {
    pdf[i] /= integral;
    cdfi += pdf[i];
    mCdf[i] = cdfi;
  }	

	if(mPdf) {
		mPdf->resize(pdf.size());
		std::copy(pdf.begin(), pdf.end(), mPdf->begin());
	}
} }

// RNG CUSTOM //

g4gen::RngCustom::RngCustom(const G4String& filename)
{
	do_file_init(filename, mXlow, mCdf, 0);
}

g4gen::RngCustom::~RngCustom()
{ }

G4double g4gen::RngCustom::DoGenerate()
{
	G4double r = g4gen::RngUniform(0,1).Generate();
	std::vector<double>::const_iterator it =
		std::lower_bound(mCdf.begin(), mCdf.end(), r);
	size_t indx = it - mCdf.begin();

	if(indx+1 < mXlow.size()) {
		return g4gen::RngUniform(mXlow[indx], mXlow[indx+1]).Generate();
	} else {
		G4cerr << "ERROR:: g4gen::RngCustom::DoGenerate:: "
					 << "Found out-of-range value in CDF. (indx, r): "
					 << indx << ", " << r << G4endl;
		throw std::range_error(std::to_string(indx));
	}
}



// RNG CUSTOM ANG DIST //

g4gen::RngCustomAngDist::RngCustomAngDist(const G4String& filename):
	g4gen::RngCustom("NULL")
{
	std::vector<G4double> pdf;
	do_file_init(filename, mXlow, mCdf, &pdf);

	G4double integral = 0;
	for(size_t i=0; i< mXlow.size() - 1; ++i) {
		G4double diff = mXlow.at(i+1) - mXlow.at(i);
		G4double Angle = mXlow.at(i) + diff/2;
		pdf.at(i) = pdf[i]*sin(Angle*CLHEP::pi/180);
		integral += pdf.at(i);
	}

  double cdfi = 0;
  mCdf.resize(pdf.size());
  for(size_t i=0; i< pdf.size(); ++i) {
    pdf[i] /= integral;
    cdfi += pdf[i];
    mCdf[i] = cdfi;
  }
}

g4gen::RngCustomAngDist::~RngCustomAngDist()
{ }






// RNG GAUS 2D //

g4gen::RngGaus2d::RngGaus2d(G4double sigma_x, G4double sigma_y, G4double rho):
	mSigmaX(sigma_x),
	mSigmaY(sigma_y),
	mRho(rho)
{  }

std::pair<G4double, G4double> g4gen::RngGaus2d::DoGenerate()
{
	// Routine taken from GSL source code (GPL licensed)
	// Available online at:
	// https://github.com/LuaDist/gsl/blob/master/randist/bigauss.c
	//
	double u, v, r2, scale;
	g4gen::RngUniform rng;
	
  do
	{
		/* choose x,y in uniform square (-1,-1) to (+1,+1) */

		u = -1 + 2 * rng.Generate();
		v = -1 + 2 * rng.Generate();

		/* see if it is in the unit circle */
		r2 = u * u + v * v;
	}
  while (r2 > 1.0 || r2 == 0);

  scale = sqrt (-2.0 * log (r2) / r2);

  double x = mSigmaX * u * scale;
  double y = mSigmaY * (mRho * u + sqrt(1 - mRho*mRho) * v) * scale;

	return std::make_pair(x, y);
}



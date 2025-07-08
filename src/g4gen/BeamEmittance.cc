#include <cmath>
#include "g4gen/BeamEmittance.hh"

g4gen::BeamEmittance::BeamEmittance():
	fEpsilon(0), fBeta(0), fAlpha(0), fGamma(0), fX0(0)
{}

g4gen::BeamEmittance::BeamEmittance(double epsilon, double alpha, double sigmaX)
{
	SetTwist(epsilon,alpha,sigmaX);
}

g4gen::BeamEmittance::BeamEmittance(double epsilon, double alpha, double sigmaX, double x0)
{
	SetTwist(epsilon,alpha,sigmaX);
	SetX0(x0);
}

void g4gen::BeamEmittance::SetTwist(double epsilon, double alpha, double sigmaX)
{
	fEpsilon = epsilon;
	fAlpha   = alpha;
	fBeta    = pow(sigmaX, 2) / fEpsilon;
	fGamma   = (1 + pow(fAlpha, 2)) / fBeta;
}

double g4gen::BeamEmittance::GetSigmaX()const 
{
	return sqrt(fBeta*fEpsilon);  
}

double g4gen::BeamEmittance::GetSigmaTX() const
{
	return sqrt(fGamma*fEpsilon);
}

double g4gen::BeamEmittance::GetRho() const
{
	return -1*fAlpha*fEpsilon / (GetSigmaX()*GetSigmaTX()); 
}


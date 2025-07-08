#include <G4SystemOfUnits.hh>
#include "g4gen/NuclearMasses.hh"
#include "g4gen/Particle.hh"


g4gen::Particle::Particle():
	fA(0), fZ(0), fM(0), fEx(0), fP(0,0,0,0), fPos(0,0,0)
{ }

g4gen::Particle::~Particle()
{ }

void g4gen::Particle::SetNucleus(G4int Z, G4int A)
{
	fZ=Z;
	fA=A;
	SetMass(g4gen::NuclearMasses::GetNuclearMass(Z, A)*MeV);
}

void g4gen::Particle::SetNucleus(const G4String& symbol)
{
	g4gen::NuclearMasses::GetZAFromSymbol(symbol, &fZ, &fA);
	SetMass(g4gen::NuclearMasses::GetNuclearMass(fZ, fA)*MeV);
}

void g4gen::Particle::SetMass(G4double mass) 
{
	fM = mass;
	fP.set(0,0,0,fM+fEx);
}

void g4gen::Particle::SetEx(G4double ex) 
{
	fEx = ex;
	fP.set(0,0,0,fM+fEx);
}

void g4gen::Particle::SetP3(const G4ThreeVector& p)
{
	G4double E = sqrt(p.mag2() + pow(fM+fEx, 2));
	fP.set(p, E);
}

void g4gen::Particle::SetP3XYZ(G4double px, G4double py, G4double pz)
{
	SetP3(G4ThreeVector(px,py,pz));
}

void g4gen::Particle::SetP3ThetaPhi(G4double p, G4double theta, G4double phi)
{
	SetP3(G4ThreeVector(p*sin(theta)*cos(phi), 
											p*sin(theta)*sin(phi), 
											p*cos(theta)));
}

void g4gen::Particle::SetEkinThetaPhi(G4double ekin, G4double theta, G4double phi)
{
	G4double E = ekin+M()+Ex();
	G4double p = sqrt(E*E - M2());
	SetP3ThetaPhi(p,theta,phi);
}

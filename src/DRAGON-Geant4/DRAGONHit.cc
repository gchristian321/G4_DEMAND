#include "DRAGONHit.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

namespace DRAGON
{

G4ThreadLocal G4Allocator<DRAGONHit>* DRAGONHitAllocator = nullptr;

DRAGONHit::DRAGONHit()
    : x_fngr(0.), y_fngr(0.), z_fngr(0.), tofg(0.), edep(0.)
{}

DRAGONHit::~DRAGONHit() {}

DRAGONHit::DRAGONHit(const DRAGONHit& right)
    : G4VHit()
{
    x_fngr = right.x_fngr;
	y_fngr = right.y_fngr;
	z_fngr = right.z_fngr;
    tofg = right.tofg;
	edep = right.edep;
}

const DRAGONHit& DRAGONHit::operator=(const DRAGONHit& right)
{
    if (this != &right) {
        x_fngr = right.x_fngr;
		y_fngr = right.y_fngr;
		z_fngr = right.z_fngr;
        tofg = right.tofg;
        edep = right.edep;
    }
    return *this;
}

G4bool DRAGONHit::operator==(const DRAGONHit& right) const
{
    return (this == &right);
}

void DRAGONHit::Draw() const
{
 
}

void DRAGONHit::Print() const
{
    G4cout << " | x_fngr: " << x_fngr / cm << " cm"
		   << " | y_fngr: " << y_fngr / cm << " cm"
		   << " | z_fngr: " << z_fngr / cm << " cm"
		   << " | tofg: " << tofg / s << " s" 
           << " | Edep: " << edep / keV << " keV" << G4endl;
}

} 


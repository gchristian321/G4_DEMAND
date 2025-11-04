#include "DRAGONDigi.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"
#include "G4SystemOfUnits.hh"


namespace DRAGON
{

G4ThreadLocal G4Allocator<DRAGONDigi>* DRAGONDigiAllocator = nullptr;

DRAGONDigi::DRAGONDigi()
 : x_fngr(0.), y_fngr(0.), z_fngr(0.), edep(0.), tofg(0.) 
{}

DRAGONDigi::~DRAGONDigi() {}

DRAGONDigi::DRAGONDigi(const DRAGONDigi& right) : G4VDigi() {
    x_fngr = right.x_fngr;
	y_fngr = right.y_fngr;
	z_fngr = right.z_fngr;
	tofg = right.tofg;
    edep = right.edep;
}

const DRAGONDigi& DRAGONDigi::operator=(const DRAGONDigi& right) {
    if (this != &right) {
         x_fngr = right.x_fngr;
	     y_fngr = right.y_fngr;
	     z_fngr = right.z_fngr;
	     tofg = right.tofg;
         edep = right.edep;
    }
    return *this;
}

G4int DRAGONDigi::operator==(const DRAGONDigi& right) const {
    return (x_fngr == right.x_fngr && y_fngr == right.y_fngr && z_fngr == right.z_fngr &&
        tofg == right.tofg &&
        edep == right.edep);
}

void DRAGONDigi::Draw() {
    // Opcional: podrías dibujar algo con G4VVisManager si querés
}

void DRAGONDigi::Print() {
    G4cout << "Digi - x_fngr: " << x_fngr / cm << " cm"
	       << "Digi - y_fngr: " << y_fngr / cm << " cm"
	       << "Digi - z_fngr: " << z_fngr / cm << " cm"
           << "tofg: " << tofg / s << " s"
           << "edep: " << edep / keV << " keV" << G4endl;
}









}

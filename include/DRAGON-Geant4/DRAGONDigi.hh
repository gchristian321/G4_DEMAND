#ifndef DRAGONDigi_h
#define DRAGONDigi_h 1

#include "G4VDigi.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

namespace DRAGON
{

class DRAGONDigi : public G4VDigi {
public:
    DRAGONDigi();
    virtual ~DRAGONDigi();
    DRAGONDigi(const DRAGONDigi& right);
    const DRAGONDigi& operator=(const DRAGONDigi& right);
    G4int operator==(const DRAGONDigi& right) const;

    virtual void Draw() override;
    virtual void Print() override;

    G4double x_fngr, y_fngr, z_fngr;
	G4double tofg;
	G4double edep;

};

extern G4ThreadLocal G4Allocator<DRAGONDigi>* DRAGONDigiAllocator;

}

#endif

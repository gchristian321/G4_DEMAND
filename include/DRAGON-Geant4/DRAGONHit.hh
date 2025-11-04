#ifndef DRAGONHIT_HH
#define DRAGONHIT_HH

#include "G4VHit.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "G4String.hh"

namespace DRAGON
{

class DRAGONHit : public G4VHit {
public:
    DRAGONHit();
    DRAGONHit(const DRAGONHit& right);
    virtual ~DRAGONHit();

    const DRAGONHit& operator=(const DRAGONHit& right);
    G4bool operator==(const DRAGONHit& right) const;

    inline void* operator new(size_t);
    inline void operator delete(void* hit);

    void Setx_fngr(G4double x_fngr_){x_fngr = x_fngr_;}
	void Sety_fngr(G4double y_fngr_){y_fngr = y_fngr_;}
	void Setz_fngr(G4double z_fngr_){z_fngr = z_fngr_;}
    void Setedep(G4double edep_){edep = edep_;}
    void Settofg(G4double tofg_){tofg = tofg_;}
    void SetPMTNumber(G4int copyNo_){copyNo = copyNo_;}
    
	G4double Getx_fngr() const {return x_fngr;}
	G4double Gety_fngr() const {return y_fngr;}
	G4double Getz_fngr() const {return z_fngr;}   
    G4double Getedep() const {return edep;}
    G4double Gettofg() const {return tofg;}   	
    G4int GetPMTNumber() const {return copyNo;}

    virtual void Draw() const;
    virtual void Print() const;

private:
    G4double x_fngr, y_fngr, z_fngr;
	G4double tofg;
	G4double edep;
    G4int copyNo;
};

// Declaración externa del allocator
extern G4ThreadLocal G4Allocator<DRAGONHit>* DRAGONHitAllocator;

// Operadores new/delete inline
inline void* DRAGONHit::operator new(size_t)
{
    if (!DRAGONHitAllocator)
        DRAGONHitAllocator = new G4Allocator<DRAGONHit>;
    return (void*)DRAGONHitAllocator->MallocSingle();
}

inline void DRAGONHit::operator delete(void* hit)
{
    DRAGONHitAllocator->FreeSingle((DRAGONHit*)hit);
}

}

#endif


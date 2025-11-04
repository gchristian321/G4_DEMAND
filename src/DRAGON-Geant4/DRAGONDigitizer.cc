#include "DRAGONHit.hh"
#include "DRAGONDigi.hh"
#include "DRAGONDigitizer.hh"
#include "DRAGONDigitizerMessenger.hh"
#include "G4DigiManager.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include "G4SDManager.hh"
#include "G4THitsCollection.hh"
#include "G4TDigiCollection.hh"
#include "DRAGONDetectorConstruction.hh"
#include "DRAGONEventAction.hh"

namespace DRAGON
{

DRAGONDigitizer::DRAGONDigitizer(const G4String& name)
 : G4VDigitizerModule(name) {
    collectionName.push_back("DRAGONDigiCollection");

fDRAGONDigitizerMessenger = new DRAGONDigitizerMessenger(this);
uvinit();
}

DRAGONDigitizer::~DRAGONDigitizer() 
{
 delete fDRAGONDigitizerMessenger;
 }

void DRAGONDigitizer::Digitize() {
	
	gudigi();
	
	/*
    auto digiCollection = new G4TDigiCollection<DRAGONDigi>(GetName(), collectionName[0]);

    G4SDManager* sdManager = G4SDManager::GetSDMpointer();

    G4int dssdHCID = sdManager->GetCollectionID("DSSDHitsCollection");
    G4int pmtHCID  = sdManager->GetCollectionID("PMTHitsCollection");
	G4int scntHCID  = sdManager->GetCollectionID("SCNTHitsCollection");	

    G4HCofThisEvent* hce = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetHCofThisEvent();
    if (!hce) return;

    if (dssdHCID >= 0) {
        auto dssdHits = dynamic_cast<G4THitsCollection<DRAGONHit>*>(hce->GetHC(dssdHCID));
        if (dssdHits) {
            for (size_t i = 0; i < dssdHits->GetSize(); ++i) {
                DRAGONHit* hit = (*dssdHits)[i];
                if (hit->Getedep() < E_threshold) continue;
                auto digi = new DRAGONDigi();
                digi->x_fngr = 0.;
				digi->y_fngr = 0.;
				digi->z_fngr = 0.;
                digi->edep = hit->Getedep();
                digi->tofg = hit->Gettofg();
                digiCollection->insert(digi);
            }
        }
    }

    if (pmtHCID >= 0) {
        auto pmtHits = dynamic_cast<G4THitsCollection<DRAGONHit>*>(hce->GetHC(pmtHCID));
        if (pmtHits) {
            for (size_t i = 0; i < pmtHits->GetSize(); ++i) {
                DRAGONHit* hit = (*pmtHits)[i];
                if (hit->Getedep() < E_threshold) continue;
                auto digi = new DRAGONDigi();
                digi->x_fngr = 0.;
				digi->y_fngr = 0.;
				digi->z_fngr = 0.;
                digi->edep = hit->Getedep(); 
                digi->tofg = hit->Gettofg();
                digiCollection->insert(digi);
            }
        }
    }
	
    if (scntHCID >= 0) {
        auto scntHits = dynamic_cast<G4THitsCollection<DRAGONHit>*>(hce->GetHC(scntHCID));
        if (scntHits) {
            for (size_t i = 0; i < scntHits->GetSize(); ++i) {
                DRAGONHit* hit = (*scntHits)[i];
                if (hit->Getedep() < E_threshold) continue;
                auto digi = new DRAGONDigi();
                digi->x_fngr = 0.;
				digi->y_fngr = 0.;
				digi->z_fngr = 0.;
                digi->edep = hit->Getedep();
                digi->tofg = hit->Gettofg();
                digiCollection->insert(digi);
            }
        }
    }
    StoreDigiCollection(digiCollection);
	*/
}

void DRAGONDigitizer::uvinit()
{
//C.
//C.    *************** USER DEFAULTS ***************
//C.
E_threshold = 0.;
pmt_thrshld = 0.;
tot_thrshld = 0.;
}


void DRAGONDigitizer::SetTHLD(const G4String& input)
{
    std::istringstream iss(input);
    iss >> tot_thrshld >> pmt_thrshld;
}



void DRAGONDigitizer::genclu(G4int nelem, G4int jelem[][500] , G4int& nclu, G4int nhit[], G4int jhit[][10])
{
    const G4int mclu = 10;
    const G4int mhit = 10;
    const G4int mstack = 30;
//C.
//C.    nelem:              Number of 'struck' modules
//C.    jelem:              Volume # of 'struck' module
//C.
//C.    nstack:             Number of neighbours still on the stack
//C.    jstack[istack]:     Volume # of stack entry 'istack'
//C.
//C.    indx_stack[istack]: Index of stack member
//C.
//C.                        index = ilevel for stack member with 'hit'
//C.                        index = index(iseed) - 1 for member without 'hit'
//C. 
//C.    nhit[nclu]:         Number of neighbours in cluster's hit list
//C.    jhit[ihit][nclu]:    Volume # of entry on hit list
//C.
    G4int nstack;
    G4int jstack[mstack];
    G4int indx_stack[mstack];

    G4int index;
    G4int i, ii, j, l, n, iseed;
    G4int jnbr;
    G4int nbr[fDetector->max_hexagon];

    static bool first = true;
    bool skip;
    static G4int ilevel = 0;

    // Inicializar vecinos
    if(first) {
        first = false;
        i = 0;
        do {
            nbr[i] = 0;
            j = 1;
            do {
                if(fDetector->n_fngr[j][i] != 0) nbr[i] = nbr[i] + 1;
                j++;
            } while(j < fDetector->Nn);
            i++;
        } while(i < fDetector->max_hexagon);
    }

    nclu = 0;

    i = 0;
    do {
        if(i >= nelem) break;
		
        iseed = jelem[0][i];

        // Verificar si iseed ya está en un cluster
        if(nclu > 0) {
            n = 0;
            do {
                l = 0;
                do {
                    if(iseed == jhit[l][n]) goto _1000_;
                    l++;
                } while(l < nhit[n]);
                n++;
            } while(n < nclu);
        }

        if(nclu >= mclu) goto _999_;
        nclu = nclu + 1;
		
		nstack = 0;
		index = ilevel;
		
        nhit[nclu] = 1;
        jhit[nhit[nclu]][nclu] = iseed;

    _100_:
//C.
//C. Loop over all nearest neighbours 
//C.                       (Note: n_fngr(1,iseed) is iseed module itself)
//C.
        j = 0;
        do {
            if(j >= nbr[iseed]) break;
			
            jnbr = fDetector->n_fngr[j+1][iseed];
//C.
//C. Loop over all members already in the hit list
//C.
            l = 0;
            do {
                if(jnbr == jhit[l][nclu]) goto _next_j_;
                l++;
            } while(l < nhit[nclu]);
//C.
//C. Loop over all members already on the stack
//C.
            l = 0;
            do {
                if(jnbr == jstack[l]) { skip = true; break; }
                l++;
            } while(l < nstack);
//C.
//C.     If the neighour is not already on the hit list or on the stack:
//C.
//C.  ** Put it on the hit list AND on the stack      if it has a  'hit' **
//C.  ** Put it                     on the stack ONLY if it has no 'hit' **
//C.  ** Put it                     on the stack ONLY if seed-index > 0  **
//C. 
            ii = 0;
            do {
                if(jnbr == jelem[0][ii]) {
                    if(nhit[nclu] >= mhit) goto _999_;
                    nhit[nclu] = nhit[nclu] + 1;
                    jhit[nhit[nclu]][nclu] = jnbr;

                    if(nstack >= mstack) goto _999_;
					nstack = nstack + 1;
                    jstack[nstack] = jnbr;
                    indx_stack[nstack] = ilevel;
                    
                    goto _next_j_;
                }
                ii++;
            } while(ii < nelem);

                if(nstack >= mstack) goto _999_;
                if(index > 0) {
					nstack = nstack + 1;
                    jstack[nstack] = jnbr;
                    indx_stack[nstack] = index - 1;               
                }
_next_j_:
            j++;
        } while(j < nbr[iseed]);

        // Tomar siguiente de stack
        if(nstack > 0) {
            iseed = jstack[nstack-1];
            index = indx_stack[nstack-1];
            nstack = nstack - 1;
            goto _100_;
        }

        i++;
    _1000_: ;
    } while(i < nelem);

    return;

_999_:
    std::cout << "*** Problem in GENCLU *** for event: " << fEventAction->ievent << std::endl;
    std::cout << "Number of clusters nclu (<=10)       " << nclu << std::endl;
    std::cout << "Number of elements in cluster (<=10) " << nhit[nclu] << std::endl;
    std::cout << "Size of search stack (<=30)          " << nstack << std::endl;

    return;
}


}

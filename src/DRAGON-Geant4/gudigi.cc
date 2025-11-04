#include "G4AnalysisManager.hh"
#include "Randomize.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"
#include "G4DigiManager.hh"
#include "DRAGONHit.hh"
#include "G4SDManager.hh"
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
#include "DRAGONEventAction.hh"
#include "DRAGONDetectorConstruction.hh"

#include "geant3functions.hh"  
#include "global_variables.hh" 

namespace DRAGON {


void DRAGONDigitizer::gudigi()
{ 

//C.
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//C                                                                      C
//C     SUBROUTINE GUDIGI is called at the end of each event             C
//C                                                                      C
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//C.

      G4int itrig;

      itrig = 0;

      anadet(itrig);

      if(itrig == 0);
//CCC        ieotri = 1
      else
         dhpmt();
	  interac();

}


void DRAGONDigitizer::anadet(G4int itrig)    
{

//C.
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//C                                                                      C
//C     SUBROUTINE ANADET reads the scintillator hits and determines     C
//C     whether or not an event has fulfilled the trigger condition      C
//C                                                                      C
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//C.

      G4SDManager* sdManager = G4SDManager::GetSDMpointer();
	  
      G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	  auto gv = GlobalVariables::GetInstance();

      G4int i, j, k, n;
	  G4int jhits;
	  G4double rnd[2];

      const static G4int nvdim = 1;
      const static G4int nhdim = 5; 
      const static G4int nhmax = 500;

      G4int numbv[nvdim][nhmax];

      G4int nhits, itra[nhmax];
      G4double hits[nhdim][nhmax];

      const static G4int mhit = 500;
      
      G4double rndm[2];
//C.
//C *** nelem:    	Number of all 'struck' scintillators
//C *** melem:            Number of 'struck' modules of detector
//C *** jelem:            Volume # of 'struck' modules
//C.
      G4int nelem, melem;
      G4int jelem[nvdim][mhit], itr[mhit];
      G4int isort[mhit];

      G4double coord[3][mhit], time[mhit], energy[mhit], fcoord[3][mhit];

      G4double E_energy;
      G4double rndvec[0], resn, thresh;

      G4int iclu;

      G4int index[mclu];
      G4double cthres = 0.0;
      
      G4double etmp[mclu];
//C JS adds e_bgos_deposit to record energy deposited in each BGO,
//C   one position in array per BGO
      G4double e_bgos_deposit[30][2];
      
      G4int num_bgos_hit_ab;
      
      G4HCofThisEvent* hce = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetHCofThisEvent();	  
      if (!hce) return;
         
      G4int scntHCID  = sdManager->GetCollectionID("SCNTHitsCollection");

      G4THitsCollection<DRAGONHit>* scntHits = dynamic_cast<G4THitsCollection<DRAGONHit>*>(hce->GetHC(scntHCID));
      if (!scntHits || scntHits->GetSize() == 0) goto _999_;

      nhits = scntHits->GetSize();
			 
      if (nhits == 0)goto _100_;
                      
      vzero(&jelem[0][0],nvdim*mhit);
      vzero(itr,mhit);
      vzero(&coord[0][0],3*mhit);
      vzero(&fcoord[0][0],3*mhit);
      vzero(time,mhit);
      vzero(energy,mhit);

      melem = 0;

      E_energy = 0.0;
//C.
//C     Threshold: testing of threshold effect (initialised in uvinit.f or 
//C     dragon_2003.ffcards).
      rndvec[0] = G4RandGauss::shoot(0., 1.);
      thresh = E_threshold*(1. + rndvec[0]*0.05);
            			 
      if (nhits > nhmax) 
         {
          G4cout << "Problem in DHPMT for event: " << fEventAction->ievent << G4endl;
          G4cout << "Number of hits " << nhits << " > nhmax " << nhmax << G4endl;
          nhits = std::min(nhmax, nhits);
          }

      j = 0;
	  do
        {
	     k = 0;
         do 
		   {
		    DRAGONHit* hit = (*scntHits)[j];
		    numbv[0][j] = hit->GetPMTNumber(); 
		    if(numbv[0][j] == jelem[0][k])
		      {
            	        coord[0][k] = coord[0][k] + hits[4][j] * hits[0][j];   //!summing hits of
            		coord[1][k] = coord[1][k] + hits[4][j] * hits[1][j];   //!same module.
            		coord[2][k] = coord[2][k] + hits[4][j] * hits[2][j];
            		fcoord[0][k] = hits[0][j];    //!AO real finger coords
            		fcoord[1][k] = hits[1][j];
            		fcoord[2][k] = hits[2][j];

            		time[k] = std::min(time[k], hits[3][j]);
            		energy[k] =  energy[k] + hits[4][j];
            		rndvec[0] = G4RandGauss::shoot(0.0, 1.0);
//C     gaussian distributed resolution dependent on energy: 
//C     Resolution based on linear fit to values of 14% for 661 keV and 
//C     4% for 10 MeV.(Temporary until I find the real resolution function).
            		resn = 0.14 - 0.001 * energy[k];
            		energy[k] = energy[k]*(1.0 + rndvec[0] * resn);

            		goto _10_;
            		k++;
                        }
                      } while (k < melem);  
		 
		  if(melem >= mhit)goto _998_;
          melem = melem + 1;

          ucopy(&numbv[0][j],&jelem[0][melem-1], 1);

          itr[melem] = itra[j];
          coord[0][melem] = hits[4][j]*hits[0][j];
          coord[1][melem] = hits[4][j]*hits[1][j];
          coord[2][melem] = hits[4][j]*hits[2][j];
          fcoord[0][melem] = hits[0][j];    //!AO real finger coords
          fcoord[1][melem] = hits[1][j]; 
          fcoord[3][melem] = hits[2][j];
          time[melem] = hits[3][j];
          energy[melem] = hits[4][j];
		  } while (j < nhits);  
      
_10_:

_100_:

//C.--> Sorts 'energy' vector in terms of decreasing energy (so that 
//C.--> 'energy(isort[0])' is highest E.
      if(melem == 0)return;

      sortzv(energy,isort,melem,1,1,0);

      k = 0;
      do {
          fEventAction->e_detect = fEventAction->e_detect + energy[k];

          coord[0][k] = coord[0][k]/energy[k];
          coord[1][k] = coord[1][k]/energy[k];
          coord[2][k] = coord[2][k]/energy[k];

          analysisManager->FillH1(gv->IDMap[33], fcoord[0][k]);
          analysisManager->FillH1(gv->IDMap[34], fcoord[1][k]);
          analysisManager->FillH1(gv->IDMap[35], fcoord[2][k]);
          analysisManager->FillH2(gv->IDMap[25], time[melem], z_react);
          analysisManager->FillH1(gv->IDMap[36], energy[k]);

         if (k < 4) 
             analysisManager->FillH1(gv->IDMap[36 + (k+1)], energy[isort[k]]);
         k++;
         } while (k < melem);

      fEventAction->x_max = coord[0][isort[0]];    //!energy weighted coords of the 
      fEventAction->y_max = coord[1][isort[0]];    //!highest energy module
      fEventAction->z_max = coord[2][isort[0]];
//C. 
//C.--> Extract time variable of highest E gamma hit
      gammatof = time[isort[0]];
//C.      print*, 'TOF Gamma = ', gammatof
//C.---> Add some timing 'resolution' (say sigma=200ps)[time in ns]
      rnd[0] = G4UniformRand();
      rnd[1] = G4UniformRand();
      gammatof = gammatof + rndm[0]*0.3;

      fEventAction->ifngr_max = jelem[0][isort[0]];
      neff[11][fEventAction->ifngr_max] = neff[11][fEventAction->ifngr_max] + 1;

      melem_gbox = melem;           //!melem = no. of modules hit
      k = 0;
      do {
          E_energy = E_energy + energy[k];

          jelem_gbox[0][k] = jelem[0][isort[k]];  
          energy_gbox[k]   = energy[isort[k]];

          k++;
          } 
      while (k < melem);

      if (fEventAction->GetRunAction()->GetIswit()[3] == 1) 
         { 
          G4cout << melem_gbox << G4endl;
          G4int k = 0;
          do 
            {
             G4cout << jelem_gbox[0][k] << " " << energy_gbox[k] << G4endl;
             k++;
             }
          while (k < melem_gbox);       //!store ordered by energy 
          } 
      else if (fEventAction->GetRunAction()->GetIswit()[3] == 2) 
              ;//{hfnt(999);}    //OJO

      nelem = 0;

      if(E_energy > 0.0)nelem = nelem + 1;

      analysisManager->FillH1( gv->IDMap[41],E_energy);

      if(nelem > 0)itrig = itrig + 1;

      if(nelem == 0) return;

      analysisManager->FillH1(gv->IDMap[31],1.*melem+0.5);
      rndvec[0] = G4RandGauss::shoot(0,1);
      
//C     Total energy.
      analysisManager->FillH1(gv->IDMap[32],fEventAction->e_detect);

      i = 0;
      do { 
          E_energy = 0.0;

          G4int j = 0;
          do { 
              if (fDetector->n_fngr[j][fEventAction->ifngr_max] != 0) 
                 {
                  G4int k = 0;
                  G4bool found = false;

                  do 
                   {
                    if (fDetector->n_fngr[j][fEventAction->ifngr_max] == jelem[0][isort[k]]) 
                       {
                        E_energy += energy[isort[k]];
                        found = true;  
                        }
                    k++;
                    } while (k < melem && !found); 

                  if (found) break; 
                  } 
               else 
                  {break;}
               j++;
               } while (j < i);  

           analysisManager->FillH2(gv->IDMap[130 + (i + 1)], 1.0*fEventAction->ifngr_max + 0.5, E_energy);

           if (E_energy > thresh) 
              {neff[i][fEventAction->ifngr_max] += 1;}

           i++;
           } while (i < fDetector->Nn);  

      genclu(melem,jelem,nclu,nhit,jhit);   //! generate clusters of modules

      vzero(eclu,mclu);
      vzero(etmp,mclu);
      vzero(xclu,mclu);
      vzero(yclu,mclu);
      vzero(zclu,mclu);

      iclu = 0;
      do { 
          if (iclu >= nclu) break;

          j = 0;
          do {
              if (j >= nhit[iclu]) break;

              k = 0;
              do 
               { 
                if (k >= melem) break;

                if (jelem[0][k] == jhit[j][iclu]) 
                   {
                    etmp[iclu] += energy[k];
                    xclu[iclu] += fDetector->x_fngr[jhit[j][iclu]] * energy[k];
                    yclu[iclu] += fDetector->y_fngr[jhit[j][iclu]] * energy[k];
                    zclu[iclu] += fDetector->z_fngr[jhit[j][iclu]] * energy[k];
                    }

                k++;
                } while (k < melem); 

             j++;
             } while (j < nhit[iclu]); 

          iclu++;
          } while (iclu < nclu); 

     n = nclu;
     nclu = 0;

     iclu = 0;
     do { 
         if (iclu >= n) break;

         if (etmp[iclu] > cthres) 
            {
             nclu++;
             etmp[nclu-1] = etmp[iclu];
             xclu[nclu-1] = xclu[iclu] / etmp[iclu];
             yclu[nclu-1] = yclu[iclu] / etmp[iclu];
             zclu[nclu-1] = zclu[iclu] / etmp[iclu];
             nhit[nclu-1] = nhit[iclu];

             //ucopy2(jhit[iclu], jhit[nclu-1], nhit[iclu]); 
             }

         iclu++;
         } while (iclu < n); 

 
 //    sortzv(etmp, index, nclu, 1, -1, 0); 

     analysisManager->FillH1(gv->IDMap[90], 1.0*nclu + 0.5);

     iclu = 0;
     do 
       { 
        if (iclu >= nclu) break;

        eclu[iclu] = etmp[index[iclu]];
        dir_clu[0][iclu] = xclu[index[iclu]];
        dir_clu[1][iclu] = yclu[index[iclu]];
        dir_clu[2][iclu] = zclu[index[iclu]];

        vunit(&dir_clu[0][iclu],3); 

        if (iclu < 3) 
           {analysisManager->FillH1(gv->IDMap[90 + (iclu+1)], eclu[iclu]);}

        iclu++;
        } while (iclu < nclu); 
//C
//C     JS code to count # of BGOs triggered
//C     and energy deposited in triggered BGOs.
//C
//C     reset variables from previous event
      i = 0;
      do 
       {
        if (i >= 30) break;

        e_bgos_deposit[i][0] = 0;
        e_bgos_deposit[i][1] = 0;

        i++; 
        } 
      while (i < 30); 

//C     loop over all hits, find out which BGO recorded the hit, numbv(1,i), and
//C     then store the energy from that hit, hits(5,i), in the array, at the
//C     index position of the hit BGO
     i = 0;
     do 
      { 
       if (i >= nhits) break;

       G4int idx = numbv[0][i]; 

       e_bgos_deposit[idx][0] = idx;
       e_bgos_deposit[idx][1] += hits[4][i];  

       i++; 
       } 
      while (i < nhits);  

//C     Sort the array with insertion sort, bringing most energetic BGO to
//C     first position, second most energetic BGO to second position, etc.
//C     Any BGO hits less than threshold (0.1 MeV) are removed.
      if (e_bgos_deposit[0][1] < 0.1)
         {
          e_bgos_deposit[0][0] = 0;
          e_bgos_deposit[0][1] = 0;
          }
// Variables supuestas:
// e_bgos_deposit[30][2], num_bgo_first_ab, e_bgo_first_ab
// num_bgos_hit, e_bgos_total

      i = 1; 
      do 
       {
        if (i >= 30) break;

        if (e_bgos_deposit[i][1] < 0.1) 
           {
            e_bgos_deposit[i][0] = 0;
            e_bgos_deposit[i][1] = 0;
            i++;
            continue; 
            }

       j = i;
       do 
        {
         if (j < 1) break;
         if (e_bgos_deposit[j][1] > e_bgos_deposit[j-1][1]) 
            {
             fEventAction->num_bgo_first_ab = e_bgos_deposit[j-1][0];
             fEventAction->e_bgo_first_ab = e_bgos_deposit[j-1][1];
             e_bgos_deposit[j-1][0] = e_bgos_deposit[j][0];
             e_bgos_deposit[j-1][1] = e_bgos_deposit[j][1];
             e_bgos_deposit[j][0] = fEventAction->num_bgo_first_ab;
             e_bgos_deposit[j][1] = fEventAction->e_bgo_first_ab;
             }
         j--;
         } 
       while (j >= 1);

       fEventAction->num_bgo_first_ab = 0;
       fEventAction->e_bgo_first_ab = 0;

       i++;
       }
      while (i < 30);

 
      i = 0; 
      do 
       {
        if (i >= 30) break;

        if (e_bgos_deposit[i][0] != 0) 
           {
            fEventAction->num_bgos_hit += 1;
            fEventAction->e_bgos_total += e_bgos_deposit[i][1];
            } 
        else 
           {break;}

        i++;
        } while (i < 30);

      fEventAction->num_bgo_first = e_bgos_deposit[0][0];
      fEventAction->e_bgo_first = e_bgos_deposit[0][1];
//C
//C	Add resolution function of form FWHM = h*sqrt(E)
//C
      rndm[0] = G4RandGauss::shoot(0.0, 1.0);
      rndm[1] = G4RandGauss::shoot(0.0, 1.0);

      fEventAction->e0_conv = fEventAction->e_bgo_first + rndm[0] * 0.1733 * std::sqrt(fEventAction->e_bgo_first) / 2.35;

      fEventAction->num_bgo_second = e_bgos_deposit[1][0]; // Fortran 2->C++ 1
      fEventAction->e_bgo_second = e_bgos_deposit[1][1];

      num_bgos_hit_ab = fEventAction->num_bgos_hit;

      if (fEventAction->num_bgos_hit != 1) 
         {
          G4int i = 0;
          do 
           {
            if (i >= fEventAction->num_bgos_hit - 1) break;
 
            G4int j = i + 1;
            do 
             { 
              if (j >= fEventAction->num_bgos_hit) break;

              if (e_bgos_deposit[i][0] != 0 && e_bgos_deposit[j][0] != 0) 
                 {
                  if (fDetector->adjacency_matrix[static_cast<G4int>(e_bgos_deposit[i][0])][static_cast<G4int>(e_bgos_deposit[j][0])] == 1) 
                     {
                      e_bgos_deposit[i][1] += e_bgos_deposit[j][1];
                      e_bgos_deposit[j][0] = 0;
                      e_bgos_deposit[j][1] = 0;
                      num_bgos_hit_ab--;
                      }
                  }

              j++;
              } 
             while (j < fEventAction->num_bgos_hit);

           i++;
           } 
          while (i < fEventAction->num_bgos_hit - 1);
          }

      i = 0;
      do 
       { 
        if (i >= fEventAction->num_bgos_hit) break;

        if (e_bgos_deposit[i][0] != 0) 
           {
            if (fEventAction->num_bgo_first_ab == 0) 
               {
                fEventAction->num_bgo_first_ab = e_bgos_deposit[i][0];
                fEventAction->e_bgo_first_ab = e_bgos_deposit[i][1];
                } 
            else if (fEventAction->num_bgo_second_ab == 0) 
                    {
                     fEventAction->num_bgo_second_ab = e_bgos_deposit[i][0];
                     fEventAction->e_bgo_second_ab = e_bgos_deposit[i][1];
                     break; 
                     }
            }

         i++;
         }
        while (i < fEventAction->num_bgos_hit);

        if (fEventAction->recoil_hit_ENDV == 1) 
           {
            analysisManager->FillH1(gv->IDMap[511], 1.0 * fEventAction->num_bgos_hit);
            analysisManager->FillH1(gv->IDMap[512], 1.0 * fEventAction->e_bgos_total);
            analysisManager->FillH2(gv->IDMap[513], 1.0 * fEventAction->num_bgos_hit, fEventAction->e_bgos_total);
            analysisManager->FillH2(gv->IDMap[514], fEventAction->e_bgo_first, fEventAction->e_bgo_second);
            }

//C     For debugging purposes
//C      std::cout << nhits << std::endl;
//C      std::cout << numbv << std::endl;
//C      std::cout << e_bgos_deposit << std::endl;
//C      std::cout << e_bgos_deposit_total << std::endl;
//C      std::cout << fEventAction->num_bgos_hit << std::endl;
//C.
//C *** nelem:  Total Number of Scintillator hits
//C.
//C *** jelem:  Array of Module indices (k=1,melem)
//C.
//C.                         jelem(1,k) -> Module    Number
//C.
//C *** itr[k]: Track number causing the kth GEANT hit
//C.
//C *** coord:  Coordinates Array of Scintillator [cm]
//C.
//C             coord[0][k]: Energy weighted x-coordinate
//C             coord[2][k]: Energy weighted y-coordinate
//C             coord[3][k]: Energy weighted z-coordinate
//C.
//C *** time:   min(tof) - abs. GEANT time of flight [ns] (from event time0)
//C.
//C *** energy: Total energy [MeV] deposited in module
//C.
//C *** E_energy: Total energy [MeV] deposited in scintillator
//C.
      return;

  _998_:

      std::cout << " *** Problem in ANADET *** for event: " << fEventAction->ievent << std::endl;
      std::cout << " Number of elements melem " << melem << " maximum " << std::endl;

      return;

_999_: 

      std::cout << " *** Problem in ANADET --- STOP!!! *** " << std::endl;
      std::cout << " GEANT HIT Bank not properly set up! " << std::endl;
	 
}




void DRAGONDigitizer::dhpmt()  
{
//C.
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//C                                                                      C
//C     SUBROUTINE DHPMT reads the PMT hits, digitizes and records ADC   C
//C     and TDC hits in the PMTs                                         C
//C                                                                      C
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//C.

      G4SDManager* sdManager = G4SDManager::GetSDMpointer();

      G4int pmtHCID  = sdManager->GetCollectionID("PMTHitsCollection");
  
	  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	  auto gv = GlobalVariables::GetInstance();

      G4int j, k;

      G4int ndet, idet;

      const static G4int nvdim = 1;
      const static G4int nhdim = 5;
      const static G4int nhmax = 10000;

      G4int numvs[nvdim] = {}, numbv[nvdim][nhmax];

      G4int nhits, itra[nhmax];
      G4double hits[nhdim][nhmax];

      G4String iudet;

      const static G4int mhit = 20;
//C.
//C *** nelem:    	Number of all 'struck' PMTs
//C *** jelem:            Volume # of 'struck' PMT
//C.
      G4int nelem;
      G4int jelem[nvdim][mhit], itr[mhit];
      G4int isort[mhit];
	  G4int jhits;

      G4double coord[3][mhit], time[mhit], energy[mhit];
	  
	  G4HCofThisEvent* hce = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetHCofThisEvent();	  
          if (!hce) return;

         auto pmtHits = dynamic_cast<G4THitsCollection<DRAGONHit>*>(hce->GetHC(pmtHCID));
        if (!pmtHits || pmtHits->GetSize() == 0) return;

      vzero(&jelem[0][0],nvdim*mhit);
      vzero(&coord[0][0],3*mhit);
      vzero(time,mhit);
      vzero(energy,mhit);

      nelem = 0;
      fEventAction->edetect = 0.0;

         if (pmtHCID >= 0) 
		    {
             auto pmtHits = dynamic_cast<G4THitsCollection<DRAGONHit>*>(hce->GetHC(pmtHCID));
			 nhits = pmtHits->GetSize();
			 
			 if (nhits == 0)goto _100_;
			 
             if (nhits > nhmax) 
                {
                 G4cout << "Problem in DHPMT for event: " << fEventAction->ievent << G4endl;
                 G4cout << "Number of hits " << nhits << " > nhmax " << nhmax << G4endl;
                 nhits = std::min(nhmax, nhits);
                 }
                      
					  j = 0;
                      do 
                        { 
                         if (j >= nhits) break;
						 
						 DRAGONHit* hit = (*pmtHits)[j];
						 numbv[0][j] = hit->GetPMTNumber(); 
						 hits[0][j] = hit->Getx_fngr();
						 hits[1][j] = hit->Gety_fngr();
						 hits[2][j] = hit->Getz_fngr();
						 hits[3][j] = hit->Getedep();
						 hits[4][j] = hit->Gettofg();
						 
                         std::cout << hit->Getedep() << std::endl;

                         G4bool found = false;
                         k = 0;
                         do 
                           { 
                            if (k >= nelem) break;

                            if (numbv[0][j] == jelem[0][k]) 
                               {
                                coord[0][k] = coord[0][k] + hits[3][j] * hits[0][j];
                                coord[1][k] = coord[1][k] + hits[3][j] * hits[1][j];
                                coord[2][k] = coord[2][k] + hits[3][j] * hits[2][j];
                                time[k] = std::min(time[k], hits[3][j]);
                                energy[k] = energy[k] + hits[4][j];
                                found = true;
                                break; 
                                }
                            k++;
                            } 
                         while (k < nelem);

                         if (!found) 
                            {
                             if (nelem >= mhit) goto _998_;
                             nelem++;

                             ucopy(&numbv[0][j], &jelem[0][nelem-1], 1);

                             coord[0][nelem-1] = hits[4][j] * hits[0][j];
                             coord[1][nelem-1] = hits[4][j] * hits[1][j];
                             coord[2][nelem-1] = hits[4][j] * hits[2][j];

                             time[nelem-1] = hits[3][j];
                             energy[nelem-1] = hits[4][j];
                             }

                         j++;
                        }
                       while (j < nhits); 
                       }
                 

_100_:

     fEventAction->ntot = nelem;
     if (nelem == 0) return;

     fEventAction->x_mean = 0.0;
     fEventAction->y_mean = 0.0;
     fEventAction->z_mean = 0.0;

     sortzv(energy, isort, nelem, 1, 1, 0);

     fEventAction->ntot = 0;

     k = 0;
     do 
      { 
       if (k >= nelem) break;

       if (energy[k] > pmt_thrshld) fEventAction->ntot++;

       fEventAction->edetect += std::max(energy[k] - pmt_thrshld, 0.0);
       fEventAction->x_mean += (coord[0][k]/energy[k]) * std::max(energy[k] - pmt_thrshld, 0.0);
       fEventAction->y_mean += (coord[1][k]/energy[k]) * std::max(energy[k] - pmt_thrshld, 0.0);
       fEventAction->z_mean += (coord[2][k]/energy[k]) * std::max(energy[k] - pmt_thrshld, 0.0);

       analysisManager->FillH1(gv->IDMap[51], energy[k]);
       if (k < 4) 
          analysisManager->FillH1(gv->IDMap[51 + (k + 1)], energy[isort[k]]);

       k++;
       } 
      while (k < nelem);

     analysisManager->FillH1(gv->IDMap[56],1.*fEventAction->ntot+0.5);
     analysisManager->FillH1(gv->IDMap[57],fEventAction->edetect);

     if(fEventAction->nloss[0] > 0)analysisManager->FillH1( gv->IDMap[72],1.*fEventAction->nloss[0]+0.5);
     if(fEventAction->nloss[1] > 0)analysisManager->FillH1( gv->IDMap[73],1.*fEventAction->nloss[1]+0.5);
     if(fEventAction->nloss[2] > 0)analysisManager->FillH1( gv->IDMap[74],1.*fEventAction->nloss[2]+0.5);
     if(fEventAction->nloss[3] > 0)analysisManager->FillH1( gv->IDMap[75],1.*fEventAction->nloss[3]+0.5);
     if(fEventAction->nloss[4] > 0)analysisManager->FillH1( gv->IDMap[76],1.*fEventAction->nloss[4]+0.5);
     if(fEventAction->nloss[5] > 0)analysisManager->FillH1( gv->IDMap[77],1.*fEventAction->nloss[5]+0.5);

     if(fEventAction->ntot == 0)return;

     fEventAction->x_mean = fEventAction->x_mean/fEventAction->edetect;
     fEventAction->y_mean = fEventAction->y_mean/fEventAction->edetect;
     fEventAction->z_mean = fEventAction->z_mean/fEventAction->edetect;
//C.
//C *** nelem:  Total Number of PMT  hits
//C.
//C *** jelem:  Array of PMT indices (k=1,nelem)
//C.
//C.                         jelem(1,k) -> PMT Number
//C.
//C *** coord:  Coordinates Array of PMT hits [cm]
//C.
//C             coord(1,k): Energy weighted x-coordinate
//C             coord(2,k): Energy weighted y-coordinate
//C             coord(3,k): Energy weighted z-coordinate
//C.
//C *** time:   min(tof) - abs. GEANT time of flight [ns] (from event time0)
//C.
      return;

  _998_: 

      std::cout << " *** Problem in DHPMT *** for event: " << fEventAction->ievent << std::endl;
      std::cout << " Number of elements nelem "<< nelem << ", maximum " << std::endl;

      return;

  _999_: 

      std::cout << " *** Problem in DHPMT --- STOP!!! *** " << std::endl;
      std::cout << " GEANT HIT Bank not properly set up! " << std::endl;
}



}

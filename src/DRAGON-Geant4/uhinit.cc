#include "DRAGONHistoManager.hh"
#include "DRAGONPrimaryGeneratorAction.hh"
#include "DRAGONDetectorConstruction.hh"
#include "DRAGONPhysicsList.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"

#include "Materials.hh"


namespace DRAGON {


void DRAGONHistoManager::uhinit()
{
 std::cout << "Creating histograms" << std::endl;

 //CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 //C                                                                      C
 //C         Defines histogram/scatterplot definitions                    C
 //C                                                                      C
 //CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 if(!macro_read)
 {
	 auto materials = Materials::Instance();
	 auto gv = GlobalVariables::GetInstance();

        G4int n, nstrip;
	G4String strip;
	const std::string chtags[] = {"GML0", "GML1", "GML2", "GML3", "GML4", "GML5", "GML6", 
		                       "GML7", "GML9", "GML10", "GML11", "GML12", "GML13", "GML14", "GML15", 
		                       "GML16", "GML18", "GML19", "GML20", "GML21", "GML22", "GML23", "GML24",
		                       "GML25", "GML27", "GML28", "GML29",
                               "HML0", "HML1", "HML2", "HML3", "HML4", "HML5", "HML6", 
                               "HML7", "HML8", "HML9"};
	const G4String num[16] = {"01", "02", "03", "04", "05", "06", "07", "08", "09", 
		                      "10", "11", "12", "13", "14", "15", "16"};
	
//C.
//C -->   Open and create a ROOT file
//C.	

  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  analysis->SetVerboseLevel(1);
  analysis->SetActivation(true); 
  //analysis->SetHistoDirectoryName("histos"); 
  analysis->SetDefaultFileType("root");

  const G4Run* run = G4RunManager::GetRunManager()->GetCurrentRun();
  G4int idrun = run->GetRunID();
  G4String rootfilename = "dragon" + std::to_string(idrun+1);
  analysis->SetFileName(rootfilename);  
  G4bool fileOpened = analysis->OpenFile();  

	if(!fileOpened)
	  {
	   G4cerr << " Error: Could not open " << rootfilename << " ROOT file " << std::endl;
	   G4cerr << " Error: Bad return from OpenFile! " << std::endl;
	   std::exit(0);
	   }

//C.
//C -->   Initialize user HBOOK histograms and scatterplots
//C.
	n = 0;

        analysis->SetFirstHistoId(1); 
	G4int id = analysis->CreateH1("h"+std::to_string(n+1)," Initial - x - ", 100, -2.0, 2.0);
	gv->IDMap[n+1] = id; 
	id = analysis->CreateH1("h"+std::to_string(n+2)," Initial - y - ", 100, -2.0, 2.0);
	gv->IDMap[n+2] = id;
	id = analysis->CreateH1("h"+std::to_string(n+3)," Initial - dx - ", 100, -100.0, 100.0);
	gv->IDMap[n+3] = id;
	id = analysis->CreateH1("h"+std::to_string(n+4)," Initial - dy - ", 100, -100.0, 100.0);
	gv->IDMap[n+4] = id;
	id = analysis->CreateH1("h"+std::to_string(n+5)," IniFin - x - Stops ", 100, -2.0, 2.0);
	gv->IDMap[n+5] = id;
	id = analysis->CreateH1("h"+std::to_string(n+6)," IniFin - y - Stops ", 100, -2.0, 2.0);
	gv->IDMap[n+6] = id;
	id = analysis->CreateH1("h"+std::to_string(n+7)," IniFin - dx - Stops ", 100, -100.0, 100.0);
	gv->IDMap[n+7] = id;
	id = analysis->CreateH1("h"+std::to_string(n+8)," IniFin - dy - Stops ", 100, -100.0, 100.0);
	gv->IDMap[n+8] = id;
	id = analysis->CreateH1("h"+std::to_string(n+9)," Initial Momentum ", 100, -5.0, 5.0);
	gv->IDMap[n+9] = id;
	id = analysis->CreateH1("h"+std::to_string(n+10)," Momentum spread (%) ", 100, -5.0, 5.0);
	gv->IDMap[n+10] = id;
	id = analysis->CreateH1("h"+std::to_string(n+11)," Final - x - ", 16, -2.4, 2.4);
	gv->IDMap[n+11] = id;
	id = analysis->CreateH1("h"+std::to_string(n+12)," Final - y - ", 16, -2.4, 2.4);
	gv->IDMap[n+12] = id;
	id = analysis->CreateH1("h"+std::to_string(n+13)," Final - dx - ", 100, -100., 100.);
	gv->IDMap[n+13] = id;
	id = analysis->CreateH1("h"+std::to_string(n+14)," Final - dy - ", 100, -50., 50.);
	gv->IDMap[n+14] = id;
	id = analysis->CreateH1("h"+std::to_string(n+15)," Final Energy ", 4000, 0., 20.);
	gv->IDMap[n+15] = id;
	id = analysis->CreateH1("h"+std::to_string(n+16)," Stop Length (cm) ", 2000, 0., 2000.);
	gv->IDMap[n+16] = id;
        //analysis->SetFirstHistoId(n+17);
        id = analysis->CreateH2("h"+std::to_string(n + 17), " X vs Stop Length (cm) ", 50, -20.0, 20.0, 200, 0.0, 2000.0);
		gv->IDMap[n+17] = id;

	n = 20;
       // analysis->SetFirstHistoId(n+1);
	id = analysis->CreateH1("h"+std::to_string(n+1)," True photon energy     ", 2000, 0., 20.);  
	gv->IDMap[n+1] = id;
	id = analysis->CreateH1("h"+std::to_string(n+2)," True photon pol. angle ", 2000, 0., 200.);
	gv->IDMap[n+2] = id;
    id = analysis->CreateH1("h"+std::to_string(n+3)," Photon conv. module   ", 29, 1., 30.);
	gv->IDMap[n+3] = id;
	id = analysis->CreateH2("h"+std::to_string(n + 4), " Photon creation time vs z_react", 600, 0.0, 300.0, 300, -15.0, 15.0);
	gv->IDMap[n+4] = id;
	id = analysis->CreateH2("h"+std::to_string(n + 5), " Photon detection time vs z_react", 600, 0.0, 300.0, 300, -15.0, 15.0);
	gv->IDMap[n+5] = id;
        //analysis->SetFirstHistoId(n+11);
	id = analysis->CreateH1("h"+std::to_string(n+11)," No. of Modules hit  ", 10, 0., 10.);
    gv->IDMap[n+11] = id;	
	id = analysis->CreateH1("h"+std::to_string(n+12)," Total energy dep.   ", 200, 0., 20.);
	gv->IDMap[n+12] = id;
	id = analysis->CreateH1("h"+std::to_string(n+13)," x-coordinates of hit", 60, -15., 15.);
	gv->IDMap[n+13] = id;
	id = analysis->CreateH1("h"+std::to_string(n+14)," y-coordinates of hit", 80, -20., 20.);   
    gv->IDMap[n+14] = id;    
	id = analysis->CreateH1("h"+std::to_string(n+15)," z-coordinates of hit", 100, -20., 20.);
	gv->IDMap[n+15] = id;
	id = analysis->CreateH1("h"+std::to_string(n+16)," Energy dep. in module     ", 200, 0., 20.);
	gv->IDMap[n+16] = id;
	id = analysis->CreateH1("h"+std::to_string(n+17)," Energy dep. in 1. module  ", 200, 0., 20.);
	gv->IDMap[n+17] = id;
	id = analysis->CreateH1("h"+std::to_string(n+18)," Energy dep. in 2. module  ", 200, 0., 20.);
    gv->IDMap[n+18] = id;	
	id = analysis->CreateH1("h"+std::to_string(n+19)," Energy dep. in 3. module  ", 200, 0., 20.);
    gv->IDMap[n+19] = id;	
	id = analysis->CreateH1("h"+std::to_string(n+20)," Energy dep. in 4. module  ", 200, 0., 20.);
	gv->IDMap[n+20] = id;
	//analysis->SetFirstHistoId(n+21);
	id = analysis->CreateH1("h"+std::to_string(n+21)," Energy dep. in crystal    ", 200, 0., 20.);
	gv->IDMap[n+21] = id;
//analysis->SetFirstHistoId(n+26);
	id = analysis->CreateH1("h"+std::to_string(n+26)," True conversion z ", 100, -20., 20.);   
gv->IDMap[n+26] = id;    
	id = analysis->CreateH1("h"+std::to_string(n+27)," Energy weighted z ", 100, -20., 20.);    
gv->IDMap[n+27] = id;	
	id = analysis->CreateH1("h"+std::to_string(n+28)," Distance: conv. and max-energy dep.    (xy)  ", 100, 0., 1.);
	gv->IDMap[n+28] = id;
	id = analysis->CreateH1("h"+std::to_string(n+29)," Distance: conv. and max-energy dep. (  xyz)  ", 100, 0., 1.);
	gv->IDMap[n+29] = id;
	id = analysis->CreateH1("h"+std::to_string(n+30)," Distance: PMT and max-energy dep. (xy)  ", 100, 0., 20.); 
	gv->IDMap[n+30] = id;
	id = analysis->CreateH1("h"+std::to_string(n+31)," Number of photons detected in PMT ", 200, 10., 10000.); 
	gv->IDMap[n+31] = id;
	id = analysis->CreateH1("h"+std::to_string(n + 32), " Photons in 1. PMT ", 2000, 0., 10000.); 
	gv->IDMap[n+32] = id;
	id = analysis->CreateH1("h"+std::to_string(n + 33), " Photons in 2. PMT ", 100, 0., 1000.);
	gv->IDMap[n+33] = id;
	id = analysis->CreateH1("h"+std::to_string(n + 34), " Photons in 3. PMT ", 100, 0., 1000.); 
	gv->IDMap[n+34] = id;
	id = analysis->CreateH1("h"+std::to_string(n + 35), " Photons in 4. PMT ", 100, 0., 1000.);  
	gv->IDMap[n+35] = id;
	id = analysis->CreateH1("h"+std::to_string(n + 36), " Number of PMTs hit above threshold", 20, 0., 20.);  
	gv->IDMap[n+36] = id;
	id = analysis->CreateH1("h"+std::to_string(n + 37), " Total Number of photons det. in PMTs > thrsld ", 2000, 0., 10000.);  
	gv->IDMap[n+37] = id;
       //analysis->SetFirstHistoId(n+40);
	id = analysis->CreateH1("h"+std::to_string(n + 40), " Reconstructed x-position ", 120, -30., 30.);  
	gv->IDMap[n+40] = id;
	id = analysis->CreateH1("h"+std::to_string(n + 41), " Reconstructed y-position ", 120, -30., 30.); 
	gv->IDMap[n+41] = id;
	id = analysis->CreateH1("h"+std::to_string(n + 42), " Reconstructed z-position ", 120, -30., 30.); 
	gv->IDMap[n+42] = id;
        //analysis->SetFirstHistoId(n+50);
	id = analysis->CreateH1("h"+std::to_string(n + 50), " Number of photons max. ", 100, 0., 15000.);  
	gv->IDMap[n+50] = id;
	id = analysis->CreateH1("h"+std::to_string(n + 51), " Number of photons generated ", 100, 0., 5000.); 
	gv->IDMap[n+51] = id;
	id = analysis->CreateH1("h"+std::to_string(n + 52), " Number of photons lost LABS ", 100, 0., 5000.); 
    gv->IDMap[n+52] = id;	
	id = analysis->CreateH1("h"+std::to_string(n + 53), " Number of photons lost REFL ", 100, 0., 5000.); 
	gv->IDMap[n+53] = id;
	id = analysis->CreateH1("h"+std::to_string(n + 54), " Number of photons lost ds < e ", 100, 0., 1000.);  
	gv->IDMap[n+54] = id;
	id = analysis->CreateH1("h"+std::to_string(n + 55), " Number of photons lost N > 1000 ", 100, 0., 1000.); 
	gv->IDMap[n+55] = id;
	id = analysis->CreateH1("h"+std::to_string(n + 56), " Number of photons unable to reflect ", 100, 0., 1000.); 
	gv->IDMap[n+56] = id;
	id = analysis->CreateH1("h"+std::to_string(n + 57), " Number of photons with error from GLISUR ", 100, 0., 1000.); 
	gv->IDMap[n+57] = id;
        //analysis->SetFirstHistoId(n+61);
	id = analysis->CreateH1("h"+std::to_string(n + 61), " Number of steps taken to PMT ", 100, 0., 200.); 
	gv->IDMap[n+61] = id;
	id = analysis->CreateH1("h"+std::to_string(n + 62), " Total track length to PMT ", 100, 0., 100.);  
	gv->IDMap[n+62] = id;
        //analysis->SetFirstHistoId(n+70);
	id = analysis->CreateH1("h"+std::to_string(n + 70), " Number of photon clusters ", 10, 0., 10.);
	gv->IDMap[n+70] = id;
	id = analysis->CreateH1("h"+std::to_string(n + 71), " Energy of Cluster #1 ", 200, 0., 20.);  
	gv->IDMap[n+71] = id;
	id = analysis->CreateH1("h"+std::to_string(n + 72), " Energy of Cluster #2 ", 200, 0., 20.); 
	gv->IDMap[n+72] = id;
	id = analysis->CreateH1("h"+std::to_string(n + 73), " Energy of Cluster #3 ", 200, 0., 20.);  
	gv->IDMap[n+73] = id;
	id = analysis->CreateH1("h"+std::to_string(n + 74), " Energy difference #1 ", 80, -2., 2.); 
	gv->IDMap[n+74] = id;
	id = analysis->CreateH1("h"+std::to_string(n + 75), " Energy difference #2 ", 80, -2., 2.); 
	gv->IDMap[n+75] = id;
	id = analysis->CreateH1("h"+std::to_string(n + 76), " Energy difference #3 ", 80, -2., 2.); 
	gv->IDMap[n+76] = id;
	id = analysis->CreateH1("h"+std::to_string(n + 77), " Dir. diff. [deg] #1 ", 90, 0., 90.);  
	gv->IDMap[n+77] = id;
	id = analysis->CreateH1("h"+std::to_string(n + 78), " Dir. diff. [deg] #2 ", 90, 0., 90.); 
	gv->IDMap[n+78] = id;
	id = analysis->CreateH1("h"+std::to_string(n + 79), " Dir. diff. [deg] #3 ", 90, 0., 90.); 
	gv->IDMap[n+79] = id;

	n = 100;
	//analysis->SetFirstHistoId(n+1);
	id = analysis->CreateH2("h"+std::to_string(n + 1), " Initial -  y - vs -  x - ", 100, -2.5, 2.5, 100, -2.5, 2.5);
	gv->IDMap[n+1] = id;
	id = analysis->CreateH2("h"+std::to_string(n + 2), " Initial - dy - vs - dx - ", 100, -100.0, 100.0, 100, -100.0, 100.0);
	gv->IDMap[n+2] = id;
	id = analysis->CreateH2("h"+std::to_string(n + 3), " IniFin -  y - vs -  x - ", 100, -2.0, 2.0, 100, -2.0, 2.0);
	gv->IDMap[n+3] = id;
	id = analysis->CreateH2("h"+std::to_string(n + 4), " IniFin - dy - vs - dx - ", 100, -100.0, 100.0, 100, -100.0, 100.0);
	gv->IDMap[n+4] = id;
	id = analysis->CreateH2("h"+std::to_string(n + 5), " Initial - dx - vs - x -", 100, -2.0, 2.0, 100, -100.0, 100.0);
	gv->IDMap[n+5] = id;
	id = analysis->CreateH2("h"+std::to_string(n + 6), " Initial - dy - vs - y -", 100, -2.0, 2.0, 100, -100.0, 100.0);
	gv->IDMap[n+6] = id;
	id = analysis->CreateH2("h"+std::to_string(n + 7), " IniFin - dx - vs - x -", 100, -2.0, 2.0, 100, -100.0, 100.0);
	gv->IDMap[n+7] = id;
	id = analysis->CreateH2("h"+std::to_string(n + 8), " IniFin - dy - vs - y -", 100, -2.0, 2.0, 100, -100.0, 100.0);
	gv->IDMap[n+8] = id;
	//analysis->SetFirstHistoId(n+11);
	id = analysis->CreateH2("h"+std::to_string(n + 11), " Final -  y - vs -  x - ", 16, -2.4, 2.4, 16, -2.4, 2.4);
	gv->IDMap[n+11] = id;
	id = analysis->CreateH2("h"+std::to_string(n + 12), " Final - dy - vs - dx - ", 100, -100.0, 100.0, 100, -50.0, 50.0);
	gv->IDMap[n+12] = id;
	id = analysis->CreateH2("h"+std::to_string(n + 13), " Final - dx - vs - x ", 100, -2.5, 2.5, 100, -100.0, 100.0);
	gv->IDMap[n+13] = id;
	id = analysis->CreateH2("h"+std::to_string(n + 14), " Final - dy - vs - y ", 100, -3.0, 3.0, 100, -50.0, 50.0);
	gv->IDMap[n+14] = id;
	id = analysis->CreateH2("h"+std::to_string(n + 15), " Final theta vs radius ", 100, 0.0, 2.0, 100, 0.0, 100.0);
	gv->IDMap[n+15] = id;
	
	n = 120;
	//analysis->SetFirstHistoId(n+1);
	id = analysis->CreateH2("h"+std::to_string(n + 1), " True conversion position ", 60, -15.0, 15.0, 80, -20.0, 20.0);
	gv->IDMap[n+1] = id;
	id = analysis->CreateH2("h"+std::to_string(n + 2), " Energy weighted position ", 60, -15.0, 15.0, 80, -20.0, 20.0);
	gv->IDMap[n+2] = id;
	id = analysis->CreateH2("h"+std::to_string(n + 3), " Reconstructed position   ", 60, -15.0, 15.0, 80, -20.0, 20.0);
	gv->IDMap[n+3] = id;
	id = analysis->CreateH2("h"+std::to_string(n + 4), " True conversion xy fngr-coordinates ", 60, -15.0, 15.0, 80, -20.0, 20.0);
	gv->IDMap[n+4] = id;
	id = analysis->CreateH2("h"+std::to_string(n + 5), " True conversion zx fngr-coordinates ", 60, -15.0, 15.0, 80, -20.0, 20.0);
	gv->IDMap[n+5] = id;
	id = analysis->CreateH2("h"+std::to_string(n + 6), " True conversion zy fngr-coordinates ", 60, -15.0, 15.0, 80, -20.0, 20.0);
	gv->IDMap[n+6] = id;
	//analysis->SetFirstHistoId(132);
	id = analysis->CreateH2("h"+std::to_string(n + 11), " Max loop =  1 ", 29, 1.0, 30.0, 100, 0.0, 10.0);
	gv->IDMap[n+11] = id;
	id = analysis->CreateH2("h"+std::to_string(n + 12), " Max loop =  2 ", 29, 1.0, 30.0, 100, 0.0, 10.0);
	gv->IDMap[n+12] = id;
	id = analysis->CreateH2("h"+std::to_string(n + 13), " Max loop =  3 ", 29, 1.0, 30.0, 100, 0.0, 10.0);
	gv->IDMap[n+13] = id;
	id = analysis->CreateH2("h"+std::to_string(n + 14), " Max loop =  4 ", 29, 1.0, 30.0, 100, 0.0, 10.0);
	gv->IDMap[n+14] = id;
	id = analysis->CreateH2("h"+std::to_string(n + 15), " Max loop =  5 ", 29, 1.0, 30.0, 100, 0.0, 10.0);
	gv->IDMap[n+15] = id;
	id = analysis->CreateH2("h"+std::to_string(n + 16), " Max loop =  6 ", 29, 1.0, 30.0, 100, 0.0, 10.0);
	gv->IDMap[n+16] = id;
	id = analysis->CreateH2("h"+std::to_string(n + 17), " Max loop =  7 ", 29, 1.0, 30.0, 100, 0.0, 10.0);
	gv->IDMap[n+17] = id;
	id = analysis->CreateH2("h"+std::to_string(n + 18), " Max loop =  8 ", 29, 1.0, 30.0, 100, 0.0, 10.0);
	gv->IDMap[n+18] = id;
	id = analysis->CreateH2("h"+std::to_string(n + 19), " Max loop =  9 ", 29, 1.0, 30.0, 100, 0.0, 10.0);
	gv->IDMap[n+19] = id;
	id = analysis->CreateH2("h"+std::to_string(n + 20), " Max loop = 10 ", 29, 1.0, 30.0, 100, 0.0, 10.0);
	gv->IDMap[n+20] = id;

	n = 200;
	//analysis->SetFirstHistoId(n+1);
	id = analysis->CreateH1("h"+std::to_string(n + 1), "Z-Stops in all col", 1000, -fDet->GetTLrms()/cm, fDet->GetTLrms()/cm); 
	gv->IDMap[n+1] = id;
	id = analysis->CreateH1("h"+std::to_string(n + 2), "R-Stops in targ entrance col", 50, 0.0, fDet->GetRrms()/cm / 2.0);
	gv->IDMap[n+2] = id;
	id = analysis->CreateH1("h"+std::to_string(n + 3), "R-Stops in targ exit     col", 50, 0.0, fDet->GetRrms()/cm / 2.0);
	gv->IDMap[n+3] = id;
	id = analysis->CreateH1("h"+std::to_string(n + 4), "TOF to TEND", 200, 0.0, 2.0);
	gv->IDMap[n+4] = id;
	id = analysis->CreateH1("h"+std::to_string(n + 5), "Reaction z-pos", 1000, -3 * materials->Gettargetl()/cm, 3 * materials->Gettargetl()/cm);    
    gv->IDMap[n+5] = id;
	id = analysis->CreateH1("h"+std::to_string(n + 6), "Beam Stops in target", 1000, -fDet->GetTLrms()/cm, fDet->GetTLrms())/cm;
	gv->IDMap[n+6] = id;
	//analysis->SetFirstHistoId(n+11);
	id = analysis->CreateH2("h"+std::to_string(n + 11), "stop/exit dist", 200, -fDet->GetTLrms()/cm, fDet->GetTLrms()/cm, 20, 0.0, fDet->GetRrms()/cm / 2.0);
	gv->IDMap[n+11] = id;
	id = analysis->CreateH2("h"+std::to_string(n + 12), "stop dist     ", 200, -fDet->GetTLrms()/cm, fDet->GetTLrms()/cm, 20, 0.0, fDet->GetRrms()/cm / 2.0);
	gv->IDMap[n+12] = id;
	id = analysis->CreateH2("h"+std::to_string(n + 13), "Reaction position R-pos vs z-pos ", 50, -3 * materials->Gettargetl()/cm, 3 * materials->Gettargetl()/cm, 50, 0.0, 1.0); 
	gv->IDMap[n+13] = id;
	id = analysis->CreateH2("h"+std::to_string(n + 14), "Exit spot", 50, -fDet->GetRrms()/cm / 2.0, fDet->GetRrms()/cm / 2.0, 50, -fDet->GetRrms()/cm / 2.0, fDet->GetRrms() / 2.0);
	gv->IDMap[n+14] = id;
	id = analysis->CreateH2("h"+std::to_string(n + 15), "cosx vs X", 100, -fDet->GetRrms()/cm / 2.0, fDet->GetRrms()/cm / 2.0, 100, -0.02, 0.02);
	gv->IDMap[n+15] = id;
	id = analysis->CreateH2("h"+std::to_string(n + 16), "cosy vs Y", 100, -fDet->GetRrms()/cm / 2.0, fDet->GetRrms()/cm / 2.0, 100, -0.02, 0.02);
	gv->IDMap[n+16] = id;
	id = analysis->CreateH2("h"+std::to_string(n + 17), "cosx vs Xtarg", 100, -1.0, 1.0, 100, -0.02, 0.02);
	gv->IDMap[n+17] = id;
	id = analysis->CreateH2("h"+std::to_string(n + 18), "cosy vs Ytarg", 100, -1.0, 1.0, 100, -0.02, 0.02);
	gv->IDMap[n+18] = id;
	
	//analysis->SetFirstHistoId(n+21);
	id = analysis->CreateH1("h"+std::to_string(n + 21), " Ini - x - Recoils ", 100, -2.0, 2.0);
	gv->IDMap[n+21] = id;
	id = analysis->CreateH1("h"+std::to_string(n + 22), " Ini - y - Recoils ", 100, -2.0, 2.0);
	gv->IDMap[n+22] = id;
	id = analysis->CreateH1("h"+std::to_string(n + 23), " Ini - dx - Recoils", 100, -100.0, 100.0);
	gv->IDMap[n+23] = id;
	id = analysis->CreateH1("h"+std::to_string(n + 24), " Ini - dy - Recoils", 100, -100.0, 100.0);
	gv->IDMap[n+24] = id;
	id = analysis->CreateH1("h"+std::to_string(n + 25), " Momentum spread (%)-Recoils", 100, -5.0, 5.0);
	gv->IDMap[n+25] = id;

    //MT adds histograms
 	n = 300;
	//analysis->SetFirstHistoId(n+0);
	id = analysis->CreateH1("h"+std::to_string(n + 0), " Recoil Stopping Length (cm) ", 2000, 0.0, 2000.0);
	gv->IDMap[n+0] = id;
	id = analysis->CreateH1("h"+std::to_string(n + 1), " Beam Particle Stopping Length (cm) ", 2000, 0.0, 2000.0);
	gv->IDMap[n+1] = id;

//    analysis->SetFirstHistoId(n+10);
//    id = analysis->CreateH2(std::to_string(n + 10), " DSSSD Energy (MeV) vs. x-strip position ",  16, -2.4, 2.4, 4000, 0.0, 20.0);
//IDMap[n+10] = id;
//    id = analysis->CreateH2(std::to_string(n + 11), " DSSSD x-strip position vs. Energy (MeV) ",  4000, 0.0, 20.0, 16, -2.4, 2.4);
//IDMap[n+11] = id;


//C.
//C.    Strip detector separate spectra

//    n = 300;
//   nstrip = 16;
//    analysis->SetFirstHistoId(n+1);
//    for(int i=1; i<=nstrip;i++) 
//       {
//        strip = "x-strip(" + num[i] + ")";
//        id = analysis->CreateH1("h"+std::to_string(n + i), strip, 16, 0., 16.0);
//        IDMap[n+i] = id;
//        strip = "y-strip(" + num[i] + ")";
//        id = analysis->CreateH1("h"+std::to_string(n + nstrip + i), strip, 16,0.,16.);
//        IDMap[n + nstrip + i] = id;
//        }

//C.    Strip detector hit patterns
//    n = 400;
//    nstrip = 16;
//    analysis->SetFirstHistoId(n+1);
//    id = analysis->CreateH1("h"+std::to_string(n + 1), "x-strip hit pattern", 16,0.,16.);
      //IDMap[n+1] = id;
//	  id = analysis->CreateH1("h"+std::to_string(n + 2), "y-strip hit pattern", 16,0.,16.);
      //IDMap[n+2] = id;

    //user defined angular distribution

    n = 250;
    //analysis->SetFirstHistoId(n);
    id = analysis->CreateH1("h"+std::to_string(n), " ang. dist.", 1000, -1.0, 1.0);
	gv->IDMap[n] = id;

    //C.    cross-section function and probability density
    n = 500;
    //analysis->SetFirstHistoId(n);
    //id = analysis->CreateH1("h"+std::to_string(n + 25), " capture cross-section", 1000, lm1, lm2); OJO 
	
	// CM energy distribution
	//analysis->SetFirstHistoId(502);
	//id = analysis->CreateH1("h502", "CM energy distribution", 1000, fPrim->Getbeamo() * (1.0 - 0.01), fPrim->Getbeamenerg() * (1.0 + 0.01)); //OJO
        id = analysis->CreateH1("h502", "CM energy distribution", 1000, 0., 1.0); 
		gv->IDMap[502] = id;

	// Beam caught after D1, scaler
	id = analysis->CreateH1("h503", "Caught beam", 1, 1.0, 2.0);
	gv->IDMap[503] = id;

	// Recoils which don't make it to ENDV
	id = analysis->CreateH2("h504", "Stopped recoil pos.", 100, -2000.0, 1000.0, 100, -1000.0, 100.0);
	gv->IDMap[504] = id;

        //analysis->SetFirstHistoId(510);
	// A quick way to keep track of how many recoils hit the end detector
	id = analysis->CreateH1("h510", "# Recoils to end detector", 1, 0.0, 1.0);
	gv->IDMap[510] = id;

	// Number of BGOs trigger in coincedence with recoil detection
	id = analysis->CreateH1("h511", "# BGOs triggered, in coin.", 10, 0.0, 10.0);
	gv->IDMap[511] = id;

	// Sum of all energy deposited in BGOs in coincedence with recoil detection
	id = analysis->CreateH1("h512", "Energy in BGOs, in coin.", 250, 0.0, 25.0);
	gv->IDMap[512] = id;

	// Sum of all energy deposited in BGOs vrs. Number of BGOs triggered,
	// in coincedence with recoil detection
	id = analysis->CreateH2("h513", "Energy in BGOs vs # BGOs, in coin.", 10, 0.0, 10.0, 250, 0.0, 25.0);
	gv->IDMap[513] = id;

	// Energy deposited in BGO with second greatest energy vrs energy deposited
	// in BGO with greatest energy, in coincedence with recoil detection
	id = analysis->CreateH2("h514", "Energy of 2nd Energetic BGO vs 1st Energergetic BGO, in coin.", 250, 0.0, 25.0, 250, 0.0, 25.0);
	gv->IDMap[514] = id;

	// # of BGO with second greatest energy vrs # of BGO with greatest
	// energy, in coincedence with recoil detection
	id = analysis->CreateH2("h515", "# 2nd Energetic BGO vrs # 1st Energetic BGO, in coin.", 30, 1.0, 30.0, 30, 1.0, 30.0);
	gv->IDMap[515] = id;

	id = analysis->CreateH1("h516", "Detected Energy (with P.H.D included)", 4000, 0.0, 20.0);
	gv->IDMap[516] = id;

	id = analysis->CreateH1("h517", "Local t.o.f (MCP-DSSSD)", 1000, 0.0, 200.0);
	gv->IDMap[517] = id;

	id = analysis->CreateH1("h518", "Local t.o.f (MCP1-MCP0)", 1000, 0.0, 200.0);
	gv->IDMap[518] = id;

    //analysis->SetFirstHistoId(520);
	id = analysis->CreateH1("h520", "Reaction x-position", 1000, -5.0, 5.0);
	gv->IDMap[520] = id;

	id = analysis->CreateH1("h521", "reaction y-position", 1000, -5.0, 5.0);
	gv->IDMap[521] = id;

	id = analysis->CreateH2("h522", "Q1 start x-y hit pattern", 16, -2.45, 2.45, 16, -2.45, 2.45);
	gv->IDMap[522] = id;

//C     Alpha source histogram(s)
	id = analysis->CreateH1("h523", "Alpha source energy", 512, 0.0, 4096.0);
	gv->IDMap[523] = id;
//	id = analysis->CreateH2("h524", "Alpha position", 100, -0.6, 0.6, 100, -0.6, 0.6);

//C.
//C.--> New ntuples 07.07.03
//C.
//C  'History' ntuple
    id = analysis->CreateNtuple("h1000","HISTORY");
    analysis->CreateNtupleDColumn("E_int");
    analysis->CreateNtupleDColumn("E_rec");
    for (int i = 0; i < 15; ++i) 
        {
         analysis->CreateNtupleDColumn("E_g(" + std::to_string(i+1) + ")");
         analysis->CreateNtupleDColumn("E_gp(" + std::to_string(i+1) + ")");
         analysis->CreateNtupleDColumn("cost_g(" + std::to_string(i+1) + ")");
         analysis->CreateNtupleDColumn("phi_g(" + std::to_string(i+1) + ")");
         analysis->CreateNtupleDColumn("cost_gp(" + std::to_string(i+1) + ")");
         }
    analysis->CreateNtupleDColumn("cost_r");
    analysis->CreateNtupleDColumn("cosp_r");
    analysis->CreateNtupleIColumn("Nodec");
    analysis->CreateNtupleIColumn("react");
    analysis->CreateNtupleIColumn("recdet");
    analysis->CreateNtupleDColumn("x_r");
    analysis->CreateNtupleDColumn("y_r");
    analysis->CreateNtupleDColumn("z_r");
    analysis->CreateNtupleDColumn("thet_r");
    analysis->CreateNtupleDColumn("xstop");
    analysis->CreateNtupleDColumn("ystop");
    analysis->CreateNtupleDColumn("zstop");
    analysis->CreateNtupleDColumn("xint");
    analysis->CreateNtupleDColumn("yint");
    analysis->CreateNtupleDColumn("zint");
    analysis->CreateNtupleDColumn("x");
    analysis->CreateNtupleDColumn("y");
    analysis->CreateNtupleDColumn("xp");
    analysis->CreateNtupleDColumn("yp");
    for (int i = 0; i < 10; ++i) 
        {
         analysis->CreateNtupleDColumn("xtest(" + std::to_string(i+1) + ")");
         analysis->CreateNtupleDColumn("ytest(" + std::to_string(i+1) + ")");
         analysis->CreateNtupleDColumn("etest(" + std::to_string(i+1) + ")");
         }    
    analysis->CreateNtupleIColumn("dsssdpos");
    analysis->CreateNtupleDColumn("beamtof");
    analysis->FinishNtuple();
  
//C.
//C.--> New nutple 20040407
//C.
//C 'Gammahit' ntuple
    //C 'Gammahit' ntuple
    id = analysis->CreateNtuple("h1001","GAMMAHITS");
    analysis->CreateNtupleIColumn("recoil_hit_ENDV");
    analysis->CreateNtupleIColumn("num_bgos_hit");
    analysis->CreateNtupleDColumn("e_bgos_total");
    analysis->CreateNtupleIColumn("num_bgo_first");
    analysis->CreateNtupleDColumn("e_bgo_first");
    analysis->CreateNtupleIColumn("num_bgo_second");
    analysis->CreateNtupleDColumn("e_bgo_second");
    analysis->CreateNtupleIColumn("pair_productions");
    analysis->CreateNtupleIColumn("num_bgos_hit_ab");
    analysis->CreateNtupleIColumn("num_bgo_first_ab");
    analysis->CreateNtupleDColumn("e_bgo_first_ab");
    analysis->CreateNtupleIColumn("num_bgo_second_ab");
    analysis->CreateNtupleDColumn("e_bgo_second_ab");
    analysis->CreateNtupleDColumn("gammatof");
    analysis->CreateNtupleDColumn("e0_conv");
    analysis->FinishNtuple();

    //C -->   Define other ntuples
    
    if(DRAGONRunAction_histo->iswit[8] == 2)
      {
       //C Gamma -HI coincidence ntuple
       id  = analysis->CreateNtuple("h100","Gamma-HI");
       const G4int nvar10 = 37;    //From higamcoinc.inc
       for (int i = 0; i<nvar10; i++)    
            analysis->CreateNtupleDColumn(chtags[i]);
       analysis->FinishNtuple();
//C.
//C.-->   Setup parameters - filled once at the end of UGINIT
//C.
       id  = analysis->CreateNtuple("h998","GBOX Geant Setup Ntuple");
       const G4int max_gbox = 29;    //From u_geom.inc
       for (int i = 0; i<max_gbox; i++) 
           {   
            analysis->CreateNtupleDColumn("x1_fngr(" + std::to_string(i+1) + ")");
            analysis->CreateNtupleDColumn("x2_fngr(" + std::to_string(i+1) + ")");
            analysis->CreateNtupleDColumn("y1_fngr(" + std::to_string(i+1) + ")");
            analysis->CreateNtupleDColumn("y2_fngr(" + std::to_string(i+1) + ")");   
            analysis->CreateNtupleDColumn("z1_fngr(" + std::to_string(i+1) + ")");
            analysis->CreateNtupleDColumn("z2_fngr(" + std::to_string(i+1) + ")");                      
            }
       analysis->FinishNtuple();

//C.
//C.-->   Event variables - filled every event at the end of GUDIGI
//C.
    id = analysis->CreateNtuple("h999","GBOX Geant Event Ntuple");
    const G4int melem_gbox = 20;    
    for (int i = 0; i < melem_gbox; ++i) 
         analysis->CreateNtupleIColumn("melem_gbox(" + std::to_string(i+1) + ")");
    const G4int jelem_gbox = 29;    
    for (int i = 0; i < jelem_gbox; ++i) 
         analysis->CreateNtupleIColumn("jelem_gbox(1," + std::to_string(i+1) + ")");

    for (int i = 0; i < melem_gbox; ++i) 
         analysis->CreateNtupleDColumn("energy_gbox(" + std::to_string(i+1) + ")");       
    analysis->FinishNtuple();
    }
 }
 else
 {
	 ;
	 
	 
 }
}

      
}




void RemoveNoneReactions(){
	
	TFile *fin = new TFile("22neag_4_9T.root");
	TTree *tin;
	fin->GetObject("h1000",tin);

	float xint, yint, zint;
	float Eint;
	float xangle, yangle;
	int react;

	tin->SetBranchAddress("xint", &xint);
	tin->SetBranchAddress("yint", &yint);
	tin->SetBranchAddress("zint", &zint);
	tin->SetBranchAddress("E_int", &Eint);
	tin->SetBranchAddress("react", &react);
	tin->SetBranchAddress("xp", &xangle);
	tin->SetBranchAddress("yp", &yangle);

	TFile *fout = new TFile("22neag_4_9T_filtered.root","RECREATE");
	TTree *tout = new TTree("h1000","h1000");

	tout->Branch("xint", &xint, "xint/F");
	tout->Branch("yint", &yint, "yint/F");
	tout->Branch("zint", &zint, "zint/F");
	tout->Branch("Eint", &Eint, "Eint/F");
	tout->Branch("xangle", &xangle, "xangle/F");
	tout->Branch("yangle", &yangle, "yangle/F");

	//cout << "i, x, y, z" << endl;

	for (int i=0; i<tin->GetEntries(); i++){
		tin->GetEntry(i);
		if (react != 0){
			//cout << i << ", " << xint << ", " << yint << ", " << zint << endl;
			tout->Fill();
		}
	}

	tout->Write();
	fout->Close();
		
}




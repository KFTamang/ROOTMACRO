
Double_t getFWHM(TH1* h, Double_t lowerbound,Double_t upperbound){

	if(upperbound<=lowerbound){
		cerr << "warning! upperbound sholbe greater than lowerbound " << endl;
		return -1;
	}
		
	Int_t upperbin = h->GetXaxis()->FindBin(upperbound);
	Int_t lowerbin = h->GetXaxis()->FindBin(lowerbound);
	Int_t bindiff = upperbin - lowerbin;

//	h->GetXaxis()->SetRangeUser(lowerbin,upperbin);
	Int_t max =0;
	Int_t value =0;
	Int_t peakBin=lowerbin;
	for(Int_t i=lowerbin;i<upperbin;i++){
		value = h->GetBinContent(i);
		if(max<value){
			max = value;
			peakBin = i;
		}
	}
	// Int_t peakBin = h->GetMaximumBin();
	// Int_t max = h->GetMaximum(peakBin);
	Int_t binup  = peakBin;
	Int_t bindown = peakBin;
//	cerr << max << " " << peakBin << endl;

	while(h->GetBinContent(binup)>max/2.){
		binup++;
	}
	while(h->GetBinContent(bindown)>max/2.){
		bindown--;
	}
	Double_t up = h->GetBinCenter(binup);
	Double_t down = h->GetBinCenter(bindown);

	Double_t fwhm = up - down;


	Int_t height =  h->GetBinContent(peakBin); 
	Double_t peakPos =  h->GetBinCenter(peakBin); 
	Double_t area =  height*(binup-bindown); 
	cerr << "Peak position: " << peakPos << endl;
	cerr << "Peak height: " << height << endl; 
	cerr << "Peak area(rough): " << area << endl;
//	cerr << "Peak area: " << h->Integral(lowerbin,upperbin) << endl;
	cerr << "Lower value:" << down << ", upper value:" << up << endl;

//	cerr << down << " " << up << " " << fwhm << endl;

	return fwhm;
}


Double_t getFWHM(TH1* h, Double_t lowerbound,Double_t upperbound,TF1* bkg){

	Int_t upperbin = h->GetXaxis()->FindBin(upperbound);
	Int_t lowerbin = h->GetXaxis()->FindBin(lowerbound);
	Int_t bindiff = upperbin - lowerbin;
	Double_t val;

	TH1* htemp = (TH1*)h->Clone("htemp");
	TH1* hdiff = (TH1*)h->Clone("hdiff");

	htemp->Reset();
	for(Int_t i=lowerbin;i<upperbin;++i){
		val = h->GetBinCenter(i);
		htemp->Fill(val,bkg->Eval(val));
		cerr << i << " " <<  val << " " << h->GetBinContent(i) - bkg->Eval(val) << " " <<  bkg->Eval(val) << endl; 
	}
	hdiff->Add(htemp,-1);

	return	getFWHM(hdiff,lowerbound,upperbound);

}






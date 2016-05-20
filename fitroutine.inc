


/*
Todo:
	D FF ableiten
	lclass fuer FF
	
*/










TVectorD calcVariatedFF(double *par){

	// Vector containing world averages of FF parameters
	static TVectorD FFval(3);
	static TMatrixD FFcov(3,3);
	// eigenvalues: elements of the diagonal form of covariance matrix
	static TVectorD lambda(3);
	// matrix with eigenvector: transformation matrix between original covariance
	// and diagonal form: C' = P^-1 C P
	static TMatrixD P(3,3);
	static TMatrixD Pinvert(3,3);


	static bool firstcall = true;
	if(firstcall){
		FFval(0) = 1.207;
		FFval(1) = 1.401;
		FFval(2) = 0.854;


		FFcov(0,0) = 7.84e-4;    FFcov(0,1) = 5.23e-4;     FFcov(0,2) = -4.52e-4;
		FFcov(1,0) = 5.23e-4;    FFcov(1,1) = 1.089e-3;    FFcov(1,2) = -5.00e-4;
		FFcov(2,0) = -4.52e-4;   FFcov(2,1) = -5.00e-4;    FFcov(2,2) = 4.00e-4;
		
		P = FFcov.EigenVectors(lambda);
		Pinvert = P; Pinvert.Invert();
		
		// Transform FF parameters in new basis
		FFval *= Pinvert;
	  firstcall = false;
	}
	
	TVectorD FFval_new(lambda);
	for(int i = 0; i <3; i++){
		FFval_new(i) = FFval(i) + par[i]*sqrt(lambda(i));

	}
	
	// Rotate back into original basis to apply event weights
	FFval_new *= P;

	return FFval_new;

}


void getNuiModMC(const double *NuiPar, TH1D** hMC_Comp){

	//decorrelate D* Nuis
	double DsPar[] = {NuiPar[DsWeight->getErrIdx(0)],NuiPar[DsWeight->getErrIdx(1)],NuiPar[DsWeight->getErrIdx(2)]};
	TVectorD decorDsPar = calcVariatedFF(DsPar);
	
	static double DssParMean[] = {-1.5,0.,0.5};
	static double DssParSigma[] = {0.5,0.5,0.5};
	static double NuiParColl[4][3];
	DWeight->SetPar(NuiPar[DWeight->getErrIdx(0)]);
	
	for(uint i = 0; i<3; i++)
	{
		DsWeight->SetPar(i,decorDsPar(i));
		//link narrow and broad  
		for(uint c= 0; c<4; c++){
			NuiParColl[c][i] = NuiPar[DssWeight[2*(c/2)]->getErrIdx(i)]; //link narrow and broad: 0 = 1, 2 =3
			DssWeight[c]->SetPar(i,DssParMean[i] + NuiParColl[c][i]*DssParSigma[i]);
			
		}
		
	}
	

	
	for(uint i = 0; i<nComponents; i++){
		hMC_Comp[i]->Reset();
	}
	
	static double BFXlnu[] = {10.8,10.1};
	static double BFXtaunu[] = {2.48,2.33};
	double dBFXlnu[2], dBFXtaunu[2];
	double b_tmp;
	for(int b =0; b<2; b++){
		b_tmp = 0;
		for(uint i = t_error::block_start[t_error::_xlnu]; i<=t_error::block_end[t_error::_xlnu]; i++)
			b_tmp += ((t_bf*)vErr[i])->BF[b]*vErr[i]->getChange(0)*NuiPar[i];
		
		dBFXlnu[b] = BFXlnu[b]/(BFXlnu[b] + b_tmp);
		
		b_tmp = 0;
		for(uint i = t_error::block_start[t_error::_xtaunu]; i<=t_error::block_end[t_error::_xtaunu]; i++)
			b_tmp += ((t_bf*)vErr[i])->BF[b]*vErr[i]->getChange(0)*NuiPar[i];
			
		dBFXtaunu[b] = BFXtaunu[b]/(BFXtaunu[b] + b_tmp);
	}
	
	//constain D BF to 100;
	/*
	double dBFDXlnu[2];
	
	int idxD = t_error::block_start[t_error::_Dxlnu];
	for(int b =0; b<2; b++){
		double BF = ((t_bf*)vErr[idxD + b])->BF[0];
		dBFDXlnu[b] = ( 100. - BF*(1.+vErr[idxD + b]->getChange(0)*NuiPar[idxD + b]) )/ (100. - BF);
	}
	*/
	
	static uint nEvents = evtlist.size();
	t_error::t_errType type;
	double c,change;
	uint nClasses, classIdx, j, iB0;
	for(uint i = 0; i<nEvents; ++i)
	{
	
		change =1.;
		nClasses = evtlist[i].classes.size();
		for(j = 0; j<nClasses; ++j)
		{
			classIdx = evtlist[i].classes[j];
			t_error &Err = *vErr[classIdx];
			c = Err.getChange(evtlist[i].bin);
			type = Err.Type;
			if(type != t_error::_FF){ 
				c = (1. + c*NuiPar[classIdx]);
				if(type == t_error::_bf){
					iB0 = evtlist[i].iB0;
					if(Err.block == t_error::_xlnu)
						c*=dBFXlnu[iB0];
					else if(Err.block == t_error::_xtaunu)
						c*=dBFXtaunu[iB0];
					//else if(Err.block == t_error::_Dxlnu && Err.subType == t_error::_dummy)
					//	c*=dBFDXlnu[classIdx-t_error::block_start[t_error::_Dxlnu]-2];
					
				}
			}
			change*= c;	
		}
		change*=1.+evtlist[i].eff_error*NuiPar[evtlist[i].eff_Idx];
		
		
		hMC_Comp[evtlist[i].component]->AddBinContent(evtlist[i].bin+1,evtlist[i].weight*change);
		
	
	}
	
}

double Chi2fcn(const double *params){
	const double *NuiPar = params + nComponents;
	
	getNuiModMC(NuiPar, hMC_Comp);
	
	hMC->Reset("ICESM");
	
	hMC->Add(hMC_Comp[0],params[0]);	
	hMC->Add(hMC_Comp[1]);	hMC->Scale(params[1]); //fit R(X)		
	for(ushort j =2; j<nComponents; j++)
		hMC->Add(hMC_Comp[j],params[j]);	
		
		
	double chi2 = 0;
	double nD;
	double nMC;
	for(uint i = 0; i< nBins; i++)
	{
		nD = hData->GetBinContent(i+1);
		nMC = hMC->GetBinContent(i+1);
		if( nD <= 0 || nMC <= 0) continue;
		double Chi_i = (nD-nMC);
		chi2 += Chi_i*Chi_i/nD; //diag covar
		
	/*	//for non diag covar
		for(int j = 0; j<nBins; j++)
		{
			if (hD[l]->GetBinContent(j+1) <= 0 || hMC[l]->GetBinContent(j+1) <= 0) continue;
			double Chi_j = (hD[l]->GetBinContent(j+1)-hMC[l]->GetBinContent(j+1));
			chi2 += Chi_i*Chi_j*(*MInvCovar[l])[i][j];
		}
	*/
	}
	
	for(uint i = 0; i<nErrors ;i++)
		chi2+=sqr(NuiPar[i]); //nui coll ausschliessen
	
	return chi2;
}



double RunChi2FitNui(double** xs, double*** errors  , int scnVar = -1, double scnVal = 0){
	ROOT::Math::Minimizer* min =    ROOT::Math::Factory::CreateMinimizer("Minuit2");
 
 
 	//get signal efficiencies
 	
 
        // set tolerance , etc...
        min->SetMaxFunctionCalls(30000); // for Minuit/Minuit2 
        min->SetMaxIterations(10000);  // for GSL 
        min->SetTolerance(1.);
        
 
        min->SetPrintLevel(4);
        // create funciton wrapper for minmizer// a IMultiGenFunction type 
        
        
        cout<<"Fit with n Pars: "<<nComponents+nErrors<<endl;
        ROOT::Math::Functor f(&Chi2fcn,nComponents+nErrors); 
        min->SetFunction(f);
 
    	// Set the free variables to be minimized!
    	for(uint i = 0; i<nComponents; i++){
    		min->SetVariable(i,vCont[i].Name, (gi_Error==0?vCont[i].startVal:xs[0][i]) ,0.01);
    		if(vCont[i].fl_const)
    			min->SetFixedVariable(i,vCont[i].Name,vCont[i].startVal);
    	}
	for(int i = 0; i<nErrors; i++){
		min->SetVariable(nComponents+i,vErr[i]->Name, (gi_Error==0?0.01:xs[0][nComponents + i]) , 0.05);
		//if(vErr[i]->Name == "D2Slnu" || vErr[i]->Name == "D*2Slnu")
			//min->SetLowerLimitedVariable(nComponents+i,vErr[i]->Name,0,0.1,-1);
		if(vErr[i]->Name.find("D2_") == 0 || vErr[i]->Name.find("D0*_") == 0 ||vErr[i]->subType == t_error::_dummy)
			min->SetFixedVariable(nComponents+i,vErr[i]->Name,0);
		if(gi_Error == i && __calcerrors==1)
			min->SetFixedVariable(nComponents+i,vErr[i]->Name,xs[0][nComponents+i]);
		if(i == 0)
			min->SetFixedVariable(nComponents+i,vErr[i]->Name,0);
	
	}
/*	for(int i = 4+nErrors -3; i<4+nErrors; i++)
		min->SetLimitedVariable(i,vErr[i-4].Name,0,0.1,-1,1);
*/	
	if(scnVar >=0)
		min->SetFixedVariable(scnVar,min->VariableName(scnVar),scnVal);
	
	
	
	
	int flag = min->Minimize();
	const double *x = min->X();
    	const double *xerr = min->Errors();
    	for(int i =0; i<nComponents+nErrors;i++){
    		xs[gi_Error][i] = x[i];
    		errors[0][gi_Error][i] = errors[1][gi_Error][i] = xerr[i];
    	}
    	
    	//min->Hesse();
    	if(0&&scnVar <0){
	    	for(int i=0; i<1; i++)
	    		min->GetMinosError(i,errors[0][gi_Error][i],errors[1][gi_Error][i]);
	 }
    	// xs contains best fit parameters for original norm values
    	
    	return 1;
    	
    	
    		
	 if(scnVar <0){
	   		
	    	cout<<"Correlation matrix"<<endl;
	    	TCanvas *cv = new TCanvas("corr","corr",1000,900);
	    	TH2F *hCorr = new TH2F("corr","corr",nComponents+nErrors,-0.5,nComponents+nErrors-0.5,nComponents+nErrors,-0.5,nComponents+nErrors-0.5);
	    	for(int i = 0; i< nComponents+nErrors; i++)
	    	{
			for(int j = 0; j< nComponents+nErrors; j++)
			{
				cout<<min->Correlation(i,j)<<"\t";
				hCorr->SetBinContent(i+1,j+1,min->Correlation(i,j));
			}
			cout<<endl;
	    	
	    	}
	    	hCorr->Draw("colz");
	    	cv->Print("correlation.pdf");
	    	delete hCorr;
	    	delete cv;
	    	
	    	cout<<"Results: "<<endl;
	    	for(int i =0; i<nComponents+nErrors;i++)
	    	{
	    		cout<<min->VariableName(i)<<" "<<xs[gi_Error][i]<<" +- "<<errors[0][gi_Error][i]<<", "<<errors[1][gi_Error][i]<<endl;
	    	}
	    	
	  } else {
	  
	  	static FILE *fp = fopen("scanChi2.dat","w");
	  	
	  	fprintf(fp,"%lf %lf %i\n", scnVal,Chi2fcn(xs[gi_Error]) ,min->Status());
	  	fflush(fp);
	  }
    delete min;

return 1;
}

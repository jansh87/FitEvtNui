
TH1D** InBinsOf(TH2F* in, int &nHists,char ax = 'x') {
	TH1D **out;
	TAxis *axisA;
	TAxis *axisB;
	if(ax == 'x')
	{
		axisA = in->GetXaxis();
		axisB = in->GetYaxis();
	}
	else if(ax == 'y')
	{
		axisA = in ->GetYaxis();
		axisB = in ->GetXaxis();
	}
	else {cout<<"Error: unknown axis"<<endl; return NULL;}

	char sName[200];

	nHists = axisB->GetNbins();
	out  = new TH1D*[nHists];
	const TArrayD *bins = axisA->GetXbins();
	for(int i =0; i<nHists; i++)
	{
		sprintf(sName,"%s_%c_%i",in->GetName(), ax, i);
		out[i] = new TH1D(sName,in->GetTitle(),axisA->GetNbins(),bins->fArray);
		out[i]->GetXaxis()->ImportAttributes(axisA);
		out[i]->GetXaxis()->SetTitle(axisA->GetTitle());
		THashList* labels=axisA->GetLabels();
		if (labels) {
			TIter iL(labels);
		        TObjString* lb;
			Int_t l = 1;
			while ((lb=(TObjString*)iL())) {
				out[i]->GetXaxis()->SetBinLabel(l,lb->String().Data());
				l++;
			}
		}

		out[i]->SetLineColor(in->GetLineColor());
		out[i]->SetFillColor(in->GetFillColor());
		out[i]->SetMarkerColor(in->GetMarkerColor());
		out[i]->SetMarkerStyle(in->GetMarkerStyle());
		//cout<<out[i]->Integral()<<endl;

		for(int j =1 ; j<=axisA->GetNbins();j++)
		{
			if(ax == 'x'){
				out[i]->SetBinContent( j, in->GetBinContent(j,i+1) );
				out[i]->SetBinError( j, in->GetBinError(j,i+1) );
			}else {
				out[i]->SetBinContent( j, in->GetBinContent(i+1,j) );
				out[i]->SetBinError( j, in->GetBinError(i+1,j) );
			}
		}

	}

	return out;
}



void PlotSlices(const char *sVarName, int axis, TH2F*** hMC, TH2F** hData, string name = ""){

			char sName[50];
			char sFileName[50];
			THStack ***stack;
			TH1D ****hTempDecomp;
			TH1D ***ToyDataDecomp;
			stack = new THStack**[2];
			hTempDecomp = new TH1D***[2];
			ToyDataDecomp = new TH1D**[2];
			
			TCanvas **cv1l = new TCanvas*[2];
			int nHists;
			
			for(int l = 0; l<2; l++)
			{
			
				cv1l[l] = new TCanvas("cv","cv",1000,1000);
				cv1l[l]->Divide(4,4);
				hTempDecomp[l] = new TH1D**[nComponents];
				//ToyDataDecomp[l] = new TH1D*[vVar[2*l+1].nBins];

				for(int k = nComponents-1; k>=0; k--)
				{
					hTempDecomp[l][k] = InBinsOf(hMC[l][k],nHists,axis>0?'x':'y');

				}

				ToyDataDecomp[l] = InBinsOf(hData[l],nHists,axis>0?'x':'y');
				stack[l] = new THStack*[nHists];
				cout<<"nh "<<nHists<<endl;
				for(int m = 0; m<nHists; m++)
				{
					sprintf(sName,"decompstack_%i_%i",l,m);
					stack[l][m] = new THStack(sName,sName);



					for(int k = nComponents-1; k>=0; k--)
					{
					//	if(k > 1) continue;
					//	cout<<l<<" "<<k<<" "<<m<<endl;
						stack[l][m]->Add(hTempDecomp[l][k][m]);
					}
					
					TPad *pad = (TPad*)cv1l[l]->cd(m+1);
					stack[l][m]->Draw("");
					stack[l][m]->GetYaxis()->SetTitle("Events");
					stack[l][m]->GetYaxis()->SetTitleOffset(0.6);
					stack[l][m]->GetXaxis()->SetTitle(!axis?hMC[l][0]-> GetXaxis() -> GetTitle():hMC[l][0]-> GetYaxis() -> GetTitle());
					stack[l][m]->Draw("hist");
					ToyDataDecomp[l][m]->Draw("e1 same");
					if(ToyDataDecomp[l][m]->GetMaximum()>stack[l][m]->GetMaximum())
						stack[l][m]->SetMaximum(1.05*(ToyDataDecomp[l][m]->GetMaximum()+sqrt(ToyDataDecomp[l][m]->GetMaximum())));
					DrawResiduals(stack[l][m],ToyDataDecomp[l][m],pad);
				}
				
			cv1l[l]->cd(0);
			
			

			sprintf(sFileName,"Dists_%s_%i_%s.pdf",sVarName,l, name.c_str());
			cv1l[l]->Print(sFileName);
			
			cv1l[l]->Delete();
			
			}


			cout<<"done"<<endl;

			for(int l = 0; l<2; l++)
			{
				for(int m = 0; m<nHists; m++)
				{
					for(int k = nComponents-1; k>=0; k--)
					{
						hTempDecomp[l][k][m]->Delete();
					}
					stack[l][m]->Delete();
					ToyDataDecomp[l][m]->Delete();
				}
				delete[] stack[l];
				delete[] ToyDataDecomp[l];
				for(int k = 3; k>=0; k--)
					delete[] hTempDecomp[l][k];
				delete[] hTempDecomp[l];

			}
				delete[] stack;
				delete[] ToyDataDecomp;
				delete[] hTempDecomp;
				cout<<"done"<<endl;
}		



void PlotProj(TH2F*** hMC, TH2F** hData, string name = ""){


	TCanvas *cv1 = new TCanvas("cv","cv",1000,1000);
	 //plot projections
	char sName[200];
	THStack*** stack = new THStack**[2];
	TH1D**** hProj = new TH1D***[2];
	TH1D*** hProjData = new TH1D**[2];
	cv1->Divide(2,2);
	for(int l = 0; l<2; l++) //e/mu
	{
		stack[l] = new THStack*[2];
		hProj[l] = new TH1D**[2];
		hProjData[l] = new TH1D*[2];
		for(int i = 0; i<2; i++) //p & mm2
		{
			sprintf(sName,"projStack_%i_%i",l,i);
			stack[l][i] = new THStack(sName,sName);

			hProjData[l][i] = !i?hData[l]->ProjectionX("_px",1,hData[l]->GetNbinsY()):hData[l]->ProjectionY("_py",1,hData[l]->GetNbinsX());
			hProj[l][i] = new TH1D*[nComponents];
			for(int m = nComponents-1; m>=0; m--)
			{
				//cout<<l<<" "<<i<<" "<<m<<endl;
				hProj[l][i][m] = !i?hMC[l][m]->ProjectionX("_px",1,hData[l]->GetNbinsY()):hMC[l][m]->ProjectionY("_py",1,hData[l]->GetNbinsX());
				hProj[l][i][m]->SetFillColor(hMC[l][m]->GetFillColor());
				hProj[l][i][m]->SetLineColor(hMC[l][m]->GetFillColor());
				
				stack[l][i] -> Add(hProj[l][i][m]);

			}

			TPad *pad = (TPad*)cv1->cd(l+2*i+1);
			stack[l][i]->Draw();
			stack[l][i] -> GetYaxis() -> SetTitle("Events");
			stack[l][i] -> GetXaxis() -> SetTitle(!i?hMC[l][0]-> GetXaxis() -> GetTitle():hMC[l][0]-> GetYaxis() -> GetTitle());
			//hProjData[l][i]->Draw("same E1");
			DrawResiduals(stack[l][i],hProjData[l][i],pad);
			

		}
	}

	cv1->cd(0);
	sprintf(sName,"ProjStack_%s.pdf",name.c_str());
	cv1->Print(sName);
	cv1->Clear();
	cv1->Delete();
	//end plot
}		
				

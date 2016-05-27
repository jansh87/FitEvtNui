
template<typename Th>
void DrawResiduals(THStack *StackMC, Th* data, TCanvas *cv1){
	TPad *pad = (TPad*)cv1 -> cd(0);
	DrawResiduals(StackMC, data,pad);
}
template<typename Th>
void DrawResiduals(THStack *StackMC, Th* data, TPad *pad){

	pad->cd();
	TPad *pad1;
	pad1 = new TPad("p1","p1",0.,0.3,1,1);
	pad1->SetTopMargin(0.1);
	pad1->SetBottomMargin(0);
	pad1->Draw();
	pad1->cd();

	StackMC->Draw("HIST");
	data->Draw("E1 X0 same");
	data->SetMarkerStyle(20);
	data->SetMarkerSize(0.5);
	data->SetMarkerColor(kBlack);
	if(data->GetMaximum()+sqrt(data->GetMaximum()) > StackMC->GetMaximum())
	{
		StackMC->SetMaximum(data->GetMaximum()+sqrt(data->GetMaximum()));
	}
	StackMC->SetMinimum(0.1);

	//cv1->Update();

	pad->cd();
	TPad *pad2;
	pad2 = new TPad("p2","p2",0.,0.,1,0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.45);
	pad2->Draw();
	pad2->cd();

	Th *sumMC;
	TList* stackList = StackMC->GetHists();
	sumMC = (Th*)stackList->At(0);
	for(int i = 1; i< stackList->GetSize(); i++)
		sumMC->Add((Th*)stackList->At(i));

	sumMC->SetStats(kFALSE);
	float val,err,dat;
	for(int i =1; i<=StackMC->GetXaxis()->GetNbins();i++)
	{
		val = sumMC->GetBinContent(i);
		dat = data->GetBinContent(i);
		if(dat != 0)
		{
		    sumMC->SetBinContent(i, (val-dat)/sqrt(dat) );
		    sumMC->SetBinError(i, 1.);
		    //sumMC->SetBinContent(i, dat/val );
		    //sumMC->SetBinError(i, sqrt(dat)/val);
		}
		else
		{
		    sumMC->SetBinContent(i, -50 );
		    sumMC->SetBinError(i, 1 );
		}

		
	}
	//gPad->SetLogy();
	sumMC->SetMarkerStyle(20);
	sumMC->SetMarkerSize(0.5);
	sumMC->SetLineColor(kBlack);
	//sumMC->SetAxisRange(0.667,1.5,"Y");
	sumMC->SetAxisRange(-5.9,5.9,"Y");
	//sumMC->SetAxisRange(-5.9,5.9,"Y");
	sumMC->GetYaxis()->SetNdivisions(505);
	//sumMC->GetYaxis()->SetTitle("#frac{Data}{MC}");
	sumMC->GetYaxis()->SetTitle("#frac{Data - MC}{#Delta}");
	sumMC->GetXaxis()->SetTitle(StackMC->GetHistogram()->GetXaxis()->GetTitle());
	sumMC->GetXaxis()->SetNdivisions(505);
	sumMC->GetYaxis()->SetTitleOffset(pad2->GetAbsHNDC()/pad->GetAbsHNDC());
	sumMC->SetLabelSize(0.05/pad2->GetAbsHNDC()*pad->GetAbsHNDC(),"x");
	sumMC->SetLabelSize(0.05/pad2->GetAbsHNDC()*pad->GetAbsHNDC(),"y");
	sumMC->SetTitleSize(0.06/pad2->GetAbsHNDC()*pad->GetAbsHNDC(),"x");
	sumMC->SetTitleSize(0.6*0.06/pad2->GetAbsHNDC()*pad->GetAbsHNDC(),"y");
	//sumMC->SetMarkerSize(pad2->GetAbsHNDC()*pad->GetAbsHNDC());
	sumMC->Draw("ep");
	TLine* line = new TLine;
	line->DrawLine(StackMC->GetHistogram()->GetXaxis()->GetXmin(),0.,StackMC->GetHistogram()->GetXaxis()->GetXmax(),0.);
	sumMC->Draw("ep same");


	//sumMC->Draw("ep");

	StackMC->GetYaxis()->SetTitleOffset(1.05*pad1->GetAbsHNDC()/pad->GetAbsHNDC());
	StackMC->GetXaxis()->SetLabelSize(0.05/pad1->GetAbsHNDC()*pad->GetAbsHNDC());
	StackMC->GetYaxis()->SetLabelSize(0.05/pad1->GetAbsHNDC()*pad->GetAbsHNDC());
	StackMC->GetXaxis()->SetTitleSize(0.06/pad1->GetAbsHNDC()*pad->GetAbsHNDC());
	StackMC->GetYaxis()->SetTitleSize(0.06/pad1->GetAbsHNDC()*pad->GetAbsHNDC());
	StackMC->GetYaxis()->Draw();
	
	pad->cd();
}

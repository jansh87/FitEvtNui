#include"TROOT.h"
void run()
{
	gROOT->ProcessLine(".L ../shared/SemiLepWeights2D_Donly.cc+O");
	gROOT->ProcessLine(".L ../shared/SemiLepDssWeights2D.cc+O");
	gROOT->ProcessLine(".L ../shared/SemiLepWeightsXc2D.cc+O");
	gROOT->ProcessLine(".L ../shared/pid/LepFake.cc+O");
	gROOT->ProcessLine(".L ../shared/pid/LepEff.cc+O");
	gROOT->ProcessLine(".L ../shared/pid/KPiEff.cc+O");
    gROOT->ProcessLine(".L ../shared/CFile.cc+O");
    gROOT->ProcessLine(".L ../../Dss/BF/DssWeightsParam.cc+O");
    gROOT->ProcessLine(".L ../../Ds/DsWeightsParam.cc+O");
	//gROOT->ProcessLine(".L ../shared/B2pilnuFF/formfactor.C+O");
	gROOT->ProcessLine(".L CreateHist.c+O");
	CreateHist();
}

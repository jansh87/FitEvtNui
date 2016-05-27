#include"../shared/include_all.h"
#include "../shared/CVar.h"
#include "../shared/CFile.h"
#include "../shared/CCont.h"
#include "../shared/tagcorrection.h"
#include "../shared/SemiLepWeights2D_Donly.h"
#include "../shared/SemiLepDssWeights2D.h"
#include "../shared/SemiLepWeightsXc2D.h"
#include "../../Dss/BF/DssWeightsParam.h"
#include "../../Ds/DsWeightsParam.h"
#include "../shared/DWeightsParam.h"
#include "../shared/branches.h"
#include "../shared/BrCorrection.h"
#include "../shared/pid/LepFake.h"
#include "../shared/pid/LepEff.h"
 
#include "TMatrixD.h"
#include "TVectorD.h"
 
#include"TEntryList.h"
 
 
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TError.h"
#include <ctime>
#include "error.h"
#include "residualPlot.h"
//globals
 
#define __pseudodata 0
#define __calcerrors 0
#define __toydata    0
 
 
 
typedef unsigned uint;
typedef unsigned short ushort;
using namespace std;
 
 
struct t_evt{
 
    float weight;
    float eff_error;
    float v1;
    float v0;
    short eff_Idx;
    short iLep;
    short bin;
    short component;
    short iB0;
    std::vector<short> classes;
 
};
struct t_evt_aux{
 
    t_evt* evt;
    //axis vars
    float v1;
    float v0;
    int iLep;
     
    t_evt_aux(t_evt* e):evt(e){};
};
 
struct t_axis{
    TAxis axis;
    vector<double>data;
    uint nBins;
 
};
 
vector<t_evt> evtlist;
 
std::vector<CVar> vVar;
std::vector<CCont> vCont;
std::vector<CFile> vFile;
std::vector<t_error*> vErr;
 
 
 
 
 
//"consts"
uint nErrors;
uint nComponents;
uint nBins;
TRandom3 rand3;
 
//vars
uint gi_Error=9999;
 
DssWeightsParam *DssWeight[4];
DsWeightsParam *DsWeight;
DWeightsParam *DWeight;
 
 
LepEff *ElecEff, *MuonEff;
LepFake *ElecFake, *MuonFake; 
 
CorrectBf CorrBf = CorrectBf();
 
TH1D *hData, *hMC;
TH1D** hMC_Comp;
TH1F **hPull;
 
clock_t *clocks;
int *counts;
uint fileCount;
t_axis axis[2][2][2];
enum {xtaunu, xlnu, other_e, other_mu};
 
 
 
t_axis makeAxis(vector<double>& v){
 
    t_axis tmp;
    tmp.axis=TAxis(v.size()-1, v.data());
    tmp.data = v;
    tmp.nBins = v.size()-1;
    return tmp;
 
}
 
void vFill(vector<double> &v, int N, double xmin, double xmax){
 
    for(int i = 0; i<=N; i++){
        v.push_back((xmax-xmin)/(N)*i + xmin);
    }
}
void makePseudoData( double *);
void DefineHists(){
 
 
//  evtlist.resize(1000000);
//components
    CCont tmp_cont;
 
    tmp_cont = CCont("Xtaunu", kRed);
    //tmp_cont.fl_const = true; tmp_cont.startVal = 0; //seitenband
    vCont.push_back(tmp_cont);
 
    tmp_cont = CCont("Xlnu", kGreen +2);
    vCont.push_back(tmp_cont);
 
    tmp_cont = CCont("other_e", kOrange -2);
    vCont.push_back(tmp_cont);
     
    tmp_cont = CCont("other_mu", kOrange -1);
    tmp_cont.fl_const = true; tmp_cont.startVal = 1;
    vCont.push_back(tmp_cont);
    nComponents = vCont.size();
    cout<<nComponents<<endl;
     
     
    vector<double> plepel = { 0.5,    0.6,    0.8,    1.,     1.2,    1.4,  2.5};//7
    vector<double> plepmu = { 0.5,    0.6,    0.8,    1.,     1.2,    1.4,     2.5};//6
     
    vector<double> plep2el = {  0.5,  0.75,   1.,       2.5};//3
    vector<double> plep2mu = {  0.7,  0.85,   1.,   2.5};//3
    vector<double> plep = {   0.3,   0.9,   2.5};//3
        //seitenband
     
    vector<double> mm2 = { 0., 2.5, 4, 5.5, 7, 10., 20.};
    vector<double> me = { -1, 1, 2, 3.5, 5.};
    vector<double> mm2SB = {0, 2.5};
    vector<double> plepSB = {1, 2.5};
    //define axis 
    axis[0][0][0] = makeAxis( plep2el );
    axis[0][0][1] = makeAxis( mm2 );
    axis[0][1][0] = makeAxis( plep2mu );
    axis[0][1][1] = makeAxis( mm2 );
     
    nBins = axis[0][0][0].axis.GetNbins()*axis[0][0][1].axis.GetNbins() + axis[0][1][0].axis.GetNbins()*axis[0][1][1].axis.GetNbins();
     
     
     
    vector<double> plepFineEl, plepFineMu, mm2Fine;
    vFill(plepFineEl,20,0.5,2.5);
    vFill(plepFineMu,20,0.7,2.5);
    vFill(mm2Fine,20,0,20);
    //fine plot binning
    axis[1][0][0] = makeAxis( plepFineEl );
    axis[1][0][1] = makeAxis( mm2Fine );
    axis[1][1][0] = makeAxis( plepFineMu );
    axis[1][1][1] = makeAxis( mm2Fine );
     
     
     
     
    hData = new TH1D("hdata","",nBins, 0, nBins);
    hMC   = new TH1D("hmc", "" ,nBins, 0, nBins);
    hMC_Comp = new TH1D*[nComponents];
     
 
     
    for(int i = 0; i<nComponents; i++){
        hMC_Comp[i] = new TH1D(Form("hmc_%i",i), "" ,nBins, 0, nBins);
    }
 
 
}
 
 
double sqr(double x){return x*x;};
#include "fitroutine.inc"
 
void Init(){
 
//files
    CFile tmp_file;
 
    string path = string("~/rootFiles/smallFiles/");
 
    float nStreams =1.;
    float nContiStreams =1.;
    
    
    
    
    //MC

//pseudodata


   

    tmp_file =CFile(path, "charged_s4.root", "tree",'o'); tmp_file.weight = 1./nStreams;
    vFile.push_back(tmp_file );
    
    tmp_file =CFile(path, "continuum_s4.root", "tree",'c'); tmp_file.weight = 1./nContiStreams;
    vFile.push_back(tmp_file );
    
    tmp_file = CFile(path,"DssMC.root","tree",'s'); tmp_file.weight = 0.405*2;// 0.388;
    vFile.push_back(tmp_file );
    
     
    tmp_file =CFile(path, "data.root", "tree",'d');
    if(!(__pseudodata || __toydata))vFile.push_back(tmp_file );
   
    tmp_file =CFile(path, "mixed_s4.root", "tree",'o'); tmp_file.weight = 1./nStreams;
    vFile.push_back(tmp_file );
    
     
 
/*  
tmp_file =CFile(path, "charged_s1.root", "tree",'o'); tmp_file.weight = 1./nStreams;
    vFile.push_back(tmp_file );
     
    tmp_file =CFile(path, "mixed_s1.root", "tree",'o'); tmp_file.weight = 1./nStreams;
    vFile.push_back(tmp_file );
     
    tmp_file =CFile(path, "continuum_s1.root", "tree",'c'); tmp_file.weight = 1./nContiStreams;
    vFile.push_back(tmp_file );

 tmp_file =CFile(path, "charged_s2.root", "tree",'o'); tmp_file.weight = 1./nStreams;
    vFile.push_back(tmp_file );
     
    tmp_file =CFile(path, "mixed_s2.root", "tree",'o'); tmp_file.weight = 1./nStreams;
    vFile.push_back(tmp_file );
     
    tmp_file =CFile(path, "continuum_s2.root", "tree",'c'); tmp_file.weight = 1./nContiStreams;
    vFile.push_back(tmp_file );
 
 
  tmp_file =CFile(path, "charged_s3.root", "tree",'o'); tmp_file.weight = 1./nStreams;
    vFile.push_back(tmp_file );
     
    tmp_file =CFile(path, "mixed_s3.root", "tree",'o'); tmp_file.weight = 1./nStreams;
    vFile.push_back(tmp_file );
     
    tmp_file =CFile(path, "continuum_s3.root", "tree",'c'); tmp_file.weight = 1./nContiStreams;
    vFile.push_back(tmp_file );
 
  tmp_file =CFile(path, "charged_s4.root", "tree",'o'); tmp_file.weight = 1./nStreams;
    vFile.push_back(tmp_file );
     
    tmp_file =CFile(path, "mixed_s4.root", "tree",'o'); tmp_file.weight = 1./nStreams;
    vFile.push_back(tmp_file );
     
    tmp_file =CFile(path, "continuum_s4.root", "tree",'c'); tmp_file.weight = 1./nContiStreams;
    vFile.push_back(tmp_file );
*/
 
    
 
    tmp_file =CFile(path, "ulnu.root", "tree",'u'); tmp_file.weight = 1./20;
    vFile.push_back(tmp_file );
    
    
    tmp_file =CFile(path, "rare.root", "tree",'r'); tmp_file.weight = 1./50;
    //vFile.push_back(tmp_file );
     
    
    
    
    
    
    
     
    
    tmp_file =CFile(path, "dilepcharged_s1.root", "tree",'o'); tmp_file.weight = 1./nStreams;
    //vFile.push_back(tmp_file );
     
    tmp_file =CFile(path, "dilepmixed_s1.root", "tree",'o'); tmp_file.weight = 1./nStreams;
    //vFile.push_back(tmp_file );
     
    tmp_file =CFile(path, "dilepcontinuum_s1.root", "tree",'c'); tmp_file.weight = 1./nContiStreams;
    //vFile.push_back(tmp_file );
 
    tmp_file = CFile(path,"dilepDssMC.root","tree",'s');    tmp_file.weight = 0.405;// 0.388;
    //vFile.push_back(tmp_file );
 
    tmp_file =CFile(path, "dilepulnu.root", "tree",'u'); tmp_file.weight = 1./20;
    //vFile.push_back(tmp_file );
    tmp_file =CFile(path, "dileprare.root", "tree",'r'); tmp_file.weight = 1./50;
    //vFile.push_back(tmp_file );
    
    //pseudodata
    tmp_file =CFile(path, "dilepdata.root", "tree",'d');
    //if(!__pseudodata)vFile.push_back(tmp_file );
 
 
//init 
 
    string sCharm[] = {"D1", "D2", "D1'", "D0*"};
    for(int i = 0; i<4 ;i++)
        DssWeight[i] = new DssWeightsParam(nBins,sCharm[i]);
    DsWeight = new DsWeightsParam(nBins, "D*");
    DWeight = new DWeightsParam(nBins, "D");
     
 
    MuonFake    = new LepFake(nBins,"muon fake1");
    ElecFake    = new LepFake(nBins,"elec fake1");
 
    MuonEff     = new LepEff(nBins,"muon eff1",0.9,1);
    ElecEff     = new LepEff(nBins,"elec eff1",0.9,0);
 
 
     
//systematics
    vErr.push_back( new t_bf("Fit",99,99,0));
    vErr.push_back( new t_bf("Dtaunu",1,1,0.1,t_error::_xtaunu));
    t_error::block_start[t_error::_xtaunu] = vErr.size()-1;
    vErr.push_back( new t_bf("D*taunu",1,3,0.15,t_error::_xtaunu));
    vErr.push_back( new t_bf("D1taunu",1,5,0.5,t_error::_xtaunu));
    vErr.push_back( new t_bf("D2*taunu",1,6,0.5,t_error::_xtaunu));
    vErr.push_back( new t_bf("D*0taunu",1,7,0.5,t_error::_xtaunu));
    vErr.push_back( new t_bf("Dp1taunu",1,8,0.5,t_error::_xtaunu));
    t_error::block_end[t_error::_xtaunu] = vErr.size()-1;
                    //pdg eval, pdg fit
                    //0.028
    vErr.push_back( new t_bf("BDlnu",2,1,0.045,t_error::_xlnu));
    t_error::block_start[t_error::_xlnu] = vErr.size()-1;
                    //0.04
    vErr.push_back( new t_bf("BD*lnu",2,3,0.045,t_error::_xlnu));
    //HFAG:
    vErr.push_back( new t_bf("D1lnu",2,5,0.5,t_error::_xlnu));
    vErr.push_back( new t_bf("D2*lnu",2,6,0.1,t_error::_xlnu));
    vErr.push_back( new t_bf("D*0lnu",2,7,0.17,t_error::_xlnu));    
    vErr.push_back( new t_bf("Dp1lnu",2,8,0.28,t_error::_xlnu));
    vErr.push_back( new t_bf("Dpilnu",2,2,1.,t_error::_xlnu)); //neu, vorher 4.5%
    vErr.push_back( new t_bf("D*pilnu",2,4,0.8,t_error::_xlnu)); //eigentlich 130%
    vErr.push_back( new t_bf("D2Slnu",2,9,1.,t_error::_xlnu));
    vErr.push_back( new t_bf("D*2Slnu",2,10,1.,t_error::_xlnu));
    t_error::block_end[t_error::_xlnu] = vErr.size()-1;
     
    /*
    vErr.push_back( new t_bf("BaryonX",8,989,0.6/6.8));
    vErr.push_back( new t_bf("XcXc",4,989,0.10));
    //vErr.push_back( new t_bf("XcXu",6,989,0.20);
    vErr.push_back( new t_bf("cc-meson",9,989,0.1));
    //vErr.push_back( new t_bf("continuum",0,0,0.1));
    vErr.push_back( new t_bf("D0->Xlnu",989,1000,0.61/13.19,t_error::_Dxlnu));
    t_error::block_start[t_error::_Dxlnu] = vErr.size()-1;
    ((t_bf*)vErr[vErr.size()-1])->BF[0] = 13.;//32.2;
    vErr.push_back( new t_bf("D+->Xlnu",989,2000,3.2/33.67,t_error::_Dxlnu));
    ((t_bf*)vErr[vErr.size()-1])->BF[0] = 32.13;
    vErr.push_back( new t_bf("D0->X",989,1100,0.,t_error::_Dxlnu,t_error::_dummy));
    vErr.push_back( new t_bf("D+->X",989,2100,0.,t_error::_Dxlnu,t_error::_dummy));
    t_error::block_end[t_error::_Dxlnu] = vErr.size()-1;
     */
    //vErr.push_back( new t_bf("f+0",989,989,(float)0.024/1.058,'a');
    
     
    vErr.push_back( new t_eff("FakeElec",0,1,ElecFake,nBins));
    vErr.push_back( new t_eff("EffElec",0,0,ElecEff,nBins));
    vErr.push_back( new t_eff("FakeMuon",1,1,MuonFake,nBins));
    vErr.push_back( new t_eff("EffMuon",1,0,MuonEff,nBins));
     
 
    vErr.push_back( new t_FF("D_rho2",2,1,DWeight,nBins));
    DWeight->setErrIdx(0,vErr.size()-1);
     
    vErr.push_back( new t_FF("D*_FF1",2,3,DsWeight,nBins));
    DsWeight->setErrIdx(0,vErr.size()-1);
    vErr.push_back( new t_FF("D*_FF2",999,999,DsWeight,nBins));
    DsWeight->setErrIdx(1,vErr.size()-1);
    vErr.push_back( new t_FF("D*_FF3",999,999,DsWeight,nBins));
    DsWeight->setErrIdx(2,vErr.size()-1);
     
     
    string sPar[] = {"_tau", "_eta", "_B12"};
    for(int i = 0; i<4; i++){
        for(int j = 0; j<3; j++){
            vErr.push_back( new t_FF(sCharm[i]+sPar[j],(j==0?2:999),5+i,DssWeight[i],nBins,i));
            DssWeight[i]->setErrIdx(j,vErr.size()-1);
        }
    }
     
    nErrors = vErr.size();
     
    clocks = new clock_t[nErrors];
    counts = new int[nErrors];
    for(int i = 1; i< nErrors; i++){
        clocks[i] = 0;
        counts[i] =0;
    }
 
 
   hPull = new TH1F*[nComponents+nErrors];
   for(int i = 0; i< nComponents+nErrors; i++)
   	hPull[i] = new TH1F((i<nComponents)? vCont[i].Name.c_str(): vErr[i-nComponents]->Name.c_str() ,(i<nComponents)? vCont[i].Name.c_str(): vErr[i-nComponents]->Name.c_str() ,21,-5,5);
 
}
 
 
int globalBin(t_evt_aux& evt){
 
     
    int b1 = axis[0][evt.iLep][1].axis.FindBin(evt.v1)-1;
    int b1n = axis[0][evt.iLep][1].axis.GetNbins();
    int b2 = axis[0][evt.iLep][0].axis.FindBin(evt.v0)-1;
    //cout<<evt.mm2<<" "<<evt.pslep<<" "<<b1<<" "<<b2<<" "<< b1 + b2*b1n<<" "<<(evt.iLep?axis[0][0].GetNbins()+axis[0][1].GetNbins():0)<<endl;
    int bin = b1 + b2*b1n + (evt.iLep?axis[0][0][0].axis.GetNbins()*axis[0][0][1].axis.GetNbins():0);
     
     
    return bin;
 
}
 
 
 
void ReadFiles(){
 
    TTree *t;
    fileCount = 0;
    for(vector<CFile>::iterator it_f = vFile.begin(); it_f != vFile.end(); ++it_f)
    {
    	fileCount++;
        if(!it_f->File)
            it_f->File = new TFile( it_f -> FileName.c_str() ,"read");
        t = (TTree*)it_f->File->Get(it_f->Tree.c_str());
        if(!t) cout<<"no tree"<<endl;
         
         
        t->SetBranchAddress("run", &run);
        t->SetBranchAddress("btag", &btag);
        t->SetBranchAddress("lep1", &lep1);
        t->SetBranchAddress("gx", &gx);
        t->SetBranchAddress("event", &event);
        t->SetBranchAddress("truth", &truth);
        t->SetBranchAddress("gmiss", &gmiss);
        t->SetBranchAddress("sBDecay", &sBDecay);
        t->SetBranchAddress("sDDecay", &sDDecay);
 
 
        int nEntries = t->GetEntries();
         
        cout<<"Read "<<it_f->FileName<<" with "<<nEntries<<" Events and weight "<<it_f->weight<<endl;
 
        //event loop
         
        for(int i =0; i<nEntries; ++i)
        {
             
            t->GetEntry(i);
        //cuts
                 
            if( btag.m_bc<5.27) continue;
            if( btag.pcode_b*lep1.q>0/* && abs(btag.pcode_b) == 521*/) continue;
            if(abs(btag.pcode_b) == 511) continue;
            if( (i%2==0) && (it_f->Type == 's')) continue; 
            
            //if( ((i+2)%4)!=0 && (it_f->Type == 's')){ continue; }
            
            if( log(btag.NB)<-4) continue;
            if(event.cos_thrAm>0.8) continue;
             
            //remove D**lnu from generic
            if( it_f->Type == 'o'){ if(  abs(event.lclass) == 2 && (event.dclass > 4/* && event.dclass != 5 && event.dclass != 7 && event.dclass != 8*/)) continue; }
           // add only true D**lnu events
            else if( it_f->Type == 's'){ if(!(abs(event.lclass) == 2 && (event.dclass > 4/* && event.dclass != 5 && event.dclass != 7 && event.dclass != 8*/)))continue;}
 
            if(event.lclass == -2 && event.dclass == 11){
                string sBDec(sBDecay);
                CorrectBf::NoCharge(sBDec);
                if(sBDec.find("511_100411") == 0){
                    event.dclass = 9;       
                } else if (sBDec.find("511_100413") == 0){
                    event.dclass = 10;
                }
                 
            } else if(event.lclass == 2 && event.dclass == 9){
                string sBDec(sBDecay);
                CorrectBf::NoCharge(sBDec);
                if(sBDec.find("521_100421") == 0){
                    event.dclass = 9;       
                } else if (sBDec.find("521_100423") == 0){
                    event.dclass = 10;
                }
            }
            
            if(event.dclass == 10) continue;
 
 
             
            //if(it_f->Type == 's') 
            gx.m2 = sqrt(gx.m2);
            //sigmal verschieben
          //  if(abs(event.lclass) == 1)
           // 	 lep1.ps += 1.1;
          //  if(gx.m2<2.4 && plep.ps<1) continue;
            //if(event.nch!=) continue;
        //vars
         
            int iLep    = lep1.fl_lep%10;
            int iFake   = (lep1.fl_lep != 10) && (lep1.fl_lep != 21);
             
            t_evt evt;
            t_evt_aux evt_aux(&evt);
         
            
         
            evt.v1 = gmiss.m2;
            evt.v0 = lep1.ps;
            evt.iLep = iLep;
             
            evt_aux.v1 = gmiss.m2; // wird nur f[r bin gebraucht
            evt_aux.v0 = lep1.ps;
            evt_aux.iLep = iLep;
             
             
             
            if(evt.v1< axis[0][iLep][1].axis.GetXmin() || evt.v1> axis[0][iLep][1].axis.GetXmax()) continue;
            if(evt.v0< axis[0][iLep][0].axis.GetXmin() || evt.v0> axis[0][iLep][0].axis.GetXmax()) continue;
         
        //component
             
             
            if(it_f -> Type != 'd'){
                if(abs(event.lclass) == 2 && !iFake)
                    evt.component = xlnu;
                else if(abs(event.lclass) == 1 && !iFake) //seitenband
                    evt.component = xtaunu;
                else if(iLep == 0)
                    evt.component = other_e;
                else if(iLep == 1)
                    evt.component = other_mu;
                else
                    exit(-1);
             
             
             
                 
             
        //weights
                double weight = 1.;
                 
                /*string tmp = string(sDDecay);
                CorrectBf::NoCharge(tmp);
                if(tmp.find("_311_") != string::npos || tmp.find("(311_") != string::npos || tmp.find("_311)") != string::npos)
                    weight*=0.9;
                */
                weight *= tagcorr(btag.b_mode, (double)btag.NB);
                if(it_f->Type == 'o' || it_f->Type == 'c')
                    weight *= genMCCorr(run.exp);
                weight *= it_f->weight;
            if(it_f->Type == 's'){ 
                if(abs(event.lclass) == 2 || abs(event.lclass) == 1){
                        weight *= CorrBf.weightB(event.lclass, event.dclass);
                    //cout<<CorrBf.weightB(event.lclass, event.dclass)<<" "<<event.lclass<<" "<<event.dclass<<endl;
                }
                else
                    weight *= CorrBf.weightB(sBDecay);
        
         
         	if(abs(event.lclass) == 2 && event.dclass >4){
				
			weight *= CorrBf.FixDss(sDDecay);
		
		}
	    }
		weight *= CorrBf.weightD(sDDecay);
         	
         	//double wtmp = weight;
            //  if(globalBin(evt_aux)<0) cout<<"minus"<<endl;
             	double eff_error;
                if(lep1.fl_lep==10){//true e
                    weight *= ElecEff->weight(lep1.fl_lep%10,run.exp,acos(lep1.th)/M_PI*180., lep1.p, globalBin(evt_aux), 0, weight);
                    eff_error = ElecEff->weightError(lep1.fl_lep%10,run.exp,acos(lep1.th)/M_PI*180., lep1.p);        
                }else if(lep1.fl_lep%10 == 0 ){//fake e
                    weight *= ElecFake->weight(lep1.fl_lep%10, lep1.fl_lep/10-1, acos(lep1.th)/M_PI*180., lep1.p,globalBin(evt_aux), 0, weight);
                     eff_error = ElecFake->weightError(lep1.fl_lep%10, lep1.fl_lep/10-1, acos(lep1.th)/M_PI*180., lep1.p);
                //muons
                }else if(lep1.fl_lep==21){//true mu
                    weight *= MuonEff->weight(lep1.fl_lep%10,run.exp,acos(lep1.th)/M_PI*180., lep1.p,globalBin(evt_aux), 0, weight);
                    eff_error = MuonEff->weightError(lep1.fl_lep%10,run.exp,acos(lep1.th)/M_PI*180., lep1.p);   
                }else if(lep1.fl_lep%10 == 1 ){//fake mu
                    weight *= MuonFake->weight(lep1.fl_lep%10, lep1.fl_lep/10-1, acos(lep1.th)/M_PI*180., lep1.p,globalBin(evt_aux),0, weight);
                    eff_error = MuonFake->weightError(lep1.fl_lep%10, lep1.fl_lep/10-1, acos(lep1.th)/M_PI*180., lep1.p);
             	} else eff_error = 0;
                
                 
                if(abs(event.lclass) == 2){
                    if(event.dclass<5){
                        if(event.dclass == 1){
                            //cout<<truth.w<<" "<<truth.costh<<" "<<globalBin(evt_aux)<<endl;
                         
                            weight *= DWeight->weight(truth.w, truth.costh,truth.q2, truth.p_l, globalBin(evt_aux)); // weight mit uebergeben
                            //cout<<truth.w<<" "<<truth.costh<<" "<<globalBin(evt_aux)<<endl;
                         
                        }else{
                            weight *= DsWeight->weight(truth.w, -truth.costh, globalBin(evt_aux), weight);
                        }
                    }
                    else
                    {
                        weight *= DssWeight[event.dclass-5]->weight(event.dclass-5,truth.w, truth.costh, globalBin(evt_aux), weight);
                    }
 
                }
 
                 
      		//weight = wtmp;
                 
                 
                 
                evt.weight = weight;
                evt.eff_error = eff_error;
                evt.bin = globalBin(evt_aux);
                evt.iB0 = abs(btag.pcode_b) == 511; //0: B+ 1: B0
                 
                
                
                
             
             
                string sDDec(sDDecay);
                CorrectBf::NoCharge(sDDec);
                vector<string> *DDec = CorrectBf::SplitDString(sDDec);
                for(int i = 0; i<DDec->size(); i++){
                    if((*DDec)[i].find("411") == 0){
                        if((*DDec)[i].find("11_12")!= string::npos || (*DDec)[i].find("13_14")!= string::npos){
                            event.dclass += 2000;
                            break;
                        } else {
                            event.dclass += 2100;
                            break;
                        }
                    } else if((*DDec)[i].find("421") == 0){
                        if((*DDec)[i].find("11_12")!= string::npos || (*DDec)[i].find("13_14")!= string::npos){
                            event.dclass += 1000;
                            break;
                        }else {
                            event.dclass += 1100;
                            break;
                        }
                    }
                    /*if((*DDec)[i].find("411_311_11_12")!= string::npos  || (*DDec)[i].find("411_311_13_14")!= string::npos){
                        event.dclass += 100;
                        break;
                    } else if((*DDec)[i].find("411_313_11_12")!= string::npos  || (*DDec)[i].find("411_313_13_14")!= string::npos){
                        event.dclass += 200;
                        break;
                    } */
                }
                delete DDec;
            //fill
         
                int lclass, dclass;
                for(int i = 0; i<nErrors; i++){
                    if(vErr[i]->Type == t_error::_eff){
                        lclass = iLep;
                        dclass = iFake;
                    } else {
                        lclass = event.lclass;
                        dclass = event.dclass%100;
                    }
                    if(vErr[i]->lclass == lclass || vErr[i]->lclass == 989 ){
                        if(vErr[i]->dclass%100 == dclass || vErr[i]->dclass == 989 ){
                         	
                            if(vErr[i]->Type == t_error::_eff)
                            	evt.eff_Idx = i;
                            else 
                            	evt.classes.push_back(i);
                         
                        }
                    }
                    if(vErr[i]->dclass>1000 && vErr[i]->dclass/100 == event.dclass/100)
                        evt.classes.push_back(i);
             
                }
             //streamtest
             //	if(fileCount < 6){
              //  	hData->AddBinContent(evt.bin+1,weight);
              //  }else{ 
                	hMC_Comp[evt.component]->AddBinContent(evt.bin+1,weight);
                	evtlist.push_back(evt);
               // }
               
                 
                //signal verschieben && _toydata
                //if(evt.component == xtaunu)
                 // hData->AddBinContent(globalBin(evt_aux)+1,weight);
         
                 
                if(__pseudodata){
                   // evt_aux.v1 += 0.04;
                    hData->AddBinContent(globalBin(evt_aux)+1,weight);
         	}
            } else {//end non-data
             if(!__toydata)
             	//shift data
             	//evt_aux.v1 += 0.04;
                hData->AddBinContent(globalBin(evt_aux)+1);
                 
            }//end data
        }//end event loop
         
        delete t;
    }// end file loop
    
    
    //toyMC
    
 
     
    THStack *stack = new THStack("","");
    for(int i = nComponents -1; i>=0; i--){
        hMC_Comp[i] -> SetFillColor(vCont[i].Color);
        hMC_Comp[i] -> SetLineColor(vCont[i].Color);
        stack->Add(hMC_Comp[i]);
    }
     
    TCanvas *cv = new TCanvas("","",1000,1000);
    DrawResiduals(stack,hData,cv);
    hData->SetMarkerStyle(20);
    hData->SetMarkerColor(20);
     
    cv->Print("hists.pdf");
     
 
}


 
void PostFitPlots(double *params);
void PerfFit(){
 	
    
   
 
 
    double **xs;
    
    double ***errors;
    double *pseudoPars = new double [nComponents + nErrors];
     if(__toydata){
     hData->Reset("ICESM");
     makePseudoData(pseudoPars);
    }
    xs = new double*[nErrors];
    errors = new double**[2];
    for(int u = 0; u<2; u++)
        errors[u] = new double*[nErrors];
    for(int e = 0; e< nErrors; e++)
    {
        xs[e] = new double[nComponents+nErrors-1];
        for(int k = 0; k < nComponents+nErrors-1; k++)
        if(k<nComponents) xs[e][k] = 1;
        else		  xs[e][k] = 0;
        for(int u = 0; u<2; u++)
        {
            errors[u][e] = new double[nComponents+nErrors-1];
        }
    }
    
    PostFitPlots(xs[0]);
   // return;
    
    //RunChi2FitNui(xs,errors);
    double **finalErr;
    int VarOfInterest =0;
    finalErr = new double*[2];
    for(int i = 0; i< 2; i++)
    {
        finalErr[i] = new double[nErrors];
    }
    FILE *fp;
    if(__toydata)
    	fp = fopen("pull_values.txt","w");
    for(gi_Error = 0; gi_Error<(__calcerrors?nErrors:1); gi_Error++)
    {
     
        if(vErr[gi_Error]->Name.find("D2_") == 0 || vErr[gi_Error]->Name.find("D0*_") == 0) continue;
        
        for(int i = 0; i<(__toydata?500:1); i++){
        	if(__toydata){
     			makePseudoData(pseudoPars);
       			cout<<"Make Pseudo Data "<<i<<endl;
     		}
		RunChi2FitNui(xs, errors);
		if(__toydata){
			for(int i = 0; i<nComponents+nErrors; i++){
				hPull[i] -> Fill( (xs[0][i] - pseudoPars[i])/fabs(errors[(xs[0][i] - pseudoPars[i]) > 0][0][i] ));
				fprintf(fp,"%s\t%i\t%lf\t%lf\t%lf\n",(i<nComponents)? vCont[i].Name.c_str(): vErr[i-nComponents]->Name.c_str(),i,pseudoPars[i],xs[0][i],fabs(errors[(xs[0][i] - pseudoPars[i]) > 0][0][i] ));
		     	}
	    	}
	}
        if(gi_Error == 0)
        {
            finalErr[0][gi_Error] = errors[0][gi_Error][VarOfInterest];
            finalErr[1][gi_Error] = errors[1][gi_Error][VarOfInterest];
        } else {
         
            finalErr[0][gi_Error] = sqrt(fabs(sqr(errors[0][gi_Error][VarOfInterest])-sqr(errors[0][0][VarOfInterest])));
            finalErr[1][gi_Error] = sqrt(fabs(sqr(errors[1][gi_Error][VarOfInterest])-sqr(errors[1][0][VarOfInterest])));
                 
        }
    }
    if(__toydata)
    	fclose(fp);
    for(gi_Error = 0; gi_Error<nErrors; gi_Error++)
        cout<<vErr[gi_Error]->Name<<"\t"<<finalErr[0][gi_Error]<<"\t"<<finalErr[1][gi_Error]<<endl;
 
    PostFitPlots(xs[0]);
}
 
#include "plotting.h"
void PostFitPlots(double *params){
 
 
    const double *NuiPar = params + nComponents;
     
     
     
    char sHistName[500];
    char sStackName[500];
     
    string sBinning[] = {"cc","cf","fc","ff"};
    TH2F ****hists = new TH2F***[4];
    TH2F ***histsData = new TH2F**[4];
    THStack ***hStack = new THStack**[4];
     
    for( int j =0; j<4;j++)//binning type: fine-fine, coarse-fine, fine-coarse- coarse-coarse 
    {
        hists[j]    = new TH2F**[2];
        histsData[j] = new TH2F*[2];
        hStack[j] = new THStack*[2];
        for(int l = 0; l<2; l++)
        {
            sprintf(sHistName,"hStack_%s_%i",sBinning[j].c_str(),l);
            hStack[j][l] = new THStack(sHistName,sHistName);
             
            sprintf(sHistName,"histData_%s_%i",sBinning[j].c_str(),l);
            histsData[j][l] = new TH2F(sHistName,sHistName,axis[j/2][l][0].nBins, axis[j/2][l][0].data.data(),axis[j%2][l][1].nBins, axis[j%2][l][1].data.data());
            histsData[j][l] -> SetLineColor(kBlack);
            histsData[j][l] -> SetFillColor(kBlack);
 
            histsData[j][l] -> SetLineWidth( 2 );
             
            hists[j][l] = new TH2F*[nComponents];
            for(int i = nComponents -1; i>=0; i--)
            {
                sprintf(sHistName,"hist_%s_%s_%i",sBinning[j].c_str(),vCont[i].Name.c_str(),l);
                hists[j][l][i] = new TH2F(sHistName,sHistName,axis[j/2][l][0].nBins, axis[j/2][l][0].data.data(),axis[j%2][l][1].nBins, axis[j%2][l][1].data.data());
                hists[j][l][i] -> SetLineColor(vCont[i].Color);
                hists[j][l][i] -> SetFillColor(vCont[i].Color);
 
                //hists[i][j] -> SetFillStyle(3354);
                hists[j][l][i] -> SetLineWidth( 2 );
     
                hStack[j][l]->Add(hists[j][l][i]);
            }
             
                 
        }
    }
     if(__toydata){
	     TCanvas *cc1 = new TCanvas("pull","pull",500,500);
	     cc1 -> Print("pulls.pdf[");
	     for(int i = 0; i< nComponents + nErrors; i++)
	     {
	     
	     	cout<<"Pull for "<<((i<nComponents)? vCont[i].Name.c_str(): vErr[i-nComponents]->Name.c_str())<<endl;
	     	hPull[i]->Draw("hist");
	     	hPull[i]->GetXaxis()->SetTitle( (i<nComponents)? vCont[i].Name.c_str(): vErr[i-nComponents]->Name.c_str() );
	     	hPull[i]->Fit("gaus","L");
	     	hPull[i]->GetFunction("gaus")->Draw("same");
	     	//hPull[i]->SetStats();
	     	cc1->Print("pulls.pdf");
	     }
	     
	     cc1->Print("pulls.pdf]");
	     delete cc1;
     }
     
    //decorrelate
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
    static double BFXtaunu[] = {2.43,2.33};
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
    uint nEvents = evtlist.size();
    t_error::t_errType type;
    double c,change, BF_del, weight;
    uint nClasses, classIdx, component, j, iB0;
    for(uint i = 0; i<nEvents; ++i)
    {
     
        change =1.;
        BF_del = 1.;
        nClasses = evtlist[i].classes.size();
        for(j = 0; j<nClasses; ++j)
        {
            classIdx = evtlist[i].classes[j];
            c = vErr[classIdx]->getChange(evtlist[i].bin);
            type = vErr[classIdx]->Type;
            if(type != t_error::_FF){ 
                c = (1. + c*NuiPar[classIdx]);
                if(type == t_error::_bf){
                    iB0 = evtlist[i].iB0;
                    if(vErr[classIdx]->block == t_error::_xlnu)
                        c*=dBFXlnu[iB0];
                    else if(vErr[classIdx]->block == t_error::_xtaunu)
                        c*=dBFXtaunu[iB0];
                }
            }
            change*= c; 
        }
         
        component = evtlist[i].component;
        weight = evtlist[i].weight;
         
        if(__pseudodata){
            for(int j = 0; j<4; j++)
                histsData[j][evtlist[i].iLep]->Fill(evtlist[i].v0,evtlist[i].v1,evtlist[i].weight);
        }
         
        evtlist[i].weight = weight*change*params[component]*(component == xtaunu?params[xlnu]:1.);
     
    }
     
     
     
     
     
    for(uint i = 0; i<nEvents; ++i)
    {
     
        for(int j = 0; j<4; j++){
                hists[j][evtlist[i].iLep][evtlist[i].component] -> Fill(evtlist[i].v0,evtlist[i].v1,evtlist[i].weight);
        }   
    }   
     
    TTree *t;
    for(vector<CFile>::iterator it_f = vFile.begin(); it_f != vFile.end(); ++it_f)
    {
        if(it_f -> Type != 'd' || __pseudodata) continue;
        if(!it_f->File)
            it_f->File = new TFile( it_f -> FileName.c_str() ,"read");
        t = (TTree*)it_f->File->Get(it_f->Tree.c_str());
        if(!t) cout<<"no tree"<<endl;
         
         
        t->SetBranchAddress("run", &run);
        t->SetBranchAddress("btag", &btag);
        t->SetBranchAddress("lep1", &lep1);
        t->SetBranchAddress("gx", &gx);
        t->SetBranchAddress("event", &event);
        t->SetBranchAddress("truth", &truth);
        t->SetBranchAddress("gmiss", &gmiss);
        t->SetBranchAddress("sBDecay", &sBDecay);
        t->SetBranchAddress("sDDecay", &sDDecay);
 
 
        int nEntries = t->GetEntries();
         
        cout<<"Read "<<it_f->FileName<<" with "<<nEntries<<" Events and weight "<<it_f->weight<<endl;
 
        //event loop
         
        for(int i =0; i<nEntries; ++i)
        {
             
            t->GetEntry(i);
        //cuts
                 
            if( btag.m_bc<5.27) continue;
            if( btag.pcode_b*lep1.q>0/* && abs(btag.pcode_b) == 521*/) continue;
            if( log(btag.NB)<-4) continue;
            if(event.cos_thrAm>0.8) continue;
           //  gx.m2 = sqrt(gx.m2);
           // if(gx.m2<2.4) continue;
        //vars
         
            int iLep    = lep1.fl_lep%10;
            int iFake   = (lep1.fl_lep != 10) && (lep1.fl_lep != 21);
             
            t_evt evt;
             
            evt.v1 = gmiss.m2;
            evt.v0 = lep1.ps;
             
             
            if(evt.v1< axis[0][iLep][1].axis.GetXmin() || evt.v1> axis[0][iLep][1].axis.GetXmax()) continue;
            if(evt.v0< axis[0][iLep][0].axis.GetXmin() || evt.v0> axis[0][iLep][0].axis.GetXmax()) continue;
         
             
             
            for(int j = 0; j<4; j++){        
                histsData[j][iLep] -> Fill(evt.v0, evt.v1);
            }
        }
        delete t;
    }
     
         
     
    for(int j = 0; j<4; j++){
        if(j == 0 || j == 3)
            PlotProj(hists[j], histsData[j], sBinning[j]);
        if(j>2) continue;
        PlotSlices("mm2",0,hists[j], histsData[j], sBinning[j]);
        PlotSlices("plep",1,hists[j], histsData[j], sBinning[j]);
    }
}
 void makePseudoData(double * par){

	for(int i = 0; i<nComponents; i++)
		par[i] = 1;
	
	for(int i = nComponents; i<nComponents+nErrors;i++)
		par[i] = rand3.Gaus(0,1);	

	getNuiModMC(par,hMC_Comp);
	
	hData->Reset("ICESM");
	for(int i = 0; i<nComponents; i++){
		double tmp = hMC_Comp[i] -> Integral();
		if(i == xlnu) tmp+=hData->Integral();
		for(int j = 1; j<=nBins; j++)
			hMC_Comp[i]->SetBinContent(j, rand3.Poisson(hMC_Comp[i]->GetBinContent(j)));
		
		hData->Add(hMC_Comp[i],par[i]);  //R(X) only if par == 1!!
		par[i] = hMC_Comp[i] -> Integral()/tmp;
		if(i == xlnu) par[i] = hData -> Integral()/tmp;
	}
		
	
}
 
void CreateHist(){
 
    gStyle -> SetHistLineWidth( 1 );
    gStyle -> SetHistFillColor( kWhite );
//  gStyle -> SetLegendFillColor ( kWhite ); // only in root 5.30
    gStyle -> SetCanvasColor(kWhite);
    gStyle -> SetFrameFillColor(kWhite);
    gStyle -> SetTitleFillColor( kWhite );
    gStyle -> SetTitleBorderSize( 0 );
    //gStyle -> SetOptStat( 0 );
    gStyle->SetHatchesSpacing(2);
    //gStyle -> SetTitleFontSize( 2 );
    //gStyle -> SetTitleAlign ( 0 );
 
 TStyle *belleStyle= new TStyle("belleStyle","Phill's  un-official plots style");
 
  // use helvetica-bold-r-normal, precision 2 (rotatable)
  Int_t belleFont = 62;
  // line thickness
  Double_t belleWidth = 1;
 
  // use plain black on white colors
  belleStyle->SetFrameBorderMode(0);
  belleStyle->SetCanvasBorderMode(0);
  belleStyle->SetPadBorderMode(0);
  belleStyle->SetPadColor(0);
  belleStyle->SetCanvasColor(0);
  belleStyle->SetStatColor(0);
  belleStyle->SetPalette(1);
  belleStyle->SetFillColor(0);
 
  // set the paper & margin sizes
  belleStyle->SetPaperSize(20,26);
  belleStyle->SetPadTopMargin(0.05);
  belleStyle->SetPadRightMargin(0.15); // increase for colz plots!!
  belleStyle->SetPadBottomMargin(0.16);
  belleStyle->SetPadLeftMargin(0.14);
 
  // use large fonts
  belleStyle->SetTextFont(belleFont);
  belleStyle->SetTextSize(0.08);
  belleStyle->SetLabelFont(belleFont,"x");
  belleStyle->SetLabelFont(belleFont,"y");
  belleStyle->SetLabelFont(belleFont,"z");
  belleStyle->SetLabelSize(0.05,"x");
  belleStyle->SetLabelSize(0.05,"y");
  belleStyle->SetLabelSize(0.05,"z");
  belleStyle->SetTitleFont(belleFont);
  belleStyle->SetTitleSize(0.06,"x");
  belleStyle->SetTitleSize(0.06,"y");
  belleStyle->SetTitleSize(0.06,"z");
    belleStyle->SetHatchesSpacing(2);
 
  // use bold lines and markers
  belleStyle->SetLineWidth(belleWidth);
  belleStyle->SetFrameLineWidth(belleWidth);
  belleStyle->SetHistLineWidth(belleWidth);
  belleStyle->SetFuncWidth(belleWidth);
  belleStyle->SetGridWidth(belleWidth);
  belleStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
 // belleStyle->SetMarkerStyle(8);
 // belleStyle->SetMarkerSize(1.0);
 
  // label offsets
  belleStyle->SetLabelOffset(0.015);
 
  // by default, do not display histogram decorations:
 // belleStyle->SetOptStat(0);
  //belleStyle->SetOptStat(1110);  // show only nent, mean, rms
  belleStyle->SetOptTitle(0);
  belleStyle->SetOptFit(0);
//  belleStyle->SetOptFit(1011); // show probability, parameters and errors
 
  // look of the statistics box:
  belleStyle->SetStatBorderSize(1);
  belleStyle->SetStatFont(belleFont);
  belleStyle->SetStatFontSize(0.05);
  belleStyle->SetStatX(0.9);
  belleStyle->SetStatY(0.9);
  belleStyle->SetStatW(0.25);
  belleStyle->SetStatH(0.15);
 
  // put tick marks on top and RHS of plots
  belleStyle->SetPadTickX(1);
  belleStyle->SetPadTickY(1);
 
  // histogram divisions: only 5 in x to avoid label overlaps
  belleStyle->SetNdivisions(505,"x");
  belleStyle->SetNdivisions(510,"y");
 
  gROOT->SetStyle("belleStyle");
  gROOT->ForceStyle();
 
  TPaveText *belleName = new TPaveText(0.65,0.8,0.9,0.9,"BRNDC");
  belleName->SetFillColor(0);
  belleName->SetTextAlign(12);
  belleName->SetBorderSize(0);
  belleName->AddText("Belle");
 
  TText *belleLabel = new TText();
  belleLabel->SetTextFont(belleFont);
  belleLabel->SetTextColor(1);
  belleLabel->SetTextSize(0.04);
  belleLabel->SetTextAlign(12);
 
  TLatex *belleLatex = new TLatex();
  belleLatex->SetTextFont(belleFont);
  belleLatex->SetTextColor(1);
  belleLatex->SetTextSize(0.04);
  belleLatex->SetTextAlign(12);
  TGaxis::SetMaxDigits(3);
 
//Add Belle Preliminary, 710 fb-1 to plots
belleLatex->SetTextSize(0.07);
 
 
 
 
 
 
 
    DefineHists();
    Init();
     
    ReadFiles();
    
    
    PerfFit();
   
 
}

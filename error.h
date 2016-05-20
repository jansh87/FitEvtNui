typedef unsigned uint;
#include <string>
#include "../../Dss/BF/FFCorrClass.h"
#include "../shared/EffCorrClass.h"
using namespace std;
class t_error{

	public:
	enum t_errType {_bf,_eff,_FF};
	enum t_ty{_active, _dummy};
	enum t_block {_none,_xtaunu, _xlnu,_Dxlnu,_nBlocks};
	string Name;
	uint lclass;
	uint dclass;
	
	t_errType Type;
	
	t_ty subType;
	
	
	t_block block;
	static uint block_start[_nBlocks];
	static uint block_end[_nBlocks];
	
	virtual double getChange(uint) = 0;
	
	
	
	

};
uint t_error::block_start[] = {0,0,0,0};
uint t_error::block_end[] = {0,0,0,0};

class t_bf: public t_error{

	public:
	
	float weight;
	
	float BF[2];
	
	t_bf(string n, uint lcl, uint dcl, float w, t_block b = _none, t_ty t = _active){
		Name =n;
		lclass = lcl;
		dclass = dcl;
		weight = w;
		Type = _bf;
		block = b;
		BF[0] = findBF(lclass,dclass);
		BF[1] = findBF(-lclass,dclass);
		subType = t;
	};
	virtual double getChange(uint){return weight;};
	
	float findBF(int lcl, int dcl){
	
		uint nBR = CorrectBf::vBBf.size();
		for(int i = 0; i<nBR; i++){
			if(CorrectBf::vBBf[i].lclass == lcl && CorrectBf::vBBf[i].dclass == dcl)
				return CorrectBf::vBBf[i].Br;
		}
		return 0;
	}

};


class t_eff: public t_error{

	public:
	EffCorrClass* corr;
	double *dict;
	
	t_eff(string n, uint ilep, uint fake, EffCorrClass *c,uint nBins = 1, t_block b = _none){
		Name = n;
		corr = c;
		lclass = ilep;
		dclass = fake;
		Type = _eff;
		subType = _active;
		block = b;
		dict = new double[nBins];
		for(int i = 0; i<nBins; i++)
			dict[i] = -10;
	}
	
	virtual double getChange(uint bin){
		return 1.;
		if(dict[bin]>-9){ 
			return dict[bin];
		} else { 
			
			dict[bin] = corr->getRelError(bin);
			return dict[bin];
		}
	};

};

class t_FF: public t_error{

	public:
	FFCorrClass *corr;
	int meson; //for D**
	
	
	t_FF(string n, uint lcl, uint dcl, FFCorrClass *c, uint nBins = 1, int m = -1, t_block b = _none){
		Name = n;
		corr = c;
		lclass = lcl;
		dclass = dcl;
		Type = _FF;
		subType = _active;
		meson = m;
		block = b;
		
	};
	
	virtual double getChange(uint bin){
		return meson<0?corr->reweight(bin):corr->reweight(meson,bin);
	};

};

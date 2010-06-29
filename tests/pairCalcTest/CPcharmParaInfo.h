#ifndef _parainfo_h_
#define _parainfo_h_

class CPcharmParaInfo {
public:

	double tol_norb;
	int cp_min_opt;
	int gen_wave;

	CPcharmParaInfo(CPcharmParaInfo &s){
		tol_norb     = s.tol_norb;
		cp_min_opt   = s.cp_min_opt;
		gen_wave     = s.gen_wave;
	}

	CPcharmParaInfo() {};
	~CPcharmParaInfo() {};

	void pup(PUP::er &p){
		p|tol_norb;
		p|cp_min_opt;
		p|gen_wave;
	}
};
#endif

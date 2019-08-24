#include "TH1F.h"
#include "TComplex.h"
#include "TMath.h"

#define POW2(x)(x*x)
#define POW3(x)(x*x*x)

//Globals

double pi_MASS  = 0.13957018; //GEV
double D_MASS   = 1.96834; //GEV

double d1_MASS  = pi_MASS;  //daughter 1 mass
double d2_MASS  = pi_MASS;
double d3_MASS  = pi_MASS;

double s12_min = (d1_MASS  + d2_MASS)*(d1_MASS + d2_MASS);
double s12_max = (D_MASS   - d2_MASS)*(D_MASS - d2_MASS);

double m12_min = s12_min;
double m12_max = s12_max;

//int slices = 60;

double  lambda (double x, double y, double z){

	double l;
	l = (x - y - z)*(x - y - z) - 4*y*z;

	return l;

}

double Form_Factor_Mother_Decay(int spin, double M, double sab, double mcsq, double mR){
    
    double s = M*M, mRsq = mR*mR;
    double fD, fD0, pstr, pstr0, q2 =0;
    double const rD2 = 25.0;
    double ret;
    
    if (spin == 0){
         ret = 1;
    }
    
    if (spin == 1) {
        pstr0 = sqrt(lambda(s,mRsq,mcsq))/(2*M);
        q2 = rD2*pstr0*pstr0;
        fD0 = sqrt(1 + q2);
        
        pstr = sqrt(lambda(s,sab,mcsq))/(2*M);
        q2 = rD2*pstr*pstr;
        fD = fD0/sqrt(1 + q2);
        ret = fD;
    }
   
    if(spin == 2){
        pstr0 = sqrt(lambda(s,mRsq,mcsq))/(2*M);
        q2 = rD2*pstr0*pstr0;
        fD0 = sqrt(9 + 3*q2 + q2*q2);
        
        pstr = sqrt(lambda(s,sab,mcsq))/(2*M);
        q2 = rD2*pstr*pstr;
        fD = fD0/sqrt(9 + 3*q2 + q2*q2);
        ret = fD;
    }

    return ret;
    
}

double Form_Factor_Resonance_Decay(int spin, double mR, double sab, double masq, double mbsq){

	double mRsq = mR*mR;
    double fR, fR0, pstr, pstr0, q2 = 0;
   
    double const rR2 = 2.25;
    double ret=-1;

	if (spin == 0){
        ret =1;
    }
    
    if (spin == 1) {

		pstr0 = sqrt(lambda(mRsq,masq,mbsq))/(2*mR);
		q2 = rR2*pstr0*pstr0;
		fR0 = sqrt(1 + q2);

		pstr = sqrt(lambda(sab,masq,mbsq))/(2*sqrt(sab));
		q2 = rR2*pstr*pstr;
		fR = fR0/sqrt(1 + q2);

		ret = fR;

	}
    
    if(spin == 2){

		pstr0 = sqrt((mRsq - masq - mbsq)*(mRsq - masq - mbsq) - 4*masq*mbsq)/(2*mR);
		q2 = rR2*pstr0*pstr0;
		fR0 = sqrt(9 + 3*q2 + q2*q2);

		//pstr = sqrt(lambda(sab,masq,mbsq))/(2*sqrt(sab));
		pstr = sqrt((sab - masq + mbsq)*(sab - masq + mbsq) - 4*sab*mbsq)/(2*sqrt(sab));
		q2 = rR2*pstr*pstr;
		fR = fR0/sqrt(9 + 3*q2 + q2*q2);

		ret = fR;
    }
    return ret;

}

double Gamma(int spin, double mR, double width, double mab, double masq, double mbsq){

    double pstr, pstr0,fR, mRsq = mR*mR, sab = mab;
    double ret=-1;

	pstr0 = sqrt(lambda(mRsq,masq,mbsq))/(2*mR);
	pstr = sqrt(lambda(sab,masq,mbsq))/(2*mab);
	if (spin == 0){
        ret =  width*(pstr/pstr0)*(mR/mab);
    }
    
    if (spin == 1){
		fR = Form_Factor_Resonance_Decay(spin, mR, sab, masq, mbsq);
		ret = width*pow((pstr/pstr0),3)*(mR/mab)*fR*fR;
    }
    
    if (spin == 2){
		fR = Form_Factor_Resonance_Decay(spin, mR, sab, masq, mbsq);
		ret = width*pow((pstr/pstr0),5)*(mR/mab)*fR*fR;
    }
    
    return ret;

}


TComplex plainBW(double *x, double *par) {

    double resmass  = par[2];
    double reswidth = par[3];
    int spin = 0;
    double mass2    = x[0];
   
    double FF_MD = Form_Factor_Mother_Decay(spin, D_MASS, mass2, POW2(pi_MASS),resmass);

    double FF_RD = Form_Factor_Resonance_Decay(spin, resmass, mass2,
		       	POW2(pi_MASS),POW2(pi_MASS));
        
	        
    double Width = Gamma(spin, resmass, reswidth, mass2,
            POW2(pi_MASS),
            POW2(pi_MASS)
        );

    //cout << FF_MD << "\t" << FF_RD << '\t' << Width << endl;
        
   
        // RBW evaluation
        double A =  mass2 - POW2(resmass) ;
        double B = resmass * Width ;
        double C = 1.0 / (POW2(A) + POW2(B));
        TComplex _BW(A * C, B * C); 

        _BW *= FF_MD*FF_RD*reswidth*resmass;

    return TComplex(par[0],par[1])*_BW;
}


TComplex flatte(double *x, double *par) {
    
    double resmass            = par[2];
    double g1                 = par[3];
    double g2                 = par[4];

    double pipmass = 0.13957018;
    double pi0mass = 0.1349766;
    double kpmass  = 0.493677;
    double k0mass  = 0.497614;

    double twopimasssq  = 4 * pipmass * pipmass;
    double twopi0masssq = 4 * pi0mass * pi0mass;
    double twokmasssq   = 4 * kpmass * kpmass;
    double twok0masssq  = 4 * k0mass * k0mass;

    TComplex ret(0., 0.);
    
        double rhopipi_real = 0, rhopipi_imag = 0;
        double rhokk_real = 0, rhokk_imag = 0;

        double s = x[0];

        if(s >= twopimasssq)
            rhopipi_real += (2. / 3) * sqrt(1 - twopimasssq / s); // Above pi+pi- threshold
        else
            rhopipi_imag += (2. / 3) * sqrt(-1 + twopimasssq / s);
        if(s >= twopi0masssq)
            rhopipi_real += (1. / 3) * sqrt(1 - twopi0masssq / s); // Above pi0pi0 threshold
        else
            rhopipi_imag += (1. / 3) * sqrt(-1 + twopi0masssq / s);
        if(s >= twokmasssq)
            rhokk_real += 0.5 * sqrt(1 - twokmasssq / s); // Above K+K- threshold
        else
            rhokk_imag += 0.5 * sqrt(-1 + twokmasssq / s);
        if(s >= twok0masssq)
            rhokk_real += 0.5 * sqrt(1 - twok0masssq / s); // Above K0K0 threshold
        else
            rhokk_imag += 0.5 * sqrt(-1 + twok0masssq / s);
            
        double A = (resmass * resmass - s) + resmass * (rhopipi_imag * g1 + rhokk_imag * g2);
        double B = resmass * (rhopipi_real * g1 + rhokk_real * g2);
        double C = 1.0 / (A * A + B * B);
        TComplex retur(A * C, B * C);
        retur *= TComplex(par[0],par[1]);
    
    return retur;
}


//Full SWave
double SWave_amp(double *x, double *par){
    return (  plainBW(x,par)  +  flatte(x,&par[4]) ).Rho2();
}

double SWave_theta(double *x, double *par){
    return (plainBW(x,par)  + flatte(x,&par[4]) ).Theta();
}

void PWACoefs(int slices,double val){

    double s = 0;
    double bin_amp_real = 0;
    double bin_amp_img = 0;

    ofstream wr("files/PWACOEFS.txt");

    double par[13] = {1.0,0.0,1.470,0.35}; 

    int temp = 0;
    int j = 1;

    for(int i = 0; i < slices; i++){


	
		s = m12_min + i*(m12_max-m12_min)/(slices - 1);
		
	


		TComplex v = (plainBW(&s,par));
    		bin_amp_real = v.Rho();
    		bin_amp_img = v.Theta();

    		printf("%lg = (%lg,%lg) \n ",s,bin_amp_real,bin_amp_img);
    		wr << sqrt(s) << " " << bin_amp_real << " "<< bin_amp_img << endl;

	
    }

    wr.close();
   
}

void isobar(){

    TF1 *f1 = new TF1("Amp",SWave_amp,m12_min,m12_max,4);
    f1->SetParameters(1.0,0.0,.480,.350 /* ,2.0,.0,.965,.165,4.21 */); 

    TF1 *f2 = new TF1("Phase",SWave_theta,m12_min,m12_max,4);
    f2->SetParameters(1.0,0.0,.480,.350 /* ,2.0,.0,.965,.165,4.21 */);

    TH1D *h1 = new TH1D("h1","Sigma(480)(BW)",100,m12_min,m12_max);
    h1->GetXaxis()->SetTitle("m^{2}(#pi^{-} #pi^{+}) [Gev]");
    h1->FillRandom("Amp",100000);

    TH1D *h0 = new TH1D("h0","Sigma(480)(BW)",100,m12_min,m12_max);
    h0->SetLineColor(kRed);
    TTree *t = new TTree("t","");
    t->ReadFile("D2PPP_toy_iso.txt","x:y:z");
    t->Draw("y>>h0");

    //TH2D *h2 = new TH2D("h2","",100,s12_min,s12_max,100,0,25000);
    //h2->SetMarkerStyle(2);
    TTree *t2 = new TTree("t2","");
    t2->ReadFile("files/PWACOEFS.txt","x:y:z");
    
    t2->Draw("(y*y + z*z):x>>h2","","*");
    //t2->Draw("TMath::ATan2(z,y):x>>h2","","*");
    //h1->DrawNormalized("L");
    //h1->Draw("L");
    //h0->DrawNormalized("Lsame"); 
    //h0->Divide(h1);
    //h0->Fit("Amp","RLL");
    f1->Draw("same");
    //f2->Draw("same");
}


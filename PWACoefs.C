#include "TH1F.h"
#include "TComplex.h"
#include "TMath.h"

#define POW2(x)(x*x)
#define POW3(x)(x*x*x)
#define torad(x)(x*M_PI/180.)
#define nKnobs 35
//Globals

double pi_MASS  = 0.13957018; //GEV
double D_MASS   = 1.96834; //GEV

double d1_MASS  = pi_MASS;  //daughter 1 mass
double d2_MASS  = pi_MASS;
double d3_MASS  = pi_MASS;

double s12_min = (d1_MASS  + d2_MASS)*(d1_MASS + d2_MASS);
double s12_max = (D_MASS   - d2_MASS)*(D_MASS - d2_MASS);

double m12_min = sqrt(s12_min);
double m12_max = sqrt(s12_max);

double *mKKlimits = nullptr;
double *real_coefs = nullptr;
double *imag_coefs = nullptr;

double Gamma(const double &m,
                        const double &m0,
                        const double &width, 
                        const double &q,    
                        const double &q0,
                        const double &BWFactor,
                        const unsigned int &spin){

                        double g = 1.;

                        if(spin==0){
                            g*= width*(q/q0)*(m0/m)*POW2(BWFactor);
                        }

                        if(spin==1){
                            g*= width*pow(q/q0 , 2*1 + 1)*(m0/m)*POW2(BWFactor);
                        }

                        if(spin==2){
                            g*= width*pow(q/q0 , 2*2 + 1)*(m0/m)*POW2(BWFactor);
                        }

                        return g;
}

double Momentum( const double &m,
                            const double &m1,
                            const double &m2
                            ) {
  
    double k1 = m*m - POW2(m1+m2);
    double k2 = m*m - POW2(m1-m2);
    double q = 0.5*sqrt(k1*k2)/m;

    return  k1*k2>0 ? q :0 ;
}

TComplex plainBW(double *x, double *par) {

    double resmass  = par[2];
    double reswidth = par[3];
    int spin = 0;
    double mass2    = x[0];
    mass2*=mass2;
    double m = sqrt(mass2);

    double m1 = 0.13957018;

    double q  = Momentum(m,m1,m1);
    double q0 = Momentum(resmass,m1,m1);
	        
    double Width = Gamma(m,resmass,reswidth,q,q0,1.,spin) ;
    //cout << resmass << endl;

    // RBW evaluation
    double A =  POW2(resmass) - mass2;
    double B = resmass * Width ;
    double C = 1.0 / (POW2(A) + POW2(B));
    TComplex _BW(A * C, B * C); 

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

    double resmass2 = POW2(resmass);
    double rho1(0.0), rho2(0.0);
   
        double s = x[0];
        s*=s;
        double dMSq = resmass2 - s;

        if (s > twopi0masssq) {
            rho1 = sqrt(1.0 - twopi0masssq/s)/3.0;
            if (s > twopimasssq) {
                rho1 += 2.0*sqrt(1.0 - twopimasssq/s)/3.0;
                if (s > twokmasssq) {
                    rho2 = 0.5*sqrt(1.0 - twokmasssq/s);
                    if (s > twok0masssq) {
                        rho2 += 0.5*sqrt(1.0 - twok0masssq/s);
                    } else {
                        // Continue analytically below higher channel thresholds
                        // This contributes to the real part of the amplitude denominator
                        dMSq += g2*resmass*0.5*sqrt(twok0masssq/s - 1.0);
                    }
                } else {
                    // Continue analytically below higher channel thresholds
                    // This contributes to the real part of the amplitude denominator
                    rho2 = 0.0;
                    dMSq += g2*resmass*(0.5*sqrt(twokmasssq/s - 1.0) + 0.5*sqrt(twok0masssq/s - 1.0));
                }
            } else {
                // Continue analytically below higher channel thresholds
                // This contributes to the real part of the amplitude denominator
                dMSq += g1*resmass*2.0*sqrt(twopimasssq/s - 1.0)/3.0;
            }
        }
    
        //the Adler-zero term fA = (m2 − sA)/(m20 − sA) can be used to suppress false 
        //kinematic singularities when m goes below threshold. For f(0)(980), sA = 0.
        
        double massFactor = s/resmass2;
        
        double width1 = g1*rho1*massFactor;
        double width2 = g2*rho2*massFactor;
        double widthTerm = width1 + width2;
    
        TComplex resAmplitude(dMSq, widthTerm);
    
        double denomFactor = dMSq*dMSq + widthTerm*widthTerm;
    
        double invDenomFactor = 0.0;
        if (denomFactor > 1e-10) {invDenomFactor = 1.0/denomFactor;}
    
        resAmplitude *= invDenomFactor;

        return resAmplitude*TComplex(par[0],par[1]);
}

TComplex MIPWA(double *x, double *par){

    double s = x[0];
    int khiAB=0,kloAB = 0;
    while(khiAB < nKnobs) {
        if(s < mKKlimits[khiAB])
            break;
        khiAB++;
    }

    kloAB                = khiAB - 1;

    double pwa_coefs_real_kloAB = real_coefs[kloAB];
    double pwa_coefs_real_khiAB = real_coefs[khiAB];
    double pwa_coefs_imag_kloAB = imag_coefs[kloAB];
    double pwa_coefs_imag_khiAB = imag_coefs[khiAB];

    double dmKK = mKKlimits[khiAB] - mKKlimits[kloAB];
    double aa   = (mKKlimits[khiAB] - s) / dmKK;
    double bb   = 1 - aa;

    TComplex ret(aa * pwa_coefs_real_kloAB + bb * pwa_coefs_real_khiAB,
                aa * pwa_coefs_imag_kloAB + bb * pwa_coefs_imag_khiAB);

    return ret;
}

//Full SWave
double SWave_amp(double *s, double *par){
    return ( flatte(s,par) + plainBW(s,&par[5]) + plainBW(s,&par[9])  ).Re();
}

double SWave_amp_MIPWA(double *x, double *par){
    return ( MIPWA(x,par) + flatte(x,par)).Rho();
}

double SWave_theta_MIPWA(double *x, double *par){
    return ( MIPWA(x,par) + flatte(x,par)).Theta();
}

double SWave_theta(double *s, double *par){
  return ( flatte(s,par) + plainBW(s,&par[5]) + plainBW(s,&par[9])  ).Im();
}

void pts(int slices,double mass,double err,double sigma){

    double s = 0;
    double bin_amp_real = 0;
    double bin_amp_img = 0;
    //printf("m12_min = %f e m12_max = %f \n",m12_min,m12_max);

    ofstream wr("files/PWACOEFS.txt");
    double par[8] = {1.,0.,0.977,0.1,0.75*cos(torad(198)),0.75*sin(torad(198)),1.434,0.172};
    //double par[13] = {1.,0.,0.965,0.165,0.693,-0.0924,-0.1413,1.505,0.109,-0.1641,-0.1018,1.434,0.172}; 
    //double par[8] = {1.,0.,0.977,0.04,0.75*cos(torad(198)),0.75*sin(torad(198)),1.434,0.172};
    //double par[9] = {-1.884286,-4.594298,0.965,0.165,0.694,-0.026994,1.889315,1.434,0.172};

    int temp = 0;
    int j = 0;
    int k = 0;
    cout << s12_max << '\n';
    for(int i = 0; i < slices; i+=1){

        s = m12_min + (i+0.5)*(m12_max-m12_min)/(slices);
        
        /* if(s<=(mass- sigma*err)){
                s = m12_min + (i+0.5)*( (mass - sigma*err) - m12_min)/(0.3*slices);
                
        }
        

        if( s>=(mass-sigma*err) && s<=(mass+sigma*err) ){
                s = (mass-sigma*err) + (j+0.5)*(2*sigma*err)/(0.4*slices);
                j++;
        }
        
        
        if(s>=(mass+sigma*err)){
                s = (mass+sigma*err) + (k+0.5)*(m12_max - mass+sigma*err)/(0.3*slices +2);
                k++;
        } */

		//TComplex v = (flatte(&s,par) + plainBW(&s,&par[5]) + plainBW(&s,&par[9]));	
		TComplex v = (plainBW(&s,par)+ plainBW(&s,&par[4]));
		//TComplex v = (flatte(&s,par)+ plainBW(&s,&par[5]));
		bin_amp_real = v.Re();
    	bin_amp_img = v.Im();

    		//printf("%d : %lg = (%lg,%lg) \n ",i,s,bin_amp_real,bin_amp_img);
    	wr << s << " " << bin_amp_real << " "<< bin_amp_img << endl;

	
    }

    wr.close();
   
}

void babarBinning(){
    ifstream rd("files/PWACOEFS_BaBar_more_pts.txt");
    ofstream wt("files/PWACOEFS.txt");

    double par[9] = {1.,0.,0.963,0.084,0.23,-0.4623,-0.6269,1.472,0.152};  
    double s , m , p;

    while(rd >> s >> m >> p){
        
        TComplex v = (flatte(&s,par) + plainBW(&s,&par[5]));	
		auto bin_amp_real = v.Rho();
    		auto bin_amp_img = v.Theta();
    	wt << s << " " << bin_amp_real << " "<< bin_amp_img << endl;
    }
    
    rd.close();
    wt.close();

}

void PWACoefs(int slices,double mass,double err,double sigma){
    //pts(slices,mass,err,sigma);
    babarBinning();
    TCanvas *foo = new TCanvas("foo","",900,600);

    TGraphErrors *a = new TGraphErrors("Fit/s_wave_fitted.txt","%lg %lg %*lg %lg %lg %*lg");
    TGraphErrors *b = new TGraphErrors("Fit/s_wave_fitted.txt","%lg %*lg %lg %lg %*lg %lg");

    TGraph *c = new TGraph("files/PWACOEFS.txt","%lg %lg %*lg");
    TGraph *d = new TGraph("files/PWACOEFS.txt","%lg %*lg %lg");

    mKKlimits = a->GetX();
    real_coefs = a->GetY();
    imag_coefs = b->GetY();

    TF1 *f1 = new TF1("AmpMIPWA",SWave_amp_MIPWA,m12_min,m12_max,5);
    TF1 *f2 = new TF1("thetaMIPWA",SWave_theta_MIPWA,m12_min,m12_max,5);
    f1->SetParameters(1.,0.,0.963,0.084,0.23);
    f2->SetParameters(1.,0.,0.963,0.084,0.23);
    
    //input onda-S isobárica
    /* 
    TGraphErrors *a = new TGraphErrors("Fit/s_wave_fitted.txt","%lg %lg %*lg %lg %lg %*lg");
    TGraphErrors *b = new TGraphErrors("Fit/s_wave_fitted.txt","%lg %*lg %lg %lg %*lg %lg");
    

    TF1 *f1 = new TF1("AmpMIPWA",SWave_amp,m12_min,m12_max,13);
    TF1 *f2 = new TF1("thetaMIPWA",SWave_theta,m12_min,m12_max,13);

    f1->SetParameters(1.,0.,0.965,0.165,0.680,-0.0924,-0.1413,1.505,0.109,-0.1641,-0.1018); //50pts n=4
    f1->SetParameter(11,1.434);
    f1->SetParameter(12,0.172);
    f2->SetParameters(1.,0.,0.965,0.165,0.680,-0.0924,-0.1413,1.505,0.109,-0.1641,-0.1018);
    f2->SetParameter(11,1.434);
    f2->SetParameter(12,0.172); */
    

    foo->Divide(2);
    foo->cd(1);
    c->Draw("AP*");
    f1->Draw("same");
    

    foo->cd(2);
    d->Draw("AP*");
    f2->Draw("same");

}


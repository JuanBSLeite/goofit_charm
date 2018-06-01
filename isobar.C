#include "TH1F.h"
#include "TComplex.h"
#include "TMath.h"

//Globals

double pi_MASS  = 0.13957018; //GEV
double D_MASS   = 1.86959; //GEV

double d1_MASS  = pi_MASS;  //daughter 1 mass
double d2_MASS  = pi_MASS;
double d3_MASS  = pi_MASS;

double s12_min = pow(d1_MASS  + d2_MASS,2);
double s12_max = pow(D_MASS   - d2_MASS,2);

int slices = 30;

double twoBodyCMmom(double rMassSq, double d1m, double d2m) {
    // For A -> B + C, calculate momentum of B and C in rest frame of A.
    // PDG 38.16.

    double kin1 = 1.0 - pow(d1m + d2m,2) / rMassSq;

    kin1 = kin1 >= 0.0 ? sqrt(kin1) : 1.0;

    double kin2 = 1.0 - pow(d1m - d2m,2) / rMassSq;
    kin2        = kin2 >= 0.0 ? sqrt(kin2) : 1.0;

    return 0.5 * sqrt(rMassSq) * kin1 * kin2;
}
double dampingFactorSquare(const double &cmmom, const int &spin, const double &mRadius) {
    double square = mRadius * mRadius * cmmom * cmmom;
    double dfsq   = 1 + square; // This accounts for spin 1
    // if (2 == spin) dfsq += 8 + 2*square + square*square; // Coefficients are 9, 3, 1.
    double dfsqres = dfsq + 8 + 2 * square + square * square;

    // Spin 3 and up not accounted for.
    // return dfsq;
    return (spin == 2) ? dfsqres : dfsq;
}

TComplex plainBW(double *x, double *par) {

    unsigned int spin = 0;
    double resmass  = par[2];
    double reswidth = par[3];

    double resmass2 = pow(resmass,2);
    
        double rMassSq    = x[0];
        double mass_daug1 = pi_MASS;
        double mass_daug2 = pi_MASS;

        double frFactor = 1;

        // Calculate momentum of the two daughters in the resonance rest frame
        // Note symmetry under interchange (dm1 <-> dm2)

        double measureDaughterMoms = twoBodyCMmom(rMassSq, mass_daug1, mass_daug2);
        double nominalDaughterMoms = twoBodyCMmom(resmass2, mass_daug1, mass_daug2);

        if(0 != spin) {
        frFactor = dampingFactorSquare(nominalDaughterMoms, spin, 1.5);
        frFactor /= dampingFactorSquare(measureDaughterMoms, spin, 1.5);
    }

        // RBW evaluation
        double A = (resmass2 - rMassSq);
        double B = resmass2 * reswidth   /* * pow(measureDaughterMoms / nominalDaughterMoms, 2.0 * spin + 1)  */ * frFactor/ sqrt(rMassSq);
        double C = 1.0 / (pow(A,2) + pow(B,2));

        TComplex ret(A * C, B * C,0.0);
        ret *= sqrt(frFactor);
        ret *= TComplex(par[0],par[1],0);

    return ret;
}

TComplex flatte(double *x, double *par) {
    // indices[1] is unused constant index, for consistency with other function types.
    double resmass            = par[2];
    double g1                 = par[3];
    double g2                 = par[4] * g1;

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
        ret = retur*TComplex(par[0],par[1],0);
    
    return ret;
}

//Full SWave
double SWave_amp(double *x, double *par){
    return (plainBW(x,par) + flatte(x,&par[4])).Rho2();
}

double SWave_theta(double *x, double *par){
    return (plainBW(x,par) + flatte(x,&par[4])).Theta();
}

void PWACoefs(){

    TF1 *amp = new TF1("amp",SWave_amp,s12_min,s12_max,9);
    amp->SetParameters(1.0,0.0,.480,.350,2.0,0.0,0.965,0.165,4.21);
    TF1 *phase = new TF1("phase",SWave_theta,s12_min,s12_max,9);
    phase->SetParameters(1.0,0.0,.480,.350,2.0,0.0,0.965,0.165,4.21);
    
    double rho = 0;
    double theta =0;
    double s = 0;

    ofstream wr("files/PWACOEFS.txt");

    for(int i = 0 ; i < slices ; i++){

    s = s12_min + i*(s12_max-s12_min)/(slices-1);
    rho = amp->Eval(s);
    theta = phase->Eval(s);

    printf("%lg = (%lg,%lg) \n ",s,rho,theta);
    wr << sqrt(s) << " " << rho << " "<< theta << endl;
    }

    wr.close();
   
}

void isobar(){

    TF1 *f1 = new TF1("Theta",SWave_theta,s12_min,s12_max,9);
    f1->SetParameters(1.0,0.0,.480,.350,2.0,0.0,0.965,0.165,4.21);

    TF1 *f3 = new TF1("Rho2",SWave_amp,s12_min,s12_max,9);
    //f3->SetParameters(1.0,0.0,.480,.350,2.0,0.0,0.965,0.165,4.21);

    //free parameters
    f3->SetParameter(0,1.0);
    f3->SetParLimits(0,-100.0,+100.0);
    f3->SetParError(0,0.01);
    
    f3->SetParameter(1,0.0);
    f3->SetParLimits(1,-100.0,+100.0);
    f3->SetParError(1,0.01);

    f3->SetParameter(4,2.0);
    f3->SetParLimits(4,-100.0,+100.0);
    f3->SetParError(4,0.01);

    f3->SetParameter(5,0.0);
    f3->SetParLimits(5,-100.0,+100.0);
    f3->SetParError(5,0.01);  

    //fixed parameters
    f3->FixParameter(2,.480);
    f3->FixParameter(3,.350);
    f3->FixParameter(6,0.965);
    f3->FixParameter(7,.165);
    f3->FixParameter(8,4.21); 

    TTree *t = new TTree("Tree","Tree");
    t->ReadFile("D2PPP_toy.txt","x:y:z");
    TH1D *h1 = new TH1D("h1","Sigma(480)(BW) + f0(980)(Flatte)",100,s12_min,s12_max);
    //h1->FillRandom("Rho2",100000);
    t->Draw("y>>h1");
    h1->GetXaxis()->SetTitle("m^{2}(#pi^{-} #pi^{+}) [Gev]");
    h1->Scale(1.0/ 91.7737,"");
    h1->Fit(f3,"RL");

    //TCanvas *c = new TCanvas("c","amp",800,500);
    //c->Divide(2,1);
    //c->cd(1);
    //f3->Draw();
    //c.SaveAs("plots/Amp.png");
    //c->Close();
    //c->cd(2);
    //TCanvas c1("c1","theta",700,600);
    //f1->Draw();
    //c->SaveAs("plots/Amp&Theta.png");
    //c1->Close();


}
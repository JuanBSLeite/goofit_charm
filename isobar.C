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

double twoBodyCMmom(double rMassSq, double d1m, double d2m) {
    // For A -> B + C, calculate momentum of B and C in rest frame of A.
    // PDG 38.16.

    double kin1 = 1.0 - pow(d1m + d2m,2) / rMassSq;

    kin1 = kin1 >= 0.0 ? sqrt(kin1) : 1.0;

    double kin2 = 1.0 - pow(d1m - d2m,2) / rMassSq;
    kin2        = kin2 >= 0.0 ? sqrt(kin2) : 1.0;

    return 0.5 * sqrt(rMassSq) * kin1 * kin2;
}


double plainBW(double *x, double *par) {

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

        // RBW evaluation
        double A = (resmass2 - rMassSq);
        double B = resmass2 * reswidth * pow(measureDaughterMoms / nominalDaughterMoms, 2.0 * spin + 1) * frFactor
                   / sqrt(rMassSq);
        double C = 1.0 / (pow(A,2) + pow(B,2));

        TComplex ret(A * C, B * C,0.0);

        TComplex c(par[0],par[1],1.0);

        ret *= c;

    return ret.Rho();
}

double flatte(double *x, double *par) {
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
        ret = retur*TComplex(par[0],par[1],1);
    
    return ret.Rho();
}


double SWave(double *x, double *par){
    return (plainBW(x,par) + flatte(x,&par[4]));
}

void isobar(){

    TF1 *f1 = new TF1("f1",plainBW,s12_min,s12_max,4);
    f1->SetParameters(1.0,0.0,1.2755,0.1867);

    TF1 *f2 = new TF1("f2",flatte,s12_min,s12_max,4);
    f2->SetParameters(2.0,0.0,0.965,0.165,4.21);

    TF1 *f3 = new TF1("f3",SWave,s12_min,s12_max,9);
    f3->SetParameters(1.0,0.0,0.480,0.350,2.0,0.0,0.965,0.165,4.21);

    f3->Draw();
}
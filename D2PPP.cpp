// ROOT stuff
#include <TApplication.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TMath.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TTree.h>
#include <TROOT.h>
#include <TMinuit.h>
#include <TNtuple.h>
#include <TComplex.h>
#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>
#include <TH2Poly.h>

//Minuit

#include <Minuit2/MnStrategy.h>

// System stuff
#include <CLI/Timer.hpp>
#include <fstream>
#include <string>
#include <time.h>

// GooFit stuff
#include <goofit/Application.h>
#include <goofit/BinnedDataSet.h>
#include <goofit/FitManager.h>
#include <goofit/fitting/FitManagerMinuit2.h>
#include <goofit/fitting/FitManagerMinuit1.h>
#include <goofit/PDFs/GooPdf.h>
#include <goofit/PDFs/basic/PolynomialPdf.h>
#include <goofit/PDFs/basic/SmoothHistogramPdf.h>
#include <goofit/PDFs/combine/AddPdf.h>
#include <goofit/PDFs/combine/ProdPdf.h>
#include <goofit/PDFs/physics/DalitzPlotPdf.h>
#include <goofit/PDFs/physics/DalitzVetoPdf.h>
#include <goofit/PDFs/physics/ResonancePdf.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/Variable.h>
#include <goofit/PDFs/physics/DalitzPlotter.h>
#include <goofit/PDFs/combine/EventWeightedAddPdf.h>
#include <goofit/PDFs/combine/CompositePdf.h>
#include <goofit/PDFs/physics/IncoherentSumPdf.h>

#include <thrust/transform_reduce.h>

using namespace std;
using namespace GooFit;
using namespace ROOT;


//Globals

double pi_MASS  = 0.13957018; //GEV
double D_MASS   = 1.96834; //GEV

double d1_MASS  = pi_MASS;  //daughter 1 mass
double d2_MASS  = pi_MASS;
double d3_MASS  = pi_MASS;

bool saveBkgPlot= false;
bool saveEffPlot= false;
bool doEffSwap  = false;
bool doBkgSwap  = false;
bool toyOn      = true;
bool bkgOn      = false;
bool effOn      = false;
bool bdt = false;

const double NevG = 1e7; 
const int bins = 500;

fptype s12_min = POW2(d1_MASS  + d2_MASS);
fptype s12_max = POW2(D_MASS   - d2_MASS);
fptype s13_min = POW2(d1_MASS  + d3_MASS);
fptype s13_max = POW2(D_MASS   - d2_MASS);

Observable s12("s12",s12_min,s12_max); //s12^{2}
Observable s13("s13",s12_min,s12_max);
Observable s23("s23",s12_min,s12_max);
EventNumber eventNumber("eventNumber");

Variable Mother_Mass("Mother_Mass",D_MASS);
Variable Daughter1_Mass("Daughter1_Mass",d1_MASS);
Variable Daughter2_Mass("Daughter2_Mass",d2_MASS);
Variable Daughter3_Mass("Daughter3_Mass",d3_MASS);

UnbinnedDataSet* Data = nullptr;

vector<fptype> HH_bin_limits;
vector<Variable> pwa_coefs_amp;
vector<Variable> pwa_coefs_phs;

TH2F *effHistogram    = nullptr;
TH2F *bkgHistogram       = nullptr;
TH2F *underlyingBins     = nullptr;

// PWA INPUT FILE NAME
const string pwa_file = "files/PWACOEFS.txt";

// Data File
const string data_name = "../../../dados/DsPPP_92.root";
const string tree_name = "DecayTree";
const string bkghisto_file = "../../../dados/bkg_histo_16.root";
const string effhisto_file = "../../../dados/eff_16.root";

//Declaring Functions 
void saveParameters(std::string file, const std::vector<Variable> &param,   bool isValid,  double fcn,  std::vector<std::vector<fptype>> ff );
void saveParameters(std::string file, const std::vector<Variable> &param);
void runtoyfit(std::string name = "D2PPP_toy.txt", int sample_number = 0);


void style(){
    TStyle *myStyle= new TStyle( "myStyle", "Josue (LHCb) official plots style" );
	Double_t lhcbWidth = 3;
	myStyle->SetPadColor(0);
	myStyle->SetCanvasColor(0);
	myStyle->SetStatColor(0); 
	myStyle->SetLineWidth(lhcbWidth);
	myStyle->SetFrameLineWidth(lhcbWidth);
	myStyle->SetHistLineWidth(lhcbWidth);
	myStyle->SetFuncWidth(lhcbWidth);
	myStyle->SetGridWidth(lhcbWidth);
	myStyle->SetMarkerStyle(8);
	myStyle->SetMarkerSize(1.5);
	myStyle->SetPadTickX(1);            
	myStyle->SetPadTickY(1);   
	myStyle->SetOptStat(0); 
	myStyle->SetOptFit(1111111);
	const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
	gROOT->SetStyle("myStyle");
}

void maketoydalitzdata(GooPdf* overallsignal,DalitzPlotPdf *signaldalitz, std::string name, size_t nEvents){
    DalitzPlotter dp(overallsignal,signaldalitz);
    Data = new UnbinnedDataSet({s12,s13,eventNumber});

    std::cout << "Toy Generation begin!" << '\n';
    dp.fillDataSetMC(*Data,nEvents);
    {
        ofstream w(name);

            for (size_t i = 0; i < Data->getNumEvents(); i++) {
                Data->loadEvent(i);
                w << i << "\t" << std::setprecision(6) << s12.getValue() << "\t" << s13.getValue() << '\n';

            }

            std::cout << "nEvents generated = " << Data->getNumEvents() << '\n';

        w.close();

    }

   
    std::cout << name << " Generation end!" << '\n';

}

ResonancePdf *loadPWAResonance(const string fname = pwa_file, bool fixAmp = false) {

    std::ifstream reader;
	//GOOFIT_INFO("LOADING FILE {}",fname);
    reader.open(fname.c_str());
    assert(reader.good());
    HH_bin_limits.clear();
    pwa_coefs_amp.clear();
    pwa_coefs_phs.clear();

    double e1, e2, e3;
    double emag, ephs;
    int i = 0;
    while(reader >> e1 >> e2 >> e3) {

        HH_bin_limits.push_back(e1*e1);

        emag = e2;//e2*cos(e3);
        ephs = e3;//e2*sin(e3);
        Variable va(fmt::format("pwa_coef_{}_real", i), emag, .0001, -20.,+20.);
        Variable vp(fmt::format("pwa_coef_{}_img", i), ephs, .0001,-20.,+20.);

	if((i<13)||(i>22)){ //aumentar
		va.setFixed(true);
		vp.setFixed(true);
	}
        pwa_coefs_amp.push_back(va);
        pwa_coefs_phs.push_back(vp);
        i++;

    }

    Variable swave_amp_real("swave_amp_real", 1.0,0.001,0,0);
    Variable swave_amp_imag("swave_amp_imag", 0.0,0.001,0,0);

    if(fixAmp) {
        swave_amp_real.setValue(1.);
        swave_amp_imag.setValue(0.);
        swave_amp_real.setFixed(true);
        swave_amp_imag.setFixed(true);
    }
    cout << "Numbers loaded: " << HH_bin_limits.size() << " / " << i << endl;

    ResonancePdf *swave_12 = new Resonances::Spline("swave_12", swave_amp_real, swave_amp_imag, HH_bin_limits, pwa_coefs_amp, pwa_coefs_phs, PAIR_12, true);

    return swave_12;
} 

GooPdf* makeEfficiencyPdf() {

    vector<Observable> lvars;
    lvars.push_back(s12);
    lvars.push_back(s13);
    BinnedDataSet *binEffData = new BinnedDataSet(lvars);
    
    TFile *f     = TFile::Open(effhisto_file.c_str());
    auto effHistogram = (TH2F *)f->Get("h_eff");
    
    int evtCounter = 0;

    for(int i = 0; i < s12.getNumBins(); ++i) {
        s12.setValue(s12.getLowerLimit() + (s12.getUpperLimit() - s12.getLowerLimit()) * (i + 0.5) / s12.getNumBins());
        for(int j = 0; j < s13.getNumBins(); ++j) {
            s13.setValue(s13.getLowerLimit() + (s13.getUpperLimit() - s13.getLowerLimit()) * (j + 0.5) / s13.getNumBins());
            if(!inDalitz(s12.getValue(), s13.getValue(), D_MASS, d1_MASS, d2_MASS, d3_MASS)){continue;}
            double weight = effHistogram->GetBinContent(effHistogram->FindBin(s12.getValue(), s13.getValue()));
            binEffData->addWeightedEvent(weight);

        }
    }


    // Smooth
    Variable *effSmoothing = new Variable("effSmoothing_eff", 0.0);
    SmoothHistogramPdf *ret = new SmoothHistogramPdf("efficiency", binEffData, *effSmoothing);
    return ret;
}

GooPdf* makeBackgroundPdf() {

   BinnedDataSet *binBkgData = new BinnedDataSet({s12, s13});
    
    TFile *f = new TFile(bkghisto_file.c_str());
    auto bkgHistogram = (TH2F*)f->Get("h_eff");
   
    int evtCounter = 0;

    for(int i = 0; i < s12.getNumBins(); ++i) {
        s12.setValue(s12.getLowerLimit() + (s12.getUpperLimit() - s12.getLowerLimit()) * (i + 0.5) / s12.getNumBins());
        for(int j = 0; j < s13.getNumBins(); ++j) {
            s13.setValue(s13.getLowerLimit() + (s13.getUpperLimit() - s13.getLowerLimit()) * (j + 0.5) / s13.getNumBins());
            if(!inDalitz(s12.getValue(), s13.getValue(), D_MASS, d1_MASS, d2_MASS, d3_MASS)){continue;}
            double weight = bkgHistogram->GetBinContent(bkgHistogram->FindBin(s12.getValue(), s13.getValue()));
            binBkgData->addWeightedEvent(weight);

        }
    }

    Variable *effSmoothing  = new Variable("effSmoothing_bkg",0.);
    SmoothHistogramPdf *ret = new SmoothHistogramPdf("background", binBkgData, *effSmoothing);

    return ret;
}
GooPdf *makeDstar_veto() {
   
    VetoInfo Dstar_veto12(Variable("Dstar_veto12_min", 2.85), Variable("Dstar_veto12_max", s12_max), PAIR_12);
    VetoInfo Dstar_veto13(Variable("Dstar_veto13_min", 2.85), Variable("Dstar_veto13_max", s12_max), PAIR_13);
    VetoInfo Fiducial_veto12(Variable("Fiducial_veto12_min", 0.), Variable("Fiducial_veto12_max", s12_min), PAIR_12);
    VetoInfo Fiducial_veto13(Variable("DFiducial_veto13_min", 0.), Variable("Fiducial_veto13_max", s12_min), PAIR_13);
    
    vector<VetoInfo> vetos;
    //vetos.push_back(Dstar_veto12);
    //vetos.push_back(Dstar_veto13);
    vetos.push_back(Fiducial_veto12);
    vetos.push_back(Fiducial_veto13);
    //vetos.push_back(kVetoInfo23);


    DalitzVetoPdf* Dstar_veto = new DalitzVetoPdf("Dstar_veto", s12, s13, Mother_Mass, Daughter1_Mass, Daughter2_Mass, Daughter3_Mass, vetos);
    return Dstar_veto;
}

DalitzPlotPdf* makesignalpdf(GooPdf* eff){

    DecayInfo3 dtoppp;
    dtoppp.motherMass   = D_MASS;
    dtoppp.daug1Mass    = d1_MASS;
    dtoppp.daug2Mass    = d2_MASS;
    dtoppp.daug3Mass    = d3_MASS;
    dtoppp.meson_radius = 1.5;

    //Mass and width

    bool fixed = true;

    TRandom3 donram(time(NULL));
    
    
    double f0_980_MASS     = 0.958;
    double f0_980_GPP     = 0.1017;
    double f0_980_GKK     = 3.31;
    double f0_980_amp     = 1.0;
    double f0_980_img    = 0.0;

    double f0_1370_MASS  = 1.468;
    double f0_1370_WIDTH = .175;
    double f0_1370_amp   = fixed ? -0.6 : 0.3*donram.Uniform(-1.,1.);
    double f0_1370_img = fixed ? -0.47 : 0.8*donram.Uniform(-1.,1.);


    double f0_1500_MASS  = 1.505;
    double f0_1500_WIDTH = .109;
    double f0_1500_amp   = fixed ? 0.3 : 0.3*donram.Uniform(-1.,1.);
    double f0_1500_img = fixed ? 0.8 : 0.8*donram.Uniform(-1.,1.);

    double f0_X_MASS  = 1.430;
    double f0_X_WIDTH = 0.320;
    double f0_X_amp   = fixed ? 0.5 : 0.5*donram.Uniform(-1.,1.);
    double f0_X_img = fixed ? 0.5 : 0.5*donram.Uniform(-1.,1.);

    double omega_MASS   = 0.78265;
    double omega_WIDTH  = 0.00849;
    double omega_amp    = fixed ? 0.8 : 0.8*donram.Uniform(-1.,1.);
    double omega_img  = fixed ? 0.8 : 0.8*donram.Uniform(-1.,1.);

    double rho770_MASS   = .77549;
    double rho770_WIDTH  = .1491;
    //E791
    //double rho770_amp    = 0.32*cos(torad(109));
    //double rho770_img  = 0.32*sin(torad(109));
    //Babar
    double rho770_amp    = -0.02;
    double rho770_img  = -0.07;
    
    
    double rho1450_MASS   = 1.465;
    double rho1450_WIDTH  = 0.4;
    //E791
    //double rho1450_amp    = 0.28*cos(torad(162));
    //double rho1450_img  = 0.28*sin(torad(162));
    //Babar
    double rho1450_amp    = 0.16;
    double rho1450_img  = -0.14;
    
    double f2_1270_MASS     = 1.2751;
    double f2_1270_WIDTH    = 0.1851;
    double f2_1270_amp      = fixed ? -0.18 : 0.3*donram.Uniform(-1.,1.);
    double f2_1270_img    = fixed ? 0.28 : 0.8*donram.Uniform(-1.,1.);

    double f2_1525_MASS     = 1.525;
    double f2_1525_WIDTH    = 0.073;
    double f2_1525_amp      = fixed ? 0.6 : 0.6*donram.Uniform(-1.,1.);
    double f2_1525_img    = fixed ? 0.3 : 0.3*donram.Uniform(-1.,1.);
    

    //omega(782)
    Variable v_omega_Mass("omega_MASS",omega_MASS,omega_MASS*0.5,omega_MASS*1.5);
    Variable v_omega_Width("omega_WIDTH",omega_WIDTH,omega_WIDTH*0.5,omega_WIDTH*1.5);
    Variable v_omega_real("omega_real",omega_amp, 0.001, -100.0, +100.0);
    Variable v_omega_img("omega_img",omega_img, 0.001, -100.0, +100.0);

//rho(770)
    Variable v_rho770_Mass("rho770_MASS",rho770_MASS,0.01,0,0);
    Variable v_rho770_Width("rho770_WIDTH",rho770_WIDTH,0.01,0,0);
    Variable v_rho770_real("rho770_real",rho770_amp, 0.01,0,0);
    Variable v_rho770_img("rho770_img",rho770_img, 0.01,0,0);
    
    //rho(1450)
    Variable v_rho1450_Mass("rho1450_MASS",rho1450_MASS,0.01,0,0);
    Variable v_rho1450_Width("rho1450_WIDTH",rho1450_WIDTH,0.01,0,0);
    Variable v_rho1450_real("rho1450_real",rho1450_amp, 0.01,0,0);
    Variable v_rho1450_img("rho1450_img",rho1450_img, 0.01,0,0);

	//f2(1270)
    Variable v_f2_1270_Mass("f2_1270_MASS",f2_1270_MASS,0.01,0,0);
    Variable v_f2_1270_Width("f2_1270_WIDTH",f2_1270_WIDTH,0.01,0,0);
    Variable v_f2_1270_real("f2_1270_real",f2_1270_amp, 0.01,0,0);
    Variable v_f2_1270_img("f2_1270_img",f2_1270_img, 0.01,0,0);


    //f2(1525)
    Variable v_f2_1525_Mass("f2_1525_MASS",f2_1525_MASS,0.0008,f2_1525_MASS*0.5,f2_1525_MASS*1.5);
    Variable v_f2_1525_Width("f2_1525_WIDTH",f2_1525_WIDTH,0.0025,f2_1525_WIDTH*0.5,f2_1525_WIDTH*1.5);
    Variable v_f2_1525_real("f2_1525_real",f2_1525_amp, 0.001, -100.0, +100.0);
    Variable v_f2_1525_img("f2_1525_img",f2_1525_img, 0.001, -100.0, +100.0);

    //f0(980)
    Variable v_f0_980_Mass("f0_980_MASS",f0_980_MASS,3.0,f0_980_MASS*0.5,f0_980_MASS*1.5);
    Variable v_f0_980_GPP("f0_980_GPP",f0_980_GPP,0.01,f0_980_GPP*0.5,f0_980_GPP*1.5);
    Variable v_f0_980_GKK("f0_980_GKK",f0_980_GKK,0.02,f0_980_GKK*0.5,f0_980_GKK*1.5);
    Variable v_f0_980_real("f0_980_real",f0_980_amp, 0.001, -100.0, +100.0);
    Variable v_f0_980_img("f0_980_img",f0_980_img, 0.001, -100.0, +100.0);

    //v_f0_980_real.setFixed(true);
    //v_f0_980_img.setFixed(true);
    //v_f2_1270_real.setFixed(true);
    //v_f2_1270_img.setFixed(true);
    
    //f0(1370)
    Variable v_f0_1370_Mass("f0_1370_MASS",f0_1370_MASS,0.4,f0_1370_MASS*0.1,f0_1370_MASS*1.5);
    Variable v_f0_1370_Width("f0_1370_Width",f0_1370_WIDTH,0.8,0.1,1.);
    Variable v_f0_1370_real("f0_1370_real",f0_1370_amp, 0.001, -100.0, +100.0);
    Variable v_f0_1370_img("f0_1370_img",f0_1370_img, 0.001, -100.0, +100.0);

    
    //f0(1500)
    Variable v_f0_1500_Mass("f0_1500_MASS",f0_1500_MASS,0.006,f0_1500_MASS*0.1,f0_1500_MASS*1.5);
    Variable v_f0_1500_Width("f0_1500_Width",f0_1500_WIDTH,0.007,0.1,1.);
    Variable v_f0_1500_real("f0_1500_real",f0_1500_amp, 0.001, -100.0, +100.0);
    Variable v_f0_1500_img("f0_1500_img",f0_1500_img, 0.001, -100.0, +100.0);


    //f0(X)
    Variable v_f0_X_Mass("f0_X_MASS",f0_X_MASS,0.006,f0_X_MASS*0.1,f0_X_MASS*1.5);
    Variable v_f0_X_Width("f0_X_Width",f0_X_WIDTH,0.007,0.1,1.);
    Variable v_f0_X_real("f0_X_real",f0_X_amp, 0.001, -100.0, +100.0);
    Variable v_f0_X_img("f0_X_img",f0_X_img, 0.001, -100.0, +100.0);


    //NR
    Variable nonr_real("nonr_real", fixed ? -0.9 : 1*donram.Uniform(-1.,1.), 0.001, -200.0, +200.0);
    Variable nonr_imag("nonr_imag", fixed ? -1.21 : 1*donram.Uniform(-1.,1.), 0.001, -200.0, +200.0);

    Variable be_real("be_real", fixed ? 1 : 1*donram.Uniform(-1.,1.), 0.001, -200.0, +200.0);
    Variable be_imag("be_imag", fixed ? 0 : 1*donram.Uniform(-1.,1.), 0.001, -200.0, +200.0);
    Variable be_coef("be_coef", fixed ? 1.9 : 1.9*donram.Uniform(-1.,1.), 0.001, -200.0, +200.0);
    //Masses and Widths fixed

    v_omega_Mass.setFixed(true);
    v_omega_Width.setFixed(true);
   
    v_rho770_Mass.setFixed(true);
    v_rho770_Width.setFixed(true);
    
    v_rho1450_Mass.setFixed(true);
    v_rho1450_Width.setFixed(true);
    
    
    v_f2_1270_Mass.setFixed(true);
    v_f2_1270_Width.setFixed(true);

    v_f2_1525_Mass.setFixed(true);
    v_f2_1525_Width.setFixed(true);

    //v_f0_980_Mass.setFixed(true);
    //v_f0_980_GKK.setFixed(true);
    //v_f0_980_GPP.setFixed(true);

    v_f0_1370_Mass.setFixed(true);
    v_f0_1370_Width.setFixed(true);


    v_f0_1500_Mass.setFixed(true);
    v_f0_1500_Width.setFixed(true);

    v_f0_X_Mass.setFixed(true);
    v_f0_X_Width.setFixed(true);

    be_coef.setFixed(true);
    be_real.setFixed(true);
    be_imag.setFixed(true);
    
   //setting resonances
ResonancePdf* rho770_12 = new Resonances::RBW("rho770",v_rho770_real,v_rho770_img,v_rho770_Mass,v_rho770_Width,1,PAIR_12,true);
    
    ResonancePdf* rho1450_12 = new Resonances::RBW("rho1450",v_rho1450_real,v_rho1450_img,v_rho1450_Mass,v_rho1450_Width,1,PAIR_12,true);


    ResonancePdf* omega_12 = new Resonances::GS("omega",v_omega_real,v_omega_img,v_omega_Mass,v_omega_Width,1,PAIR_12,true);
    
    ResonancePdf* f2_1270_12 = new Resonances::RBW("f2",v_f2_1270_real,v_f2_1270_img,v_f2_1270_Mass,v_f2_1270_Width,2,PAIR_12,true);

    ResonancePdf* f2_1525_12 = new Resonances::RBW("f2",v_f2_1525_real,v_f2_1525_img,v_f2_1525_Mass,v_f2_1525_Width,2,PAIR_12,true);
  
    ResonancePdf* f0_980_12 = new Resonances::FLATTE("f0_980",v_f0_980_real,v_f0_980_img,v_f0_980_Mass,v_f0_980_GPP,v_f0_980_GKK,PAIR_12,true);

    ResonancePdf* f0_1370_12 = new Resonances::RBW("f0_1370_12",v_f0_1370_real,v_f0_1370_img,v_f0_1370_Mass,v_f0_1370_Width,(unsigned int)0,PAIR_12,true);

    ResonancePdf* f0_1500_12 = new Resonances::RBW("f0_1500_12",v_f0_1500_real,v_f0_1500_img,v_f0_1500_Mass,v_f0_1500_Width,(unsigned int)0,PAIR_12,true);  

    ResonancePdf* f0_X_12 = new Resonances::RBW("f0_X_12",v_f0_X_real,v_f0_X_img,v_f0_X_Mass,v_f0_X_Width,(unsigned int)0,PAIR_12,true);  
     
    ResonancePdf *nonr = new Resonances::NonRes("nonr", nonr_real, nonr_imag);

    ResonancePdf *be   = new Resonances::BoseEinstein("be",be_real,be_imag,be_coef);

    //MIPWA
    ResonancePdf *swave_12 = loadPWAResonance(pwa_file, true);
    //Pushing Resonances 

    //dtoppp.resonances.push_back(omega_12);
     dtoppp.resonances.push_back(rho770_12); 
    dtoppp.resonances.push_back(rho1450_12);
     dtoppp.resonances.push_back(f2_1270_12);
    //dtoppp.resonances.push_back(f2_1525_12);
    dtoppp.resonances.push_back(f0_980_12);
    //dtoppp.resonances.push_back(f0_1370_12);
    //dtoppp.resonances.push_back(f0_1500_12);
    //dtoppp.resonances.push_back(f0_X_12);
    dtoppp.resonances.push_back(nonr);
    //dtoppp.resonances.push_back(be);
    dtoppp.resonances.push_back(swave_12);

    if(!eff) {
        // By default create a constant efficiency.
        vector<Variable> offsets;
        vector<Observable> observables;
        vector<Variable> coefficients;
        Variable constantOne("constantOne", 1);
        Variable constantZero("constantZero", 0);
        observables.push_back(s12);
        observables.push_back(s13);
        offsets.push_back(constantZero);
        offsets.push_back(constantZero);
        coefficients.push_back(constantOne);
        eff = new PolynomialPdf("constantEff", observables, coefficients, offsets, 0); //No efficiency
    }

    return new DalitzPlotPdf("signalPDF", s12, s13, eventNumber, dtoppp, eff);
}

std::pair<GooPdf*,DalitzPlotPdf*> TotalPdf(){ 

    if(effOn){
        std::vector<PdfBase *> comps;
        GooPdf* Dstar_veto = makeDstar_veto();
        GooPdf* effdalitz = makeEfficiencyPdf();
        comps.push_back(Dstar_veto);
        comps.push_back(effdalitz);   
        ProdPdf *effWithVeto = new ProdPdf("effWithVeto", comps);

        comps.clear();
        DalitzPlotPdf *signaldalitz = makesignalpdf(effWithVeto);
        comps.push_back(signaldalitz);

        
        
        if(bkgOn){
            GooPdf* bkgdalitz = makeBackgroundPdf();
            bkgdalitz->setParameterConstantness(true);
        
            comps.push_back(bkgdalitz);
            Variable constant("constant",0.95);
            std::vector<Variable> weights;
            weights.push_back(constant);

            AddPdf* overallPdf = new AddPdf("overallPdf",weights,comps);

            return std::pair<GooPdf*,DalitzPlotPdf*>(overallPdf,signaldalitz);
        }else{
            ProdPdf* overallPdf = new ProdPdf("overallPdf",comps);
            return std::pair<GooPdf*,DalitzPlotPdf*>(overallPdf,signaldalitz);
        }

        

             
    }else{
        std::vector<PdfBase *> comps;
        GooPdf* Dstar_veto = makeDstar_veto();
        comps.push_back(Dstar_veto);
        ProdPdf *effWithVeto = new ProdPdf("effWithVeto", comps);
        

        comps.clear();
        DalitzPlotPdf *signaldalitz = makesignalpdf(effWithVeto);
       
        comps.push_back(signaldalitz);

        
        
        if(bkgOn){
            GooPdf* bkgdalitz = makeBackgroundPdf();
            bkgdalitz->setParameterConstantness(true);
      
            comps.push_back(bkgdalitz);
            Variable constant("constant",0.95);
            std::vector<Variable> weights;
            weights.push_back(constant);

            AddPdf* overallPdf = new AddPdf("overallPdf",weights,comps);
           
            return std::pair<GooPdf*,DalitzPlotPdf*>(overallPdf,signaldalitz);
        }else{
            ProdPdf* overallPdf = new ProdPdf("overallPdf",comps);
            return std::pair<GooPdf*,DalitzPlotPdf*>(overallPdf,signaldalitz);
        }


    }
    

    
}


void getdata(std::string name){

    std::cout << "get data begin!" << '\n';

    Data = new UnbinnedDataSet({s12,s13,eventNumber});

    if(toyOn){
        std::ifstream reader(name.c_str());

        if(reader.is_open()){

        while(reader >> eventNumber >> s12 >> s13){

                Data->addEvent();

        }

        reader.close();
        }else{
            std::cout << "file not open" << std::endl;
        }

    }else{

    TFile *f = TFile::Open(data_name.c_str());
    TTree *t = (TTree *)f->Get(tree_name.c_str());

    double _s12, _s13;

    t->SetBranchAddress("s12_pipi_DTF",&_s12);
    t->SetBranchAddress("s13_pipi_DTF",&_s13);
  
    
   
    int j = 0;
    
    for(size_t i = 0; i < t->GetEntries() ; i++){

         t->GetEntry(i);

       
            if((_s12>s12_min)&&(_s12<s12_max)&&(_s13>s12_min)&&(_s13<s13_max)){
                s12.setValue(_s12);
                s13.setValue(_s13);
                eventNumber.setValue(i);
                Data->addEvent();

                 
                j++;

                if(j==250000){
                    break;
                }
                
            }    
               
             
          
        
            
            
    }


    f->Close();

    }
    
    std::cout << "get data end!" << '\n';
}

void runtoygen(std::string name, size_t events){

    s12.setNumBins(1500);
    s13.setNumBins(1500);

    auto overallPdf = TotalPdf();  

    
    {
        maketoydalitzdata(overallPdf.first,overallPdf.second,name,events);
    }
}

void saveParameters(std::string file, const std::vector<Variable> &param){

   // std::cout << param[0].GetName() << std::endl;

    ofstream wr(file.c_str());

    for(int i = 0; i < param.size(); i++){
        if( !( param[i].IsFixed()) ){
            wr << param[i].getName() <<'\t'<< param[i].getValue()<< '\t' << param[i].getError() << endl;
        }
    }


    wr.close();

}

void saveParameters(std::string file, const std::vector<Variable> &param,   bool isValid, double fcn, std::vector<std::vector<fptype>> ff ){

    DalitzPlotPdf* signaldalitz = makesignalpdf(0);
    size_t n_res = signaldalitz->getDecayInfo().resonances.size();
    
    ofstream wr(file.c_str());

    for(int i = 0; i < param.size(); i++){
        if( !(param[i].IsFixed()) ){
            wr << param[i].getName() <<'\t'<< param[i].getValue()<< '\t' << param[i].getError() << endl;
        }
    }

    for(int i = 0; i < n_res; i++){
        wr << ("FF_"+ to_string(i) ).c_str() << '\t' << ff[i][i]  << '\t' << 0 << endl;
    }

    wr << "FCN" << '\t' << fcn << '\t' << 0.0 << '\t' << 0.0 << endl;
    wr << "Status" << '\t' << isValid << '\t' << 0.0 << '\t' << 0.0 << endl;
    
    wr.close();

}

void runtoyfit(std::string name, int sample_number) {

    s12.setNumBins(bins);
    s13.setNumBins(bins);

    getdata(name);

    if(toyOn){
        GOOFIT_INFO("Using "+ name);
    }else{
        GOOFIT_INFO("Using DATA: {}",data_name );
    }

    GOOFIT_INFO("Number of Events in dataset: {}", Data->getNumEvents());

    auto overallPdf = TotalPdf();

    overallPdf.first->setData(Data);
    overallPdf.second->setDataSize(Data->getNumEvents());
 
    
    GooFit::FitManagerMinuit1 fitter(overallPdf.first);
    fitter.setVerbosity(1);
    fitter.setMaxCalls(4000);

    auto var = fitter.getMinuitObject()->getVaraibles();
    
    std::string command = "mkdir -p Fit";
    if(system(command.c_str()) != 0)
        throw GooFit::GeneralError("Making `Fit` directory failed");

    auto params = fitter.getMinuitObject()->getVaraibles();
    string input_name = fmt::format("Fit/fit_parameters_inicial.txt");
    saveParameters(input_name,params);
 
    fitter.getMinuitObject()->fIstrat = 1;
    fitter.getMinuitObject()->SetErrorDef(0.5);
    fitter.setMaxCalls(200000);
    fitter.useHesse(false);
    fitter.useHesseBefore(false);
    //fitter.useImprove(true);
    fitter.fit();
    
    auto ff = overallPdf.second->fit_fractions();
    params  =  fitter.getMinuitObject()->getVaraibles();
    fptype foo, fcn = -1;
    int foo2, status = -1;
    fitter.getMinuitObject()->mnstat(fcn, foo,foo,foo2,foo2,status);
    string output_name = fmt::format("Fit/fit_parameters_{0}.txt",sample_number);
    saveParameters(output_name , params ,status, fcn , ff );

//
    // remove comment for plotting 
//    DalitzPlotter dp(overallPdf.first,overallPdf.second);
 // dp.Plot(D_MASS,d1_MASS,d2_MASS,d3_MASS,"#pi^{-} #pi^{+}","#pi^{-} #pi^{+}","#pi^{-} #pi^{+}","plots",*Data);

std::ofstream wr("plots/PWACoefs_fitted.txt");

	for(int i = 0; i < HH_bin_limits.size(); i++){

		wr << sqrt( HH_bin_limits[i] ) << "\t" << pwa_coefs_amp[i].getValue() << "\t" << pwa_coefs_phs[i].getValue() << std::endl;
	}

wr.close();

}


void plot(string name){

    getdata(name);

    auto overallPdf = TotalPdf();
    overallPdf.first->setData(Data);
    overallPdf.second->setDataSize(Data->getNumEvents());

    DalitzPlotter dp(overallPdf.first,overallPdf.second);
    dp.Plot(D_MASS,d1_MASS,d2_MASS,d3_MASS,"#pi^{-} #pi^{+}","#pi^{-} #pi^{+}","#pi^{-} #pi^{+}","plots",*Data);
}

void genfitplot(int nsamples,int nvar){

    double vec_inicial[nvar];
       
    DalitzPlotPdf* signaldalitz = makesignalpdf(0);

   //getting inicial parameters

    int index = 0;
    double inicial_var = 0;
    double inicial_error =0;
    
    size_t n_res = signaldalitz->getDecayInfo().resonances.size();
    std::string inicial_name;
    ifstream r_inicial("Fit/fit_parameters_inicial.txt");
	
    TFile *f = new TFile("GenFit_Results.root","CREATE");
    TTree *t = new TTree("tree","Variables"); 

    double val[nvar];
    double val_error[nvar];
    double ff_val[n_res];
      
    while(r_inicial >> inicial_name >> inicial_var >> inicial_error){
            vec_inicial[index] = inicial_var;
	    t->Branch(inicial_name.c_str(),val+index);
	    t->Branch((inicial_name+"_error").c_str(),val_error+index);
            index++;
    }

    for(int i=0; i < n_res ; i++){
        t->Branch( ("FF_"+to_string(i)).c_str(), ff_val + i);
    }

    

    double FCN, Status;
    t->Branch("FCN",&FCN);
    t->Branch("Status",&Status);
    

    //getting fit parameters
   std::string name;
   std::string line;
 
       for(int i = 0; i < nsamples; i++){
        
            name = fmt::format("Fit/fit_parameters_{0}.txt",i);
            ifstream r_fit(name.c_str());

            if(r_fit.is_open()){

                int index_fit = 0;
                double fit_var = 0;
		        double fit_error=0;
		        std::string fit_name;
                int count = 0;
                int index_ff = 0;

                while(std::getline(r_fit,line)){

                    std::stringstream ss(line);

                    if(ss >> name >> fit_var >> fit_error){

                        if(count<nvar){
                            val[index_fit] = fit_var;
		                    val_error[index_fit] = fit_error; 
                            index_fit++;  
                        }

                         if( (count >= nvar)&&(count < (nvar+n_res)) ){
                            ff_val[index_ff] = fit_var;
                            index_ff++;
                        }

                        if(count==(nvar+n_res)){
                            FCN = fit_var;
                                                   
                        }
                    
                        if(count==(nvar+n_res+1)){
                            Status = fit_var;
                        }

                        count++;

                    }
                    
                }

		        t->Fill();

                r_fit.close();

            }else{
                printf("Error \n");
            }
       }

   t->Write();
   f->Close();
}


int main(int argc, char **argv){

    auto sample_number = 0;
    GooFit::Application app{"D2PPP",argc,argv};
    app.add_option("-i,--int", sample_number, "sample number", true);

    size_t  nevents = 100000;
    auto gen = app.add_subcommand("gen","generate toy data");
    gen->add_option("-e,--events",nevents,"The number of events to generate",true);

    auto toyfit = app.add_subcommand("fit","fit data/toyMC");
    auto plotpdf = app.add_subcommand("plot","plot");

    int ns = 10;
    int nv = 10;
    auto gfplot = app.add_subcommand("gfplot","genfit plot");
    gfplot->add_option("-n,--ns",ns,"n samples",true);
    gfplot->add_option("-v,--nv",nv,"n variables",true);
    app.require_subcommand();


    GOOFIT_PARSE(app);

    /// Make the plot directory if it does not exist
    std::string command = "mkdir -p plots";
    if(system(command.c_str()) != 0)
        throw GooFit::GeneralError("Making `plots` directory failed");

    /// Make the MC directory if it does not exist
    std::string command2 = "mkdir -p MC";
    if(system(command2.c_str()) != 0)
        throw GooFit::GeneralError("Making `MC` directory failed");

    std::string name = fmt::format("MC/MC_Toy_{0}",sample_number);

    if(*gen){
        CLI::AutoTimer timer("MC Generation");
        runtoygen(name,nevents);
    }

    if(*toyfit){
        CLI::AutoTimer timer("FIT");
        runtoyfit(name,sample_number);
    }

    if(*plotpdf){
        CLI::AutoTimer timer("plot");
        plot(name);
    }

    if(*gfplot){
        CLI::AutoTimer timer("plot genfit result");
        genfitplot(ns,nv);
    }



}

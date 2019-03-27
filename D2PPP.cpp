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
#include <TGraphErrors.h>
#include <TVector.h>

// System stuff
#include <CLI/Timer.hpp>
#include <fstream>
#include <string>
#include <fmt/format.h>

// GooFit stuff
#include <goofit/Application.h>
#include <goofit/BinnedDataSet.h>
#include <goofit/FitManager.h>
#include <goofit/fitting/FitManagerMinuit2.h>
#include <goofit/fitting/Params.h>
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
double D_MASS   = 1.86959; //GEV

double d1_MASS  = pi_MASS;  //daughter 1 mass
double d2_MASS  = pi_MASS;
double d3_MASS  = pi_MASS;

bool saveBkgPlot= true;
bool saveEffPlot= true;
bool doEffSwap  = true;
bool toyOn      = true;
bool bkgOn      = false;

const double NevG = 1e7; 

fptype s12_min = POW2(d1_MASS  + d2_MASS);
fptype s12_max = POW2(D_MASS   - d2_MASS);
fptype s13_min = POW2(d1_MASS  + d3_MASS);
fptype s13_max = POW2(D_MASS   - d2_MASS);

Observable s12("s12",s12_min,s12_max); //s12^{2}
Observable s13("s13",s13_min,s13_max);
EventNumber eventNumber("eventNumber");

DalitzPlotPdf* signaldalitz = nullptr;
SmoothHistogramPdf* bkgdalitz = nullptr;

UnbinnedDataSet* Data = nullptr;

std::vector<PdfBase *> comps;
vector<fptype> HH_bin_limits;
vector<Variable> pwa_coefs_amp;
vector<Variable> pwa_coefs_phs;

Variable massSum("massSum", POW2(D_MASS) + POW2(d1_MASS) + POW2(d2_MASS) + POW2(d3_MASS));


// PWA INPUT FILE NAME
const string pwa_file = "files/PWACOEFS.txt";

// Data File
const string data_name = "/mnt/DATA/Dropbox/D2PPP/juan/ntuples/dados_finais/DsPPP_TIS_PosTS_par_sPloted_large__MVA.root";
const string tree_name = "sTree";
const string bkghisto_file = "files/bkghisto.root";

//functions
fptype cpuGetM23(fptype massPZ, fptype massPM) { return (massSum.getValue() - massPZ - massPM); }

DalitzPlotPdf *makesignalpdf(GooPdf *eff = 0);
void getdata(std::string name);

TH2F *weightHistogram    = nullptr;
TH2F *bkgHistogram       = nullptr;
TH2F *underlyingBins     = nullptr;


void maketoydalitzdata(GooPdf* overallsignal,std::string name, size_t nEvents){

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

    std::cout << "toy Generation end!" << '\n';

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

        HH_bin_limits.push_back(e1);

        emag = e2;
        ephs = e3;

        Variable va(fmt::format("pwa_coef_{}_real", i), emag,0.001,-100.0,+100.0);
        Variable vp(fmt::format("pwa_coef_{}_img", i), ephs,0.001,-100.0,+100.0);

        pwa_coefs_amp.push_back(va);
        pwa_coefs_phs.push_back(vp);
        i++;

    }

    Variable swave_real("swave_real", 1.0,-100.0,+100.0);
    Variable swave_imag("swave_imag", 0.0,-100.0,+100.0);

    if(fixAmp) {
        swave_real.setValue(1.);
        swave_imag.setValue(0.);
        swave_real.setFixed(true);
        swave_imag.setFixed(true);
    }
    cout << "Numbers loaded: " << HH_bin_limits.size() << " / " << i << endl;

    ResonancePdf *swave_12 = new Resonances::Spline("swave_12", swave_real, swave_imag, HH_bin_limits, pwa_coefs_amp, pwa_coefs_phs, PAIR_12, true);

    return swave_12;
} 


SmoothHistogramPdf* makeEfficiencyPdf() {

    vector<Observable> lvars;
    lvars.push_back(s12);
    lvars.push_back(s13);
    BinnedDataSet *binEffData = new BinnedDataSet(lvars);
    
    TFile *f     = TFile::Open(data_name.c_str());
    TH2F *bkgHistogram = (TH2F *)f->Get("h0");
    bkgHistogram->SetStats(false);

    TRandom3 donram(0);
    for(int i = 0; i < NevG; i++) {
        do {
            s12.setValue(donram.Uniform(s12.getLowerLimit(), s12.getUpperLimit()));
            s13.setValue(donram.Uniform(s13.getLowerLimit(), s13.getUpperLimit()));
        } while(!inDalitz(s12.getValue(), s13.getValue(), D_MASS, d1_MASS, d2_MASS, d3_MASS));

        double weight = weightHistogram->GetBinContent(weightHistogram->FindBin(s12.getValue(), s13.getValue()));
        binEffData->addWeightedEvent(weight);

        if(doEffSwap) {
            double swapmass = s12.getValue();
            s12.setValue(s13.getValue());
            s13.setValue(swapmass);
            weight = weightHistogram->GetBinContent(weightHistogram->FindBin(s12.getValue(), s13.getValue()));
            binEffData->addWeightedEvent(weight);
        }
    }
    if(saveEffPlot) {
        TCanvas foo;
        foo.cd();
        weightHistogram->Draw("colz");
        foo.SaveAs("plots/efficiency_bins.png");
        foo.SetLogz(true);
        foo.SaveAs("plots/efficiency_bins_log.png");
    }
    // Smooth
    Variable effSmoothing("effSmoothing", 1);
    SmoothHistogramPdf *ret = new SmoothHistogramPdf("efficiency", binEffData, effSmoothing);
    return ret;
}

SmoothHistogramPdf* makeBackgroundPdf() {


    s12.setNumBins(120);
    s13.setNumBins(120);

   BinnedDataSet *binBkgData = new BinnedDataSet({s12, s13});
    
    TFile *f = new TFile(bkghisto_file.c_str());
    
    
    TH2F* bkgHistogram = (TH2F*)f->Get("h0");

       TRandom3 donram(50);
    for(int i = 0; i < NevG; i++) {
        do {
            s12.setValue(donram.Uniform(s12.getLowerLimit(), s12.getUpperLimit()));
            s13.setValue(donram.Uniform(s13.getLowerLimit(), s13.getUpperLimit()));
        } while(!inDalitz(s12.getValue(), s13.getValue(), D_MASS, d1_MASS, d2_MASS, d3_MASS));

        double weight = bkgHistogram->GetBinContent(bkgHistogram->FindBin(s12.getValue(), s13.getValue()));

        binBkgData->addWeightedEvent(weight);

        if(doEffSwap) {
            double swapmass = s12.getValue();
            s12.setValue(s13.getValue());
            s13.setValue(swapmass);
            weight = bkgHistogram->GetBinContent(bkgHistogram->FindBin(s12.getValue(), s13.getValue()));
            binBkgData->addWeightedEvent(weight);
           
        }
    }

    for(int i = 0; i < 10 ; i++){
       cout << i << "\t" <<  s12.getValue() << '\t' << s13.getValue() << '\t' << binBkgData->getBinContent(i) << endl;
    }


    if(saveBkgPlot) {
        TCanvas foo;
        bkgHistogram->Draw("colz");
        foo.SetLogz(false);
        foo.SaveAs("plots/background_bins.png");
        foo.SetLogz(true);
        foo.SaveAs("plots/background_bins_log.png");
    }

    

    Variable *effSmoothing  = new Variable("effSmoothing", 1.0);
    SmoothHistogramPdf *ret = new SmoothHistogramPdf("efficiency", binBkgData, *effSmoothing);

    s12.setNumBins(1500);
    s13.setNumBins(1500);
    return ret;
}


DalitzPlotPdf* makesignalpdf(GooPdf* eff){

    DecayInfo3 dtoppp;
    dtoppp.motherMass   = D_MASS;
    dtoppp.daug1Mass    = d1_MASS;
    dtoppp.daug2Mass    = d2_MASS;
    dtoppp.daug3Mass    = d3_MASS;
    dtoppp.meson_radius = 1.5;

    //Mass and width
    
    double f0_980_MASS     = 0.965;
    double f0_980_GPP     = 0.165;
    double f0_980_GKK     = 0.695;
    double f0_980_amp     = 1.0;
    double f0_980_img    = 0.0;

    double f0_1500_MASS  = 1.505;
    double f0_1500_WIDTH = .109;
    double f0_1500_amp   = 0.3;
    double f0_1500_img = 0.8;

    double f0_X_MASS  = 1.430;
    double f0_X_WIDTH = 0.320;
    double f0_X_amp   = 0.5;
    double f0_X_img = 0.5;

    double omega_MASS   = 0.78265;
    double omega_WIDTH  = 0.00849;
    double omega_amp    = .8;
    double omega_img  = .8;

    double f2_1270_MASS     = 1.2755;
    double f2_1270_WIDTH    = 0.1867;
    double f2_1270_amp      = .3;
    double f2_1270_img    = .8;

    double f2_1525_MASS     = 1.525;
    double f2_1525_WIDTH    = 0.073;
    double f2_1525_amp      = .6;
    double f2_1525_img    = .3;
    

    //omega(782)
    Variable v_omega_Mass("omega_MASS",omega_MASS);
    Variable v_omega_Width("omega_WIDTH",omega_WIDTH);
    Variable v_omega_real("omega_real",omega_amp, 0.001, -100.0, +100.0);
    Variable v_omega_img("omega_img",omega_img, 0.001, -100.0, +100.0);

    //f2(1270)
    Variable v_f2_1270_Mass("f2_1270_MASS",f2_1270_MASS,0.0008,f2_1270_MASS*0.5,f2_1270_MASS*1.5);
    Variable v_f2_1270_Width("f2_1270_WIDTH",f2_1270_WIDTH,0.0025,f2_1270_WIDTH*0.5,f2_1270_WIDTH*1.5);
    Variable v_f2_1270_real("f2_1270_real",f2_1270_amp, 0.001, -100.0, +100.0);
    Variable v_f2_1270_img("f2_1270_img",f2_1270_img, 0.001, -100.0, +100.0);

    //f2(1525)
    Variable v_f2_1525_Mass("f2_1525_MASS",f2_1525_MASS,0.0008,f2_1525_MASS*0.5,f2_1525_MASS*1.5);
    Variable v_f2_1525_Width("f2_1525_WIDTH",f2_1525_WIDTH,0.0025,f2_1525_WIDTH*0.5,f2_1525_WIDTH*1.5);
    Variable v_f2_1525_real("f2_1525_real",f2_1525_amp, 0.001, -100.0, +100.0);
    Variable v_f2_1525_img("f2_1525_img",f2_1525_img, 0.001, -100.0, +100.0);

    //f0(980)
    Variable v_f0_980_Mass("f0_980_MASS",f0_980_MASS,3.0,f0_980_MASS*0.5,f0_980_MASS*1.5);
    Variable v_f0_980_GPP("f0_980_GPP",f0_980_GPP,0.01,f0_980_GPP*0,f0_980_GPP*20.0);
    Variable v_f0_980_GKK("f0_980_GKK",f0_980_GKK,0.02,f0_980_GKK*0,f0_980_GKK*3.0);
    Variable v_f0_980_real("f0_980_real",f0_980_amp, 0.001, -100.0, +100.0);
    Variable v_f0_980_img("f0_980_img",f0_980_img, 0.001, -100.0, +100.0);

    v_f0_980_real.setFixed(true);
    v_f0_980_img.setFixed(true);

    
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
    Variable nonr_real("nonr_real", 1.0, 0.001, -200.0, +200.0);
    Variable nonr_imag("nonr_imag", 1.0, 0.001, -200.0, +200.0);

    Variable be_real("be_real", 1.0, 0.001, -200.0, +200.0);
    Variable be_imag("be_imag", 0., 0.001, -200.0, +200.0);
    Variable be_coef("be_coef", 1.9, 0.001, -200.0, +200.0);
    //Masses and Widths fixed

    v_omega_Mass.setFixed(true);
    v_omega_Width.setFixed(true);
   
    v_f2_1270_Mass.setFixed(true);
    v_f2_1270_Width.setFixed(true);

    v_f2_1525_Mass.setFixed(true);
    v_f2_1525_Width.setFixed(true);

    v_f0_980_Mass.setFixed(true);
    v_f0_980_GKK.setFixed(true);
    v_f0_980_GPP.setFixed(true);

    v_f0_1500_Mass.setFixed(true);
    v_f0_1500_Width.setFixed(true);

    v_f0_X_Mass.setFixed(true);
    v_f0_X_Width.setFixed(true);

    be_coef.setFixed(true);
    
   //setting resonances
   
    ResonancePdf* omega_12 = new Resonances::GS("omega",v_omega_real,v_omega_img,v_omega_Mass,v_omega_Width,1,PAIR_12,true);
    
    ResonancePdf* f2_1270_12 = new Resonances::RBW("f2",v_f2_1270_real,v_f2_1270_img,v_f2_1270_Mass,v_f2_1270_Width,2,PAIR_12,true);

    ResonancePdf* f2_1525_12 = new Resonances::RBW("f2",v_f2_1525_real,v_f2_1525_img,v_f2_1525_Mass,v_f2_1525_Width,2,PAIR_12,true);
  
    ResonancePdf* f0_980_12 = new Resonances::FLATTE("f0_980",v_f0_980_real,v_f0_980_img,v_f0_980_Mass,v_f0_980_GPP,v_f0_980_GKK,PAIR_12,true);

    ResonancePdf* f0_1500_12 = new Resonances::RBW("f0_1500_12",v_f0_1500_real,v_f0_1500_img,v_f0_1500_Mass,v_f0_1500_Width,(unsigned int)0,PAIR_12,true);  

    ResonancePdf* f0_X_12 = new Resonances::RBW("f0_X_12",v_f0_X_real,v_f0_X_img,v_f0_X_Mass,v_f0_X_Width,(unsigned int)0,PAIR_12,true);  
     
    ResonancePdf *nonr = new Resonances::NonRes("nonr", nonr_real, nonr_imag);

    ResonancePdf *be   = new Resonances::BoseEinstein("be",be_real,be_imag,be_coef);

    //MIPWA
   //ResonancePdf *swave_12 = loadPWAResonance(pwa_file, true);

    //Pushing Resonances 

    dtoppp.resonances.push_back(omega_12);
    dtoppp.resonances.push_back(f2_1270_12);
    dtoppp.resonances.push_back(f2_1525_12);
    dtoppp.resonances.push_back(f0_980_12);
    dtoppp.resonances.push_back(f0_1500_12);
    dtoppp.resonances.push_back(f0_X_12);
    //dtoppp.resonances.push_back(nonr);
    //dtoppp.resonances.push_back(be);
    //dtoppp.resonances.push_back(swave_12);

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





void getdata(std::string name){

    std::cout << "get data begin!" << '\n';

    Data = new UnbinnedDataSet({s12,s13,eventNumber});

if(toyOn){
    std::ifstream reader(name.c_str());

    while(reader >> eventNumber >> s12 >> s13){

        	Data->addEvent();

    }

    reader.close();

}else{

    TFile *f = TFile::Open(data_name.c_str());
    TTree *t = (TTree *)f->Get("sTree");

    double _s12, _s13;

    t->SetBranchAddress("s12_pipi_DTF",&_s12);
    t->SetBranchAddress("s13_pipi_DTF",&_s13);
   

    for(size_t i = 0; i < 100000 ; i++){

         t->GetEntry(i);

       if( (_s12<s12.getUpperLimit()) && (_s13<s13.getUpperLimit()) ){
            s12.setValue(_s12);
            s13.setValue(_s13);
            eventNumber.setValue(i);
            Data->addEvent();
       }
            
    }


    f->Close();

}
    
    std::cout << "get data end!" << '\n';
}

void runtoygen(std::string name, size_t events){

    s12.setNumBins(1500);
    s13.setNumBins(1500);

    signaldalitz = makesignalpdf(0);
    
    Variable constant("constant",1);
    std::vector<Variable> weights;
    weights.push_back(constant);
    
    comps.clear();
    comps.push_back(signaldalitz);

    if(bkgOn){
        bkgdalitz = makeBackgroundPdf();
        bkgdalitz->setParameterConstantness(true);
 	    comps.push_back(bkgdalitz);
        weights.push_back(constant);
    }
    
    AddPdf* overallPdf = new AddPdf("overallPdf",weights,comps);
    
    
    {
        maketoydalitzdata(overallPdf,name,events);
    }
}

void PrintFF(std::vector<std::vector<fptype>> ff){

    size_t nEntries = signaldalitz->getCachedWave(0).size();
    size_t n_res = signaldalitz->getDecayInfo().resonances.size();
    fptype sum = 0;

    std::cout << "nEntries= " << nEntries << '\n';
    for(size_t i = 0; i < n_res ; i++){

        for(size_t j = 0; j< n_res ; j++){
            std::cout << "FF[" << i << "," << j <<"]= " << ff[i][j] << std::endl;

        }

        sum+=ff[i][i];
    }

    std::cout << "Sum[i,i]= " << sum << std::endl;
}


void drawFitPlotsWithPulls(TH1 *hd, TH1 *ht, string plotdir) {
    const char *hname = hd->GetName();
    char obsname[10];
    for(int i = 0;; i++) {
        if(hname[i] == '_')
            obsname[i] = '\0';
        else
            obsname[i] = hname[i];
        if(obsname[i] == '\0')
            break;
    }
    ht->Scale(hd->Integral() / ht->Integral()*5);
    ht->SetLineColor(kRed);
    ht->SetLineWidth(3);
    ht->SetMarkerStyle(0);

    hd->SetMarkerColor(kBlack);
    hd->Rebin(5);


    TCanvas foo;

    hd->Draw("E");
    ht->Draw("HIST C same");


    foo.SaveAs(TString::Format("plots/%s_fit.png",obsname));


}


void makeToyDalitzPdfPlots(GooPdf *overallSignal, string plotdir = "plots") {
    TH1F s12_dat_hist("s12_dat_hist", "", s12.getNumBins(), s12.getLowerLimit(), s12.getUpperLimit());
    s12_dat_hist.GetXaxis()->SetTitle("m^{2}(K^{-} K^{+}) [GeV]");
    s12_dat_hist.GetYaxis()->SetTitle(TString::Format("Events / %.1f MeV", 1e3 * s12_dat_hist.GetBinWidth(1)));

    TH1F s12_pdf_hist("s12_pdf_hist", "", s12.getNumBins(), s12.getLowerLimit(), s12.getUpperLimit());

    TH1F s13_dat_hist("s13_dat_hist", "", s13.getNumBins(), s13.getLowerLimit(), s13.getUpperLimit());
    s13_dat_hist.GetXaxis()->SetTitle("m^{2}(K^{-} K^{+}) [GeV]");
    s13_dat_hist.GetYaxis()->SetTitle(TString::Format("Events / %.1f MeV", 1e3 * s13_dat_hist.GetBinWidth(1)));

    TH1F s13_pdf_hist("s13_pdf_hist", "", s13.getNumBins(), s13.getLowerLimit(), s13.getUpperLimit());

    TH1F s23_dat_hist("s23_dat_hist", "", s13.getNumBins(), s13.getLowerLimit(), s13.getUpperLimit());
    s23_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{+} #pi^{-}) [GeV]");
    s23_dat_hist.GetYaxis()->SetTitle(TString::Format("Events / %.1f MeV", 1e3 * s13_dat_hist.GetBinWidth(1)));

    TH1F s23_pdf_hist("s23_pdf_hist", "", s13.getNumBins(), s13.getLowerLimit(), s13.getUpperLimit());

    double totalPdf = 0;
    double totalDat = 0;
    TH2F dalitzpp0_dat_hist("dalitzpp0_dat_hist",
                            "",
                            s12.getNumBins(),
                            s12.getLowerLimit(),
                            s12.getUpperLimit(),
                            s13.getNumBins(),
                            s13.getLowerLimit(),
                            s13.getUpperLimit());
    dalitzpp0_dat_hist.SetStats(false);
    dalitzpp0_dat_hist.GetXaxis()->SetTitle("m^{2}(K^{-} K^{+}) [GeV]");
    dalitzpp0_dat_hist.GetYaxis()->SetTitle("m^{2}(K^{-} K^{+}) [GeV^{2}]");
    TH2F dalitzpp0_pdf_hist("dalitzpp0_pdf_hist",
                            "",
                            s12.getNumBins(),
                            s12.getLowerLimit(),
                            s12.getUpperLimit(),
                            s13.getNumBins(),
                            s13.getLowerLimit(),
                            s13.getUpperLimit());

    dalitzpp0_pdf_hist.GetXaxis()->SetTitle("m^{2}(K^{-} K^{+}) [GeV^{2}]");
    dalitzpp0_pdf_hist.GetYaxis()->SetTitle("m^{2}(K^{-} K^{+}) [GeV^{2}]");
    dalitzpp0_pdf_hist.SetStats(false);
    std::vector<Observable> vars;
    vars.push_back(s12);
    vars.push_back(s13);
    vars.push_back(eventNumber);
    UnbinnedDataSet currData(vars);
    int evtCounter = 0;

    for(int i = 0; i < s12.getNumBins(); ++i) {
        s12.setValue(s12.getLowerLimit() + (s12.getUpperLimit() - s12.getLowerLimit()) * (i + 0.5) / s12.getNumBins());
        for(int j = 0; j < s13.getNumBins(); ++j) {
            s13.setValue(s13.getLowerLimit()
                         + (s13.getUpperLimit() - s13.getLowerLimit()) * (j + 0.5) / s13.getNumBins());
            if(!inDalitz(s12.getValue(), s13.getValue(), D_MASS, d1_MASS, d2_MASS, d3_MASS)){
                continue;}
            eventNumber.setValue(evtCounter);
            evtCounter++;
            currData.addEvent();
        }
    }
    overallSignal->setData(&currData);
    signaldalitz->setDataSize(currData.getNumEvents());
    std::vector<std::vector<double>> pdfValues = overallSignal->getCompProbsAtDataPoints();
    for(unsigned int j = 0; j < pdfValues[0].size(); ++j) {
        double currs12 = currData.getValue(s12, j);
        double currs13 = currData.getValue(s13, j);

        dalitzpp0_pdf_hist.Fill(currs12, currs13, pdfValues[0][j]);
        s12_pdf_hist.Fill(currs12, pdfValues[0][j]);
        s13_pdf_hist.Fill(currs13, pdfValues[0][j]);
        s23_pdf_hist.Fill(cpuGetM23(currs12, currs13), pdfValues[0][j]);
        totalPdf += pdfValues[0][j];
    }

    TCanvas foo;
    foo.SetLogz(false);
    dalitzpp0_pdf_hist.Draw("colz");

    foo.SaveAs("plots/dalitzpp0_pdf.png");

    for(unsigned int evt = 0; evt < Data->getNumEvents(); ++evt) {
        double data_s12 = Data->getValue(s12, evt);
        s12_dat_hist.Fill(data_s12);
        double data_s13 = Data->getValue(s13, evt);
        s13_dat_hist.Fill(data_s13);
        dalitzpp0_dat_hist.Fill(data_s12, data_s13);
        s23_dat_hist.Fill(cpuGetM23(data_s12, data_s13));
        totalDat++;
    }
    dalitzpp0_dat_hist.Draw("colz");
    foo.SaveAs("plots/dalitzpp0_dat.png");

    drawFitPlotsWithPulls(&s12_dat_hist, &s12_pdf_hist, plotdir);
    drawFitPlotsWithPulls(&s13_dat_hist, &s13_pdf_hist, plotdir);
    drawFitPlotsWithPulls(&s23_dat_hist, &s23_pdf_hist, plotdir);
}

void runMakeToyDalitzPdfPlots(std::string name){

    s12.setNumBins(1500);
    s13.setNumBins(1500);

    getdata(name);

    signaldalitz = makesignalpdf(0);
    
    if(bkgOn){

        bkgdalitz = makeBackgroundPdf();
        bkgdalitz->setParameterConstantness(true);

    }
    Variable constant("constant",1);
    std::vector<Variable> weights;
    weights.push_back(constant);
    
    comps.clear();
    comps.push_back(signaldalitz);
    if(bkgOn){
    comps.push_back(bkgdalitz);
    }

    AddPdf* overallPdf = new AddPdf("overallPdf",weights,comps);
    overallPdf->setData(Data);
    signaldalitz->setDataSize(Data->getNumEvents());

    makeToyDalitzPdfPlots(overallPdf);

}

void saveParameters(const std::vector<ROOT::Minuit2::MinuitParameter> &param, std::string file){

   // std::cout << param[0].GetName() << std::endl;

    ofstream wr(file.c_str());

    for(int i = 0; i < param.size(); i++){
        if( !(param[i].IsConst() || param[i].IsFixed()) ){
            wr << param[i].GetName() <<'\t'<< param[i].Value()<< '\t' << param[i].Error() << endl;
        }
    }
    
    wr.close();

}




double runtoyfit(std::string name, int sample_number,int bins){

    s12.setNumBins(bins);
    s13.setNumBins(bins);

    getdata(name);

    if(toyOn){
        GOOFIT_INFO("Using toyDATA");
    }else{
        GOOFIT_INFO("Using DATA: {}",data_name );
    }

    GOOFIT_INFO("Number of Events in dataset: {}", Data->getNumEvents());
   
    if(signaldalitz == nullptr){ 
    	signaldalitz = makesignalpdf(0);
    }

    Variable constant("constant",1);
    std::vector<Variable> weights;
    weights.push_back(constant);
    
    comps.clear();
    comps.push_back(signaldalitz);
    
    if(bkgOn){
	if(bkgdalitz == nullptr){
        bkgdalitz = makeBackgroundPdf();
        }
        bkgdalitz->setParameterConstantness(true);
        comps.push_back(bkgdalitz);
       
    }

    AddPdf* overallPdf = new AddPdf("overallPdf",weights,comps);
    overallPdf->setData(Data);
    signaldalitz->setDataSize(Data->getNumEvents());

    //ROOT::Math::MinimizerOptions opt;
    //opt.SetTolerance(0.01);
    FitManagerMinuit2 fitter(overallPdf);
    fitter.setVerbosity(2);
    //fitter.setPrint(0);

    std::string command = "mkdir -p Fit";
    if(system(command.c_str()) != 0)
        throw GooFit::GeneralError("Making `Fit` directory failed");


    auto params = fitter.getParams()->Parameters();

    string input_name = fmt::format("Fit/fit_parameters_inicial.txt");
    
    saveParameters(params, input_name); 

    
    auto func_min = fitter.fit();

    params.clear();
    params  = fitter.getParams()->Parameters();

    string output_name = fmt::format("Fit/fit_parameters_{0}.txt",sample_number);
    
    saveParameters(params, output_name);



    return overallPdf->normalize();

    }

//normalization study

void normStudy(std::string name, int sample_number,int nbins_i,int nbins_f){

  std::ofstream wr("norms.txt");

  while(nbins_i != nbins_f){
	
         wr << nbins_i << '\t' << runtoyfit(name,sample_number,nbins_i) << std::endl;
        
        nbins_i = nbins_i + 100.0;		
   
  }	

  wr.close();

  TGraph h1("norms.txt","%lg %lg");
  TCanvas foo;
  h1.Draw("AP*"); 
  foo.SaveAs("normStudies.jpg");

}




/*
void genfitplot(const int nsamples,const int nvar){

    ifstream r_i("Fit/fit_parameters_inicial.txt");

    double val_i,error_i;
    std::string foo_i;
    double vec_val_i[nvar];
    double vec_error_i[nvar];

    size_t l =0;
    while(r_i >> foo_i >> val_i >> error_i){

        vec_val_i[l] = val_i;
        vec_error_i[l] = error_i;
        l++;
    }

    r_i.close();

    std::string name;
    double val;
    double error;
    double vec_val[nvar][nsamples];
    double vec_error[nvar][nsamples];
    std::string vec_name[nvar];

    size_t val_index, var_index = 0;

    for(size_t i = 0; i < nsamples; i++){
       
        name = fmt::format("Fit/fit_parameters_{0}.txt",i);
        ifstream rd(name.c_str());
    
        while(rd >> name >> val >> error){

            vec_val[var_index][val_index] = val;
            vec_error[var_index][val_index] = error;
            vec_name[var_index] = name;
            var_index++;
   
         }

        rd.close();
        val_index++;
        var_index = 0;

    }

    TH1F *hist_var[nvar];
    TH1F *hist_error[nvar];
    double graph_values[nvar];
    double graph_errors[nvar];
    double index[nvar];
    double index_error[nvar] ;

    TFile f("histo.root","RECREATE");

    for(size_t i = 0; i < nvar; i++){

        index[i] = i;
        std::sort(vec_val[i], vec_val[i] + nsamples);
        std::sort(vec_error[i], vec_error[i] + nsamples);

        hist_var[i] = new TH1F( (vec_name[i] + "_value").c_str(),(vec_name[i] + "_value").c_str(),1,
                               vec_val[i][0],vec_val[i][nsamples]);
        hist_error[i] = new TH1F( (vec_name[i] + "_error").c_str(),(vec_name[i] + "_error").c_str(),1,
                               vec_error[i][0],vec_error[i][nsamples]);

        for(size_t j=0; j < nsamples; j++){
            printf("(%d,%d) = [%lg,%lg] \n",i,j,vec_val[i][j],vec_error[i][j]);
            hist_var[i]->Fill(vec_val[i][j]);
            hist_error[i]->Fill(vec_error[i][j]);
        }

        hist_var[i]->Write();
        hist_error[i]->Write();
       
        graph_values[i] = hist_var[i]->GetMean();
        graph_errors[i] = hist_error[i]->GetMean();
    }

    f.Close();
    
    TGraphErrors *graph_result = new TGraphErrors(nvar,index,graph_values,index_error,graph_errors);
    TGraphErrors *graph_initial = new TGraphErrors(nvar,index,vec_val_i,index_error,vec_error_i);
    TCanvas foo;
    graph_result->SetMarkerStyle(kCircle);
    graph_result->SetMarkerColor(kRed);
    graph_result->Draw("AP");
    graph_initial->SetMarkerStyle(20);
    graph_result->SetMarkerColor(kBlue);
    graph_initial->Draw("Psame");

    foo.SaveAs("test.png");

}*/


void genfitplot(int nsamples,int nvar){

    double vec_inicial[nvar];
    
   //getting inicial parameters

    printf("Initial Parameters");
    int index = 0;
    double inicial_var = 0;
    double inicial_error =0;
    std::string inicial_name;
    ifstream r_inicial("Fit/fit_parameters_inicial.txt");
	
    TFile f("Fit/results.root","RECREATE");
    TTree *t = new TTree("t","Variables"); 

    double val[nvar];
    double val_error[nvar];
    //std::fill(val,val+nvar,0);
    //std::fill(val_error,val_error+nvar,0);
  
    while(r_inicial >> inicial_name >> inicial_var >> inicial_error){
            vec_inicial[index] = inicial_var;
	    t->Branch(inicial_name.c_str(),val+index);
	    t->Branch((inicial_name+"_error").c_str(),val_error+index);
            index++;
    }

    //getting fit parameters
   std::string name;
 
       for(int i = 0; i < nsamples; i++){
        
            name = fmt::format("Fit/fit_parameters_{0}.txt",i);
            ifstream r_fit(name.c_str());

            if(r_fit.is_open()){

                int index_fit = 0;
                double fit_var = 0;
		double fit_error=0;
		std::string fit_name;

                while(r_fit >> fit_name >> fit_var >> fit_error){

                    val[index_fit] = fit_var;
		    val_error[index_fit] = fit_error;

		    cout << val[index_fit] << " " << val_error[index_fit] << endl;	     
		    index_fit++;
                }

		 t->Fill();

            

                r_fit.close();
            }else{
                printf("Error \n");
            }
       }

   t->Write();
   f.Close();
}



int main(int argc, char **argv){

    //TROOT* groot = ROOT::GetROOT();
    //groot->SetBatch(1);
    //groot->ProcessLine( "gErrorIgnoreLevel = 1001;");
    int sample_number = 0 ;

    GooFit::Application app{"D2PPP",argc,argv};
    app.add_option("-i,--int", sample_number, "sample number", true);

    
    size_t  nevents = 100000;

    auto gen = app.add_subcommand("gen","generate toy data");
    gen->add_option("-e,--events",nevents,"The number of events to generate");

    int nbins = 1500;
    auto toyfit = app.add_subcommand("fit","fit toy data/toyMC");
    toyfit->add_option("-b,--bins",nbins,"The number of bins for normalization");

    auto plot = app.add_subcommand("plot","plot signal");

    int ns = 10;
    int nv = 10;

    auto gfplot = app.add_subcommand("gfplot","genfit plot");
    gfplot->add_option("-n,--ns",ns,"n samples");
    gfplot->add_option("-v,--nv",nv,"n variables");
    app.require_subcommand();

    int bins_i = 1000;
    int bins_f = 3000;
    auto normStudies = app.add_subcommand("normstudies","norm studies");
    normStudies->add_option("-i,--bi",bins_i,"initial number of bins");
    normStudies->add_option("-f,--bf",bins_f,"final number of bins");


    GOOFIT_PARSE(app);

    std::string name = fmt::format("MC/MC_Toy_{0}",sample_number);

    /// Make the plot directory if it does not exist
    std::string command = "mkdir -p plots";
    if(system(command.c_str()) != 0)
        throw GooFit::GeneralError("Making `plots` directory failed");

    
    std::string command2 = "mkdir -p MC";
    if(system(command2.c_str()) != 0)
        throw GooFit::GeneralError("Making `MC` directory failed");


    if(*gen){
        CLI::AutoTimer timer("MC Generation");
        cout << name << endl;

        runtoygen(name,nevents);
    }

    if(*toyfit){
        CLI::AutoTimer timer("FIT");
        runtoyfit(name,sample_number,nbins);
    }

    
    if(*plot){
        CLI::AutoTimer timer("Plot model and data");
        runMakeToyDalitzPdfPlots("MC/MC_Toy_0.txt");
    }
	
     if(*gfplot){
        CLI::AutoTimer timer("plot genfit result");
        genfitplot(ns,nv);
    }

     if(*normStudies){
        CLI::AutoTimer timer("norm studies");
        normStudy(name,sample_number,bins_i,bins_f);
    }


}

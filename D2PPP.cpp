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

Observable s12("s12",s12_min,2.85); //s12^{2}
Observable s13("s13",s13_min,2.85);
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
    {
        TCanvas foo;
        auto th1 = dp.make2D();
        th1->Rebin2D(5,5);
        th1->GetXaxis()->SetTitle("#pi^{-}#pi^{+} [Gev/c^{2}]");
        th1->GetYaxis()->SetTitle("#pi^{-}#pi^{+} [Gev/c^{2}]");
        th1->SetStats(0);
        th1->Draw("COLZ");
        foo.SaveAs("plots/PDF.png");
        std::cout << "PDF plotted" << '\n';
    }

        dp.fillDataSetMC(*Data,nEvents);
        TH2F th2("toyData", "", 200, s12.getLowerLimit(), s12.getUpperLimit(), 200, s13.getLowerLimit(),
                         s13.getUpperLimit());
    th2.GetXaxis()->SetTitle("#pi^{-}#pi^{+} [Gev/c^{2}]");
    th2.GetYaxis()->SetTitle("#pi^{-}#pi^{+} [Gev/c^{2}]");

    {
        ofstream w(name);

            for (size_t i = 0; i < Data->getNumEvents(); i++) {
                Data->loadEvent(i);
                th2.Fill(s12, s13);
                w << i << "\t" << std::setprecision(6) << s12.getValue() << "\t" << s13.getValue() << '\n';

            }

            std::cout << "nEvents generated = " << Data->getNumEvents() << '\n';

        w.close();

    }

    TCanvas foo;
    th2.Draw("COLZ");
    th2.SetStats(0);
    foo.SaveAs("plots/toyData.png");

    std::cout << "toy data plotted" << '\n';
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

    Variable swave_amp_real("swave_amp_real", 1.0,-100.0,+100.0);
    Variable swave_amp_imag("swave_amp_imag", 0.0,-100.0,+100.0);

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

     /*
     TTree *s = (TTree*)f->Get(tree_name.c_str());

     double _s12, _s13;

    s->SetBranchAddress("s12_pipi_DTF",&_s12);
    s->SetBranchAddress("s13_pipi_DTF",&_s13);
   

    for(size_t i = 0; i < 100000 ; i++){

            s->GetEntry(i);
            s12.setValue(_s12);
            s13.setValue(_s13);
            eventNumber.setValue(i);
            binBkgData->addEvent();

    }

for(int i = 0; i < 100000 ; i++){
       cout << i << "\t" <<  s12.getValue() << '\t' << s13.getValue() << '\t' << binBkgData->getBinContent(i) << endl;
    } */

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
    double f0_980_GKK     = 4.21*f0_980_GPP;
    double f0_980_amp     = 4.0;
    double f0_980_phase    = 0.0;

    double f0_1370_MASS     = 1.400;
    double f0_1370_WIDTH    = 0.3;
    double f0_1370_amp      = 1.0;
    double f0_1370_phase    = 0;

    double rho_MASS     = 0.77526;
    double rho_WIDTH    = 0.1478;
    double rho_amp      = 0.7;
    double rho_phase    = 0;

    double rho_1450_MASS     = 1.32635;
    double rho_1450_WIDTH    = 0.32413;
    double rho_1450_amp      = 0.5;
    double rho_1450_phase    = 0;

    double omega_MASS   = 0.78265;
    double omega_WIDTH  = 0.00849;
    double omega_amp    = 0.5;
    double omega_phase  = 0.0;

    double f2_MASS     = 1.2755;
    double f2_WIDTH    = 0.1867;
    double f2_amp      = 0.5;
    double f2_phase    = 0;

    double sigma_MASS  = 0.480;
    double sigma_WIDTH = 0.350;
    double sigma_amp   = 0.5;
    double sigma_phase = 0.0;



    //rho(770)
    Variable v_rho_Mass("rho_MASS",rho_MASS,0.01,rho_MASS*0.95,rho_MASS*1.05);
    Variable v_rho_Width("rho_WIDTH",rho_WIDTH,0.01,rho_WIDTH*0.95,rho_WIDTH*1.05);
    Variable v_rho_amp_real("rho_amp_real",rho_amp*cos(rho_phase), 0.001, -100.0, +100.0);
    Variable v_rho_amp_img("rho_amp_img",rho_amp*sin(rho_phase), 0.001, -100.0, +100.0);

    //rho(1450)
    Variable v_rho_1450_Mass("rho_1450_MASS",rho_1450_MASS,0.01,rho_1450_MASS*0.95,rho_1450_MASS*1.05);
    Variable v_rho_1450_Width("rho_1450_WIDTH",rho_1450_WIDTH,0.01,rho_1450_WIDTH*0.95,rho_1450_WIDTH*1.05);
    Variable v_rho_1450_amp_real("rho_1450_amp_real",rho_1450_amp*cos(rho_1450_phase), 0.001, -100.0, +100.0);
    Variable v_rho_1450_amp_img("rho_1450_amp_img",rho_1450_amp*sin(rho_1450_phase), 0.001, -100.0, +100.0);


    //omega(782)
    Variable v_omega_Mass("omega_MASS",omega_MASS);
    Variable v_omega_Width("omega_WIDTH",omega_WIDTH);
    Variable v_omega_amp_real("omega_amp_real",omega_amp*cos(omega_phase), 0.001, -100.0, +100.0);
    Variable v_omega_amp_img("omega_amp_img",omega_amp*sin(omega_phase), 0.001, -100.0, +100.0);

    //f2(1270)
    Variable v_f2_Mass("f2_MASS",f2_MASS,0.01,f2_MASS*0.95,f2_MASS*1.05);
    Variable v_f2_Width("f2_WIDTH",f2_WIDTH,0.01,f2_WIDTH*0.95,f2_WIDTH*1.05);
    Variable v_f2_amp_real("f2_amp_real",f2_amp*cos(f2_phase), 0.001, -100.0, +100.0);
    Variable v_f2_amp_img("f2_amp_img",f2_amp*sin(f2_phase), 0.001, -100.0, +100.0);

    //sigma(480)
    Variable v_sigma_Mass("sigma_MASS",sigma_MASS,0.01,sigma_MASS*0.95,sigma_MASS*1.05);
    Variable v_sigma_Width("sigma_WIDTH",sigma_WIDTH,0.01,sigma_WIDTH*0.95,sigma_WIDTH*1.05);
    Variable v_sigma_amp_real("sigma_amp_real",sigma_amp*cos(sigma_phase), 0.001, -100.0, +100.0);
    Variable v_sigma_amp_img("sigma_amp_img",sigma_amp*sin(sigma_phase), 0.001, -100.0, +100.0);

    //f0(980)
    Variable v_f0_980_Mass("f0_980_MASS",f0_980_MASS,0.01,f0_980_MASS*0.95,f0_980_MASS*1.05);
    Variable v_f0_980_GPP("f0_980_GPP",f0_980_GPP,0.01,f0_980_GPP*0.95,f0_980_GPP*1.05);
    Variable v_f0_980_GKK("f0_980_GKK",f0_980_GKK,0.01,f0_980_GKK*0.95,f0_980_GKK*1.05);
    Variable v_f0_980_amp_real("f0_980_amp_real",f0_980_amp*cos(f0_980_phase), 0.001, -100.0, +100.0);
    Variable v_f0_980_amp_img("f0_980_amp_img",f0_980_amp*sin(f0_980_phase), 0.001, -100.0, +100.0);

    
    //f0(1370)
    Variable v_f0_1370_Mass("f0_1370_MASS",f0_1370_MASS,0.01,f0_1370_MASS*0.95,f0_1370_MASS*1.05);
    Variable v_f0_1370_Width("f0_1370_Width",f0_1370_WIDTH,0.01,f0_1370_WIDTH*0.95,f0_1370_WIDTH*1.05);
    Variable v_f0_1370_amp_real("f0_amp_1500_real",f0_1370_amp*cos(f0_1370_phase), 0.001, -100.0, +100.0);
    Variable v_f0_1370_amp_img("f0_1370_amp_img",f0_1370_amp*sin(f0_1370_phase), 0.001, -100.0, +100.0);


    //NR
    Variable nonr_amp_real("nonr_amp_real", 1.0, 0.001, -100.0, +100.0);
    Variable nonr_amp_imag("nonr_amp_imag", 0.0, 0.001, -100.0, +100.0);

    //Masses and Widths fixed

    v_rho_Mass.setFixed(true);
    v_rho_Width.setFixed(true);

    v_rho_1450_Mass.setFixed(true);
    v_rho_1450_Width.setFixed(true);

    v_omega_Mass.setFixed(true);
    v_omega_Width.setFixed(true);

    v_f2_Mass.setFixed(true);
    v_f2_Width.setFixed(true);

    v_sigma_Mass.setFixed(true);
    v_sigma_Width.setFixed(true);

    v_f0_980_Mass.setFixed(true);
    v_f0_980_GKK.setFixed(true);
    v_f0_980_GPP.setFixed(true);


    v_f0_1370_Mass.setFixed(true);
    v_f0_1370_Width.setFixed(true);


    //setting resonances
    ResonancePdf* rho_12 = new Resonances::GS("rho",v_rho_amp_real,v_rho_amp_img,v_rho_Mass,v_rho_Width,(unsigned int)1,PAIR_12,true);
       
    ResonancePdf* rho_1450_12 = new Resonances::GS("rho_1450",v_rho_1450_amp_real,v_rho_1450_amp_img,v_rho_1450_Mass,v_rho_1450_Width,(unsigned int)1,PAIR_12,true);
    
    ResonancePdf* omega_12 = new Resonances::RBW("omega",v_omega_amp_real,v_omega_amp_img,v_omega_Mass,v_omega_Width,1,PAIR_12,true);
    
    ResonancePdf* f2_12 = new Resonances::RBW("f2",v_f2_amp_real,v_f2_amp_img,v_f2_Mass,v_f2_Width,2,PAIR_12,true);
       
    ResonancePdf* sigma_12 = new Resonances::RBW("sigma",v_sigma_amp_real,v_sigma_amp_img,v_sigma_Mass,v_sigma_Width,(unsigned int)0,PAIR_12,true);
   
    ResonancePdf* f0_1370_12 = new Resonances::RBW("f0_1370_12",v_f0_1370_amp_real,v_f0_1370_amp_img,v_f0_1370_Mass,v_f0_1370_Width,(unsigned int)0,PAIR_12,true);
      
    ResonancePdf* f0_980_12 = new Resonances::FLATTE("f0_980",v_f0_980_amp_real,v_f0_980_amp_img,v_f0_980_Mass,v_f0_980_GPP,v_f0_980_GKK,PAIR_12,true);
      
    ResonancePdf *nonr = new Resonances::NonRes("nonr", nonr_amp_real, nonr_amp_imag);

    //MIPWA
    ResonancePdf *swave_12 = loadPWAResonance(pwa_file, true);

    //Pushing Resonances 

    dtoppp.resonances.push_back(rho_12);
    //dtoppp.resonances.push_back(rho_1450_12);
    //dtoppp.resonances.push_back(omega_12);
    //dtoppp.resonances.push_back(f2_12);
    //dtoppp.resonances.push_back(sigma_12);
    dtoppp.resonances.push_back(f0_980_12);
    //dtoppp.resonances.push_back(f0_1370_12);
    //dtoppp.resonances.push_back(nonr);
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

void saveParameters(const std::vector<double> &param, std::string file){

    ofstream wr(file.c_str());

    for(int i = 0; i < param.size(); i++){
        wr << i <<'\t'<< param[i] << endl;
    }
    
    wr.close();

}



void runtoyfit(std::string name, int sample_number) {

    s12.setNumBins(1500);
    s13.setNumBins(1500);

    getdata(name);

    if(toyOn){
        GOOFIT_INFO("Using toyDATA");
    }else{
        GOOFIT_INFO("Using DATA: {}",data_name );
    }

    GOOFIT_INFO("Number of Events in dataset: {}", Data->getNumEvents());
 
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
    }

    AddPdf* overallPdf = new AddPdf("overallPdf",weights,comps);
    overallPdf->setData(Data);
    signaldalitz->setDataSize(Data->getNumEvents());

    FitManagerMinuit2 fitter(overallPdf);
    fitter.setVerbosity(3);

    std::string command = "mkdir -p Fit";
    if(system(command.c_str()) != 0)
        throw GooFit::GeneralError("Making `Fit` directory failed");


    const std::vector<double> params_iniciais  = fitter.getParams()->make_minuit_vector();

    string input_name = fmt::format("Fit/fit_parameters_inicial.txt");
    
    saveParameters(params_iniciais, input_name);

    
    auto func_min = fitter.fit();

    const std::vector<double> params  = fitter.getParams()->make_minuit_vector();

    string output_name = fmt::format("Fit/fit_parameters_{0}.txt",sample_number);
    
    saveParameters(params, output_name);

    /*const std::vector<std::vector<double>> params2  = fitter.getParams()->get_recorded();
    for(int i = 0; i < params.size(); i++){
        for(int k = 0; k < params[0].size(); k++){
            printf("[%d][%d] = %lg \n",i,k,params2[i][k]);
        }
    }*/

    }

void genfitplot(int nsamples,int nvar){

    double vec_inicial[nvar];
    
    //creating histograms
    printf("creating histograms");
    TH1D* hist[nvar];
    string var_name;
    for(int i=0; i < nvar; i++){
        var_name = fmt::format("variable_{0}",i);
        hist[i] = new TH1D(var_name.c_str(),var_name.c_str(),20,-1,1);
    }

    //getting inicial parameters
    int index = 0;
    double inicial_var = 0;
    ifstream r_inicial("Fit/fit_parameters_inicial.txt");
    while(r_inicial >> index >> inicial_var){

        vec_inicial[index] = inicial_var;
        printf("%d --> %lg \n",index,inicial_var);

    }

    //getting fit parameters
    double test;
    string name;

       for(int i = 0; i < nsamples; i++){
        
            name = fmt::format("Fit/fit_parameters_{0}.txt",i);
            ifstream r_fit(name.c_str());

            if(r_fit.is_open()){

                printf("file %d is open!",i);

                int index_fit = 0;
                double fit_var = 0;

                while(r_fit >> index_fit >> fit_var){

                    test = vec_inicial[index_fit] - fit_var;
                    printf("[%d] = %lg",index_fit,test);
                    hist[index_fit]->Fill(test);

                }
            

                r_fit.close();
            }else{
                printf("Error \n");
            }
       }

    TFile f("Fit/results.root","RECREATE");

       for(int i = 0; i < nvar ; i++){
            hist[i]->Write();       
       } 
    f.Close();
}



int main(int argc, char **argv){

    int sample_number = 0 ;

    GooFit::Application app{"D2PPP",argc,argv};
    app.add_option("-i,--int", sample_number, "sample number", true);

    
    size_t  nevents = 100000;

    auto gen = app.add_subcommand("gen","generate toy data");
    gen->add_option("-e,--events",nevents,"The number of events to generate");

    auto toyfit = app.add_subcommand("fit","fit toy data/toyMC");

    auto plot = app.add_subcommand("plot","plot signal");

    int ns = 10;
    auto gfplot = app.add_subcommand("gfplot","genfit plot");
    gfplot->add_option("-n,--ns",ns,"n samples");

    app.require_subcommand();
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
        runtoyfit(name,sample_number);
    }

    
    if(*plot){
        CLI::AutoTimer timer("FIT");
        runMakeToyDalitzPdfPlots("D2PPP_toy.txt");
    }

     if(*gfplot){
        CLI::AutoTimer timer("FIT");
        genfitplot(ns,12);
    }

}

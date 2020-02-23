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

#include <TGraphErrors.h>



#define torad(x)(x*M_PI/180.)

using namespace std;
using namespace GooFit;
using namespace ROOT;


//Globals

double pi_MASS  = 0.13957018; 
double D_MASS   = 1.96834; 

double Mother_MASS = D_MASS;
double d1_MASS  = pi_MASS; 
double d2_MASS  = pi_MASS;
double d3_MASS  = pi_MASS;

const int bins = 400;

fptype s12_min = POW2(d1_MASS  + d2_MASS)*0.9;
fptype s12_max = POW2(D_MASS   - d2_MASS)*1.1;
fptype s13_min = POW2(d1_MASS  + d3_MASS)*0.9;
fptype s13_max = POW2(D_MASS   - d2_MASS)*1.1;


Observable s12("s12",s12_min,s12_max); 
Observable s13("s13",s13_min,s13_max);
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

// PWA INPUT FILE NAME
const string pwa_file = "files/PWACOEFS.txt";

//Declaring Functions 
void saveParameters(std::string file, const std::vector<Variable> &param,   bool isValid,  double fcn,  std::vector<std::vector<fptype>> ff );
    
ResonancePdf *loadPWAResonance(const string fname = pwa_file, bool fixAmp = true) {

    std::ifstream reader;
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

        emag = e2*cos(e3);
        ephs = e2*sin(e3);
	    
        Variable va(fmt::format("pwa_coef_{}_mag", i), emag,0.01,0,0);//-100.,+100.);
        Variable vp(fmt::format("pwa_coef_{}_phase", i), ephs,0.01,0,0);//-100.,+100.);
        //va.setFixed(true);
        //vp.setFixed(true);		
        pwa_coefs_amp.push_back(va);
        pwa_coefs_phs.push_back(vp);
        i++;
    }


    Variable swave_amp_real("swave_real_coef", 1.0,0.01,-100.,+100.);
    Variable swave_amp_imag("swave_imag_coef", 0.0,0.01,-100.,+100.);

    if(fixAmp) {
        swave_amp_real.setValue(1.);
        swave_amp_imag.setValue(0.);
        swave_amp_real.setFixed(true);
        swave_amp_imag.setFixed(true);
    }

    ResonancePdf *swave_12 = new Resonances::Spline("swave_12", swave_amp_real, swave_amp_imag, HH_bin_limits, pwa_coefs_amp, pwa_coefs_phs, PAIR_12, true);

    return swave_12;
} 

GooPdf* makeEfficiencyPdf(std::string efffile, std::string effhist, Observable s12, Observable s13) {

    vector<Observable> lvars;
    lvars.push_back(s12);
    lvars.push_back(s13);
    BinnedDataSet *binEffData = new BinnedDataSet(lvars);
    
    TFile *f     = TFile::Open(efffile.c_str());
    auto effHistogram = (TH2F *)f->Get(effhist.c_str());
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

    Variable *effSmoothing = new Variable("effSmoothing_eff", 0.0);
    SmoothHistogramPdf *ret = new SmoothHistogramPdf("efficiency", binEffData, *effSmoothing);
    return ret;
}

GooPdf* makeBackgroundPdf(std::string bkgfile, std::string bkghist, Observable s12, Observable s13) {

   BinnedDataSet *binBkgData = new BinnedDataSet({s12, s13});
    
    TFile *f = new TFile(bkgfile.c_str());
    auto bkgHistogram = (TH2F*)f->Get(bkghist.c_str());
   
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

    Variable *effSmoothing  = new Variable("Smoothing_bkg",0.);
    SmoothHistogramPdf *ret = new SmoothHistogramPdf("background", binBkgData, *effSmoothing);

    return ret;
}

GooPdf *makeDstar_veto() {
   
    VetoInfo Dstar_veto12(Variable("Dstar_veto12_min", 2.85), Variable("Dstar_veto12_max", s12_max), PAIR_12);
    VetoInfo Dstar_veto13(Variable("Dstar_veto13_min", 2.85), Variable("Dstar_veto13_max", s12_max), PAIR_13);
    VetoInfo Fiducial_veto12(Variable("Fiducial_veto12_min", 0.), Variable("Fiducial_veto12_max", s12.getLowerLimit()), PAIR_12);
    VetoInfo Fiducial_veto13(Variable("DFiducial_veto13_min", 0.), Variable("Fiducial_veto13_max", s13.getLowerLimit()), PAIR_13);
    
    vector<VetoInfo> vetos;
    //vetos.push_back(Dstar_veto12);
    //vetos.push_back(Dstar_veto13);
    vetos.push_back(Fiducial_veto12);
    vetos.push_back(Fiducial_veto13);

    DalitzVetoPdf* Dstar_veto = new DalitzVetoPdf("Dstar_veto", s12, s13, 
        Variable("Mother_Mass",D_MASS), Variable("Daughter1_Mass",d1_MASS), 
        Variable("Daughter2_Mass",d2_MASS), Variable("Daughter3_Mass",d3_MASS), 
        vetos);

    return Dstar_veto;
}

DalitzPlotPdf* makesignalpdf( Observable s12, Observable s13, EventNumber eventNumber, GooPdf* eff = 0){

    DecayInfo3 dtoppp;
    dtoppp.motherMass   = D_MASS;
    dtoppp.daug1Mass    = d1_MASS;
    dtoppp.daug2Mass    = d2_MASS;
    dtoppp.daug3Mass    = d3_MASS;
    dtoppp.meson_radius = 1.5;

    //Mass and width

    double f0_980_MASS    = 0.960;
    double f0_980_GPP     = 0.02;
    double f0_980_GKK     = 4.5*0.02;
    double f0_980_WIDTH   = 0.04;
    double f0_980_amp     = 1.0;
    double f0_980_img     = 0.0;
  
    
    double f0_1370_MASS  = 1.370;
    double f0_1370_WIDTH = .3;
    double f0_1370_amp   = 0.75*cos(torad(198));
    double f0_1370_img = 0.75*sin(torad(198));

    double f0_1500_MASS  = 1.505;
    double f0_1500_WIDTH = .109;
    double f0_1500_amp   = 1.;
    double f0_1500_img = 0.;

    double omega_MASS   = 0.78265;
    double omega_WIDTH  = 0.00849;
    double omega_amp    = 1.;
    double omega_img  =   0.;

    double rho770_MASS   = .77549;
    double rho770_WIDTH  = .1491;
    //E791
    //double rho770_amp    = 0.32*cos(torad(109));
    //double rho770_img  = 0.32*sin(torad(109));
    //Babar
    double rho770_amp    = 0.19*cos(1.1);
    double rho770_img  = 0.19*sin(1.1);
  


    double rho1450_MASS   = 1.465;//+ 0*0.025;
    double rho1450_WIDTH  = 0.4 ;//+ 0*0.06;//0.4;
    //E791
    //double rho1450_amp    = 0.28*cos(torad(162));
    //double rho1450_img  = 0.28*sin(torad(162));
    //Babar
    double rho1450_amp    = 1.2*cos(4.1);
    double rho1450_img  = 1.2*sin(4.1);
 

    
    
    double f2_1270_MASS     = 1.2751;
    double f2_1270_WIDTH    = 0.1851;
    //E791
    //double f2_1270_amp      = 0.59*cos(torad(133));
    //double f2_1270_img    = 0.59*sin(torad(133));
    //Babar
    double f2_1270_amp      = 1.;
    double f2_1270_img    = 0.;

    double f2_1525_MASS     = 1.525;
    double f2_1525_WIDTH    = 0.073;
    double f2_1525_amp      = 1.;
    double f2_1525_img    = 0.;
    

    //omega(782)
    Variable v_omega_Mass("omega_MASS",omega_MASS,0.01,0,0);
    Variable v_omega_Width("omega_WIDTH",omega_WIDTH,0.01,0,0);
    Variable v_omega_real("omega_real",omega_amp, 0.01,0,0);
    Variable v_omega_img("omega_img",omega_img, 0.01,0,0);

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
    Variable v_f2_1525_Mass("f2_1525_MASS",f2_1525_MASS,0.01,0,0);
    Variable v_f2_1525_Width("f2_1525_WIDTH",f2_1525_WIDTH,0.01,0,0);
    Variable v_f2_1525_real("f2_1525_real",f2_1525_amp, 0.01,0,0);
    Variable v_f2_1525_img("f2_1525_img",f2_1525_img, 0.01,0,0);

    //f0(980)
    Variable v_f0_980_Mass("f0_980_MASS",f0_980_MASS,.02,0,0);
    Variable v_f0_980_GPP("f0_980_GPP",f0_980_GPP,0.01,0.,0);
    Variable v_f0_980_GKK("f0_980_GKK",f0_980_GKK,0.01,0.,0.);
    Variable v_f0_980_Width("f0_980_WIDTH",f0_980_WIDTH,0.01,0,0);
    Variable v_f0_980_real("f0_980_real",f0_980_amp, 0.01,0,0);
    Variable v_f0_980_img("f0_980_img",f0_980_img, 0.01,0,0);

    //v_f0_980_real.setFixed(true);
    //v_f0_980_img.setFixed(true);
    v_f2_1270_real.setFixed(true);
    v_f2_1270_img.setFixed(true);

    
    //f0(1370)
    Variable v_f0_1370_Mass("f0_1370_MASS",f0_1370_MASS,0.025,0,0);
    Variable v_f0_1370_Width("f0_1370_Width",f0_1370_WIDTH,0.035,0,0);
    Variable v_f0_1370_real("f0_1370_real",f0_1370_amp, 0.01,0,0);
    Variable v_f0_1370_img("f0_1370_img",f0_1370_img, 0.01, 0,0);

    
    //f0(1500)
    Variable v_f0_1500_Mass("f0_1500_MASS",f0_1500_MASS,0.06,0,0);
    Variable v_f0_1500_Width("f0_1500_Width",f0_1500_WIDTH,0.07,0,0);
    Variable v_f0_1500_real("f0_1500_real",f0_1500_amp, 0.01,0,0);
    Variable v_f0_1500_img("f0_1500_img",f0_1500_img, 0.01, 0,0);

    //NR
    Variable nonr_real("nonr_real",0.09*cos(torad(181)), 0.01,0,0);
    Variable nonr_imag("nonr_imag",0.09*sin(torad(181)), 0.01,0,0);
    //Variable nonr_real("nonr_real",1., 0.01,0,0);
     // Variable nonr_imag("nonr_imag",0., 0.01,0,0);


    Variable be_real("be_real",1., 0.01,0,0);
    Variable be_imag("be_imag", 0., 0.01,0,0);
    Variable be_coef("be_coef", 1.9,0.01,0,0);
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

    v_f0_980_Mass.setFixed(true);
    v_f0_980_GKK.setFixed(true);
    v_f0_980_GPP.setFixed(true);
    v_f0_980_Width.setFixed(true);

    v_f0_1370_Mass.setFixed(true);
    v_f0_1370_Width.setFixed(true);

    v_f0_1500_Mass.setFixed(true);
    v_f0_1500_Width.setFixed(true);

    be_coef.setFixed(true);
        
   //setting resonances
    //ResonancePdf* sigma_12 = new Resonances::POLE("sigma",Variable("v_sigma_real",1.),Variable("v_sigma_img",0.),Variable("v_sigma__pole_real",0.47),Variable("v_sigma_pole_img",0.22),0,PAIR_12,true);
    
    ResonancePdf* omega_12 = new Resonances::RBW("omega",v_omega_real,v_omega_img,v_omega_Mass,v_omega_Width,1,PAIR_12,true);
    
    ResonancePdf* rho770_12 = new Resonances::RBW("rho770",v_rho770_real,v_rho770_img,v_rho770_Mass,v_rho770_Width,1,PAIR_12,true);
    
    ResonancePdf* rho1450_12 = new Resonances::RBW("rho1450",v_rho1450_real,v_rho1450_img,v_rho1450_Mass,v_rho1450_Width,1,PAIR_12,true);
    
    ResonancePdf* f2_1270_12 = new Resonances::RBW("f2",v_f2_1270_real,v_f2_1270_img,v_f2_1270_Mass,v_f2_1270_Width,2,PAIR_12,true);

    ResonancePdf* f2_1525_12 = new Resonances::RBW("f2",v_f2_1525_real,v_f2_1525_img,v_f2_1525_Mass,v_f2_1525_Width,2,PAIR_12,true);
  
    ResonancePdf* f0_980_12 = new Resonances::FLATTE("f0_980",v_f0_980_real,v_f0_980_img,v_f0_980_Mass,v_f0_980_GPP,v_f0_980_GKK,PAIR_12,true);

    //ResonancePdf* f0_980_12 = new Resonances::RBW("f0_980",v_f0_980_real,v_f0_980_img,v_f0_980_Mass,v_f0_980_Width,(unsigned int)0,PAIR_12,true);
	
    ResonancePdf* f0_1370_12 = new Resonances::RBW("f0_1370_12",v_f0_1370_real,v_f0_1370_img,v_f0_1370_Mass,v_f0_1370_Width,(unsigned int)0,PAIR_12,true);

    ResonancePdf* f0_1500_12 = new Resonances::RBW("f0_1500_12",v_f0_1500_real,v_f0_1500_img,v_f0_1500_Mass,v_f0_1500_Width,(unsigned int)0,PAIR_12,true);  

    ResonancePdf *nonr = new Resonances::NonRes("nonr", nonr_real, nonr_imag);

    ResonancePdf *be   = new Resonances::BoseEinstein("be",be_real,be_imag,be_coef);

    //MIPWA
    ResonancePdf *swave_12 = loadPWAResonance(pwa_file, true);

    //swave_12->recalculateCache();

    //Pushing Resonances 
    //dtoppp.resonances.push_back(sigma_12);
    dtoppp.resonances.push_back(omega_12); 
    dtoppp.resonances.push_back(rho770_12); 
    dtoppp.resonances.push_back(rho1450_12);
    dtoppp.resonances.push_back(f2_1270_12);
    //dtoppp.resonances.push_back(f2_1525_12);
    //dtoppp.resonances.push_back(f0_980_12);
    //dtoppp.resonances.push_back(f0_1500_12);
    //dtoppp.resonances.push_back(f0_1370_12);
    //dtoppp.resonances.push_back(nonr);
    //sdtoppp.resonances.push_back(be);
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

void saveParameters(const std::vector<ROOT::Minuit2::MinuitParameter> &param, bool status, double fcn, std::string file){

    std::vector<std::string> v1;
    std::vector<fptype> v2;
    std::vector<fptype> v3;

    std::ofstream output_file(file.c_str(),std::ofstream::out);

    for(size_t i = 0 ; i < param.size() ; i++){
	        v1.push_back(param[i].GetName());
            v2.push_back(param[i].Value());
            v3.push_back(param[i].Error());
     }

    for(size_t i = 0; i < v1.size(); i++) {
        output_file << v1[i] << "\t" << std::fixed << std::setprecision(6) << v2[i] << "\t" << v3[i] << std::endl;
    }


   output_file << "Min FCN \t" << std::fixed << std::setprecision(6) << fcn << std:: endl;
   output_file << "Status \t" << status << std::endl;


   output_file.close();
}

void getToyData(std::string toyFileName, GooFit::Application &app, DataSet &data, bool toy) {

    toyFileName = app.get_filename(toyFileName, "MC/");

    auto obs               = data.getObservables();
    Observable s12         = obs.at(0);
    Observable s13         = obs.at(1);
    Observable eventNumber = obs.at(2);

    auto openRoot   = new TFile(toyFileName.c_str());
    auto tree       = (TTree*)openRoot->Get("DecayTree");
    auto s12_val(0.);
    auto s13_val(0.);

    if(toy){
        tree->SetBranchAddress("s12",&s12_val);
        tree->SetBranchAddress("s13",&s13_val);
    }else{
        
        tree->SetBranchAddress("s12_pipi_DTF",&s12_val);
        tree->SetBranchAddress("s13_pipi_DTF",&s13_val);
    }

    size_t j = 0;
    for(size_t i = 0 ; i < tree->GetEntries(); i++){
        tree->GetEntry(i);
        s12.setValue(s12_val);
        s13.setValue(s13_val);
        eventNumber.setValue(data.getNumEvents());
        if((s12.getValue()<s12.getUpperLimit())
            &&(s13.getValue()<s13.getUpperLimit())
            &&(s12.getValue()>s12.getLowerLimit())
            &&(s13.getValue()>s13.getLowerLimit()))
        {
            data.addEvent();
            if(j<10) printf("[%d] = (%.2f , %.2f)\n",i,s12.getValue(),s13.getValue());
            j++;
            if(!toy && data.getNumEvents()==200000) break;
        }
    }


}

void to_root(UnbinnedDataSet& toyMC , std::string name ){

    auto obs               = toyMC.getObservables();
    Observable s12         = obs.at(0);
    Observable s13         = obs.at(1);
    Observable eventNumber = obs.at(2);

    double _s12, _s13,_s23;
    auto f = new TFile(name.c_str(),"recreate");
    auto t = new TTree("DecayTree","toyMC");
    auto b_s12 = t->Branch("s12",&_s12,"s12/D");
    auto b_s13 = t->Branch("s13",&_s13,"s13/D");
    auto b_s23 = t->Branch("s23",&_s23,"s23/D");

   
    for(int i = 0; i < toyMC.getNumEvents(); i++){
		toyMC.loadEvent(i);
		t->GetEntry(i);
		_s12 = s12.getValue();
		_s13 = s13.getValue();
		_s23 = POW2(Mother_MASS) + POW2(d1_MASS) + POW2(d2_MASS) + POW2(d3_MASS) - s12.getValue() - s13.getValue() ;
		t->Fill();
    }
	t->Write("",TObject::kOverwrite);
	f->Write();
	f->Close();
    std::cout << "------------------------------------------" << std::endl;
    std::cout << "toyMC --> " << name.c_str() << " was saved!" << std::endl;
    std::cout << "------------------------------------------" << std::endl;

}

int genFit(GooPdf *totalPdf,DalitzPlotPdf *signal, UnbinnedDataSet *data, int rank) {

    totalPdf->setData(data);
    signal->setDataSize(data->getNumEvents());
    FitManager datapdf(totalPdf);
    datapdf.setVerbosity(2);
    auto func_min = datapdf.fit();
   
    auto param = datapdf.getParams()->Parameters();
    
    auto output = fmt::format("Fit/fit_parameters_{0}.txt",rank);

    saveParameters(param, func_min.IsValid(), func_min.Fval(), output);	
    
    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "Sample " << rank << " fitted ! "		         << std::endl;
    std::cout << "nEvents --> " << data->getNumEvents()          << std::endl;
    std::cout << "minFCN --> "   << func_min.Fval()              << std::endl;
    std::cout << "AmpNorm --> "   << totalPdf->normalize()       << std::endl;
    std::cout << "---------------------------------------------" << std::endl;

    return datapdf;
}


int runFit(GooPdf *totalPdf,DalitzPlotPdf *signal, UnbinnedDataSet *data, std::string name) {

    auto obs               = data->getObservables();
    Observable s12         = obs.at(0);
    Observable s13         = obs.at(1);
    Observable eventNumber = obs.at(2);

    totalPdf->setData(data);
    signal->setDataSize(data->getNumEvents());
    FitManager datapdf(totalPdf);
    datapdf.setVerbosity(2);
    auto func_min = datapdf.fit();
   
    auto param = datapdf.getParams()->Parameters();

    auto output = name;

    saveParameters(param, func_min.IsValid(), func_min.Fval(), output);	
    
    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "Sample --> fitted ! "		                     << std::endl;
    std::cout << "nEvents --> " << data->getNumEvents()          << std::endl;
    std::cout << "minFCN --> "   << func_min.Fval()              << std::endl;
    std::cout << "AmpNorm --> "   << totalPdf->normalize()       << std::endl;
    std::cout << "Output fit saved in --> "   << output          << std::endl;
    std::cout << "---------------------------------------------" << std::endl;

    DalitzPlotter dplotter{totalPdf, signal};

    TCanvas foo;
    dplotter.Plot(D_MASS,d1_MASS,d2_MASS,d3_MASS,"s_{#pi^{-} #pi^{+}}","s_{#pi^{-} #pi^{+}}","s_{#pi^{-} #pi^{+}}","plots",*data);

    auto toyMC = new UnbinnedDataSet({s12,s13,eventNumber});
    dplotter.fillDataSetMC(*toyMC,10000000);

    size_t npar = 0; //number of free parameters
    for(size_t i = 0; i < param.size(); i++){
        if((!param[i].IsFixed()) && (!param[i].IsConst())){
            npar++;	
        }
    }
    dplotter.chi2(npar,"Ds3pi_bins.txt",0.05,1.95,0.3,3.4,10000000,*data,*toyMC);

    //saving s-wave qmpiwa
    size_t N = HH_bin_limits.size();
    double amp[N];
    double amp_er[N];
    double phase[N];
    double phase_er[N];
    double s[N];
    double s_er[N];
    
    std::ofstream wr("Fit/s_wave_fitted.txt");
    for(size_t i = 0; i < N ; i++){
        	s[i] = sqrt(HH_bin_limits[i]);
       		amp[i] = pwa_coefs_amp[i].getValue();
        	amp_er[i] = pwa_coefs_amp[i].getError();
        	phase[i] = pwa_coefs_phs[i].getValue();
        	phase_er[i] = pwa_coefs_phs[i].getError();
        	s_er[i] = 0.0;
            wr << std::fixed << std::setprecision(4) << sqrt(HH_bin_limits[i])<< '\t' << amp[i] << '\t' << phase[i] 
	 	<< '\t' << s_er[i] << '\t' << amp_er[i] << '\t' << phase_er[i] << std::endl;
    }
    wr.close();

    TGraphErrors *gr_real = new TGraphErrors(N,s,amp,s_er,amp_er);
    TGraphErrors *gr_img = new TGraphErrors(N,s,phase,s_er,phase_er);

    foo.Divide(2);
    foo.cd(1);
    gr_real->Draw("ALP");
    foo.cd(2);
    gr_img->Draw("ALP");
    foo.SaveAs("plots/Swave.png");

    return datapdf;
}

void genfitplot(int nsamples,int nvar,int n_res){
	printf("Not implemented yet! \n");
}


int main(int argc, char **argv){
    
    GooFit::Application app{"Genfit",argc,argv};
    
    std::string filename = "toyMC";
    app.add_option("-f,--filename,filename", filename, "File to read in", true)->check(GooFit::ExistingFile);
    
    bool make_toy;
    std::string Path = "MC";
    app.add_flag("-m,--make-toy", make_toy, "Make a toy instead of reading a file in");

    bool save_toy;
    app.add_flag("-s,--save-toy", save_toy, "Save toy in a ROOT file") ;

    GOOFIT_PARSE(app);

    // Make the mc directory if it does not exist
    std::string command = "mkdir -p MC";
    if(system(command.c_str()) != 0)
        throw GooFit::GeneralError("Making `MC` directory failed");

    // Make the Fit directory if it does not exist
    command = "mkdir -p Fit";
    if(system(command.c_str()) != 0)
        throw GooFit::GeneralError("Making `Fit` directory failed");

    // Make the Fit directory if it does not exist
    command = "mkdir -p plots";
    if(system(command.c_str()) != 0)
        throw GooFit::GeneralError("Making `plots` directory failed");

    Observable s12("s12", s12_min, s12_max);
    Observable s13("s13", s13_min, s13_max);
    EventNumber eventNumber("eventNumber");

    s12.setNumBins(bins);
    s13.setNumBins(bins);


    const string bkgfile = "../../../dados/bkg_histo_16.root";
    const string efffile = "../../../dados/eff_16.root";
    const string bkghist = "h_eff";
    const string effhist = "h_eff";

    UnbinnedDataSet data({s12, s13, eventNumber});

    auto efficiency = makeEfficiencyPdf(efffile,effhist,s12,s13);
    auto background = makeBackgroundPdf(bkgfile,bkghist,s12,s13);
    auto signal = makesignalpdf(s12, s13, eventNumber,efficiency);

    AddPdf *prodpdf = new AddPdf("prodpdf", Variable("frac",0.925), signal, background) ;

    auto localRank = app.get_localRank();

    DalitzPlotter dplotter{prodpdf, signal};

    try{

        if(make_toy) {
            dplotter.fillDataSetMC(data, 200000);
            std::cout << "----------------------------------------------------------" << std::endl;
            std::cout <<   data.getNumEvents() << " events in DataSet rank = " << localRank << " was generated!"  << std::endl;
            std::cout << "----------------------------------------------------------" << std::endl;
            auto fullName= fmt::format("MC/{0}_{1}.root",filename,localRank);
            if(save_toy) to_root(data,fullName); 
            return genFit(prodpdf,signal, &data, localRank);
        }else{
            auto fullName= fmt::format("{0}",filename);
            getToyData(fullName, app, data,false);
            auto output = "Fit/fit_result.txt";	
            return runFit(prodpdf,signal, &data, output);
        }

    }catch(const std::runtime_error &e) {
        std::cerr << e.what() << std::endl;
        return 7;
    }
       
}

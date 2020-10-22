// ROOT stuff
#include <TApplication.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TTree.h>
#include <TROOT.h>
#include <TMinuit.h>
#include <TNtuple.h>
#include <TComplex.h>
#include <TFile.h>
#include <TStyle.h>
#include <TH2Poly.h>
#include <TGraphErrors.h>

//Minuit
#include <Minuit2/MnStrategy.h>
#include <Minuit2/Minuit2Minimizer.h>

// System stuff
#include <CLI/Timer.hpp>
#include <fstream>
#include <string>
#include <time.h>

// GooFit stuff
#include <goofit/Application.h>
#include <goofit/BinnedDataSet.h>
#include <goofit/FitManager.h>
#include <goofit/PDFs/GooPdf.h>
#include <goofit/PDFs/basic/SmoothHistogramPdf.h>
#include <goofit/PDFs/basic/PolynomialPdf.h>
#include <goofit/PDFs/combine/AddPdf.h>
#include <goofit/PDFs/combine/ProdPdf.h>
#include <goofit/PDFs/physics/DalitzPlotPdf.h>
#include <goofit/PDFs/physics/DalitzVetoPdf.h>
#include <goofit/PDFs/physics/ResonancePdf.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/Variable.h>
#include <goofit/PDFs/physics/DalitzPlotter.h>
#include <goofit/FunctorWriter.h>

//Matrix 
#include <Eigen/Core>
#include <Eigen/LU>


#define torad(x)(x*M_PI/180.)

using namespace std;
using namespace GooFit;
using namespace ROOT;


//Initial and final states parameters
double pi_MASS  = 0.13957018; 
double D_MASS   = 1.96834; 
double k_MASS = 0.493677;

double Mother_MASS = D_MASS;
double d1_MASS  = pi_MASS; 
double d2_MASS  = pi_MASS;
double d3_MASS  = pi_MASS;

Variable Mother_Mass("Mother_Mass",D_MASS);
Variable Daughter1_Mass("Daughter1_Mass",d1_MASS);
Variable Daughter2_Mass("Daughter2_Mass",d2_MASS);
Variable Daughter3_Mass("Daughter3_Mass",d3_MASS);

//Bins for grid normalization
const int bins = 400;

//N Bins for eff and bkg scanning
const int bins_eff_bkg = bins;

//Dalitz Limits
fptype s12_min = POW2(d1_MASS  + d2_MASS);
fptype s12_max = POW2(D_MASS   - d3_MASS);
fptype s13_min = POW2(d1_MASS  + d3_MASS);
fptype s13_max = POW2(D_MASS   - d2_MASS);
fptype s23_min = POW2(d2_MASS  + d3_MASS);
fptype s23_max = POW2(D_MASS   - d1_MASS);

//Observables
Observable s12("s12",s12_min,s12_max); 
Observable s13("s13",s13_min,s13_max);
Observable s23("s23",s23_min,s23_max);
EventNumber eventNumber("eventNumber"); 

// For MIPWA
vector<fptype> HH_bin_limits; // bins over s_pipi spectrum
vector<Variable> pwa_coefs_amp; // mag(real)_coefs
vector<Variable> pwa_coefs_phs; // phase(imag)_coefs
const string pwa_file = "files/PWACOEFS_50pts.txt"; // input points for MIPWA

ResonancePdf *loadPWAResonance(const std::string fname = pwa_file, bool polar=true) {
//load the MIPWA resonance function

    std::ifstream reader;
    reader.open(fname.c_str());
    assert(reader.good());
    HH_bin_limits.clear();
    pwa_coefs_amp.clear();
    pwa_coefs_phs.clear();

    double e1, e2, e3;
    double emag, ephs;
    const int n = 10;
    const int m = 14;
    double min = 0.990 - n*0.01;
    double max = 0.990 + m*0.01;

    //reading input file
    size_t i = 0;
    while(reader >> e1 >> e2 >> e3) {

			
           HH_bin_limits.push_back(e1*e1); //MIPWA first input

           if(!polar){
	    	emag = e2*cos(e3);
                ephs = e2*sin(e3);
		//Instantiation of fit parameters for MIPWA
		Variable va(fmt::format("pwa_coef_{}_real", i), emag,0.01,-100.0,+100.0);
            	Variable vp(fmt::format("pwa_coef_{}_imag", i), ephs,0.01,-100.0,+100.0);
		pwa_coefs_amp.push_back(va);
            	pwa_coefs_phs.push_back(vp);
	    }else{
		emag = e2;
            	ephs = e3;
		Variable va(fmt::format("pwa_coef_{}_mag", i), emag,0.01,0.,+50.);
            	Variable vp(fmt::format("pwa_coef_{}_phase", i), ephs,0.01,0,0);//-2.*M_PI,+2.*M_PI);
		pwa_coefs_amp.push_back(va);
            	pwa_coefs_phs.push_back(vp);
	    } 
    		
           
            i++;

	

    }

    std::cout << "------------------------------------------" << std::endl;
    std::cout << pwa_coefs_amp.size() << " QMIPWA points loaded!" << std::endl;
    std::cout << "------------------------------------------" << std::endl;

    //global phase for S-wave (ever fixed)
    Variable swave_amp_real("swave_real_coef", 1.0);
    Variable swave_amp_imag("swave_imag_coef", 0.0);

    if(polar){

    	auto swave_12 = new Resonances::SplinePolar("MIPWA-Polar", swave_amp_real, swave_amp_imag, HH_bin_limits, pwa_coefs_amp, pwa_coefs_phs, PAIR_12, true);
	return swave_12;

    }else{

	auto swave_12 = new Resonances::Spline("MIPWA-Real", swave_amp_real, swave_amp_imag, HH_bin_limits, pwa_coefs_amp, pwa_coefs_phs, PAIR_12, true);
	return swave_12;

    }

    
} 

GooPdf* makeHistogramPdf(std::string efffile, std::string effhist, Observable s12, Observable s13, std::string name) {
//make a loop over DP Histogram and return a HistogramPDF

    std::vector<Observable> lvars = {s12,s13};
    auto binEffData = new BinnedDataSet(lvars);
    
    auto n = 0;

    //Open root file and get TH2 pointer
    auto f     = TFile::Open(efffile.c_str());
    auto effHistogram = (TH2F *)f->Get(effhist.c_str());
    
    //loop over DP
    for(int i = 0; i < bins_eff_bkg; ++i) {
        s12.setValue(s12.getLowerLimit() + (s12.getUpperLimit() - s12.getLowerLimit()) * (i + 0.5) / bins_eff_bkg);
        for(int j = 0; j < bins_eff_bkg; ++j) {
            s13.setValue(s13.getLowerLimit() + (s13.getUpperLimit() - s13.getLowerLimit()) * (j + 0.5) / bins_eff_bkg);
	    //check if pair (s12,s13) is inside DP
            if(!inDalitz(s12.getValue(), s13.getValue(), D_MASS, d1_MASS, d2_MASS, d3_MASS)){continue;}
            double weight = effHistogram->GetBinContent(effHistogram->FindBin(s12.getValue(), s13.getValue()));
            binEffData->addWeightedEvent(weight);
            n++;
        }
    }

    printf("Hist avaliated in %d events inside dalitz \n",n);
    Variable *effSmoothing = new Variable(name.c_str(), 0.0);
    effSmoothing->setBlind(0.0);
    auto ret = new SmoothHistogramPdf("efficiency", binEffData, *effSmoothing);
    return ret;

}

GooPdf *polyEff( Observable s12 , Observable s13){
//Polynomial default PDF if no eff provided
    vector<Variable> offsets;
    vector<Observable> observables;
    vector<Variable> coefficients;
    Variable constantOne("c1", 1);
    Variable constantZero("c0", 0);
    observables.push_back(s12);
    observables.push_back(s13);
    offsets.push_back(constantZero);
    offsets.push_back(constantZero);
    coefficients.push_back(constantOne);
    PolynomialPdf* eff = new PolynomialPdf("constantEff", observables, coefficients, offsets, 0); //No efficiency

    return eff;
}


GooPdf *make_veto(Observable s12, Observable s13) {
//Make veto if needed. For Ds we only apply fiducial 
   
    VetoInfo Fiducial_veto12(Variable("Fiducial_veto12_min", 0.), Variable("Fiducial_veto12_max", s12.getLowerLimit()), PAIR_12);
    VetoInfo Fiducial_veto13(Variable("Fiducial_veto13_min", 0.), Variable("Fiducial_veto13_max", s13.getLowerLimit()), PAIR_13);
    
    VetoInfo D0_veto12(Variable("D0_veto12_min", 3.31), Variable("D0_veto12_max", s12.getUpperLimit()+0.1), PAIR_12);
    VetoInfo D0_veto13(Variable("D0_veto13_min", 3.31), Variable("D0_veto13_max", s13.getUpperLimit()+0.1), PAIR_13);
    
    std::vector<VetoInfo> vetos;
    vetos.push_back(Fiducial_veto12);
    vetos.push_back(Fiducial_veto13);
    vetos.push_back(D0_veto12);
    vetos.push_back(D0_veto13);

    auto Dstar_veto = new DalitzVetoPdf("Dstar_veto", s12, s13, 
        Variable("Mother_Mass",D_MASS), Variable("Daughter1_Mass",d1_MASS), 
        Variable("Daughter2_Mass",d2_MASS), Variable("Daughter3_Mass",d3_MASS), 
        vetos);

    return Dstar_veto;
}

DalitzPlotPdf* makesignalpdf( Observable s12, Observable s13, EventNumber eventNumber, GooPdf* eff = 0){
// This function create our signal model 

    //set up the decay channel
    DecayInfo3 dtoppp;
    dtoppp.motherMass   = D_MASS;
    dtoppp.daug1Mass    = d1_MASS;
    dtoppp.daug2Mass    = d2_MASS;
    dtoppp.daug3Mass    = d3_MASS;
    dtoppp.meson_radius = 1.5;
	
    printf("Setting Signal model For Ds->3pi channel. \n");

    //Mass and width

    //parameters from Laura++
    double f0_980_MASS    = 0.990;
    double f0_980_GPP     = 0.165;
    double f0_980_GKK     = 4.2*0.165;
    double f0_980_WIDTH   = 0.4;
    double f0_980_amp     = 1.;
    double f0_980_img     = 0.0;

    //from PDG
    double a0_980_MASS    = 0.980 ;
    double a0_980_WIDTH   = 0.05;
    double a0_980_amp     = 1.0;
    double a0_980_img     = 0.0;
    
    //from E791 paper
    double f0_1370_MASS  = 1.434;
    double f0_1370_WIDTH = 0.172;
    double f0_1370_amp   = -0.8357;
    double f0_1370_img = -0.5730;

    //from PDG
    double f0_1500_MASS  = 1.505;
    double f0_1500_WIDTH = .109;
    double f0_1500_amp   = 1.;
    double f0_1500_img = 0.;

    //from PDG
    double omega_MASS   = 0.78265;
    double omega_WIDTH  = 0.00849;
    double omega_amp    = -0.016251;
    double omega_img  = 0.00207564;

    //From PDG and BABAR 
    double rho770_MASS   = 0.77549;
    double rho770_WIDTH  = 0.1491;
    double rho770_amp    = -0.0102842;//10.19*cos(1.1);
    double rho770_img  =  0.171685;//0.19*sin(1.1);

    //From PDG and BABAR 
    double rho1450_MASS   = 1.465 ;
    double rho1450_WIDTH  = 0.4  ;
    double rho1450_amp    = -0.167428;//1.2*cos(4.1);
    double rho1450_img  =  -0.769214;//1.2*sin(4.1);

    //From PDG - Ref resonance
    double f2_1270_MASS     = 1.2751;
    double f2_1270_WIDTH    = 0.1851;
    double f2_1270_amp      = 1.;
    double f2_1270_img    = 0.;
    
    // Setting fit parameters
    // Variable(name,value) for fixed parameters
    // Variable(name,value,error,lower_limit,upper_limit) for free parametes with limits
    // if you want without limits just do lower_limit=upper_limit=0

    // P-Wave
    //omega(782)
    Variable v_omega_Mass("omega_MASS",omega_MASS);
    Variable v_omega_Width("omega_WIDTH",omega_WIDTH);
    Variable v_omega_real("omega_REAL",omega_amp,0.01,0,0);
    Variable v_omega_img("omega_IMAG",omega_img,0.01,0,0);

    //rho(770)
    Variable v_rho770_Mass("rho770_MASS",rho770_MASS);
    Variable v_rho770_Width("rho770_WIDTH",rho770_WIDTH);
    Variable v_rho770_real("rho770_REAL",rho770_amp,0.01,0,0);
    Variable v_rho770_img("rho770_IMAG",rho770_img,0.01,0,0);
    
    //rho(1450)
    Variable v_rho1450_Mass("rho1450_MASS",rho1450_MASS);
    Variable v_rho1450_Width("rho1450_WIDTH",rho1450_WIDTH);
    Variable v_rho1450_real("rho1450_REAL",rho1450_amp,0.01,0,0);
    Variable v_rho1450_img("rho1450_IMAG",rho1450_img,0.01,0,0);

    //rho(1700)
    Variable v_rho1700_Mass("rho1700_MASS",1.720);
    Variable v_rho1700_Width("rho1700_WIDTH",0.25);
    Variable v_rho1700_real("rho1700_REAL",1.80595,0.01,0,0);
    Variable v_rho1700_img("rho1700_IMAG",-0.574122,0.01,0,0);


    //D-wave
    //f2(1270) Reference
    Variable v_f2_1270_Mass("f2_1270_MASS",f2_1270_MASS);
    Variable v_f2_1270_Width("f2_1270_WIDTH",f2_1270_WIDTH);
    Variable v_f2_1270_real("f2_1270_REAL",f2_1270_amp);
    Variable v_f2_1270_img("f2_1270_IMAG",f2_1270_img);

    //S-wave
    //f0(980)
    Variable v_f0_980_Mass("f0_980_MASS",f0_980_MASS);
    Variable v_f0_980_GPP("f0_980_GPP",f0_980_GPP,0.01,0,0);//0.05*f0_980_GPP,2.5*f0_980_GPP);
    Variable v_f0_980_GKK("f0_980_GKK",f0_980_GKK,0.01,0,0);//0.5*f0_980_GKK,2.);
    Variable v_f0_980_Width("f0_980_WIDTH",f0_980_WIDTH);
    Variable v_f0_980_real("f0_980_REAL",f0_980_amp, 0.01,0,0);
    Variable v_f0_980_img("f0_980_IMAG",f0_980_img, 0.01,0,0);

    //a0(980)
    Variable v_a0_980_Mass("a0_980_MASS",a0_980_MASS);
    Variable v_a0_980_Width("a0_980_WIDTH",a0_980_WIDTH);
    Variable v_a0_980_real("a0_980_REAL",a0_980_amp,0.01,0,0);
    Variable v_a0_980_img("a0_980_REAL",a0_980_img, 0.01,0,0);
    
    //f0(1370)
    Variable v_f0_1370_Mass("f0_1370_MASS",f0_1370_MASS);
    Variable v_f0_1370_Width("f0_1370_WIDTH",f0_1370_WIDTH);
    Variable v_f0_1370_real("f0_1370_REAL",f0_1370_amp,0.01,0,0);
    Variable v_f0_1370_img("f0_1370_IMAG",f0_1370_img,0.01,0,0);

    //f0(1500)
    Variable v_f0_1500_Mass("f0_1500_MASS",f0_1500_MASS);
    Variable v_f0_1500_Width("f0_1500_WIDTH",f0_1500_WIDTH);
    Variable v_f0_1500_real("f0_1500_REAL",f0_1500_amp, 0.01,0,0);
    Variable v_f0_1500_img("f0_1500_IMAG",f0_1500_img, 0.01, 0,0);

    //NR
    Variable nonr_real("nonr_REAL",1./*0.09*cos(torad(181))*/, 0.01,0,0);
    Variable nonr_imag("nonr_IMAG",0./*0.09*sin(torad(181))*/, 0.01,0,0);

    //Bose-Einstein - Parameter R from CMS paper
    Variable be_real("be_REAL",1.33775/*0.9*/,0.01,0,0);
    Variable be_imag("be_IMAG",8.58481/*-2.9*/,0.01,0,0);
    Variable be_coef("be_RCOEF",1.);
    Variable be_delta("be_RDELTA",73.e-3);
    //it is possible to initial variables above with random values in a range
    //e.g. v_omega_real.setRandomValue(-0.0160 - 5*0.0009,-0.0160 + 5*0.0009)
   
    //Instatiation of resonances
    auto omega = new Resonances::RBW("omega",v_omega_real,v_omega_img,v_omega_Mass,v_omega_Width,1,PAIR_12,true);
    
    auto rho770 = new Resonances::RBW("rho770",v_rho770_real,v_rho770_img,v_rho770_Mass,v_rho770_Width,1,PAIR_12,true);
    
    auto rho1450 = new Resonances::RBW("rho1450",v_rho1450_real,v_rho1450_img,v_rho1450_Mass,v_rho1450_Width,1,PAIR_12,true);

    auto rho1700 = new Resonances::RBW("rho1700",v_rho1700_real,v_rho1700_img,v_rho1700_Mass,v_rho1700_Width,1,PAIR_12,true);    

    auto f2_1270 = new Resonances::RBW("f2",v_f2_1270_real,v_f2_1270_img,v_f2_1270_Mass,v_f2_1270_Width,2,PAIR_12,true);

    auto f0_980 = new Resonances::FLATTE("f0_980",v_f0_980_real,v_f0_980_img,v_f0_980_Mass,
	v_f0_980_GPP,v_f0_980_GKK,PAIR_12,true);

    auto f0_980_RBW = new Resonances::RBW("f0_980_RBW",v_f0_980_real,v_f0_980_img,v_f0_980_Mass,v_f0_980_Width,(unsigned int)0,PAIR_12,true);
	
    //still implementing
    //auto f0_Mix = new Resonances::f0_MIXING("f0_a0_mixing",Variable("ar",1.,0.0001,0,0),Variable("ai",0.0,0.0001,0,0),Variable("a0_gkk",1.03*0.105,0.0001,0,0),Variable("f0_Gkk",0.165*4.2),Variable("a0_getapi",0.105,0.0001,0,0),Variable("f0_gpipi",0.165),Variable("a0_mass",0.982,0.0001,0,0),Variable("f0_mass",0.9773),PAIR_12,true);
 
    auto a0_980 = new Resonances::RBW("a0_980",v_a0_980_real,v_a0_980_img,v_a0_980_Mass,v_a0_980_Width,(unsigned int)0,PAIR_12,true);

    auto f0_1370 = new Resonances::RBW("f0_1370_12",v_f0_1370_real,v_f0_1370_img,v_f0_1370_Mass,v_f0_1370_Width,(unsigned int)0,PAIR_12,true);

    auto f0_1500 = new Resonances::RBW("f0_1500_12",v_f0_1500_real,v_f0_1500_img,v_f0_1500_Mass,v_f0_1500_Width,(unsigned int)0,PAIR_12,true);  

    auto nonr = new Resonances::NonRes("nonr", nonr_real, nonr_imag);

    auto BEC   = new Resonances::BoseEinstein("be",be_real,be_imag,be_coef,be_delta);

    //MIPWA
    ResonancePdf *MIPWA = loadPWAResonance(pwa_file, true);

    //If you want include a resonance in your model, just push into the vector 'vec_resonances'
    std::vector<ResonancePdf *> vec_resonances;
   
    vec_resonances.push_back(omega); 
    vec_resonances.push_back(rho770); 
    vec_resonances.push_back(rho1450);
    vec_resonances.push_back(rho1700);
    vec_resonances.push_back(f2_1270);
    vec_resonances.push_back(BEC);
    vec_resonances.push_back(MIPWA);

    //not included
    //vec_resonances.push_back(a0_980);
    //vec_resonances.push_back(f0_980);
    //vec_resonances.push_back(f0_Mix);
    //vec_resonances.push_back(f0_1500);
    //vec_resonances.push_back(f0_1370);
    //vec_resonances.push_back(nonr);

    dtoppp.resonances = vec_resonances;

    if(!eff)
        eff = polyEff(s12,s13);
   
    return new DalitzPlotPdf("signalPDF", s12, s13, eventNumber, dtoppp, eff);
}

void getData(std::string toyFileName, GooFit::Application &app, DataSet &data, bool toy) {
    //load data in a GooFit::dataset

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
        if((s12.getValue()<3.31)
            &&(s13.getValue()<3.31)
            &&(s12.getValue()>s12.getLowerLimit())
            &&(s13.getValue()>s13.getLowerLimit()))
        {
            data.addEvent();
            if(j<10) printf("[%d] = (%f , %f)\n",i,s12.getValue(),s13.getValue());
            j++;
            if(!toy  &&  data.getNumEvents()==200000) break; 
        }
    }


}

void to_root(UnbinnedDataSet& toyMC , std::string name ){
//save GooFit::Dataset in a root file
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



TH1F* plotSigComps(DalitzPlotPdf* signal,const unsigned int index){
//plot Signal componentes -- not working yet
    Observable s12("s12", s12_min, s12_max);
    Observable s13("s13", s13_min, s13_max);
    s12.setNumBins(bins);
    s13.setNumBins(bins);

    DecayInfo3 dtoppp;
    dtoppp.motherMass   = D_MASS;
    dtoppp.daug1Mass    = d1_MASS;
    dtoppp.daug2Mass    = d2_MASS;
    dtoppp.daug3Mass    = d3_MASS;
    dtoppp.meson_radius = 1.5;

    auto eff = polyEff(s12,s13);
    
    ResonancePdf* res = signal->getDecayInfo().resonances[index];
    dtoppp.resonances.push_back(res);
    DalitzPlotPdf *compPdf = new DalitzPlotPdf((fmt::format("compPdf_{0}",index)).c_str(), s12, s13, eventNumber, dtoppp, eff);
    
    ProdPdf ProdRes((fmt::format("ProdcompPdf_{0}",index)).c_str(),{compPdf});
    DalitzPlotter dplotter{&ProdRes, compPdf};
    auto hist = dplotter.make2D();
    auto proj_s12 = (TH1F*)hist->ProjectionX((fmt::format("Projs12_{0}",index)).c_str());

    return proj_s12;
    
}

std::vector<std::vector<fptype>> fractions(DalitzPlotPdf* signal,UnbinnedDataSet* flatMC){
//calc the fit fractions - a FlatMC need to be given as input for integration

    cout << "Fit Fractions" << endl;
    auto obs               = flatMC->getObservables();
    Observable s12         = obs.at(0);
    Observable s13         = obs.at(1);
    Observable eventNumber = obs.at(2);

    TFile *f = new TFile("MC/toyMC_norm.root","read");
    TTree *s = (TTree*)f->Get("DecayTree");

    double _s12 = -1, _s13 = -1;
    s->SetBranchAddress("s12",&_s12);
    s->SetBranchAddress("s13",&_s13);

    for(int i =0; i< s->GetEntries();i++){
        s->GetEntry(i);
        s12.setValue(_s12);
        s13.setValue(_s13);
        flatMC->addEvent();
	
    }

    signal->setParameterConstantness(true); 
    ProdPdf overallsignal{"overallsignal",{signal}};
    overallsignal.setData(flatMC);
    signal->setDataSize(flatMC->getNumEvents());

    FitManagerMinuit2 fitter(&overallsignal);
    fitter.setVerbosity(0);
    fitter.fit();
    
    return signal->fit_fractions();
   
}

int genFit(GooPdf *totalPdf,DalitzPlotPdf *signal, UnbinnedDataSet *data, std::string rank) {
//genfit 
    totalPdf->setData(data);
    signal->setDataSize(data->getNumEvents());
    
    FitManager datapdf(totalPdf);
    datapdf.setVerbosity(0);
    datapdf.setMaxCalls(20000);
    datapdf.setTolerance(0.1);
    auto func_min = datapdf.fit();
    
    auto output = fmt::format("GenFit/{0}_fitResult.txt",rank);
    writeToFile(signal, output.c_str());

    std::ofstream open(output,std::ios_base::app);
    open << "FCN" << "\t" << func_min.Fval() << "\t" << 0 << "\t" << 0 << "\t" << 0 << std::endl;
    open << "Status" << "\t" << func_min.IsValid() << "\t" << 0 << "\t" << 0 << "\t" << 0 << std::endl;
    open.close();

    auto norm = signal->normalize() ;
    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "Sample " << rank << " fitted ! "		         << std::endl;
    std::cout << "Status " << func_min.IsValid() 		<< std::endl;
    std::cout << "nEvents --> " << data->getNumEvents()          << std::endl;
    std::cout << "minFCN --> "   << func_min.Fval()              << std::endl;
    std::cout << "AmpNorm --> "   << norm       << std::endl;
    std::cout << "---------------------------------------------" << std::endl;

    return datapdf;
}

void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
{
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.bottomRows(numRows-rowToRemove);

    matrix.conservativeResize(numRows,numCols);
}

void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove)
{
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.rightCols(numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
}

//convert real,imag fit result to mag,phase
//still implementing
double dfMag_x(double x, double y){
	return x/sqrt(x*x + y*y) ;
}

double dfMag_y(double x, double y){
        return y/sqrt(x*x + y*y) ;
}

double dfPhs_x(double x, double y){
	return -y/(x*x + y*y);
}

double dfPhs_y(double x, double y){
        return x/(x*x + y*y);
}

void convertToMagPhase(Eigen::MatrixXd &m,std::vector<ROOT::Minuit2::MinuitParameter> &params,size_t totalUserPar, size_t totalFitPar,std::string name){

   auto output = fmt::format("Fit/{0}/ConvertedToMagPhaseResults.txt",name.c_str());
   std::ofstream wt(output);
   
   std::vector<double> vec_values;
   std::vector<double> vec_error;

   for(int i=0; i<totalUserPar; i++){

	if(params[i].IsConst()){
		continue;
	}else{
		vec_values.push_back(params[i].Value());
		vec_error.push_back(params[i].Error());
	}

   }


   for(int row=0; row<totalFitPar; row+=2){
	double mag = 0;
	double mag_error = 0;
	double phs = 0;
	double phs_error =0;

	double real = vec_values[row];
	double imag = vec_values[row+1];
	double real_error = vec_error[row];
	double imag_error = vec_error[row+1];

	mag = sqrt(POW2(real) + POW2(imag));
	phs = atan2(imag,real);
	
	mag_error = POW2(dfMag_x(real,imag)*real_error) + POW2(dfMag_y(real,imag)*imag_error);
	phs_error = POW2(dfPhs_x(real,imag)*real_error) + POW2(dfPhs_y(real,imag)*imag_error);


	//for(int col=0; col<totalFitPar; col++){

			mag_error += 2*dfMag_x(real,imag)*dfMag_y(real,imag)*m(row+1,row);
			phs_error += 2*dfPhs_x(real,imag)*dfPhs_y(real,imag)*m(row+1,row);
			
	//}

	mag_error = sqrt(abs(mag_error));
	phs_error = sqrt(abs(phs_error));

	wt << fmt::format("par[{0}] = ( {1} +/- {2} , {3} +/- {4} )",row/2,mag,mag_error,phs*180./M_PI,phs_error*180./M_PI) << std::endl;
	//printf(" par[%d] = ( %f +/- %f , %f +/- %f ) \n", row/2, mag,mag_error, phs*180./M_PI, phs_error*180./M_PI);
		
   }

	wt.close();

}


DalitzPlotPdf* runFit(GooPdf *totalPdf,DalitzPlotPdf *signal, UnbinnedDataSet *data, std::string name, bool save_toy) {
//This function run the data fit

    //Setting data and EventNumber
    totalPdf->setData(data);
    signal->setDataSize(data->getNumEvents());

    //Fitter (it uses ROOT::FunctionMinimum API)
    FitManager datapdf(totalPdf);
    datapdf.setVerbosity(2);
    datapdf.setMaxCalls(200000);
    datapdf.setTolerance(0.1);

    //run the fit
    auto func_min = datapdf.fit();

    auto output = fmt::format("Fit/{0}/fit_result.txt",name.c_str());
    //save external state user parameters after fit
	//here include Const and Free parameters
    writeToFile(totalPdf, output.c_str());

    std::ofstream open(output,std::ios_base::app);
    open << "FCN" << "\t" << func_min.Fval() << "\t" << 0 << "\t" << 0 << "\t" << 0 << std::endl;
    open << "Status" << "\t" << func_min.IsValid() << "\t" << 0 << "\t" << 0 << "\t" << 0 << std::endl;
	open.close();

	//getting the number of free parameters for covMatrix and chi2
	size_t npar = 0; 
	auto param = datapdf.getParams()->Parameters();
    for(size_t i = 0; i < param.size(); i++){
        if(!param[i].IsConst()){
            npar++;	
        }
    }

	//User CovMatrix (include rows and cols of non free parameters)
	output = fmt::format("Fit/{0}/covMatrix.txt",name.c_str());
	std::ofstream open2(output);
    auto covMatrix = func_min.UserState().Covariance().Data();
    Eigen::MatrixXd m(param.size(),param.size());

	//std::cout << covMatrix[1 + 0*(0+1)/2] << '\t' << covMatrix[0 + 1*(1+1)/2] << std::endl;
	for(int i=0; i<param.size(); i++){
		for(int j=0; j<param.size(); j++){
			if(i>j){
				m(i,j) = covMatrix[j + i*(i+1)/2];
			}else{
				m(i,j) = covMatrix[i + j*(j+1)/2];
			}
		}
	}

	//we don't want non free parameters rows and columns, so we must to remove them
	for(int i =0; i< param.size(); i++){
		if(param[i].IsConst()){
			removeRow(m,i);
			removeColumn(m,i);
		}
	}

   //saving sub CovMatrix in a file
   open2.precision(5);
   open2 <<  m << std::endl;
   open2.close();

   //not working properly yet
//   convertToMagPhase(m,param,param.size(),npar,name);

    //TIP!
    //If you want to free some paramters after previous fit
    //signalpdf->getParameterByName("parameter")->setFixed(false)
    //FitManager datapdf2(totalPdf);
    //auto func_min2 = datapdf2.fit();

    //Get the value of normalization
    auto norm = signal->normalize() ;
    
	//Status
    std::cout << "---------------------------------------------" << std::endl;
    if(func_min.IsValid()){
        std::cout << "Sample --> fitted successfully! "		     << std::endl;
    }else{
        std::cout << "Sample --> Not Converged! "                << std::endl;
    }
    std::cout << "HasAccurateCovar --> " << func_min.HasAccurateCovar() << std::endl;
    std::cout << "nEvents --> " << data->getNumEvents()          << std::endl;
    std::cout << "minFCN --> "   << func_min.Fval()              << std::endl;
    std::cout << "SignalNorm --> "   <<   norm    << std::endl;
    std::cout << fmt::format("Output of fit saved in Fit folder Fit/{0}",name)  << std::endl;
    std::cout << "---------------------------------------------" << std::endl; 

    //Save a toy of the model if you want
    DalitzPlotter dplotter{totalPdf, signal};
    {
        if(save_toy){

            auto obs               = data->getObservables();
            Observable s12         = obs.at(0);
            Observable s13         = obs.at(1);
            Observable eventNumber = obs.at(2);
            
            UnbinnedDataSet toyMC({s12,s13,eventNumber});
            dplotter.fillDataSetMC(toyMC,10000000);
            to_root(toyMC,fmt::format("Fit/{0}/toyMC.root",name.c_str()));
        }

    }

    //Plot model and Data projection
    TCanvas foo;
    dplotter.Plot("s_{#pi^{-} #pi^{+}}","s_{#pi^{-} #pi^{+}}","s_{#pi^{+} #pi^{+}}",fmt::format("Fit/{0}",name.c_str()),*data);
   
    //Calc chi2
    dplotter.chi2(npar,"Ds3pi_bins.txt",0.05,1.95,0.3,3.4,*data,fmt::format("Fit/{0}",name.c_str())); 

    
    //saving s-wave qmpiwa
    size_t N = HH_bin_limits.size();
    double amp[N];
    double amp_er[N];
    double phase[N];
    double phase_er[N];
    double s[N];
    double s_er[N];
    
    std::ofstream wr(fmt::format("Fit/{0}/s_wave_fitted.txt",name.c_str()));
    for(size_t i = 0; i < N ; i++){
        	s[i] = HH_bin_limits[i];
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

/*
    //plot Signal Components --- only plotting the first one. Why???
    auto nRes = signal->getDecayInfo().resonances.size();
    for(int i =0; i<nRes; i++){
        auto h0 = plotSigComps(signal,i);
        h0->Draw("HIST");
        auto n = fmt::format("plots/{0}.png",signal->getDecayInfo().resonances[i]->getName());
        foo.SaveAs(n.c_str());
        delete h0;
    }*/
   
    foo.Divide(2);
    foo.cd(1);
    gr_real->Draw("ALP");
    foo.cd(2);
    gr_img->Draw("ALP");
    output = fmt::format("Fit/{0}/Swave.png",name.c_str());
    foo.SaveAs(output.c_str());

    return signal;
}


int main(int argc, char **argv){
    
    GooFit::Application app{"Genfit",argc,argv};
    
    std::string input_data_name = "input.root";
    std::string fit_name = "Fit";
    bool save_toy = false;
    bool is_toy = false;
    auto fit = app.add_subcommand("fit","fit data");
    fit->add_option("-f,--file",input_data_name,"name_of_file.root");
    fit->add_option("-t,--isToy", is_toy, "Get toyData for fit") ;
    fit->add_option("-s,--saveToy",save_toy,"save toy in root file");
    fit->add_option("-n,--fitName",fit_name,"name of this fit(useful to save results)")->required(true);
   
    size_t Nevents=1000000;
    std::string toyName = "MC.root";
    auto makeToy = app.add_subcommand("makeToy","make a toy");
    makeToy->add_option("-e,--nevents",Nevents,"number of events");
    makeToy->add_option("-n,--name",toyName,"output_toy_name.root");
    makeToy->add_option("-s,--saveToy",save_toy,"save toy in root file");	
    
	size_t sampleNumber = 0;
    auto genfit = app.add_subcommand("genfit","perform genfit");
	genfit->add_option("-n,--sampleNumber",sampleNumber,"sample number");

    auto NormStudies = app.add_subcommand("NormStudies","perform norm studies");

	size_t fitNumber = 0;
	auto fitsys = app.add_subcommand("fitsys","perform fits with random initial parameters.");
	fitsys->add_option("-n,--fitNumber",fitNumber,"number of fits to perform");
    fitsys->add_option("-f,--file",input_data_name,"name_of_file.root");
    fitsys->add_option("-t,--isToy", is_toy, "Get toyData for fit") ;

    app.require_subcommand();

    GOOFIT_PARSE(app);

    // Make the MC directory if it does not exist
    std::string command = "mkdir -p MC";
    if(system(command.c_str()) != 0)
        throw GooFit::GeneralError("Making `MC` directory failed");

    // Make the Fit directory if it does not exist
    command = "mkdir -p Fit";
    if(system(command.c_str()) != 0)
        throw GooFit::GeneralError("Making `Fit` directory failed");

    // Make the Plots directory if it does not exist
    command = "mkdir -p plots";
    if(system(command.c_str()) != 0)
        throw GooFit::GeneralError("Making `plots` directory failed");
    
    // Make the NormStudies directory if it does not exist
    command = "mkdir -p NormStudies";
    if(system(command.c_str()) != 0)
        throw GooFit::GeneralError("Making `NormStudies` directory failed");

    command = "export CUDA_VISIBLE_DEVICES=1";
    if(system(command.c_str()) != 0)
        throw GooFit::GeneralError("Set GPU 1 faliled");

    Observable s12("s12", s12_min, s12_max);
    Observable s13("s13", s13_min, s13_max);
    EventNumber eventNumber("eventNumber");

    s12.setNumBins(bins);
    s13.setNumBins(bins);

    const string bkgfile = "../../../dados/bkgBW_16_BDT0.18_Smoothed.root";
    const string efffile = "../../../dados/acc_15_MC_TIS_RW_BDT0.18_SigRegion_Smoothed.root";
    const string bkghist = "h_eff";
    const string effhist = "h_eff";

    auto vetos = make_veto(s12,s13);
    auto efficiency = makeHistogramPdf(efffile,effhist,s12,s13,"eff_coef");
    auto effWithVeto = new ProdPdf("effWithVeto",{efficiency,vetos});

    auto background = makeHistogramPdf(bkgfile,bkghist,s12,s13,"bkg_coef");
    auto bkgWithVeto = new ProdPdf("bkgWithVeto",{background,vetos});

    auto signal = makesignalpdf(s12, s13, eventNumber,effWithVeto);
    AddPdf *prodpdf = new AddPdf("prodpdf", Variable("frac",0.926), signal, bkgWithVeto) ;


    if(*makeToy) {
	DalitzPlotter dplotter{prodpdf, signal};
	UnbinnedDataSet data({s12, s13, eventNumber});
	dplotter.fillDataSetMC(data, Nevents);
        std::cout << "----------------------------------------------------------" << std::endl;
        std::cout <<   data.getNumEvents() << " events was generated!"  << std::endl;
        std::cout << "----------------------------------------------------------" << std::endl;
        if(save_toy) {
		auto fullName= fmt::format("MC/{0}",toyName);
        	to_root(data,fullName);
		std::cout <<   toyName << " root file was saved in MC folder" << std::endl;
        	std::cout << "----------------------------------------------------------" << std::endl;
        }
	return 0;
    }    

    if(*fit){

	auto command = fmt::format("mkdir -p Fit/{0}",fit_name);
    	if(system(command.c_str()) != 0)
        	throw GooFit::GeneralError("Making directory failed");

        std::cout << "------------------------------------------" << std::endl;
        std::cout << "Reading file --> " << input_data_name       << std::endl; 
        std::cout << "------------------------------------------" << std::endl;
	UnbinnedDataSet data({s12, s13, eventNumber});
        getData(input_data_name, app, data,is_toy);	
        auto output_signal = runFit(prodpdf,signal, &data, fit_name, save_toy);
        auto flatMC = new UnbinnedDataSet({s12,s13,eventNumber});
        auto frac = fractions(output_signal,flatMC);
        ofstream wt(fmt::format("Fit/{0}/fractions.txt",fit_name));
        std::cout << "Fit Fractions Interference" << '\n';
		auto nres = signal->getDecayInfo().resonances.size();
		Eigen::MatrixXd m(nres,nres);
		for(int i=0; i < nres; i++)
			m.row(i) = Eigen::Map<Eigen::VectorXd>(&frac[i][0],nres);

	std::cout << m << std::endl;
	wt << m << std::endl;

    }
        
            
    if(*genfit){
		writeToFile(signal, "GenFit/InitialParameters.txt");
		UnbinnedDataSet data({s12, s13, eventNumber});
		DalitzPlotter dplotter{prodpdf, signal};
		dplotter.fillDataSetMC(data, 100000);	
		std::cout << "------------------------------------------" << std::endl;
		std::cout << "Fitting Sample --> " << sampleNumber		          << std::endl; 
		std::cout << "------------------------------------------" << std::endl;
		auto sample_name = fmt::format("Sample_{0}",sampleNumber);	
		auto output_signal = genFit(prodpdf,signal, &data, sample_name);
    }

    if(*fitsys){
		UnbinnedDataSet data({s12, s13, eventNumber});
		getData(input_data_name, app, data,is_toy);	
		std::cout << "------------------------------------------" << std::endl;
		std::cout << "Fit Number --> " << fitNumber	          	  << std::endl; 
		std::cout << "------------------------------------------" << std::endl;
		auto sample_name = fmt::format("Fit_{0}",fitNumber);	
		auto output_signal = genFit(prodpdf,signal, &data, sample_name);
    }

            
    if(*NormStudies){
        auto fullName= fmt::format("{0}",input_data_name);
		UnbinnedDataSet data({s12, s13, eventNumber});
        std::cout << "------------------------------------------" << std::endl;
        std::cout << "Reading file --> " << fullName              << std::endl; 
        std::cout << "------------------------------------------" << std::endl;
        getData(fullName, app, data,is_toy);
                
        for(int i=1200; i<=6000; i+=100){
                    auto output = fmt::format("NormStudies/fit_bins_{0}.txt",i);
                    s12.setNumBins(i);
                    s13.setNumBins(i);
                    auto output_signal = runFit(prodpdf,signal, &data, output, save_toy);
        }
        
     }
                
}

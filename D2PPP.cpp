	// ROOT stuff
	#include <TApplication.h>
	#include <TCanvas.h>
	#include <TFile.h>
	#include <TRandom3.h>
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

	double rho1450_MASS   = 1.465 ;
	double rho1450_WIDTH  = 0.4 ;

	double rho1700_MASS   = 1.720 ;
	double rho1700_WIDTH  = 0.250 ;


	//Initial and final states parameters
	double pi_MASS  = 0.13957018; 
	//D_MASS from inv Mass fit
	double D_MASS   = 1.96834;//MeV/c²  
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
	const int bins = 1000;

	//N Bins for eff and bkg scanning
	const int bins_eff_bkg = 3000;

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
	const string pwa_file = "files/pwa_coefs.txt";//"files/PWACOEFS_50pts.txt"; // input points for MIPWA

	const int NevG = 1e7; 

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
				Variable vp(fmt::format("pwa_coef_{}_phase", i), ephs,0.01,-2.*M_PI,+2.*M_PI);
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


	GooPdf* makeHistogramPdf(std::string efffile, std::string effhist, Observable s12, Observable s13,bool eff, bool doEffSwap, bool saveEffPlot) {

	    std::vector<Observable> lvars = {s12,s13};
	    BinnedDataSet *binEffData = new BinnedDataSet(lvars);
	    
	    TFile *f     = TFile::Open(efffile.c_str());
	    auto bkgHistogram = (TH2F *)f->Get(effhist.c_str());
	    bkgHistogram->SetStats(false);  

	    TRandom3 donram(16);

	    if(eff){
		std::cout << "---------------------------------------------" << std::endl;
		std::cout << "Loading Efficiency PDF" << "\n";
		std::cout << "---------------------------------------------" << std::endl;    
	    }else{
		std::cout << "---------------------------------------------" << std::endl;
		std::cout << "Loading Background PDF" << "\n";
		std::cout << "---------------------------------------------" << std::endl;
	    }

	    for(int i = 0; i < bins_eff_bkg; ++i) {
		s12.setValue(s12.getLowerLimit() + (s12.getUpperLimit() - s12.getLowerLimit()) * (i + 0.5) / bins_eff_bkg);
		for(int j = 0; j < bins_eff_bkg; ++j) {
		    s13.setValue(s13.getLowerLimit() + (s13.getUpperLimit() - s13.getLowerLimit()) * (j + 0.5) / bins_eff_bkg);
		    //check if pair (s12,s13) is inside DP
		    if(!inDalitz(s12.getValue(), s13.getValue(), D_MASS, d1_MASS, d2_MASS, d3_MASS)){continue;}
		    double weight = bkgHistogram->GetBinContent(bkgHistogram->FindBin(s12.getValue(), s13.getValue()));
		    binEffData->addWeightedEvent(weight);

		    if(doEffSwap) {
			double swapmass = s12.getValue();
			s12.setValue(s13.getValue());
			s13.setValue(swapmass);
			weight = bkgHistogram->GetBinContent(bkgHistogram->FindBin(s12.getValue(), s13.getValue()));
			binEffData->addWeightedEvent(weight);
		    }
		}
	    } 
	    
	    if(saveEffPlot) {
		TCanvas foo;
		if(eff){
		    foo.cd();
		    bkgHistogram->Draw("colz");
		    foo.SaveAs("plots/efficiency_bins.png");
		    foo.SetLogz(true);
		    foo.SaveAs("plots/efficiency_bins_log.png");
		}else{
		    foo.cd();
		    bkgHistogram->Draw("colz");
		    foo.SaveAs("plots/background_bins.png");
		    foo.SetLogz(true);
		    foo.SaveAs("plots/background_bins_log.png");
		}
		
	    }
	  
	   
	    Variable *effSmoothing= nullptr;
	    
	    if(eff){
		effSmoothing = new Variable("effSmoothing", 0.);
		return new SmoothHistogramPdf("eff_pdf", binEffData, *effSmoothing);
	    }else{
		effSmoothing = new Variable("bkgSmoothing", 0.);
		return new SmoothHistogramPdf("bkg_pdf", binEffData, *effSmoothing);
	    }

	    
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
	    double f0_980_MASS    = 0.965;
	    double f0_980_GPP     = 0.342;
	    double f0_980_GKK     = 1.44;
	    double f0_980_WIDTH   = 0.4;
	    double f0_980_amp     = 1.;
	    double f0_980_img     = 0.0;

	    //from PDG 2020
	    double omega_MASS   = 0.78265;
	    double omega_WIDTH  = 0.00849;
	    double omega_amp    = -0.0173295;
	    double omega_img  = 0.00528296;

	    //From PDG 2020 CHARGED ONLY, HADROPRODUCED
	    double rho770_MASS   = 0.77549;
	    double rho770_WIDTH  = 0.1491;
	    double rho770_amp    = 0.0215937;
	    double rho770_img  =  0.131332;
	    double rho770_MASS_lower    = rho770_MASS - 2*0.01;
	    double rho770_MASS_upper  =  rho770_MASS + 2*0.01;
	    double rho770_WIDTH_lower    = rho770_WIDTH - 2*0.06;
	    double rho770_WIDTH_upper  =  rho770_WIDTH + 2*0.06;

	    //From PDG 2020
	    double rho1450_amp    = -0.304873;
	    double rho1450_img  =  -1.09686;
	    double rho1450_MASS_lower    = 0.0;//rho1450_MASS - 5*0.025;
	    double rho1450_MASS_upper  =  0.0;//1.8;//rho1450_MASS + 5*0.025;
	    double rho1450_WIDTH_lower    = 0.0;//rho1450_WIDTH - 5*0.06;
	    double rho1450_WIDTH_upper  =  0.0;//rho1450_WIDTH + 5*0.06;

	    //From PDG 2020 
	    double rho1700_amp    = 0.590158;
	    double rho1700_img  = -1.12964;
	    double rho1700_MASS_lower    = 0.0;//1.4;//rho1700_MASS - 5*0.02;
	    double rho1700_MASS_upper  =  0.0;//2.0;//rho1700_MASS + 5*0.02;
	    double rho1700_WIDTH_lower    = 0.0;//rho1700_WIDTH - 1*0.1;
	    double rho1700_WIDTH_upper  =  0.0;//rho1700_WIDTH + 5*0.1;

	    //From PDG 2020 - ABLIKIM
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
	    Variable v_rho770_Mass("rho770_MASS",rho770_MASS,0.01,rho770_MASS_lower,rho770_MASS_upper);
	    Variable v_rho770_Width("rho770_WIDTH",rho770_WIDTH,0.01,rho770_WIDTH_lower,rho770_WIDTH_upper);
	    Variable v_rho770_real("rho770_REAL",rho770_amp,0.01,0,0);
	    Variable v_rho770_img("rho770_IMAG",rho770_img,0.01,0,0);
	    
	    //rho(1450)
	    Variable v_rho1450_Mass("rho1450_MASS",rho1450_MASS,0.01,rho1450_MASS_lower,rho1450_MASS_upper);
	    Variable v_rho1450_Width("rho1450_WIDTH",rho1450_WIDTH,0.01,rho1450_WIDTH_lower,rho1450_WIDTH_upper);
	    Variable v_rho1450_real("rho1450_REAL",rho1450_amp,0.01,0,0);
	    Variable v_rho1450_img("rho1450_IMAG",rho1450_img,0.01,0,0);

	    //rho(1700)
	    Variable v_rho1700_Mass("rho1700_MASS",rho1700_MASS,0.01,rho1700_MASS_lower,rho1700_MASS_upper);
	    Variable v_rho1700_Width("rho1700_WIDTH",rho1700_WIDTH,0.01,rho1700_WIDTH_lower,rho1700_WIDTH_upper);
	    Variable v_rho1700_real("rho1700_REAL",rho1700_amp,0.01,0,0);
	    Variable v_rho1700_img("rho1700_IMAG",rho1700_img,0.01,0,0);

	    //D-wave
	    //f2(1270) Reference
	    Variable v_f2_1270_Mass("f2_1270_MASS",f2_1270_MASS);
	    Variable v_f2_1270_Width("f2_1270_WIDTH",f2_1270_WIDTH);
	    Variable v_f2_1270_real("f2_1270_REAL",f2_1270_amp);
	    Variable v_f2_1270_img("f2_1270_IMAG",f2_1270_img);

	    //S-wave
	    //f0(980)
	    Variable v_f0_980_Mass("f0_980_MASS",f0_980_MASS);
	    Variable v_f0_980_GPP("f0_980_GPP",f0_980_GPP);
	    Variable v_f0_980_GKK("f0_980_GKK",f0_980_GKK);
	    Variable v_f0_980_Width("f0_980_WIDTH",f0_980_WIDTH);
	    Variable v_f0_980_real("f0_980_REAL",f0_980_amp,0.01,0,0);
	    Variable v_f0_980_img("f0_980_IMAG",f0_980_img,0.01,0,0);

	    //Bose-Einstein - Parameter R from CMS paper
	    Variable be_real("be_REAL",0.,0.01,0,0);
	    Variable be_imag("be_IMAG",0.,0.01,0,0);
	    Variable be_coef("be_RCOEF",1.0);
	    Variable be_delta("be_RDELTA",0.);//73.e-3);

        //v_rho1450_Mass.setGaussianRandomValue(rho1450_MASS, 0.025);
        //v_rho1450_Width.setGaussianRandomValue(rho1450_WIDTH,0.06);
        //v_rho1700_Mass.setGaussianRandomValue(rho1700_MASS, 0.02);
        //v_rho1700_Width.setGaussianRandomValue(rho1700_WIDTH, 0.1);

        std::cout << fmt::format("Rho770: mass={} width={}",v_rho770_Mass.getValue(),v_rho770_Width.getValue()) << '\n';
        std::cout << fmt::format("Rho1450: mass={} width={}",v_rho1450_Mass.getValue(),v_rho1450_Width.getValue()) << '\n';
        std::cout << fmt::format("Rho1700: mass={} width={}",v_rho1700_Mass.getValue(),v_rho1700_Width.getValue()) << '\n';

	v_rho770_Mass.setFixed(true);
	v_rho770_Width.setFixed(true);
	v_rho1450_Mass.setFixed(true);
	v_rho1450_Width.setFixed(true);
	v_rho1700_Mass.setFixed(true);
	v_rho1700_Width.setFixed(true);

	    //it is possible to initial variables above with random values in a range
	    //e.g. v_omega_real.setRandomValue(-0.0160 - 5*0.0009,-0.0160 + 5*0.0009)
	   
	    //Instatiation of resonances
	    
	    auto omega = new Resonances::RBW("omega",v_omega_real,v_omega_img,v_omega_Mass,v_omega_Width,1,PAIR_12,true,true);
	    auto rho770 = new Resonances::GS("rho770",v_rho770_real,v_rho770_img,v_rho770_Mass,v_rho770_Width,1,PAIR_12,true,true);
	    auto rho1450 = new Resonances::GS("rho1450",v_rho1450_real,v_rho1450_img,v_rho1450_Mass,v_rho1450_Width,1,PAIR_12,true,true);
	    auto rho1700 = new Resonances::GS("rho1700",v_rho1700_real,v_rho1700_img,v_rho1700_Mass,v_rho1700_Width,1,PAIR_12,true,true);    
	    auto f2_1270 = new Resonances::RBW("f2",v_f2_1270_real,v_f2_1270_img,v_f2_1270_Mass,v_f2_1270_Width,2,PAIR_12,true,true);
	    auto BEC   = new Resonances::BoseEinstein("be",be_real,be_imag,be_coef,be_delta);
	    
	   /* 
	    auto omega = new Resonances::RBW("omega",v_omega_real,v_omega_img,v_omega_Mass,v_omega_Width,1,PAIR_12,true,false);
	    auto rho770 = new Resonances::GS("rho770",v_rho770_real,v_rho770_img,v_rho770_Mass,v_rho770_Width,1,PAIR_12,true,false);
	    auto rho1450 = new Resonances::GS("rho1450",v_rho1450_real,v_rho1450_img,v_rho1450_Mass,v_rho1450_Width,1,PAIR_12,true,false);
	    auto rho1700 = new Resonances::GS("rho1700",v_rho1700_real,v_rho1700_img,v_rho1700_Mass,v_rho1700_Width,1,PAIR_12,true,false);    
	    auto f2_1270 = new Resonances::RBW("f2",v_f2_1270_real,v_f2_1270_img,v_f2_1270_Mass,v_f2_1270_Width,2,PAIR_12,true,false);
	    auto BEC   = new Resonances::BoseEinstein("be",be_real,be_imag,be_coef,be_delta);*/
	    //MIPWA
	    ResonancePdf *MIPWA = loadPWAResonance(pwa_file, false);

	    //If you want include a resonance in your model, just push into the vector 'vec_resonances'
	    std::vector<ResonancePdf *> vec_resonances;
	   
	    vec_resonances.push_back(omega); 
	    vec_resonances.push_back(rho770); 
	    vec_resonances.push_back(rho1450);
	    vec_resonances.push_back(rho1700);
	    vec_resonances.push_back(f2_1270);
	    //vec_resonances.push_back(BEC);
	    vec_resonances.push_back(MIPWA);

	    //not included
	    //vec_resonances.push_back(f0_980);

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
        if((s12.getValue()<s12.getUpperLimit())
            &&(s13.getValue()<s13.getUpperLimit())
            &&(s12.getValue()>s12.getLowerLimit())
            &&(s13.getValue()>s13.getLowerLimit()))
        {
            data.addEvent();
            if(j<10) printf("[%d] = (%f , %f)\n",i,s12.getValue(),s13.getValue());
            j++;
//            if(!toy  &&  data.getNumEvents()==100000) break; 
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


DalitzPlotPdf* runFit(GooPdf *totalPdf,DalitzPlotPdf *signal, UnbinnedDataSet *data, std::string name, bool save_toy) {
    //This function run the data fit
    //TIP!
    //If you want to free some paramters after previous fit
    //signalpdf->getParameterByName("parameter")->setFixed(false)
    //FitManager datapdf2(totalPdf);
    //auto func_min2 = datapdf2.fit();
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
    datapdf.printParams(fmt::format("Fit/{0}/fit_result_mag_phase.txt",name.c_str()));

    auto output = fmt::format("Fit/{0}/fit_result.txt",name.c_str());
    writeToFile(totalPdf, output.c_str());

    std::ofstream open(output,std::ios_base::app);
    open << "FCN" << "\t" << func_min.Fval() << "\t" << 0 << "\t" << 0 << "\t" << 0 << std::endl;
    open << "Status" << "\t" << func_min.IsValid() << "\t" << 0 << "\t" << 0 << "\t" << 0 << std::endl;
    open.close();

    size_t npar = 0; 
    auto param = datapdf.getParams()->Parameters();
    for(size_t i = 0; i < param.size(); i++){
        if(!param[i].IsConst()){
            npar++;	
        }
    }


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
    //dplotter.chi2(npar,"D2PPP_DpAdpBinning.txt",0.05,1.95,0.3,3.4,*data,fmt::format("Fit/{0}",name.c_str())); 
    
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
    output = fmt::format("Fit/{0}/Swave.png",name.c_str());
    foo.SaveAs(output.c_str());

    return signal;
}


int main(int argc, char **argv){
    
    GooFit::Application app{"Genfit",argc,argv};
    
    std::string input_data_name = "input.root";
    std::string fit_name = "Fit";
	std::string acc_file = "acc_hist_0_Smoothed.root";
    bool save_toy = false;
    bool is_toy = false;

    auto fit = app.add_subcommand("fit","fit data");
    fit->add_option("-f,--file",input_data_name,"name_of_file.root");
    fit->add_option("-t,--isToy", is_toy, "Get toyData for fit") ;
    fit->add_option("-s,--saveToy",save_toy,"save toy in root file");
    fit->add_option("-n,--fitName",fit_name,"name of this fit(useful to save results)")->required(true);
    fit->add_option("--r1-Mass",rho1450_MASS,"rho 1450 Mass (default= 1.466 GeV)");
    fit->add_option("--r1-Width",rho1450_WIDTH,"rho 1450 Width (default= 0.4 GeV)");
    fit->add_option("--r2-Mass",rho1700_MASS,"rho 1700 Mass (default= 1.720 GeV)");
    fit->add_option("--r2-Width",rho1700_WIDTH,"rho 1700 Width (default= 0.25 GeV)");
    fit->add_option("-a,--acc",acc_file,"name of acc file");
   
    size_t Nevents=1000000;
    std::string toyName = "MC.root";
    auto makeToy = app.add_subcommand("makeToy","make a toy");
    makeToy->add_option("-e,--nevents",Nevents,"number of events");
    makeToy->add_option("-n,--name",toyName,"output_toy_name.root");
    makeToy->add_option("-s,--saveToy",save_toy,"save toy in root file");	

    std::string model_name = "model";
    int n_samples = 100;
    auto runRandom = app.add_subcommand("rdmFit","fit the same data set starting with diferent parameters states");
    runRandom->add_option("-m,--model",model_name,"model's name");
    runRandom->add_option("-n,--samples",n_samples,"number of samples");
    runRandom->add_option("-f,--file",input_data_name,"name_of_file.root");

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

    //const string bkgfile = "/home/juan/juan/work/DsPPP_Analysis/7-AccAndBkg/bkg_hist_15_95_uBoost_Smoothed.root";
    //const string efffile = "/home/juan/juan/work/DsPPP_Analysis/7-AccAndBkg/acc_hist_15_95_uBoost_Smoothed.root";
    const string bkgfile = "/home/juan/juan/work/DsPPP_Analysis/7-AccAndBkg/bkg_hist_15_95_Smoothed.root";
    const string efffile = "/home/juan/juan/work/DsPPP_Analysis/7-AccAndBkg/acc_hist_15_95_Smoothed.root";

    std::cout << "using acc " << efffile << '\n';

    const string bkghist = "h_eff";
    const string effhist = "h_eff";

    auto efficiency = makeHistogramPdf(efffile,effhist,s12,s13,true,false,false);
    auto background = makeHistogramPdf(bkgfile,bkghist,s12,s13,false,false,false);
    auto signal = makesignalpdf(s12, s13, eventNumber,efficiency); 
    AddPdf *prodpdf = new AddPdf("prodpdf", Variable("frac",0.95),signal, background) ;
    //readFromFile(prodpdf,"Fit/model1-parent/fit_results.txt");

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
        std::cout << "------------------------------------------" << std::endl;
        std::cout << "Num Entries Loaded =  " << data.getNumEvents()       << std::endl; 
        std::cout << "------------------------------------------" << std::endl;
        auto output_signal = runFit(prodpdf,signal, &data, fit_name, save_toy);
        auto flatMC = new UnbinnedDataSet({s12,s13,eventNumber});
        auto frac = signal->fit_fractions(4000);
        ofstream wt(fmt::format("Fit/{0}/fractions.txt",fit_name));
        std::cout << "Fit Fractions Interference" << '\n';
		auto nres = signal->getDecayInfo().resonances.size();
		Eigen::MatrixXd m(nres,nres);
		for(int i=0; i < nres; i++)
			m.row(i) = Eigen::Map<Eigen::VectorXd>(&frac[i][0],nres);

        std::cout << m << std::endl;
        wt << m << std::endl;

    }

    if(*runRandom){
	auto command = fmt::format("mkdir -p Fit_Random/{0}",fit_name);
	if(system(command.c_str()) != 0)
                throw GooFit::GeneralError("Making directory failed");

        UnbinnedDataSet data({s12, s13, eventNumber});
        getData(input_data_name, app, data,is_toy); 
        prodpdf->setData(&data);
	signal->setDataSize(data.getNumEvents());

 	FitManager fitter(prodpdf);
	fitter.setVerbosity(2);
        fitter.setMaxCalls(200000); 
	fitter.setRandMinuitValues(n_samples);
        
	for(int i=0;i<n_samples;i++){
		std::cout<< "################################################" << '\n';
		std::cout<< "Fitting model " << fit_name << " sample " << i << '\n';
		std::cout<< "################################################" << '\n';
		
		fitter.loadSample(i);
		fitter.printOriginalParams();
                auto min = fitter.fit();

		if(min.IsValid()){

			command = fmt::format("mkdir -p Fit_Random/{0}/Sample_{1}",model_name,i);
        			if(system(command.c_str()) != 0)

                	throw GooFit::GeneralError("Making directory failed");
			auto output = fmt::format("Fit_Random/{0}/Sample_{1}/fit_result.txt",model_name,i);
    			writeToFile(prodpdf, output.c_str());
    			std::ofstream open(output,std::ios_base::app);
    			open << "FCN" << "\t" << min.Fval() << "\t" << 0 << "\t" << 0 << "\t" << 0 << std::endl;
    			open << "Status" << "\t" << min.IsValid() << "\t" << 0 << "\t" << 0 << "\t" << 0 << std::endl;
    			open.close();
        		auto frac = signal->fit_fractions(4000);
        		ofstream wt(fmt::format("Fit_Random/{0}/Sample_{1}/fractions.txt",model_name,i));
        		std::cout << "Fit Fractions Interference" << '\n';
        		auto nres = signal->getDecayInfo().resonances.size();
        		Eigen::MatrixXd m(nres,nres);
            		for(int i=0; i < nres; i++)
                		m.row(i) = Eigen::Map<Eigen::VectorXd>(&frac[i][0],nres);

        		std::cout << m << std::endl;
        		wt << m << std::endl;
			
			std::cout<< "################################################" << '\n';
                	std::cout<< "Fitting model " << fit_name << " converged " << i << '\n';
                	std::cout<< "################################################" << '\n';
			std::cout << '\n';
		}else{
			std::cout<< "################################################" << '\n';
                        std::cout<< "Fitting model " << fit_name << "not converged " << i << '\n';
                        std::cout<< "################################################" << '\n';
			std::cout << '\n';
		continue;}
	}      
	 

    }
        
}

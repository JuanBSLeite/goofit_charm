#include <math.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TVector.h>


void compare_results(const int N,const int n){

std::ifstream r("files/PWACOEFS.txt");
double s[N];
double pr[N];
double pi[N];
double pr_e[N];
double pi_e[N];
double foo_e[N];

int foo;
double _s,_pr,_pi;
int j = 0;

printf("Initial Points \n");
while(r >> _s >> _pr >> _pi){
	s[j]=_s;
	pr[j]=_pr;
	pi[j]=_pi;
	pr_e[j]=0;
	pi_e[j]=0;
	foo_e[j]=0;
	
	if(j<10)
		printf("%d = (%lg,%lg,%lg,%lg,%lg) \n",j,s[j],pr[j],pi[j],pr_e[j],pi_e[j]); 

	 j++;
}

r.close();
string _name = "Fit/fit_parameters_"+to_string(n)+".txt";
std::ifstream fr(_name.c_str());
if(!fr.is_open()){
	cout << "nao abriu" << endl;
}else{
	cout << "abriu" << endl;
}

double fpr[N];
double fpi[N];
double fpr_e[N];
double fpi_e[N];
double _fpr,_fpi,_fpr_e,_fpi_e;
string _foo;

j=0;
int d,f=0;
printf("Fit Points \n");
while(fr >> _foo >> _fpr >> _fpr_e ){
        string n_r = "pwa_coef_"+to_string(d)+"_real";
	string n_i = "pwa_coef_"+to_string(f)+"_img";
	
	if(_foo.compare(n_r)==0){
		fpr[d]=_fpr;
		fpr_e[d]=_fpr_e;
		printf("%d = (%lg,%lg) \n",d,_fpr,_fpr_e);
		d++;
	}

	if(_foo.compare(n_i)==0){
		fpi[f]=_fpr;
		fpi_e[f]=_fpr_e;
		printf("%d = (%lg,%lg) \n",f,_fpr,_fpr_e);
		f++;
	}

		
}

fr.close();

TCanvas c;

TGraphErrors *gr1 = new TGraphErrors(N,s,pr,foo_e,pr_e);
gr1->SetMarkerStyle(4);
gr1->SetMinimum(-30);
gr1->SetMaximum(30);

TGraphErrors *gr2 = new TGraphErrors(N,s,fpr,foo_e,fpr_e);
gr2->SetLineColor(kRed);
gr2->SetMarkerColor(kRed);
gr2->SetMarkerStyle(6);

gr1->Draw("AP");
gr2->Draw("Psame");
TString name;
name.Form("plots/parte_real_%d.png",N);
c.SaveAs(name);

TGraphErrors *gr3 = new TGraphErrors(N,s,pi,foo_e,pi_e);
gr3->SetMarkerStyle(4);
gr3->SetMinimum(-30);
gr3->SetMaximum(30);

TGraphErrors *gr4 = new TGraphErrors(N,s,fpi,foo_e,fpi_e);
gr4->SetLineColor(kRed);
gr4->SetMarkerColor(kRed);
gr4->SetMarkerStyle(6);

gr3->Draw("AP");
gr4->Draw("Psame");
name.Form("plots/parte_imaginaria_%d.png",N);
c.SaveAs(name);



double cp,cpi;
int count_cp  = 0;
int count_icp = 0;
for(int i = 0 ; i < N ; i++){

	cp = abs(pr[i]-fpr[i])/fpr_e[i];
	if(cp <=2 ){
		printf("real_coef(%d) = compatible(%lg) \n",i,cp);
		count_cp++;
	}else{
		printf("real_coef(%d) = incompatible(%lg) \n",i,cp);	
		count_icp++;
	}

	cpi = abs(pi[i]-fpi[i])/fpi_e[i];
	if(cpi <=2 ){
		printf("img_coef(%d) = compatible(%lg) \n",i,cpi);
		count_cp++;
	}else{
		printf("img_coef(%d) = incompatible(%lg) \n",i,cpi);	
		count_icp++;
	}
}

printf("\n\n %d compatibles points and %d incompatible points \n\n",count_cp,count_icp);

/*
TCanvas foo;

TGraphErrors gr1("Parametros_iniciais.txt","%lg%lg%*lg%*lg%*lg");
gr1.Draw("AP*");
TGraphErrors gr2("Parametros_fit.txt","%lg%lg%*lg%lg%*lg");
gr2.SetMarkerColor(kRed);
gr2.SetLineColor(kRed);
gr2.Draw("Psame");
foo.SaveAs("plots/parte_imaginaria.png");


TGraphErrors gr3("Parametros_iniciais.txt","%lg%*lg%lg%*lg%lg");
gr3.Draw("AP*");
TGraphErrors gr4("Parametros_fit.txt","%lg%*lg%lg%*lg%lg");
gr4.SetMarkerColor(kRed);
gr4.SetLineColor(kRed);
gr4.Draw("Psame");
foo.SaveAs("plots/parte_real.png");


*/

}

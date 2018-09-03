{

#include <math.h>
#include <TMath.h>

const int N = 100;

double ev[N];
double amp[N];
double amp_e[N];
double phase[N];
double phase_e[N];
double ev_e[N];
ifstream r("Parametros_iniciais.txt");

double x,y,z,w,t = 0;

int i = 0;

while(r >> x >> y >> z >> w >> t){

	ev[i] = i;
	amp[i] = sqrt(z*z + y*y);
	amp_e[i] = sqrt( pow(w*y/sqrt(z*z+y*y),2) + pow(t*z/sqrt(z*z+y*y),2) );
	phase[i] = TMath::ATan(z/y);
	phase_e[i] = sqrt(pow(-w*z/(z*z + y*y),2) + pow(t*y/(z*z + y*y),2)) ;
	ev_e[i] = 0.0;

	i++;
	
}

r.close();

double fev[N];
double famp[N];
double famp_e[N];
double fphase[N];
double fphase_e[N];
double fev_e[N];

ifstream fr("Parametros_fit.txt");

double fx,fy,fz,fw,ft = 0;

int j = 0;

while(fr >> fx >> fy >> fz >> fw >> ft){
	
	fev[j] = j;
	famp[j] = sqrt(fz*fz + fy*fy);
	famp_e[j] = sqrt( pow(fw*fy/sqrt(fz*fz+fy*fy),2) + pow(ft*fz/sqrt(fz*fz+fy*fy),2) );
	fphase[j] = TMath::ATan(fz/fy);
	fphase_e[j] = sqrt(pow(-fw*fz/(fz*fz + fy*fy),2) + pow(ft*fy/(fz*fz + fy*fy),2)) ;
	fev_e[j] = 0.0;

	j++;

}

fr.close();


double shift_amp[N];
double shift_amp_e[N];
double shift_phase[N];
double shift_phase_e[N];

for(int k = 0; k < N ; k++){

    shift_amp[k] = amp[k]-famp[k];
    shift_amp_e[k] = sqrt(famp_e[k]*famp_e[k] + amp_e[k]*amp_e[k]);
    shift_phase[k] = phase[k] - fphase[k];
    shift_phase_e[k] = sqrt( phase_e[k]*phase_e[k] + fphase_e[k]*fphase_e[k]);
}

TCanvas *c = new TCanvas("c","",900,500);


c->Divide(2,2);

c->cd(1);

TGraphErrors *gr = new TGraphErrors(N,ev,amp,ev_e,amp_e);
gr->SetMarkerStyle(20);
gr->SetMarkerColor(kRed);
gr->SetTitle("Magnetude");
gr->GetYaxis()->SetTitle("");
gr->GetXaxis()->SetTitle("Event");
gr->SetMarkerSize(0.7);
TGraphErrors *fgr = new TGraphErrors(N,fev,famp,ev_e,famp_e);
fgr->SetMarkerStyle(20);
fgr->SetMarkerColor(kViolet);
fgr->SetMarkerSize(0.7);

gr->Draw("AP");
fgr->Draw("Psame");

c->cd(2);

TGraphErrors *gi = new TGraphErrors(N,ev,phase,ev_e,phase_e);
gi->SetMarkerStyle(20);
gi->SetMarkerColor(kRed);
gi->SetTitle("Phase");
gi->GetYaxis()->SetTitle("");
gi->GetXaxis()->SetTitle("Event");
gi->SetMinimum(-1);
gi->SetMarkerSize(0.7);

TGraphErrors *fgi = new TGraphErrors(N,fev,fphase,ev_e,fphase_e);
fgi->SetMarkerStyle(20);
fgi->SetMarkerColor(kViolet);
fgi->SetMarkerSize(0.7);

gi->Draw("AP");
fgi->Draw("Psame");

c->cd(3);

TGraphErrors *gsamp = new TGraphErrors(N,ev,shift_amp,ev_e,shift_amp_e);
gsamp->SetMarkerStyle(20);
gsamp->SetMarkerColor(kRed);
gsamp->SetTitle("Shift of S-Wave Magnetude");
gsamp->GetYaxis()->SetTitle("");
gsamp->GetXaxis()->SetTitle("Event");
gsamp->SetMarkerSize(0.7);

gsamp->Draw("AP");

c->cd(4);

TGraphErrors *gsphase = new TGraphErrors(N,ev,shift_phase,ev_e,shift_phase_e);
gsphase->SetMarkerStyle(20);
gsphase->SetMarkerColor(kRed);
gsphase->SetTitle("Shift of S-Wave Phase");
gsphase->GetYaxis()->SetTitle("");
gsphase->GetXaxis()->SetTitle("Event");
gsphase->SetMarkerSize(0.7);
gsphase->Draw("AP");

c->SaveAs("Parameters.png");
}

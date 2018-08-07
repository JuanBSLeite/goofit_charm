{
const int N = 50;

double ev[N];
double real[N];
double real_e[N];
double img[N];
double img_e[N];
double ev_e[N];

ifstream r("Parametros_iniciais.txt");

double x,y,z,w,t = 0;

int i = 0;

while(r >> x >> y >> z >> w >> t){

	ev[i] = i;
	real[i] = y;
	real_e[i] = w;
	img[i] = z;
	img_e[i] = t;
	ev_e[i] = 0.0;

	i++;
	
}

r.close();

double fev[N];
double freal[N];
double freal_e[N];
double fimg[N];
double fimg_e[N];

ifstream fr("Parametros_fit.txt");

double fx,fy,fz,fw,ft = 0;

int j = 0;

while(fr >> fx >> fy >> fz >> fw >> ft){
	
	fev[j] = j;
	freal[j] = fy;
	freal_e[j] = fw;
	fimg[j] = fz;
	fimg_e[j] = ft;

	j++;

}

fr.close();

TCanvas *c = new TCanvas("c","",900,500);

c->Divide(2,1);

c->cd(1);

TGraphErrors *gr = new TGraphErrors(N,ev,real,ev_e,real_e);
gr->SetMarkerStyle(20);
gr->SetMarkerColor(kRed);
gr->SetTitle("Parte Real");
gr->GetYaxis()->SetTitle("Valor");
gr->GetXaxis()->SetTitle("Evento");

TGraphErrors *fgr = new TGraphErrors(N,fev,freal,ev_e,freal_e);
fgr->SetMarkerStyle(20);
fgr->SetMarkerColor(kViolet);

gr->Draw("AP");
fgr->Draw("Psame");

c->cd(2);

TGraphErrors *gi = new TGraphErrors(N,ev,img,ev_e,img_e);
gi->SetMarkerStyle(20);
gi->SetMarkerColor(kRed);
gi->SetTitle("Parte Imagin#acute{a}ria");
gi->GetYaxis()->SetTitle("Valor");
gi->GetXaxis()->SetTitle("Evento");
gi->SetMinimum(-1);

TGraphErrors *fgi = new TGraphErrors(N,fev,fimg,ev_e,fimg_e);
fgi->SetMarkerStyle(20);
fgi->SetMarkerColor(kViolet);

gi->Draw("AP");
fgi->Draw("Psame");


c->SaveAs("Parameters.png");
}

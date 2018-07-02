{

double pi_MASS  = 0.13957018; //GEV
double D_MASS   = 1.86959; //GEV

double d1_MASS  = pi_MASS;  //daughter 1 mass
double d2_MASS  = pi_MASS;
double d3_MASS  = pi_MASS;

double s12_min = (d1_MASS  + d2_MASS)*(d1_MASS + d2_MASS);
double s12_max = (D_MASS   - d2_MASS)*(D_MASS - d2_MASS);

double m12_min = s12_min;
double m12_max = s12_max;

TH1F *h0 = new TH1F("h0","",100,s12_min,s12_max);
h0->SetLineColor(kRed);
TH1F *h1 = new TH1F("h1","",100,s12_min,s12_max);

TTree t("t","");
t.ReadFile("D2PPP_toy_iso.txt","x:y:z");

TTree t2("t2","");
t2.ReadFile("D2PPP_toy_pwa.txt","x:y:z");


t.Draw("y>>h0");
t2.Draw("(y)>>h1","","same");

}

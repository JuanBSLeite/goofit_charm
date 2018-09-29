{
gStyle->SetOptStat(0);

double pi_MASS  = 0.13957018; //GEV
double D_MASS   = 1.86959; //GEV

double d1_MASS  = pi_MASS;  //daughter 1 mass
double d2_MASS  = pi_MASS;
double d3_MASS  = pi_MASS;

double s12_min = (d1_MASS  + d2_MASS)*(d1_MASS + d2_MASS);
double s12_max = (D_MASS   - d2_MASS)*(D_MASS - d2_MASS);

double m12_min = s12_min;
double m12_max = s12_max;

TCanvas *c1= new TCanvas("c1","",700,500);
c1->SetGridy();
TH1F *h0 = new TH1F("h0","",120,s12_min,s12_max);
h0->SetLineColor(kBlue);
h0->SetFillColor(kBlue);
h0->SetTitle("Ratio");
//h0->SetFillStyle(3001);

TH1F *h1 = new TH1F("h1","",120,s12_min,s12_max);
h1->SetLineColor(kRed);
h1->SetLineWidth(2);
//h1->SetFillColorAlpha(kBlue,0.60);
//h1->SetFillStyle(3001);

TTree t("t","");
t.ReadFile("D2PPP_toy_iso.txt","x:y:z");

TTree t2("t2","");
t2.ReadFile("D2PPP_toy_pwa.txt","x:y:z");

t.Draw("y>>h0");
t2.Draw("(y)>>h1","","same");

h0->GetXaxis()->SetTitle("m^{2}(#pi^{-}#pi^{+}) GeV");
h0->GetYaxis()->SetTitle("");
h0->GetYaxis()->SetTitleOffset(1.5);


TLegend *lg = new TLegend(0.8,0.4,1.0,0.6);
lg->AddEntry(h0,"Isobar","l");
lg->AddEntry(h1,"MIPWA","l");
lg->Draw(); 

h0->Sumw2();
h0->Divide(h1);
h0->Draw("E");


c1->SaveAs("toyDataComparacao_ratio.png");


}

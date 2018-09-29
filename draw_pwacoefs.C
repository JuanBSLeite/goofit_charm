{

gStyle->SetOptStat(0);

double pi_MASS  = 0.13957018; //GEV
double D_MASS   = 1.86959; //GEV

double d1_MASS  = pi_MASS;  //daughter 1 mass
double d2_MASS  = pi_MASS;
double d3_MASS  = pi_MASS;

double s12_min = (d1_MASS  + d2_MASS)*(d1_MASS + d2_MASS);
double s12_max = (D_MASS   - d2_MASS)*(D_MASS - d2_MASS);

TCanvas c;
c.SetGridy();
c.Divide(2,1);

c.cd(1);

TH2F *h1 = new TH2F("h1","Magnitude",100,s12_min,s12_max,100,0,160);
h1->SetTitle("Magnitude");
h1->SetMarkerStyle(20);
h1->SetMarkerColor(kBlue);
h1->SetMarkerSize(0.7);
h1->GetXaxis()->SetTitle("m^{2}(#pi^{-}#pi^{+}) GeV");
h1->GetYaxis()->SetTitle("Magnitude");
h1->GetYaxis()->SetTitleOffset(1.5);

TTree t("t","");
t.ReadFile("files/PWACOEFS.txt","EventNumber:s12:s13");
t.Draw("(s13*s13+s12*s12):EventNumber>>h1","","");

c.cd(2);

TH2F *h2 = new TH2F("h2","Phase",100,s12_min,s12_max,100,0.5,3.5);
h2->SetTitle("Phase");
h2->SetMarkerStyle(20);
h2->SetMarkerSize(0.7);
h2->SetMarkerColor(kBlue);
h2->GetXaxis()->SetTitle("m^{2}(#pi^{-}#pi^{+}) GeV");
h2->GetYaxis()->SetTitle("Phase");
h2->GetYaxis()->SetTitleOffset(1.5);

t.Draw("(TMath::ATan2(s13,s12)):EventNumber>>h2","","");

c.SaveAs("PWACOEFs.png");

}

{

gStyle->SetOptStat(0);

double pi_MASS  = 0.13957018; //GEV
double D_MASS   = 1.86959; //GEV

double d1_MASS  = pi_MASS;  //daughter 1 mass
double d2_MASS  = pi_MASS;
double d3_MASS  = pi_MASS;

double s12_min = (d1_MASS  + d2_MASS)*(d1_MASS + d2_MASS)*(1+0.001);
double s12_max = (D_MASS   - d2_MASS)*(D_MASS - d2_MASS)*(1+0.001);

TCanvas c;
c.SetGridy();
c.Divide(2,1);

c.cd(1);

TGraph t("files/PWACOEFS.txt","%lg%lg%*lg");
t.Draw("AP*");

c.cd(2);

TGraph t1("files/PWACOEFS.txt","%lg%*lg%lg");
t1.Draw("AP*");

c.SaveAs("PWACOEFs.png");

}

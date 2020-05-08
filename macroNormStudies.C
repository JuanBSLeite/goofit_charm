void macroNormStudies(int bins_i = 400, int bins_f = 4000 , int step=100)
{
	TFile output("NormStudies/NormStudies.root","recreate");
	TTree *tree = new TTree("tree","");
 
	auto file = Form("NormStudies/fit_bins_%d.txt",bins_i);
	std::cout << file << '\n';
	std::ifstream read(file);

	if(read.is_open())
		std::cout << "First file opened!" << '\n';

	const size_t N = 116;
	double value[N] = {0};
	double error[N] = {0};
	
	double foo ;
	std::string sfoo;
	size_t index = 0;	
	int Bins;

	tree->Branch("Bins",&Bins);
	
	while(read >> sfoo >> foo >> foo)
	{
		if(sfoo=="frac" || sfoo=="effSmoothing_eff" || sfoo=="Smoothing_bkg" || sfoo=="Min_FCN" || sfoo=="Norm" || sfoo=="Status")
		{
			tree->Branch(sfoo.c_str(), value + index);
			index++;
		}else
		{	
			tree->Branch(sfoo.c_str(), value + index);
			tree->Branch(Form("%s_error",sfoo.c_str()), error + index);
			index++;
		}
	}

	double _value;
	double _error;

	for(int i = bins_i; i <= bins_f; i+=step)
	{
		file = Form("NormStudies/fit_bins_%d.txt",i);
		std::ifstream read(file);
		index = 0;

		while(read >> sfoo >> _value >> _error)
		{
			if(sfoo=="frac" || sfoo=="effSmoothing_eff" || sfoo=="Smoothing_bkg" || sfoo=="Min_FCN" || sfoo=="Norm" || sfoo=="Status")
			{
				Bins = i;
				value[index] = _value;
				index++;
			}else
			{
				Bins = i;
				value[index] = _value;
				error[index] = _error;
				index++;
			}
		}		
		
		tree->Fill();
		read.close();
	} 

	tree->Print();
	tree->Write("",TObject::kOverwrite);
	output.Close();

}

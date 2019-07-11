// A simple fit in energy for signal and a background
#include <string>
#include <vector>
#include <math.h>
#include <Rand.h>
#include <fstream>
#include <iostream>

#include <TH1D.h>
#include <TH2D.h>
#include <TPaveText.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <ROOTNtuple.h>
#include <TFile.h>
#include <TMath.h>
#include <string>
#include <TNtuple.h>
#include <iostream>
#include <TObject.h>
#include <TRandom3.h>
#include <TVector3.h>

#include <BinnedED.h>
#include <BinnedEDGenerator.h>
#include <SystematicManager.h>
#include <BinnedNLLH.h>
#include <FitResult.h>
#include <Minuit.h>
#include <DistTools.h>
#include <Minuit.h>
#include <Convolution.h>
#include <Scale.h>
#include <BoolCut.h>
#include <BoxCut.h>
#include <Gaussian.h>
#include <ParameterDict.h>
#include <ContainerTools.hpp>
#include <NuOsc.h>
#include <SurvProb.h>
#include <TH1D.h>
#include <ROOTMultiPlot.h>
#include <TRandom3.h>

#include <TVectorD.h>

#include <locale>
#include <boost/algorithm/string.hpp>
int unoscPdfint;
int oscPdfint;

//Surv Probablity function for antinu of KE nuE and distance traveled of baseline.
double probability(double nuE, double baseline, double delmsqr21, double sinsqrtheta12, double sinsqrtheta13) {
  double fSSqr2Theta12 = pow(sin(2.0 * TMath::ASin(sqrt(sinsqrtheta12))), 2.0);
  double fS4 = pow(sinsqrtheta13, 2.0);
  double fC4 = pow(1.0-sinsqrtheta13, 2.0);
  double scale = 1.267e3; // for nuE in [MeV] and baseline in [km]
  double sSqrDmBE = pow(sin(scale * delmsqr21 * baseline / nuE), 2.0);
  double fOscProb = (fC4 * (1.0 - fSSqr2Theta12 * sSqrDmBE) + fS4);
  return fOscProb;
}
  
//Coincidence tagging - function to apply cuts to isolate positron and neutron events
// At present, looking at the MC event number and EV index number to ensure no pairs are missed
void OscPromptE (const std::string infile, const std::string outfile,  double delmsqr21, double sinsqrtheta12, double sinsqrtheta13, double Distance) {
  
  unoscPdfint = 0;
  oscPdfint = 0;

  TFile *fin = TFile::Open(infile.c_str());
  TNtuple *T = (TNtuple *) fin->Get("nt");

  TFile *fout = new TFile(outfile.c_str(), "RECREATE");
  TNtuple* ntout = new TNtuple("nt","nt","E1");
  float parke, Energy;
  
  T->SetBranchAddress("ParKE", &parke);
  T->SetBranchAddress("E1", &Energy);
  //T->SetBranchAddress("ReactorDistance", &ReactorDistance); // Average distance taken from info file, not event by event (though it can be)

  double survprob;
  TRandom3 *r1 = new TRandom3();
  r1->SetSeed(0);    
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    
    survprob = probability(parke, Distance, delmsqr21, sinsqrtheta12, sinsqrtheta13);
    unoscPdfint +=1;
    if (survprob > r1->Rndm()){
      oscPdfint +=1;
      ntout->Fill(Energy);
    }
  }

  ntout->Write();
  fin->Close();
  fout->Close();
  delete fin;
  delete fout;
}

void readInfoFile(const std::string &runInfoFileName, std::vector<std::string> &reactorNames, std::vector<double> &distances, std::vector<std::string> &reactorTypes, std::vector<int> &nCores, std::vector<double> &powers ) {
  // Read couchDB run-info text file
  std::ifstream in;
  in.open(runInfoFileName.c_str());
  std::cout << "opening file: " << runInfoFileName.c_str() << std::endl;

  std::fill(reactorNames.begin(), reactorNames.end(), "");
  std::fill(distances.begin(), distances.end(), 0.);
  std::fill(reactorTypes.begin(), reactorTypes.end(), "");
  std::fill(nCores.begin(), nCores.end(), 0);
  std::fill(powers.begin(), powers.end(), 0.);

  std::string reactorName,distance,reactorType,nCore,power;
  int lineNo = 0;

  // read until end of file.
  while(in.good()){
    std::getline(in,reactorName,',');
    std::getline(in,distance,',');
    std::getline(in,reactorType,',');
    std::getline(in,nCore,',');
    std::getline(in,power,'\n');

    if (lineNo>0){ //skip csv header
      if (strcmp(reactorName.c_str(),"")!=0) {
	reactorNames.push_back(reactorName);
	distances.push_back(atof(distance.c_str()));
	reactorTypes.push_back(reactorType.c_str());
	nCores.push_back(atoi(nCore.c_str()));
	powers.push_back(atof(power.c_str()));

	std::cout << "v: reactorName: " << reactorNames[lineNo-1] << ", distance: " << distances[lineNo-1] << ", reactorType: " << reactorTypes[lineNo-1] << ", nCore: " << nCores[lineNo-1] << ", power: " << powers[lineNo-1] << std::endl; //debug check ('-1' for header)
      }
    }
    lineNo++;
  }
  in.close();
}

int main(int argc, char *argv[]) {
  int FakeDataArgc = 9 + 1;
  int NoFakeDataArgc = 7 + 1;

  if (argc < NoFakeDataArgc) {
    std::cout<<"7 (9) arguments expected: \n \n 1: Data Emin \n 2: DataEmax \n 3: numbins \n \n 4: location/filename for reactor info file \n \n 5: filename of rootfile containting constraints and Data (no .root!) \n 6: temporary root file to fill ntuples (can delete afterwards) \n \n 7: location/filename for data spectrum to be fit (pruned ntuple with just one entry, EPrompt) \n (7): Fake data at d21 \n (8):Fake data at s12 \n (9): location of UnOsc pruned Ntuple file "<<std::endl;
  }
  else{
    const std::string &infoFile = argv[4];
    
    std::vector<std::string> reactorNames;
    std::vector<double> distances;
    std::vector<std::string> reactorTypes;
    std::vector<int> nCores;
    std::vector<double> powers;
    std::cout<<"starting filling info....... \n"<<std::endl;
    readInfoFile(infoFile, reactorNames, distances, reactorTypes, nCores, powers); // get reactor information
    std::cout<<"\n finished filling info"<<std::endl;
    int numPdfs = reactorNames.size();

    std::cout<<"Number of reactors: "<<numPdfs<<std::endl;
    
    for (int i = 0; i< numPdfs; i++){
      std::cout<<"Reactor name: "<<reactorNames[i]<<" Distance:"<<distances[i]<<std::endl;
    }

    //axes range and number of bins, must not include any bins with zero entries!!!
    double Emin = atof(argv[1]);
    double Emax = atof(argv[2]);
    int numbins = atoi(argv[3]);
    
    //SET UP DATA PDFS
    // Only interested in first bit of data ntuple
    ObsSet dataRep(0);
    // Set up binning
    AxisCollection axes;
    axes.AddAxis(BinAxis("EPrompt (MeV)", Emin, Emax, numbins));

    BinnedED dataSetPdf("dataSetPdf",axes);
    dataSetPdf.SetObservables(dataRep);

    const std::string &dataconstraintoutFile = argv[5];
    std::stringstream ssout;

    std::vector<double> normconstraints;
    std::vector<double> normconstraintsuncert;

    // Reading constraints from csv file
    ifstream ip;
    ip.open("passes_kamland_rat6169_flux1_penergy_cleanround1_plots.csv");
    //ip.open("3CAD_200x1yrrat6169_un.csv");
    if (!ip.is_open()) std::cout<<"Error in opening csv file"<<std::endl;

    std::string reactorname_csv;
    std::string fit_norm_csv;
    std::string fit_norm_err_csv;
    std::string fit_norm_sigma_csv;
    std::string fit_norm_sigma_err_csv;

    std::vector<std::string> reactornamevecUnSorted;
    std::vector<double> fitnormmeanvecUnSorted;
    std::vector<double> fitnormsigmavecUnSorted;

    double totalnormmean;
    double totalsigma;
    
    std::string line;
    std::string cell;
    while (getline(ip,line,'\n')){
      int i = 0;
      std::stringstream linestream(line);
      while(getline(linestream,cell,',')){
	if (i == 0)
	  reactorname_csv = cell;
	if (i == 1)
	  fit_norm_csv = cell;
	if (i == 2)
	  fit_norm_err_csv = cell;
	if (i == 3)
	  fit_norm_sigma_csv = cell;
	if (i == 4)
	  fit_norm_sigma_err_csv = cell;
	i += 1;
      }
      if (reactorname_csv != "reactor"){
	if (reactorname_csv == "combinedall"){
	  totalnormmean = atof(fit_norm_csv.c_str());
	  totalsigma = atof(fit_norm_sigma_csv.c_str());
	}else{
	  reactornamevecUnSorted.push_back(reactorname_csv);
	  fitnormmeanvecUnSorted.push_back(atof(fit_norm_csv.c_str()));
	  fitnormsigmavecUnSorted.push_back(atof(fit_norm_sigma_csv.c_str()));
	  std::cout<<reactorname_csv<<"  "<<fit_norm_csv<<"  "<<fit_norm_err_csv<<"  "<<fit_norm_sigma_csv<<" "<<fit_norm_sigma_err_csv<<std::endl;
	}
      }
    }
    ip.close();
    
    std::cout<<"Total Combined reactors fit norm: "<<totalnormmean<<" total sigma: "<<totalsigma<<std::endl;

    if (reactornamevecUnSorted.size() != numPdfs){
      std::cout<<"number data norm constraints dont match number of Pdfs!!"<<std::endl;
      return 0;
    }

    for (int i = 0; i < reactorNames.size(); i++){
      std::string NAME = reactorNames[i];
      std::string name = boost::algorithm::to_lower_copy(NAME);

      std::vector<std::string>::iterator it = std::find(reactornamevecUnSorted.begin(), reactornamevecUnSorted.end(), name);
      
      if (it == reactornamevecUnSorted.end())
	std::cout << "Element Not Found" << std::endl;
      
      int index = std::distance(reactornamevecUnSorted.begin(), it);

      normconstraints.push_back(fitnormmeanvecUnSorted[index]);
      //double sigma_over_norm_sqrd = pow(fitnormsigmavecUnSorted[index]/fitnormmeanvecUnSorted[index],2);
      //double total_sigma_over_norm_sqrd = pow(totalsigma/totalnormmean,2);
   
      //normconstraintsuncert.push_back(sqrt(total_sigma_over_norm_sqrd+sigma_over_norm_sqrd));
      normconstraintsuncert.push_back(fitnormsigmavecUnSorted[index]/fitnormmeanvecUnSorted[index]);
      std::cout<<"reactorname "<<reactorNames[i]<<" norm "<<normconstraints[i]<<" sigma/norm "<<normconstraintsuncert[i]<<std::endl;
    }
    

    double antinueffic = 1.;//0.85;

    double FitConstraintsTotalInt = 0.;
    std::cout<<"data vector length: "<<normconstraints.size()<<std::endl;
    for (int i = 0; i < normconstraints.size(); i++){
      normconstraints[i] = (normconstraints[i]*antinueffic);
      FitConstraintsTotalInt += normconstraints[i];
    }
    std::cout<<"Fit Constraints Int: "<<FitConstraintsTotalInt<<std::endl;

    if (normconstraints.size() != numPdfs){
      std::cout<<"number data norm constraints dont match number of Pdfs!!"<<std::endl;
      return 0;
    }
    
    const std::string &tempFile = argv[6];
    
    std::vector<double> databincontents;

    if (argc == NoFakeDataArgc){
      const std::string &dataFile = argv[7];
      ROOTNtuple dataNtp(dataFile, "nt");
 
      for(size_t i = 0; i < dataNtp.GetNEntries(); i++)
	dataSetPdf.Fill(dataNtp.GetEntry(i));
	
      databincontents = dataSetPdf.GetBinContents();
      ssout<<dataconstraintoutFile<<".root";
    }
    else if(argc == FakeDataArgc){
      double FakeDatad21 = atof(argv[7]);
      double FakeDatas12 = atof(argv[8]);
      const std::string &UnOscfile = argv[9];

      ssout<<dataconstraintoutFile<<"ds21_"<<FakeDatad21<<"s12_"<<FakeDatas12<<".root";

      
      //Using data from paper (digitised)
      double OscDatacontentarray[15] = {12.5397, 89.0563, 118.848, 120.102, 100.032, 144.249, 158.047, 130.764, 139.858, 104.109, 63.6553, 41.3902, 29.4737, 10.6581, 8.46297};

      if (numbins != 15){
	std::cout<<"wrong number of bins!!"<<std::endl;
	return 0;
      }
      
      double DataInt = 0.;
      
      std::cout<<"\n  Digitised KL Paper binc contents: "<<std::endl;
      for (int i = 0; i < 15; i++){
	databincontents.push_back(OscDatacontentarray[i]);
	std::cout<<"i="<<i<<"  "<<databincontents[i]<<std::endl;
	DataInt += OscDatacontentarray[i];
      }
       
      std::cout<<"Integral of measured (oscillated) KL Data: "<<DataInt<<std::endl;
      //end of Kl digitised data
      
      /*
      // Generating fake data for 3CAD (given chosen osc params that are not told to the fitter)
      for (int i = 0; i< numPdfs; i++){
	
	OscPromptE(UnOscfile, tempFile,FakeDatad21,FakeDatas12,0.,distances[i]);   //oscillate each pdf -> outputs a pruned ntuple with an oscillated EPrompt positron spectrum called tempFile, to be used below to fill a reactorPdf

	BinnedED *dataPdf = new BinnedED(reactorNames[i],axes);
	dataPdf->SetObservables(0);
    
	ROOTNtuple dataNtp(tempFile, "nt");
	for(size_t j = 0; j < dataNtp.GetNEntries(); j++)
	  dataPdf->Fill(dataNtp.GetEntry(j));

	std::cout<<"pre osc Int: "<<unoscPdfint<<std::endl;   //integral of input unoscillated pdfs
	//int postoscint = (int)(reactorPdf->Integral());
	std::cout<<"post osc Int: "<<oscPdfint<<std::endl;  //integral of pdf after KE->EPrompt + oscillation
	double osc_loss = oscPdfint/(double)unoscPdfint;      // calculate loss in pdf integral
	std::cout<<" osc loss factor: "<<osc_loss<<std::endl;
      
	dataPdf->Normalise();
	dataPdf->Scale(normconstraints[i]*osc_loss);
  
	dataSetPdf.Add(*dataPdf,1);
      }
      
      databincontents = dataSetPdf.GetBinContents();
      // End of fake data
      */      

      
    }
    else{
      std::cout<<"Wrong number of arguments!!"<<std::endl;
      return 0;
    }	
    
    //TH1D DataHist;
    //DataHist = DistTools::ToTH1D(dataSetPdf);
    //DataHist.SetName("DataHist");
    /*
    double s12min = atof(argv[5]);
    double s12max = atof(argv[6]);
    double d21min = atof(argv[8]);//1e-5;
    double d21max = atof(argv[9]);//20e-5;
    
    int num_s12 = atoi(argv[7]);
    int num_d21 = atoi(argv[10]);
    
    TH2D *h2 = new TH2D("lh2d", "lh2d", num_s12,s12min,s12max,num_d21,d21min,d21max);
    
    h2->GetXaxis()->SetTitle("(sintheta12)^2");
    h2->GetXaxis()->SetTitleSize(0.05);
    h2->GetXaxis()->SetTitleOffset(0.9);
    h2->GetXaxis()->CenterTitle();
    h2->GetYaxis()->SetTitle("(delm21)^2 (eV^2)");
    h2->GetYaxis()->SetTitleSize(0.05);
    h2->GetYaxis()->SetTitleOffset(1);
    h2->GetYaxis()->CenterTitle();
    h2->SetLabelSize(0.05,"xyz");
    h2->SetTitle("Best Fit LH value vs Osc Parameters");
    std::cout<<"\n \n here \n \n"<<std::endl;
    */
    
    TFile *fileOut = new TFile(ssout.str().c_str(),"RECREATE");
    TNtuple* datacontentNT = new TNtuple("databincontents","databincontents","counts");
    for (int i = 0; i< databincontents.size(); i++)
      datacontentNT->Fill(databincontents[i]);
    
    TNtuple* normconstraintsNT = new TNtuple("normconstraints","normconstraints","norm:normuncert");
    for (int i = 0; i< numPdfs; i++)
      normconstraintsNT->Fill(normconstraints[i],normconstraintsuncert[i]);
    
    datacontentNT->Write();
    normconstraintsNT->Write();

    fileOut->Close();
  }
}


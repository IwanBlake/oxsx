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

bool badfit = true;

int unoscPdfint;
int oscPdfint;

ParameterDict bestFit;

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
    //if (survprob > r1->Rndm()){
    if (survprob > -1.){
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

double LHFit(const std::string UnOscfile, const std::string tempFile, AxisCollection axes, BinnedED dataSetPdf, int numPdfs, std::vector<std::string> reactorNames, std::vector<double> reactorDistances, double d21fix, double s12fix,std::vector<double> normconstraints,int Normmin, int Normmax, std::vector<double> normconstraintsuncert){
  
  char name[100];
  TRandom3 *r1 = new TRandom3();
  r1->SetSeed(0);

  ////////////////////
  // 1. Set Up PDFs //
  ////////////////////
  
  ParameterDict minima;
  ParameterDict maxima;
  ParameterDict initialval;
  ParameterDict initialerr;
 
  BinnedNLLH lhFunction;
  lhFunction.SetBufferAsOverflow(true);
  int Buff = 0;
  lhFunction.SetBuffer(0,Buff,Buff);  
  lhFunction.SetDataDist(dataSetPdf); // initialise withe the data set
  double rand = r1->Rndm();

  for (int i = 0; i < numPdfs; i++){
    OscPromptE(UnOscfile,tempFile ,d21fix,s12fix,0.0215,reactorDistances[i]);   //oscillate each input reactor KE pdfs -> outputs a pruned ntuple with an oscillated EPrompt positron spectrum called tempFile, to be used below to fill a reactorPdf
    
    BinnedED *reactorPdf = new BinnedED(reactorNames[i],axes);
    reactorPdf->SetObservables(0);
    
    ROOTNtuple reactorNtp(tempFile, "nt");
    for(size_t j = 0; j < reactorNtp.GetNEntries(); j++)
      reactorPdf->Fill(reactorNtp.GetEntry(j));
    //std::cout<<"pre osc Int: "<<unoscPdfint<<std::endl;   //integral of input unoscillated pdfs
    //std::cout<<"post osc Int: "<<oscPdfint<<std::endl;  //integral of pdf after KE->EPrompt + oscillation
    double osc_loss = oscPdfint/(double)unoscPdfint;      // calculate loss in pdf integral
    //std::cout<<" osc loss factor: "<<osc_loss<<std::endl;
    reactorPdf->Normalise();

    // Setting optimisation limits
    sprintf(name,"%s_norm",reactorNames[i].c_str());
    minima[name] = Normmin;
    maxima[name] = Normmax;
    initialval[name] = (rand*(maxima[name]-minima[name]))+minima[name];
    initialerr[name] = 0.5*initialval[name];
    
    lhFunction.AddDist(*reactorPdf);
    //double constraint = normconstraints[i] * osc_loss; //vector of unoscillated norm contraints made in main() Jeff you can probs improve this
    double constraint = normconstraints[i]; //vector of unoscillated norm contraints made in main() Jeff you can probs improve this
    std::cout<<reactorNames[i]<<" constraint: "<<constraint<<" uncert: "<<(normconstraintsuncert[i]*constraint)<<std::endl;
    lhFunction.SetConstraint(name, constraint, (normconstraintsuncert[i]*constraint));  // I just randomly picked 10% for contraint uncertainty 
  }
  
  ////////////
  // 4. Fit //
  ////////////
  Minuit min;
  min.SetMethod("Migrad");
  min.SetMaxCalls(10000000);
  min.SetTolerance(0.001);
  min.SetMinima(minima);
  min.SetMaxima(maxima);
  min.SetInitialValues(initialval);
  min.SetInitialErrors(initialerr);
  
  /////////////////////////////////////////////
  ////////        Fit Result        ///////////
  /////////////////////////////////////////////
  FitResult fitResult = min.Optimise(&lhFunction);
  bestFit = fitResult.GetBestFit();

  bool fitValid = fitResult.GetValid();
  if (fitValid){
    fitResult.Print();
    badfit = false;
    std::cout<<"goodfit!"<<std::endl;
  }else{
    std::cout<<"INVALID FIT!! \n Trying Again:"<<std::endl;
    badfit = true;
  }
  
  lhFunction.SetParameters(bestFit);
  double lhval =(-1)*lhFunction.Evaluate();
  
  printf("\n-------------------------------------------------------------------------------------\n");
  std::cout<<"LH val at best fit: "<<lhval<<std::endl;
  printf("-------------------------------------------------------------------------------------\n");
	
  return lhval;
}

int main(int argc, char *argv[]) {
  int FakeDataArgc = 9 + 1;
  //int NoFakeDataArgc = 6 + 1;

  if (argc < FakeDataArgc) {
    std::cout<<"9 arguments expected: \n \n 1: Data Emin \n 2: DataEmax \n 3: numbins \n \n 4: location/filename for reactor info file \n \n 5: temporary root file to fill ntuples (can delete afterwards) \n \n 6: location/filename for data spectrum to be fit (pruned ntuple with just one entry, EPrompt) \n 7: Fake data at d21 \n 8:Fake data at s12 \n 9: location of UnOsc pruned Ntuple file "<<std::endl;
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

    const std::string &tempFile = argv[5];
    const std::string &UnOscfile = argv[9];
    double FakeDatad21= atof(argv[7]);
    double FakeDatas12= atof(argv[8]);

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
    //axes.AddAxis(BinAxis("EPrompt (MeV)",0.425*2,8.075,17));

    BinnedED dataSetPdf("dataSetPdf",axes);
    dataSetPdf.SetObservables(dataRep);
    
    
    const std::string &dataFile = argv[6];

    //ROOTNtuple dataNtp(dataFile, "nt");
    //for(size_t i = 0; i < dataNtp.GetNEntries(); i++)
    //dataSetPdf.Fill(dataNtp.GetEntry(i));

    //TFile *fin = TFile::Open(dataFile.c_str(),"READ");
    //TNtuple *datacontentNT = (TNtuple*) fin->Get("data");
    //float Count;
    //datacontentNT->SetBranchAddress("counts", &Count);
    //std::vector<double> Datacontent;
    
    //for(int i = 0; i < datacontentNT->GetEntries(); i++){
    //datacontentNT->GetEntry(i);
    //Datacontent.push_back(Count);
    //}

    //PAPER 1
    //std::vector<double> Datacontent;
    /*//Datacontent.push_back(9.);
    //Datacontent.push_back(8.);
    //Datacontent.push_back(6.);
    //Datacontent.push_back(9.);
    Datacontent.push_back(7.);
    Datacontent.push_back(11.);
    Datacontent.push_back(9.);
    Datacontent.push_back(8.);
    Datacontent.push_back(8.);
    Datacontent.push_back(4.);
    Datacontent.push_back(5.);
    Datacontent.push_back(2.);
    Datacontent.push_back(0.001);
    Datacontent.push_back(0.001);
    Datacontent.push_back(0.001);
    Datacontent.push_back(0.001);
    Datacontent.push_back(0.001);
    */
    // PAPER 2
    /*Datacontent.push_back(12.5);
    Datacontent.push_back(90.);
    Datacontent.push_back(119.);
    Datacontent.push_back(120.);
    Datacontent.push_back(100.);
    Datacontent.push_back(114.);
    Datacontent.push_back(158.);
    Datacontent.push_back(131.);
    Datacontent.push_back(140.);
    Datacontent.push_back(104.);
    Datacontent.push_back(64.);
    Datacontent.push_back(41.);
    Datacontent.push_back(29.5);
    Datacontent.push_back(11.);
    Datacontent.push_back(8.5);
    Datacontent.push_back(2.5);
    Datacontent.push_back(5.);
    */

    //paper 2
    //double Datacontentarray[17] = {12.5397, 89.0563, 118.848, 120.102, 100.032, 144.249, 158.047, 130.764, 139.858, 104.109, 63.6553, 41.3902, 29.4737, 10.6581, 8.46297, 2.50471, 4.69986};
    //unosc paper 2
    double Datacontentarray[17] = {29.1601, 125.747, 206.34, 258.71, 278.78, 278.466, 244.912, 203.204, 165.259, 127.942, 92.1922, 63.9689, 41.3902, 23.2018, 11.9125, 4.69986, 2.19112};
    std::vector<double> Datacontent;
    double Dataintegral = 0;
    for (int i = 0; i < 17; i++){
      Datacontent.push_back(Datacontentarray[i]);
      Dataintegral += Datacontentarray[i];
    }
    
    dataSetPdf.SetBinContents(Datacontent);

    
    double dataintegral = dataSetPdf.Integral();
    
    std::cout<<"\n \n Data Pdf Integral: "<<dataintegral<<" Data content sum: "<<Dataintegral<<"\n "<<std::endl;
    
    AxisCollection AxisCollection = dataSetPdf.GetAxes();
    BinAxis BinAxis = AxisCollection.GetAxis("EPrompt (MeV)");
    std::vector<double> bincentres = BinAxis.GetBinCentres();
    
    std::cout<<"bin centres size: "<<bincentres.size()<<" Datacontentsize: "<<Datacontent.size()<<std::endl;
    for (int i = 0; i < Datacontent.size(); i++){
      std::cout<<"E: "<<bincentres[i]<<"   "<<Datacontent[i]<<std::endl;
    }

    if (Datacontent.size() != bincentres.size()){
      std::cout<<"data content size != bincentres size"<<std::endl;
      return 0;
    }
    /*std::vector<double> normconstraints;
    std::vector<double> normconstraintsuncert;

    //3CAD (only 1 year)
    //normconstraints.push_back(93.);
    //normconstraints.push_back(21.);
    //normconstraints.push_back(20.);
    //normconstraintsuncert.push_back(0.1);
    //normconstraintsuncert.push_back(0.1);
    //normconstraintsuncert.push_back(0.1);
    
    //98% constribution kamland
    normconstraints.push_back(6.46945701333); //shika
    normconstraints.push_back(24.0056564207); //tsuruga
    normconstraints.push_back(25.2114015133); //mihama
    normconstraints.push_back(30.4530774217);// kashiwazaki
    normconstraints.push_back(38.3548046866);// ohi
    normconstraints.push_back(25.7183490858);// takahama
    normconstraints.push_back(10.2448833212);// hamaoka
    normconstraints.push_back(5.33996987795);// tokai
    normconstraints.push_back(5.82295624693);//fukushima
    normconstraints.push_back(2.477021545);//shimane
    normconstraints.push_back(3.16155975112);//onogame
    normconstraints.push_back(2.72220900133);//hanul
    normconstraints.push_back(2.54944943237);//kori
    normconstraints.push_back(2.31381716777);//genkai
    
    normconstraintsuncert.push_back(2.3383735463/normconstraints[0]); //shika
    normconstraintsuncert.push_back(3.96083389604/normconstraints[1]); //tsuruga
    normconstraintsuncert.push_back(4.46522707503/normconstraints[2]); //mihama
    normconstraintsuncert.push_back(5.2986370078/normconstraints[3]);// kashiwazaki
    normconstraintsuncert.push_back(5.47494204929/normconstraints[4]);// ohi
    normconstraintsuncert.push_back(5.83647201958/normconstraints[5]);// takahama
    normconstraintsuncert.push_back(2.82538208432/normconstraints[6]);// hamaoka
    normconstraintsuncert.push_back(2.18901624773/normconstraints[7]);// tokai
    normconstraintsuncert.push_back(2.33036484295/normconstraints[8]);//fukushima
    normconstraintsuncert.push_back(1.72800472684/normconstraints[9]);//shimane
    normconstraintsuncert.push_back(1.91823012974/normconstraints[10]);//onogame
    normconstraintsuncert.push_back(1.35248778871/normconstraints[11]);//hanul
    normconstraintsuncert.push_back(1.47061082114/normconstraints[12]);//kori
    normconstraintsuncert.push_back(1.0901559326/normconstraints[13]);//genkai
    //double antinueffic = 0.85;
    */


    double normconstraintsarray[14] = {6.46945701333,24.0056564207,25.2114015133,30.4530774217,38.3548046866,25.7183490858,10.2448833212,5.33996987795,5.82295624693,2.477021545,3.16155975112,2.72220900133,2.54944943237,2.31381716777};
    std::vector<double> normconstraints;
    for (int i = 0; i < 14; i++)
      normconstraints.push_back(normconstraintsarray[i]);
    
    double normconstraintsuncertarray[14] = {2.3383735463,3.96083389604,4.46522707503,5.2986370078,5.47494204929,5.83647201958,2.82538208432,2.18901624773,2.33036484295,1.72800472684,1.91823012974,1.35248778871,1.47061082114,1.0901559326};
    std::vector<double> normconstraintsuncert;
    for (int i = 0; i < 14; i++)
      normconstraintsuncert.push_back(normconstraintsuncertarray[i]);
    
    for (int i = 0; i < normconstraintsuncert.size(); i++){
      normconstraintsuncert[i] = (normconstraintsuncert[i]/normconstraints[i]);
    }


    double constrainttotal = 0;
    //std::cout<<"data vector length: "<<normconstraints.size()<<std::endl;
    //for (int i = 0; i < normconstraints.size(); i++){
    for (int i = 0; i < normconstraints.size(); i++){
      //normconstraints[i] = (normconstraints[i]*antinueffic);
      constrainttotal += normconstraints[i];
    }

    //for (int i = 0; i < normconstraints.size(); i++){
    for (int i = 0; i < normconstraints.size(); i++){
      std::cout<<"reactor: "<<reactorNames[i]<<" OG constraint: "<<normconstraints[i]<<" new constraint: "<<((normconstraints[i]/constrainttotal)*dataintegral)<<std::endl;
      normconstraints[i] = ((normconstraints[i]/constrainttotal)*dataintegral);
    }
    
    if (normconstraints.size() != numPdfs){
      std::cout<<"number data norm constraints dont match number of Pdfs!!"<<std::endl;
      return 0;
    }
  
    // Generating fake data for 3CAD, norms put in by hand below, taken from number of event expected from the 3CAD reactors using flux=1, 1yr data
      
    /*for (int i = 0; i< numPdfs; i++){
      OscPromptE(UnOscfile, tempFile,FakeDatad21,FakeDatas12,0.0215,distances[i]);   //oscillate each pdf -> outputs a pruned ntuple with an oscillated EPrompt positron spectrum called tempFile, to be used below to fill a reactorPdf
      
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
	std::cout<<reactorNames[i]<<"norm: "<<normconstraints[i]*osc_loss<<std::endl;
	dataSetPdf.Add(*dataPdf,1);
	}
    */
    
    int Normmin = 0;
    int Normmax = 200;

    double LHval;
    int fitattempts = 0;
    badfit = true;
    
    
    while (badfit){
      fitattempts += 1;
      std::cout<<"fit attempt: "<<fitattempts<<std::endl;
      LHval = LHFit(UnOscfile,tempFile,axes,dataSetPdf,numPdfs,reactorNames,distances,FakeDatad21,FakeDatas12,normconstraints,Normmin,Normmax,normconstraintsuncert);  
    }
    


    BinnedED * Result = new BinnedED("Result",axes);
    Result->SetObservables(0);
    BinnedED TotalResult("TotalResult",axes);
    TotalResult.SetObservables(0);
    ROOTMultiPlot* Plot = new ROOTMultiPlot;
    
    TPaveText pt(0.75,0.35,1.0,0.65,"NDC");
    
    double resultintegral = 0;
    for (int i = 0; i< numPdfs; i++){    
      char name[100];
      OscPromptE(UnOscfile,tempFile ,FakeDatad21,FakeDatas12,0.0215,distances[i]);   //oscillate each input reactor KE pdfs -> outputs a pruned ntuple with an oscillated EPrompt positron spectrum called tempFile, to be used below to fill a reactorPdf
    
      BinnedED *reactorPdf = new BinnedED(reactorNames[i],axes);
      reactorPdf->SetObservables(0);
	
      ROOTNtuple reactorNtp(tempFile, "nt");
      for(size_t j = 0; j < reactorNtp.GetNEntries(); j++)
	reactorPdf->Fill(reactorNtp.GetEntry(j));
      //std::cout<<"pre osc Int: "<<unoscPdfint<<std::endl;   //integral of input unoscillated pdfs
      //std::cout<<"post osc Int: "<<oscPdfint<<std::endl;  //integral of pdf after KE->EPrompt + oscillation
      double osc_loss = oscPdfint/(double)unoscPdfint;      // calculate loss in pdf integral
      //std::cout<<" osc loss factor: "<<osc_loss<<std::endl;
      reactorPdf->Normalise();
      //sprintf(name,"ReactorPdf%d_norm",i);
      sprintf(name,"%s_norm",reactorNames[i].c_str());
      
      std::cout<<"Best Fit Norm "<<name<<": "<<bestFit.at(name)<<std::endl;
      resultintegral += bestFit.at(name);
      reactorPdf->Scale(bestFit.at(name));
      //std::cout<<"Norm constraint: "<<normconstraints[i]<<std::endl;
      //resultintegral += normconstraints[i];
      //reactorPdf->Scale(normconstraints[i]);

      TotalResult.Add(*reactorPdf,1);
    
      pt.AddText(Form("norm = %.5f" ,bestFit[name]));
      std::stringstream Name;
      Name<<reactorNames[i]<<"_norm = "<<bestFit[name];
      pt.AddText(Name.str().c_str());
      Result->Empty();
    }
    std::cout<<"\n \n Result Integral: "<<resultintegral<<"\n "<<std::endl;


    pt.AddText(Form("#Delta m_{21} = %.6f",bestFit["d21"]));
    pt.AddText(Form("#theta_{12} = %.3f",bestFit["s12"]));
    pt.SetFillColor(kWhite);
    pt.SetShadowColor(kWhite);
 
    TH1D DataHist;
    TH1D FitHist;
  
    DataHist = DistTools::ToTH1D(dataSetPdf);
    FitHist = DistTools::ToTH1D(TotalResult);
  
    TH1D FullFit("FullFit","",FitHist.GetNbinsX(),FitHist.GetXaxis()->GetXmin(),FitHist.GetXaxis()->GetXmax());
    FullFit.Add(&FitHist);
    DataHist.Sumw2();

    TLegend* leg = new TLegend(0.75,0.8,1.0,1.0);
    leg->AddEntry(&DataHist,"Data","lf");
    leg->AddEntry(&FitHist,"Fit Result","lf");

    TCanvas* c1 = new TCanvas("c1");
    c1->cd();  
    DataHist.SetTitle("Data to Fit");  
    DataHist.GetYaxis()->SetTitle(Form("Counts"));
    DataHist.Draw();
    FitHist.SetLineColor(kRed);
    FitHist.SetLineWidth(3);
    FitHist.Draw("same e");
    leg->Draw();
  
    pt.Draw();
    c1->cd();
    TFile *fileOut = new TFile("~/oxsx/examples/SinglePointLH2D.root","RECREATE");

    FitHist.Write();
    DataHist.Write();
    c1->Write();

    fileOut->Close();
  }
}


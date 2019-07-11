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
#include <ParameterDict.h>
#include <ContainerTools.hpp>
#include <ROOTMultiPlot.h>
#include <TRandom3.h>

// flag which is used for repeating fits that weren't good
bool badfit = true;

double LHFit(const std::string UnOscfile, AxisCollection axes, BinnedED dataSetPdf, int numPdfs, std::vector<std::string> reactorNames, std::vector<double> reactorDistances,int Normmin, int Normmax, std::vector<double> normconstraints, std::vector<double> normconstraintsuncert){
  
    //Fit, Initialisation of lhfunction

    double s12 = atof(argv[4]);
    double d21 = atof(argv[5]);

    const std::string &tempFile = argv[12];

    std::cout<<"s12: "<<s12<<"      d21: "<<d21<<std::endl;
    double LHval;
    int fitattempts = 0;
    badfit = true;
    while (badfit){
      fitattempts += 1;
      std::cout<<"fit attempt: "<<fitattempts<<std::endl;
      LHval = LHFit(UnOscfile,tempFile,axes,dataSetPdf,numPdfs,reactorNames,distances,d21,s12,normconstraints,Normmin,Normmax,normconstraintsuncert);
    }
    
    //TH2D *h_2d = (TH2D*) fin->Get("h_2d");
    int bin_i = atoi(argv[6]);
    int bin_j = atoi(argv[7]);

    std::stringstream sstxt;
    const std::string &textfilename = argv[11];
    //sstxt<<textfilename<<"binx_"<<bin_i<<"biny_"<<bin_j<<".txt"
    ofstream textfile;
    textfile.open (textfilename.c_str(),std::ios_base::app);
    textfile <<"binx_"<<bin_i<<"biny_"<<bin_j<<"LHval_"<<LHval<<"\n";
    textfile.close();
    //h_2d->SetBinContent(bin_i,bin_j,LHval);

    //TH1D DataHist;
    //DataHist = DistTools::ToTH1D(dataSetPdf);
    //DataHist.SetName("DataHist");

    //std::cout<<" saving h2 file: "<<data_constraintsstream.str()<<"_.root"<<std::endl;
    
    //TFile *fileOut = new TFile((data_constraintsstream.str()+"_.root").c_str(),"RECREATE");

    //h_2d->Write();
    //datacontentNT->Write();
    //DataHist.Write();
    //normconstraintsNT->Write();
    fin->Close();
    //fileOut->Close();
    
    //const char *command1 = ("rm " + data_constraintsstream.str() + ".root").c_str();
    //const char *command2 = ("mv " + data_constraintsstream.str() + "_.root " + data_constraintsstream.str() + ".root").c_str();
    //system(command1);
    //system(command2);











  
  //////////////////////////
  // 1. Set Up Likelihood //
  //////////////////////////

  BinnedNLLH lhFunction;

  // A buffer region is used when systematics are used that involve smearing across bins.
  // These the buffer number set here dictates the bins (from the left and right-hand side of the pdf)
  // which won't be used the likelihood calculation.
  lhFunction.SetBufferAsOverflow(true);
  int Buff = 0; // zero buffer region, no systematics used in this example
  lhFunction.SetBuffer(0,Buff,Buff);  

  //////////////////////////////
  // 2. Set Data Dist and pdfs//
  //////////////////////////////

  // Only interested in first bit of data ntuple
  ObsSet dataRep(0);

  //set up axes, binning
  double Emin = 1.;
  double Emax = 8.;
  int numbins = 30;
  AxisCollection axes;
  axes.AddAxis(BinAxis("EPrompt (MeV)", Emin, Emax, numbins));


  /// Data distribution ///
  // #1 //
  const std::string dataFile = "";
  const std::string dataTreeName = "";
  ROOTNtuple dataNt(dataFile, dataTreeName);

  // #2 //
  BinnedED dataSetPdf("dataSetPdf",axes);
  dataSetPdf.SetObservables(dataRep);
  ROOTNtuple dataNtp(dataFile, "nt");
  for(size_t i = 0; i < dataNtp.GetNEntries(); i++)
    dataSetPdf.Fill(dataNtp.GetEntry(i));

  lhFunction.SetDataDist(dataSetPdf); // initialise withe the data set
  

  ////////////////////
  // 1. Set Up PDFs //
  ////////////////////

  // for random intial values of parameters to be fit
  TRandom3 *r1 = new TRandom3();
  r1->SetSeed(0);
  double rand = r1->Rndm();
  
  // objects which contain the min,max,intial parameter values for each parameter in the fit.
  ParameterDict minima;
  ParameterDict maxima;
  ParameterDict initialval;
  ParameterDict initialerr;
 
  ///   REACTOR PDFS   ///
  
  // use prepared info file containing reactor names, distances etc
  // creates vectors for name and distances which are used for oscillating each reactors pdf
  const std::string &UnOscfile = argv[1];
  const std::string &infoFile = argv[3];
  
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

  // looping over list of reactors and creating pdfs for each
  int Normmin = 0;
  int Normmax = 1000;

  char name[100];
  for (int i = 0; i < numPdfs; i++){
    OscPromptE(UnOscfile,tempFile ,d21fix,s12fix,0.,reactorDistances[i]); //or s13 = 0.0215
    double osc_loss = oscPdfint/(double)unoscPdfint;      // calculate loss in pdf integral
    std::cout<<" osc loss factor: "<<osc_loss<<std::endl;
    normconstraints[i] = normconstraints[i]*osc_loss;
    oscconstraintint += normconstraints[i];
  }
  
  //renomalising oscillated constraints, to KL Oscillated spectrum integral
  double DataInt = dataSetPdf.Integral();
  std::cout<<" Data Integral: "<<DataInt<<std::endl;
  double oscrenormInt = 0.;
  for (int i = 0; i < numPdfs; i++){
    normconstraints[i] = ((normconstraints[i]/oscconstraintint)*DataInt); 
    oscrenormInt += normconstraints[i];
    std::cout<<reactorNames[i]<<"   mean: "<<normconstraints[i]<<"  sigma: "<<(normconstraintsuncert[i]*normconstraints[i])<<std::endl;
  }
  
  std::cout<<"Integral of measured (oscillated) KL Data: "<<DataInt<<std::endl;
  std::cout<<"Integral of oscillated renormalised constraints: "<<oscrenormInt<<std::endl;
  
  for (int i = 0; i < numPdfs; i++){
    OscPromptE(UnOscfile,tempFile ,d21fix,s12fix,0.,reactorDistances[i]);   //oscillate each input reactor KE pdfs -> outputs a pruned ntuple with an oscillated EPrompt positron spectrum called tempFile, to be used below to fill a reactorPdf
    
    BinnedED *reactorPdf = new BinnedED(reactorNames[i],axes);
    reactorPdf->SetObservables(0);
    
    ROOTNtuple reactorNtp(tempFile, "nt");
    for(size_t j = 0; j < reactorNtp.GetNEntries(); j++)
      reactorPdf->Fill(reactorNtp.GetEntry(j));
    //std::cout<<"pre osc Int: "<<unoscPdfint<<std::endl;   //integral of input unoscillated pdfs
    //int postoscint = (int)(reactorPdf->Integral());
    //std::cout<<"post osc Int: "<<oscPdfint<<std::endl;  //integral of pdf after KE->EPrompt + oscillation
    double osc_loss = oscPdfint/(double)unoscPdfint;      // calculate loss in pdf integral
    std::cout<<" osc loss factor: "<<osc_loss<<std::endl;
    reactorPdf->Normalise();

    // Setting optimisation limits
    sprintf(name,"%s_norm",reactorNames[i].c_str());
    minima[name] = Normmin;
    maxima[name] = Normmax;
    
    double sigma = (normconstraintsuncert[i]*normconstraints[i]);
    
    initialval[name] = normconstraints[i];//(rand*(maxima[name]-minima[name]))+minima[name];
    initialerr[name] = sigma;//0.5*initialval[name];
    
    lhFunction.AddDist(*reactorPdf);
    //double constraint = normconstraints[i] * osc_loss; //vector of unoscillated norm contraints made in main() Jeff you can probs improve this
    //lhFunction.SetConstraint(name, constraint, (normconstraintsuncert[i]*constraint));  // I just randomly picked 10% for contraint uncertainty 
    lhFunction.SetConstraint(name, normconstraints[i], sigma);
    std::cout<<reactorNames[i]<<"   m: "<<normconstraints[i]<<"    sigma: "<<sigma<<"\n"<<std::endl;
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
  ParameterDict bestFit = fitResult.GetBestFit();

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
  double lhval = (-1.)*lhFunction.Evaluate();//(-1)*lhFunction.Evaluate();
  
  printf("\n-------------------------------------------------------------------------------------\n");
  std::cout<<"LH val at best fit: "<<lhval<<std::endl;
  printf("-------------------------------------------------------------------------------------\n");


  //If Want to plot for testing purposes:  
  BinnedED TotalResult("TotalResult",axes);
  BinnedED ConstrainedTotalResult("ConstrainedTotalResult",axes);
  TotalResult.SetObservables(0);
  ConstrainedTotalResult.SetObservables(0);
  ROOTMultiPlot* Plot = new ROOTMultiPlot;
  
  TPaveText pt(0.75,0.35,1.0,0.65,"NDC");
  for (int i = 0; i< numPdfs; i++){    
    OscPromptE(UnOscfile,tempFile ,d21fix,s12fix,0.0215,reactorDistances[i]);
    
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
    
    sprintf(name,"%s_norm",reactorNames[i].c_str());
    std::cout<<"Best Fit Norm: "<<bestFit.at(name)<<std::endl;
    reactorPdf->Scale(bestFit.at(name));
    
    TotalResult.Add(*reactorPdf,1);
    
    //pt.AddText(Form("norm = %.5f" ,bestFit[name]));
    std::stringstream Name;
    Name<<reactorNames[i]<<"_norm = "<<bestFit[name];
    pt.AddText(Name.str().c_str());

    OscPromptE(UnOscfile,tempFile ,d21fix,s12fix,0.0215,reactorDistances[i]);
    
    BinnedED *ConstrainedreactorPdf = new BinnedED(reactorNames[i],axes);
    ConstrainedreactorPdf->SetObservables(0);
    
    ROOTNtuple ConstrainedreactorNtp(tempFile, "nt");
    for(size_t j = 0; j < ConstrainedreactorNtp.GetNEntries(); j++)
      ConstrainedreactorPdf->Fill(ConstrainedreactorNtp.GetEntry(j));
    ConstrainedreactorPdf->Normalise();

    std::cout<<reactorNames[i]<<"                "<<normconstraints[i]<<std::endl;
    ConstrainedreactorPdf->Scale(normconstraints[i]);

    ConstrainedTotalResult.Add(*ConstrainedreactorPdf,1);
  }
  
  pt.SetFillColor(kWhite);
  pt.SetShadowColor(kWhite);
 
  TH1D DataHist;
  TH1D FitHist;
  TH1D ConstrainedHist;
  
  DataHist = DistTools::ToTH1D(dataSetPdf);
  FitHist = DistTools::ToTH1D(TotalResult);
  ConstrainedHist = DistTools::ToTH1D(ConstrainedTotalResult);

  std::cout<<"Data Int: "<<DataHist.Integral()<<std::endl;
  std::cout<<"Best Fit Pdf Int: "<<FitHist.Integral()<<std::endl;
  std::cout<<"Constrained No Fit Pdf Int: "<<FitHist.Integral()<<std::endl;
  
  TH1D FullFit("FullFit","",FitHist.GetNbinsX(),FitHist.GetXaxis()->GetXmin(),FitHist.GetXaxis()->GetXmax());
  FullFit.Add(&FitHist);
  TH1D ConstrainedFullFit("ConstrainedFullFit","",ConstrainedHist.GetNbinsX(),ConstrainedHist.GetXaxis()->GetXmin(),ConstrainedHist.GetXaxis()->GetXmax());
  ConstrainedFullFit.Add(&ConstrainedHist);
  DataHist.Sumw2();

  TLegend* leg = new TLegend(0.75,0.8,1.0,1.0);
  leg->AddEntry(&DataHist,"Data","lf");
  leg->AddEntry(&FitHist,"Fit Result","lf");
  leg->AddEntry(&ConstrainedHist,"Constrained NoFit Result","lf");

  TCanvas* c1 = new TCanvas("c1");
  c1->cd();  
  DataHist.SetTitle("KL Data(blue) best fit(red)");  
  DataHist.GetYaxis()->SetTitle(Form("Counts"));
  DataHist.Draw("same");
  FitHist.SetLineColor(kRed);
  FitHist.SetLineWidth(3);
  FitHist.Draw("same");//FitHist.Draw("same e");
  ConstrainedHist.SetLineColor(kGreen);
  ConstrainedHist.SetLineWidth(3);
  //ConstrainedHist.Draw("same");
  leg->Draw();
  
  pt.Draw();
  c1->cd();
  
  TFile * fitout = new TFile("TestFit.root","RECREATE");
  c1->Write();
  
  fitout->Close();
    
  return lhval;
}

int main(int argc, char *argv[]) {

  int FakeDataArgc = 14 + 1;
  int NoFakeDataArgc = 12 + 1;

  if (argc < NoFakeDataArgc) {
    std::cout<<"12 (14) arguments expected: \n 1: location of UnOsc pruned Ntuple file (needed for oscillating of prompt energy event spectrum) \n 2: filename of rootfile containting constraints and Data (no .root!) \n 3: location/filename for reactor info file \n \n 4: s12 \n 5: d21 \n 6: s12bin \n 7: d21bin \n \n 8: Emin \n 9: Emax \n 10: numbins \n 11: result text file name/loc \n 12: temporary root file to fill ntuples (can delete afterwards) \n \n (13): Fake data at d21 \n (14):Fake data at s12"<<std::endl;
  }
  else{
    int Normmin = 0;
    int Normmax = 1000;
    
    std::stringstream data_constraintsstream;
    const std::string &data_constraintsFile= argv[2];
    
    if (argc == NoFakeDataArgc)
      data_constraintsstream<<data_constraintsFile;
    else if(argc == FakeDataArgc){
      double FakeDatad21 = atof(argv[13]);
      double FakeDatas12 = atof(argv[14]);
      
      data_constraintsstream<<data_constraintsFile<<"ds21_"<<FakeDatad21<<"s12_"<<FakeDatas12;
    }
    else{
      std::cout<<"Wrong number of arguments!!"<<std::endl;
      return 0;
    }	

    TFile *fin = TFile::Open((data_constraintsstream.str()+".root").c_str(),"READ");
    
    //Data
    TNtuple *datacontentNT = (TNtuple*) fin->Get("databincontents");
    
    float Count;
    datacontentNT->SetBranchAddress("counts", &Count);
    //const std::vector<double> Datacontent;
    std::vector<double> Datacontent;
    
    std::cout<<"\n  Digitised KL Paper binc contents (2nd time): "<<std::endl;
    for(int i = 0; i < datacontentNT->GetEntries(); i++){
      datacontentNT->GetEntry(i);
      std::cout<<"i="<<i<<"  "<<Count<<std::endl;
      Datacontent.push_back(Count);
    }
    
    double Emin = atof(argv[8]);
    double Emax = atof(argv[9]);
    int numbins = atoi(argv[10]);
    ObsSet dataRep(0);
    AxisCollection axes;
    axes.AddAxis(BinAxis("EPrompt (MeV)", Emin, Emax, numbins));
    
    BinnedED dataSetPdf("dataSetPdf",axes);
    dataSetPdf.SetObservables(dataRep);

    dataSetPdf.SetBinContents(Datacontent);

    //norm constraints
    TNtuple* normconstraintsNT = (TNtuple*) fin->Get("normconstraints");
    std::vector<double> normconstraints;
    std::vector<double> normconstraintsuncert;
    float norm;
    float normuncert;
    normconstraintsNT->SetBranchAddress("norm", &norm);
    normconstraintsNT->SetBranchAddress("normuncert", &normuncert);
    for(int i = 0; i < normconstraintsNT->GetEntries(); i++){
      normconstraintsNT->GetEntry(i);
      normconstraints.push_back(norm);
      normconstraintsuncert.push_back(normuncert);
    }
    
    
    //Fit, Initialisation of lhfunction
    const std::string &UnOscfile = argv[1];
    const std::string &infoFile = argv[3];
    
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

    double s12 = atof(argv[4]);
    double d21 = atof(argv[5]);

    const std::string &tempFile = argv[12];

    std::cout<<"s12: "<<s12<<"      d21: "<<d21<<std::endl;
    double LHval;
    int fitattempts = 0;
    badfit = true;
    while (badfit){
      fitattempts += 1;
      std::cout<<"fit attempt: "<<fitattempts<<std::endl;
      LHval = LHFit(UnOscfile,tempFile,axes,dataSetPdf,numPdfs,reactorNames,distances,d21,s12,normconstraints,Normmin,Normmax,normconstraintsuncert);
    }
    
    //TH2D *h_2d = (TH2D*) fin->Get("h_2d");
    int bin_i = atoi(argv[6]);
    int bin_j = atoi(argv[7]);

    std::stringstream sstxt;
    const std::string &textfilename = argv[11];
    //sstxt<<textfilename<<"binx_"<<bin_i<<"biny_"<<bin_j<<".txt"
    ofstream textfile;
    textfile.open (textfilename.c_str(),std::ios_base::app);
    textfile <<"binx_"<<bin_i<<"biny_"<<bin_j<<"LHval_"<<LHval<<"\n";
    textfile.close();
    //h_2d->SetBinContent(bin_i,bin_j,LHval);

    //TH1D DataHist;
    //DataHist = DistTools::ToTH1D(dataSetPdf);
    //DataHist.SetName("DataHist");

    //std::cout<<" saving h2 file: "<<data_constraintsstream.str()<<"_.root"<<std::endl;
    
    //TFile *fileOut = new TFile((data_constraintsstream.str()+"_.root").c_str(),"RECREATE");

    //h_2d->Write();
    //datacontentNT->Write();
    //DataHist.Write();
    //normconstraintsNT->Write();
    fin->Close();
    //fileOut->Close();
    
    //const char *command1 = ("rm " + data_constraintsstream.str() + ".root").c_str();
    //const char *command2 = ("mv " + data_constraintsstream.str() + "_.root " + data_constraintsstream.str() + ".root").c_str();
    //system(command1);
    //system(command2);

  }
}

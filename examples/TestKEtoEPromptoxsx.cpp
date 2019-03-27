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
#include <TCutG.h>

int main(int argc, char *argv[]) {
  TFile * f = new TFile("/data/snoplus/blakei/antinu/test/h2Bruce5yr1000flux.root");
  TH2D * KEvsE = (TH2D*)f->Get("antinuKEvsEPrompt");
  
  int h2Nbins = KEvsE->GetNbinsY();
  double h2min = KEvsE->GetYaxis()->GetXmin();
  double h2max = KEvsE->GetYaxis()->GetXmax();
  
  std::cout<<"KEvsE:     xnbins: "<<h2Nbins<<" xmin: "<<h2min<<" xmaxs: "<<h2max<<std::endl;

  //KE axes and ERecon axes -> need same energy range per bin?
  AxisCollection KEaxes;
  KEaxes.AddAxis(BinAxis("ParKE", 2., 8., 200));
  //KEaxes.AddAxis(BinAxis("EPrompt", 1., 8., 35));
  AxisCollection Dataaxes;
  Dataaxes.AddAxis(BinAxis("E1", 0., 9., 60));
  
  BinnedED *reactorPdfKE = new BinnedED("BruceKE",KEaxes);
  reactorPdfKE->SetObservables(0);
  BinnedED *reactorPdfE = new BinnedED("BruceE",Dataaxes);
  reactorPdfE->SetObservables(0);

  ROOTNtuple reactorKENtp("/data/snoplus/blakei/antinu/test/Bruce5yr1000fluxKEds21_7.4e-05_ss12_0.297_ss13_0.0215_oxsx.root", "nt");
  for(size_t j = 0; j < reactorKENtp.GetNEntries(); j++)
    reactorPdfKE->Fill(reactorKENtp.GetEntry(j));
  reactorPdfKE->Normalise();
  
  BinnedED *reactorPdfEev = new BinnedED("BruceEev",Dataaxes);
  reactorPdfEev->SetObservables(0);
  ROOTNtuple reactorEevNtp("/data/snoplus/blakei/antinu/test/Bruce5yr1000fluxE1ds21_7.4e-05_ss12_0.297_ss13_0.0215_oxsx.root", "nt");
  for(size_t j = 0; j < reactorEevNtp.GetNEntries(); j++)
    reactorPdfEev->Fill(reactorEevNtp.GetEntry(j));
  reactorPdfEev->Normalise();
  
  Histogram histKE = reactorPdfKE->GetHistogram();
  AxisCollection AxisCollectionKE = histKE.GetAxes();
  BinAxis BinAxisKE = AxisCollectionKE.GetAxis("ParKE");
  //BinAxis BinAxisKE = AxisCollectionKE.GetAxis("EPrompt");
  double minKE = BinAxisKE.GetMin();
  double maxKE = BinAxisKE.GetMax();
  size_t NumBinsKE = BinAxisKE.GetNBins();
  
  BinAxis BinAxisE = Dataaxes.GetAxis("E1");
  double minE = BinAxisE.GetMin();
  double maxE = BinAxisE.GetMax();
  size_t NumBinsE = BinAxisE.GetNBins();
  
  std::cout<<"minKE: "<<minKE<<" maxKE: "<<maxKE<<" numbinsKE: "<<NumBinsKE<<std::endl;
  
  std::vector<double> binlowedgesKE = BinAxisKE.GetBinLowEdges();
  std::vector<double> bincentresKE = BinAxisKE.GetBinCentres();
  std::vector<double> binhighedgesKE = BinAxisKE.GetBinHighEdges();
  
  std::vector<double> BincontentsKE = reactorPdfKE->GetBinContents();
  
  std::vector<double> binlowedgesE = BinAxisE.GetBinLowEdges();
  std::vector<double> bincentresE = BinAxisE.GetBinCentres();
  std::vector<double> binhighedgesE = BinAxisE.GetBinHighEdges();
  
  std::cout<<binlowedgesKE.size()<<":"<<bincentresKE.size()<<":"<<binhighedgesKE.size()<<":"<<BincontentsKE.size()<<":"<<NumBinsKE<<std::endl;
  std::cout<<binlowedgesE.size()<<":"<<bincentresE.size()<<":"<<binhighedgesE.size()<<":"<<NumBinsE<<"\n"<<std::endl;
  
  std::cout<<"bin : "<<"low : "<<"cen : "<<"high : "<<"cont : "<<std::endl;
  std::cout<<"---------------------------------------------------------------------"<<std::endl;
  for (int i = 0; i<BincontentsKE.size() ; i++){
    std::cout<<" "<<i<<"  :  "<<binlowedgesKE[i]<<"  :  "<<bincentresKE[i]<<"  :  "<<binhighedgesKE[i]<<"  :  "<<BincontentsKE[i]<<std::endl;
  }
  
  std::cout<<"\n \n minE: "<<minE<<" maxE: "<<maxE<<" numbinsE: "<<NumBinsE<<std::endl;
  
  std::cout<<"blah "<<binlowedgesKE[NumBinsKE-1]<<std::endl;
  //for ranges and numbinsE/KE, num of entries in slicecontentsE
  std::vector<double> BincontentsE (NumBinsE,0.);
  double slicecontentsE [NumBinsE];
  
  //std::cout<<"INT TEST "<<KEvsE->Integral(0,499,121,121)<<std::endl;

  
  for (int i = 1; i< NumBinsKE-1; i++){
  //for (int i = 100; i< 102; i++){
    int ybinmin = KEvsE->GetYaxis()->FindBin(binlowedgesKE[i]);
    int ybinmax = (KEvsE->GetYaxis()->FindBin(binhighedgesKE[i]))-1;
    //int ybinminnext = KEvsE->GetYaxis()->FindBin(binlowedgesKE[i+1]);
    //int ybinmaxb4 = KEvsE->GetYaxis()->FindBin(binhighedgesKE[i-1]);
    //int ybin = KEvsE->GetYaxis()->FindBin(bincentresKE[i]);

    //int ybinmin = KEvsE->GetYaxis()->FindBin(bincentresKE[i]);
    //int ybinmax = KEvsE->GetYaxis()->FindBin(bincentresKE[i+1])-1;
    //int ybinminnext = KEvsE->GetYaxis()->FindBin(binlowedgesKE[i+1]);
    //int ybinmaxb4 = KEvsE->GetYaxis()->FindBin(binhighedgesKE[i-1]);
    
    std::cout<<"minbin: "<<ybinmin<<" maxbin: "<<ybinmax<<std::endl;
    //ybinmaxb4 = ybinmax;
    double sliceintegral = 0.;
    
    for (int j = 0; j < NumBinsE; j++){
      int xbinmin = KEvsE->GetXaxis()->FindBin(binlowedgesE[j]);
      int xbinmax = (KEvsE->GetXaxis()->FindBin(binhighedgesE[j]))-1;
      //int xbinminnext = KEvsE->GetXaxis()->FindBin(binlowedgesE[j+1]);
      //int xbinmaxb4 = KEvsE->GetXaxis()->FindBin(binhighedgesE[j-1]);
      
      //int xbinmin = KEvsE->GetXaxis()->FindBin(bincentresE[j]);
      //int xbinmax = KEvsE->GetXaxis()->FindBin(bincentresE[j+1])-1;
      double weight = ((KEvsE->Integral(xbinmin,xbinmax,ybinmin,ybinmax))/(double)(ybinmax-ybinmin+1))/(double)(xbinmax-xbinmin+1);
      //double weight = KEvsE->Integral(xbinmin,xbinmax,ybinmin,ybinmax);
      //std::cout<<weight<<std::endl;
      //std::cout<<"j: "<<j<<" xbinmin: "<<xbinmin<<" xbinmax: "<<xbinmax<<" weight "<<weight<<std::endl;
      /*if (i < NumBinsKE){
	if (ybinmax == ybinminnext)
	  weight -= (KEvsE->Integral(xbinmin,xbinmax,ybinmax,ybinmax)/2.);
      }
      if (i > 1){
	if (ybinmin == ybinmaxb4 && xbinmin != xbinmaxb4)
	  weight -= (KEvsE->Integral(xbinmin,xbinmax,ybinmin,ybinmin)/2.);
	  }*/
      sliceintegral += weight;
      slicecontentsE[j] = weight;
      //std::cout<<"KE: "<<binlowedgesKE[i]<<"_"<<binhighedgesKE[i]<<" -> E: "<<binlowedgesE[j]<<"_"<<binhighedgesE[j]<<" -> weight: "<<weight<<std::endl;
    }
    
    for (int j = 0; j < NumBinsE; j++){
      double origEbincontent = BincontentsE[j];
      BincontentsE[j] = origEbincontent + ((slicecontentsE[j]/sliceintegral)*BincontentsKE[i]);
    }
    
    //std::cout<<"-----------------------------------------------------------------------"<<std::endl;
    //}
  }
  
  std::cout<<"bin : "<<"low : "<<"cen : "<<"high : "<<"cont : "<<std::endl;
  std::cout<<"---------------------------------------------------------------------"<<std::endl;
  for (int i = 0; i<BincontentsE.size() ; i++){
    std::cout<<" "<<i<<"  :  "<<binlowedgesE[i]<<"  :  "<<bincentresE[i]<<"  :  "<<binhighedgesE[i]<<"  :  "<<BincontentsE[i]<<std::endl;
  }

  reactorPdfE->SetBinContents(BincontentsE);
  
  if (binlowedgesKE.size() != bincentresKE.size() || binhighedgesKE.size() != bincentresKE.size()){
    std::cout<<" sizes don't match"<<std::endl;
    return 0;
  }
 
  TH1D ReactorHistKE;
  ReactorHistKE = DistTools::ToTH1D(*reactorPdfKE);
  TH1D ReactorHistE;
  ReactorHistE = DistTools::ToTH1D(*reactorPdfE);
  TH1D ReactorHistEev;
  ReactorHistEev = DistTools::ToTH1D(*reactorPdfEev);
  
  TCanvas* c1 = new TCanvas("c1");
  c1->Divide(1,2);
  c1->cd(1);
  ReactorHistKE.Draw();
  c1->cd(2);
  ReactorHistEev.SetLineColor(2);
  ReactorHistEev.Draw("same");
  ReactorHistE.SetLineColor(3);
  ReactorHistE.Draw("same");

  TCanvas* c2 = new TCanvas("c2");
  c2->cd();
  ReactorHistKE.Draw();

  TFile * fitout = new TFile("~/oxsx/examples/Test.root","RECREATE");
  c1->Write();
  c2->Write();
  fitout->Close();
}



  //TH1D * slice = new TH1D("slice","slice",KEvsE->GetNbinsX(),KEvsE->GetXaxis()->GetXmin(),KEvsE->GetXaxis()->GetXmax());
  //TH1D * slice = new TH1D("slice","slice",NumBinsE,minE,maxE);
  /*
  int ybinmin = KEvsE->GetYaxis()->FindBin(4.);
  int ybinmax = KEvsE->GetYaxis()->FindBin(5.);
    
  for (int i = 0; i < NumBins2; i++){
    int xbinmin = KEvsE->GetXaxis()->FindBin(binlowedges2[i]);
    int xbinmax = KEvsE->GetXaxis()->FindBin(binhighedges2[i]);
    double weight = KEvsE->Integral(xbinmin,xbinmax,ybinmin,ybinmax);
    slice->SetBinContent(i+1,weight);
    std::cout<<"bin num: "<<i<<" xlow: "<<binlowedges2[i]<<" xhigh: "<<binhighedges2[i]<<" h2xlow: "<<xbinmin<<" h2xhigh: "<<xbinmax<<" weight: "<<weight<<std::endl;
  }
  */

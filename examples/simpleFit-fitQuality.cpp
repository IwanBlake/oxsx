#include <ROOTNtuple.h>
#include <BinAxis.h>
#include <AxisCollection.h>
#include <BinnedED.h>
#include <DataSet.h>
#include <vector>
#include <string>
#include <iostream>
#include <Event.h>
#include <DistFiller.h>
#include <BoolCut.h>
#include <CutCollection.h>
#include <DataSetGenerator.h>
#include <OXSXDataSet.h>
#include <IO.h>
#include <Rand.h>
#include <BinnedNLLH.h>
#include <BinnedEDManager.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <Minuit.h>

ParameterDict maximumlikelihood(BinnedED bgPdf, BinnedED signalPdf, BinnedED dataSetPdf ){
	
    //Comment this if I want to check I get the same data.
    Rand::SetSeed(0);

    BinnedNLLH lh;
    lh.SetDataDist(dataSetPdf);
    lh.AddPdf(bgPdf);
	lh.AddPdf(signalPdf);

	ParameterDict values;
	values["minima"]= 0;
	values["maxima"]= 10000;
	values["steps"]= 2;

	ParameterDict initVals;
	initVals["bgHist_norm"]= 20000;
	initVals["sigHist_norm"]= 40000;

	ParameterDict initErrs;
	initErrs["bgHist_norm"]= 1480;
	initErrs["sigHist_norm"]= 2730;

	//ParameterDict expectedRates;
	//expectedRates["vals"]= 14889;
	//expectedRates["vals"]= 27323;

    Minuit min;
    min.SetInitialValues(initVals);
    min.SetInitialErrors(initErrs);
	//min.SetMinima(values);
    //min.SetMaxima(values);

    min.Optimise(&lh);

    FitResult fitResult = min.GetFitResult();

    ParameterDict bestFit = fitResult.GetBestFit();
    fitResult.Print();


    /////***** Checking what the likelihood looks like *****/////

    lh.SetParameters(bestFit);
    std::cout << "Likelihood value at best fit point" << lh.Evaluate() << std::endl;

    //lh.SetParameters(expectedRates);
    //std::cout << "Likelihood value at expected rate point" << lh.Evaluate() << std::endl;

    return(bestFit);
}

int main()
{
	Rand::SetSeed(0);

	// Load data and make PDF's
	//const Histogram bgHist = IO::LoadHistogram("/data/snoplus/turnere/OXSX_HighStats/PDFs/JustRopes_SixMonthsBi214AvAndBothDustsPdf_ValidFit_2E3_058r077.hdf5");
	//const Histogram sigHist = IO::LoadHistogram("/data/snoplus/turnere/OXSX_HighStats/PDFs/JustRopes_SixMonthsTl208AvAndBothDustsPdf_ValidFit_2E3_058r077.hdf5");
	//BinnedED bgPdf("bgHist",bgHist);
	//BinnedED sigPdf("sigHist",sigHist);
	//const Histogram dataSetHist = IO::LoadHistogram("/data/snoplus/turnere/OXSX_HighStats/FakeDataSets/JustRopes_SixMonthsFakeDataSetPdf_ValidFit_2E3_058r077.hdf5");
	//BinnedED dataSetPdf("dataHist",dataSetHist);
	////
	
	TH1D bgRateHist("bgRateHist", "bgRateHist", 1000, -10000.0, 50000.0);
	TH1D sigRateHist("sigRateHist", "sigRateHist", 1000, -10000.0, 50000.0);

	const std::string bgMCfile    = "/data/snoplus/lidgard/OXSX/ntp/TeLoadedTe130_2n2b.root";
	const std::string sigMCfile   = "/data/snoplus/lidgard/OXSX/ntp/TeLoadedTe130_0n2b-partB.root";
	const std::string bgTreeName  = "output";
	const std::string sigTreeName = "output";

	const std::string dataFile = "/data/snoplus/lidgard/OXSX/ntp/TeLoadedTe130_0n2b-partA.root";
	const std::string dataTreeName = "output";

	// Set up binning
    AxisCollection axes;
    axes.AddAxis(BinAxis("energy", 2, 3, 10, "Energy"));

    // Only interested in first bit of data ntuple
    ObsSet dataRep(0);

    // Set up pdf with these bins in this observable
    BinnedED bgPdf("bgHist",axes);
    bgPdf.SetObservables(dataRep);
    BinnedED  sigPdf("sigHist",axes);
    sigPdf.SetObservables(dataRep);

	ROOTNtuple bgMCNtp(bgMCfile, bgTreeName);
    ROOTNtuple sigMCNtp(sigMCfile, sigTreeName);

    for(size_t i = 0; i < bgMCNtp.GetNEntries(); i++){
        bgPdf.Fill(bgMCNtp.GetEntry(i));
    }
    for(size_t i = 0; i < sigMCNtp.GetNEntries(); i++){
        sigPdf.Fill(sigMCNtp.GetEntry(i));
    }

	bgPdf.Normalise();
    sigPdf.Normalise();
	
	BinnedED  dataSetPdf("dataSetPdf",axes);
    dataSetPdf.SetObservables(dataRep);
	
	ROOTNtuple dataNt(dataFile, dataTreeName);
	for(size_t i = 0; i < dataNt.GetNEntries(); i++){
        dataSetPdf.Fill(dataNt.GetEntry(i));
    }

	int noOfExperiments = 50;
	int noOfBins = dataSetPdf.GetNBins();

	for (int j = 0; j < noOfExperiments; j++){
		for (int i = 0; i < noOfBins; i++){
		  dataSetPdf.SetBinContent(i, Rand::Poisson(dataSetPdf.GetBinContent(i)));
		}

		std::cout << "Total number of events in fake data PDF " << dataSetPdf.Integral() << std::endl;
		ParameterDict bestFit = maximumlikelihood(bgPdf, sigPdf, dataSetPdf);

		std::cout << "Best Fit Values: " << std::endl;
		for(ParameterDict::const_iterator it = bestFit.begin(); it != bestFit.end(); ++it){
			std::cout << "Parameter name: " << it->first << "\tParameter Value: " << it->second << std::endl;

			if (it->first=="bgHist_norm")
				bgRateHist.Fill(it->second);
			if (it->first=="sigHist_norm")
				sigRateHist.Fill(it->second);
		}
	}

  //bgRateHist.SaveAs("/data/snoplus/lidgard/OXSX/JustRopes_Bi214AllExceptRopesRates_ValidFit_2E3_058r077.root");
  //sigRateHist.SaveAs("/data/snoplus/lidgard/OXSX/JustRopes_Tl208AllExceptRopesRates_ValidFit_2E3_058r077.root");
  bgRateHist.SaveAs("/data/snoplus/lidgard/OXSX/0n.root");
  sigRateHist.SaveAs("/data/snoplus/lidgard/OXSX/2n.root");
}

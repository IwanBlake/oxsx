#include <BinnedED.h>
#include <ROOTNtuple.h>
#include <BinnedNLLH.h>
#include <ParameterDict.h>
#include <Rand.h>
#include <TH1D.h>
#include <Minuit.h>

ParameterDict maximumlikelihood(BinnedED bgPdf, BinnedED signalPdf, BinnedED dataSetPdf ){

    Rand::SetSeed(0); //Comment this if I want to check I get the same data.

	// Set Up LH function  
    BinnedNLLH lhFunction;
    lhFunction.SetDataDist(dataSetPdf);
    lhFunction.AddPdf(bgPdf);
	lhFunction.AddPdf(signalPdf);

	// set fit parameter limits
	ParameterDict minima;
	minima["bgPdf_norm"]= 0;
	minima["signalPdf_norm"]= 0;

	ParameterDict maxima;
	maxima["bgPdf_norm"]= 100000;
	maxima["signalPdf_norm"]= 100000;

	ParameterDict initVals;
	initVals["bgPdf_norm"]= 20000;
	initVals["signalPdf_norm"]= 40000;

	ParameterDict initErrs;
	initErrs["bgPdf_norm"]= 1480;
	initErrs["signalPdf_norm"]= 2730;

	ParameterDict expectedRates;
	expectedRates["bgPdf_norm"]= 14889;
	expectedRates["signalPdf_norm"]= 27323;

	// Setup and run the minimiser
    Minuit min;
    min.SetInitialValues(initVals);
    min.SetInitialErrors(initErrs);
	min.SetMinima(minima);
    min.SetMaxima(maxima);

    min.Optimise(&lhFunction);

    FitResult fitResult = min.GetFitResult();
    ParameterDict bestFit = fitResult.GetBestFit();
    fitResult.Print();

    // Check what the likelihood looks like
    lhFunction.SetParameters(bestFit);
    std::cout << "Likelihood value at best fit point: " << lhFunction.Evaluate() << std::endl;
    lhFunction.SetParameters(expectedRates);
    std::cout << "Likelihood value at expected rate point: " << lhFunction.Evaluate() << std::endl;

    return(bestFit);
}

int main(){
	Rand::SetSeed(0);

	// data (ntuples) to load
	const std::string bgMCfile    = "/data/snoplus/lidgard/OXSX/ntp/TeLoadedTe130_2n2b.root";
	const std::string signalMCfile   = "/data/snoplus/lidgard/OXSX/ntp/TeLoadedTe130_0n2b-partB.root";
	const std::string bgTreeName  = "output";
	const std::string signalTreeName = "output";

	const std::string dataFile = "/data/snoplus/lidgard/OXSX/ntp/TeLoadedTe130_0n2b-partA.root";
	const std::string dataTreeName = "output";

	// Set up binning
    AxisCollection axes;
    axes.AddAxis(BinAxis("energy", 2, 3, 10, "Energy"));

    // Only interested in first bit of data ntuple
    ObsSet dataRep(0);

    // Set up pdf with these bins in this observable
    BinnedED bgPdf("bgPdf",axes);
	BinnedED  signalPdf("signalPdf",axes);
	BinnedED  dataSetPdf("dataSetPdf",axes);
    
	bgPdf.SetObservables(dataRep);
    signalPdf.SetObservables(dataRep);
    dataSetPdf.SetObservables(dataRep);

	ROOTNtuple bgMCNtp(bgMCfile, bgTreeName);
    ROOTNtuple signalMCNtp(signalMCfile, signalTreeName);
	ROOTNtuple dataNtp(dataFile, dataTreeName);

	// Fill with data
    for(size_t i = 0; i < bgMCNtp.GetNEntries(); i++)
        bgPdf.Fill(bgMCNtp.GetEntry(i));
    for(size_t i = 0; i < signalMCNtp.GetNEntries(); i++)
        signalPdf.Fill(signalMCNtp.GetEntry(i));
	for(size_t i = 0; i < dataNtp.GetNEntries(); i++)
        dataSetPdf.Fill(dataNtp.GetEntry(i));

	bgPdf.Normalise();
    signalPdf.Normalise();

	// run fitting a number of times to check the fit results are normally distributed
	int noOfExperiments = 100;
	int noOfBins = dataSetPdf.GetNBins();
	TH1D bgRateHist("bgRateHist", "bgRateHist", 1000, -10000.0, 50000.0);
	TH1D signalRateHist("signalRateHist", "signalRateHist", 1000, -10000.0, 50000.0);

	for (int j = 0; j < noOfExperiments; j++){
		for (int i = 0; i < noOfBins; i++)
		  dataSetPdf.SetBinContent(i, Rand::Poisson(dataSetPdf.GetBinContent(i))); // adds some noise

		std::cout << "Total number of events in fake data PDF " << dataSetPdf.Integral() << std::endl;
		ParameterDict bestFit = maximumlikelihood(bgPdf, signalPdf, dataSetPdf); // do the fitting

		std::cout << "Best Fit Values: " << std::endl;
		for(ParameterDict::const_iterator it = bestFit.begin(); it != bestFit.end(); ++it){ // print out results
			std::cout << "Parameter name: " << it->first << "\tParameter Value: " << it->second << std::endl;

			if (it->first=="bgPdf_norm")
				bgRateHist.Fill(it->second); // histogram results
			if (it->first=="signalPdf_norm")
				signalRateHist.Fill(it->second);
		}
	}

	//save results
	bgRateHist.SaveAs("/data/snoplus/lidgard/OXSX/0n.root");
	signalRateHist.SaveAs("/data/snoplus/lidgard/OXSX/2n.root");
}

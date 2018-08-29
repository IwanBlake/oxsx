// A simple fit in energy for signal and a background
#include <BinnedED.h>
#include <ROOTNtuple.h>
#include <BinnedNLLH.h>
#include <BinnedOscNLLH.h>
#include <Minuit.h>
#include <ParameterDict.h>

 // data (ntuples) to load
 const std::string bgMCfile    = "/data/snoplus/blakei/antinu/mc/ntuples/Oscbg.root";
 const std::string signalMCfile   = "/data/snoplus/blakei/antinu/mc/ntuples/Oscsig.root";
 const std::string bgTreeName  = "nt";
 const std::string signalTreeName = "nt";

 const std::string dataFile = "/data/snoplus/blakei/antinu/mc/ntuples/Oscdata.root";
 const std::string dataTreeName = "nt";

int main(){
    ////////////////////
    // 1. Set Up PDFs //
    ////////////////////

    // Set up binning
    AxisCollection axes;
    axes.AddAxis(BinAxis("ParKE", 2, 8, 30, "ParKE"));

    // Only interested in first bit of data ntuple
    ObsSet dataRep(0);

    // Set up pdf with these bins in this observable
    BinnedED bgPdf("bgPdf",axes);
    bgPdf.SetObservables(dataRep);
    BinnedED  signalPdf("signalPdf",axes);
    signalPdf.SetObservables(dataRep);
	BinnedED  dataSetPdf("dataSetPdf",axes);
    dataSetPdf.SetObservables(dataRep);

    std::cout << "Initialised Pdfs" << std::endl;

    ///////////////////////
    // 2. Fill with data //
    ///////////////////////
    ROOTNtuple bgMCNtp(bgMCfile, bgTreeName);
    ROOTNtuple signalMCNtp(signalMCfile, signalTreeName);
	ROOTNtuple dataNtp(dataFile, dataTreeName);

    for(size_t i = 0; i < bgMCNtp.GetNEntries(); i++)
        bgPdf.Fill(bgMCNtp.GetEntry(i));
    for(size_t i = 0; i < signalMCNtp.GetNEntries(); i++)
        signalPdf.Fill(signalMCNtp.GetEntry(i));
	for(size_t i = 0; i < dataNtp.GetNEntries(); i++)
        dataSetPdf.Fill(dataNtp.GetEntry(i));

    bgPdf.Normalise();
    signalPdf.Normalise();

    std::cout << "Filled pdfs " << std::endl;

    /////////////////////////////////////////////
    // 3. Set Up LH function & fit parameters  //
    /////////////////////////////////////////////
    BinnedNLLH lhFunction;
	//BinnedOscNLLH lhFunction;
    lhFunction.SetDataDist(dataSetPdf); // initialise withe the data set
    lhFunction.AddPdf(bgPdf);  //const std::string& name_,
    lhFunction.AddPdf("norm", signalPdf);
	lhFunction.AddPdf("osc", signalPdf);
    std::cout << "Built LH function " << std::endl;

	// set fit parameter limits
	ParameterDict minima;
	minima["bgPdf_norm"]= 0;
	minima["signalPdf_norm"]= 0;
	minima["signalPdf_delmsqr21"]  = 0;
    minima["signalPdf_sinsqrtheta12"]= 0;
    minima["signalPdf_sinsqrtheta13"]= 0;

	ParameterDict maxima;
	maxima["bgPdf_norm"]= 100000;
	maxima["signalPdf_norm"]= 100000;
	maxima["signalPdf_delmsqr21"]  = 1;
    maxima["signalPdf_sinsqrtheta12"]= 1;
    maxima["signalPdf_sinsqrtheta13"]= 1;

	ParameterDict initVals;
	initVals["bgPdf_norm"]= 20000;
	initVals["signalPdf_norm"]= 40000;
	initVals["signalPdf_delmsqr21"]  = 7.4e-5;
    initVals["signalPdf_sinsqrtheta12"]= 0.297;
    initVals["signalPdf_sinsqrtheta13"]= 0.0215;

	ParameterDict initErrs;
	initErrs["bgPdf_norm"]= 1480;
	initErrs["signalPdf_norm"]= 2730;
	initErrs["signalPdf_delmsqr21"]  = 7.4e-6;
    initErrs["signalPdf_sinsqrtheta12"]= 0.0297;
    initErrs["signalPdf_sinsqrtheta13"]= 0.00215;

	////////////
    // 4. Fit //
    ////////////
    Minuit min;
	min.SetMethod("Migrad");
	min.SetMaxCalls(1000000);
    min.SetInitialValues(initVals);
    min.SetInitialErrors(initErrs);
	min.SetMinima(minima);
    min.SetMaxima(maxima);

    min.Optimise(&lhFunction);

    FitResult fitResult = min.GetFitResult();

    ParameterDict bestFit = fitResult.GetBestFit();
    fitResult.Print();

	// save results
    //result.SaveAs("simpleFit_result_%d.txt");

    return 0;
}

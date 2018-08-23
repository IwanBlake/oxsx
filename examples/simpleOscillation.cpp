// A simple fit in energy for signal and a background
#include <BinnedED.h>
#include <ROOTNtuple.h>
#include <BinnedNLLH.h>
#include <Minuit.h>
#include <ParameterDict.h>

// data (ntuples) to load
// const std::string bgMCfile    = "/data/snoplus/lidgard/OXSX/ntp/TeLoadedTe130_0n2b.root";
// const std::string signalMCfile   = "/data/snoplus/lidgard/OXSX/OXSXntp/3CADReactorScintFit10Day.root";
// const std::string bgTreeName  = "output";
// const std::string signalTreeName = "nt";

// const std::string dataFile = "/data/snoplus/lidgard/OXSX/OXSXntp/3CADReactorScintFit10DayOscillated.root";
// const std::string dataTreeName = "nt";

const std::string bgMCfile    = "/data/snoplus/lidgard/OXSX/ntp/TeLoadedTe130_2n2b.root";
const std::string signalMCfile   = "/data/snoplus/lidgard/OXSX/ntp/TeLoadedTe130_0n2b-partB.root";
const std::string bgTreeName  = "output";
const std::string signalTreeName = "output";

const std::string dataFile = "/data/snoplus/lidgard/OXSX/ntp/TeLoadedTe130_0n2b-partA.root";
const std::string dataTreeName = "output";

int main(){
    ////////////////////
    // 1. Set Up PDFs //
    ////////////////////

    // Set up binning
    AxisCollection axes;
    axes.AddAxis(BinAxis("energy", 2, 3, 10, "Energy"));

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
    lhFunction.SetDataDist(dataSetPdf); // initialise withe the data set
    lhFunction.AddPdf(bgPdf);
    lhFunction.AddPdf(signalPdf);
    std::cout << "Built LH function " << std::endl;

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
	
    // apply oscillation to events
    // ParameterDict oscillation_values;
    // oscillation_parameters["delta_msqr_21"]  = 7.4e-5;
    // oscillation_parameters["sinsqr_theta_12"]= 0.297;
    // oscillation_parameters["sinsqr_theta_13"]= 0.0215;
    // lhFunction.SetParameters(oscillation_values);
	

	////////////
    // 4. Fit //
    ////////////
    Minuit min;
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

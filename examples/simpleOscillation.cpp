// A simple fit in energy for signal and a background
#include <BinnedED.h>
#include <ROOTNtuple.h>
#include <BinnedNLLH.h>
#include <GridSearch.h>
#include <ParameterDict.h>

const std::string bgMCfile    = "/data/snoplus/lidgard/OXSX/ntp/TeLoadedTe130_0n2b.root";
const std::string sigMCfile   = "/data/snoplus/lidgard/OXSX/OXSXntp/3CADReactorScintFit10Day.root";
const std::string bgTreeName  = "output";
const std::string sigTreeName = "nt";

const std::string dataFile = "/data/snoplus/lidgard/OXSX/OXSXntp/3CADReactorScintFit10DayOscillated.root";
const std::string dataTreeName = "nt";

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
    BinnedED bgPdf("bgPDF",axes);
    bgPdf.SetObservables(dataRep);
    BinnedED  signalPdf("signalPDF",axes);
    signalPdf.SetObservables(dataRep);

    std::cout << "Initialised Pdfs" << std::endl;

    /////////////////////////////////////
    // 2. Fill with data and normalise //
    /////////////////////////////////////

    ROOTNtuple bgMC(bgMCfile, bgTreeName);
    ROOTNtuple signalMC(sigMCfile, sigTreeName);

    for(size_t i = 0; i < bgMC.GetNEntries(); i++){
        bgPdf.Fill(bgMC.GetEntry(i));
    }

    for(size_t i = 0; i < signalMC.GetNEntries(); i++){
        signalPdf.Fill(signalMC.GetEntry(i));
    }

    bgPdf.Normalise();
    signalPdf.Normalise();

    std::cout << "Filled pdfs " << std::endl;

    ////////////////////////////
    // 3. Set Up LH function  //
    ////////////////////////////
    ROOTNtuple dataNt(dataFile, dataTreeName);
    BinnedOscNLLH lhFunction;
    lhFunction.SetDataSet(&dataNt); // initialise withe the data set
    lhFunction.AddPdf(bgPdf);
    lhFunction.AddPdf(signalPdf);
    std::cout << "Built LH function " << std::endl;


    // apply oscillation to events
    ParameterDict oscillation_values;
    oscillation_values["delta_msqr_21"]  = 7.4e-5;
    oscillation_values["sinsqr_theta_12"]= 0.297;
    oscillation_values["sinsqr_theta_13"]= 0.0215;
    lhFunction.SetParameters(oscillation_values);


    // Set up the optimisation
    GridSearch gSearch;

    ParameterDict values;
    values["minima"]= 0;
    values["maxima"]= 1000;
    values["steps"]= 2;

    gSearch.SetMaxima(values);
    gSearch.SetMinima(values);
    gSearch.SetStepSizes(values);

    ////////////
    // 4. Fit //
    ////////////
    FitResult result = gSearch.Optimise(&lhFunction);

    ParameterDict fit = result.GetBestFit();
    result.Print();
    sprintf(name, "simpleFit_result_%d.txt",parameter_i);
    result.SaveAs(name);

    return 0;
}

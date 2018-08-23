// A simple fit in energy for signal and a background
#include <BinnedED.h>
#include <ROOTNtuple.h>
#include <BinnedNLLH.h>
#include <GridSearch.h>
#include <ParameterDict.h>

const std::string bgMCfile    = "/data/snoplus/lidgard/OXSX/ntp/TeLoadedTe130_2n2b.root";
const std::string sigMCfile   = "/data/snoplus/lidgard/OXSX/ntp/TeLoadedTe130_0n2b-partB.root";
const std::string bgTreeName  = "output";
const std::string sigTreeName = "output";

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
    BinnedNLLH lhFunction;
    lhFunction.SetDataSet(&dataNt); // initialise withe the data set
    lhFunction.AddPdf(bgPdf);
    lhFunction.AddPdf(signalPdf);

    std::cout << "Built LH function " << std::endl;

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
    result.SaveAs("simpleFit_result.txt");
    return 0;
}

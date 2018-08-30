#include <Minuit.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinimize.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/MnPlot.h>
#include <Minuit2/MnSimplex.h>
#include <Minuit2/MnUserParameters.h>
#include <Minuit2/MnUserCovariance.h>
#include <Minuit2/MnContours.h>
#include <FitResult.h>
#include <Exceptions.h>
#include <Formatter.hpp>
#include <math.h> //sqrt

using ContainerTools::GetValues;
using ContainerTools::HasSameKeys;
using ContainerTools::GetKeys;
using ContainerTools::ToString;


using ROOT::Minuit2::MnMigrad;
using ROOT::Minuit2::MnMinimize;
using ROOT::Minuit2::MnSimplex;
using ROOT::Minuit2::MnUserParameters;

// FIXME:: CHeck the length of the vectors as you set them match fNparams

Minuit::~Minuit(){
    delete fMinimiser;
}

void
Minuit::SetMethod(const std::string& method_){
    fMethod = method_;   
}

std::string
Minuit::GetMethod() const{
    return fMethod;
}

void
Minuit::SetInitialErrors(const ParameterDict& errs_){
    fInitialErrors = errs_;
}

void
Minuit::SetInitialValues(const ParameterDict& vals_){
    fInitialValues = vals_;
}

void
Minuit::SetUpperContourEdge(double val_){
    fMinuitFCN.SetUp(val_);
}

double
Minuit::GetUpperContourEdge() const{
    return fMinuitFCN.GetUp();
}

void
Minuit::SetMaxima(const ParameterDict& maxima_) {fMaxima = maxima_;}

ParameterDict 
Minuit::GetMaxima() const {return fMaxima;}

void
Minuit::SetMinima(const ParameterDict& minima_) {fMinima = minima_;}

ParameterDict 
Minuit::GetMinima() const {return fMinima;}

void 
Minuit::Initialise(TestStatistic * testStat_){
    std::cout << "Start of changes." << std::endl;
    delete fMinimiser;
    fMinimiser = NULL;

    // check everything makes sense
    if( !HasSameKeys(fMinima, fMaxima))
        throw LogicError(Formatter()
                         << "Minuit initialisation error "
                         << " minima/maxima parameters don't match:\n"
                         << "Minima for :\n" << ToString(GetKeys(fMinima)) << "\n"
                         << "Maxima for :\n" << ToString(GetKeys(fMaxima)) << "\n"
                         );

    if( !HasSameKeys(fInitialErrors, fInitialValues))
        throw LogicError(Formatter()
                         << "Minuit initialisation error "
                         << "Initial value/error parameters don't match:\n"
                         << "Initial Values for :\n" << ToString(GetKeys(fInitialValues)) << "\n"
                         << "Initial Errors for :\n" << ToString(GetKeys(fInitialErrors)) << "\n"
                         );

    for(ParameterDict::iterator it = fMinima.begin(); it != fMinima.end();
        ++it){
        if(fMinima[it->first] >= fMaxima[it->first])
            throw LogicError(Formatter()
                             << "Minuit initialisation error "
                             << "Invalid fit range for " << it->first
                             << " : " << it->second << " - " << fMaxima[it->first]
                             << "\n"
                             );
    }


    // take a copy of these here to make sure we hit the same order each time - minuit uses vector we use set/map...
    fParameterNames = GetKeys(fInitialErrors);

    // pass this ordering on to the called function so it knows the right order too
    fMinuitFCN = MinuitFCN(testStat_, fParameterNames);

    // max or min?
    if(fMaximising)
        fMinuitFCN.SetSignFlip(true);

	// std::cout << "fParameterNames: " << std::endl;
	// for (int i=0; i<fParameterNames.size(); ++i){
		// std::cout << "param names" << i << ":" << ToString(fParameterNames) << std::endl;
	// }
	// std::cout << "---fInitialValues: " << ToString(GetKeys(fInitialValues)) << std::endl;
	// std::cout << "---fInitialErrors: " << ToString(GetKeys(fInitialErrors)) << std::endl;

    // Create parameters and set limits
    // MnUserParameters params(GetValues(fInitialValues, fParameterNames), GetValues(fInitialErrors, fParameterNames));
    MnUserParameters params;

    // std::vector<double> initValues = GetValues(fInitialValues, fParameterNames);
    // std::vector<double> initErrors = GetValues(fInitialErrors, fParameterNames);
    // std::vector<std::string> names  = GetValues(fParameterNames);
    for (ParameterDict::const_iterator i = fInitialValues.begin(); i != fInitialValues.end(); ++i) {
      std::cout << i->first <<" "<<fInitialValues.at(i->first)<<" "<<fInitialErrors.at(i->first) << std::endl;
      params.Add(i->first.c_str(),fInitialValues.at(i->first),fInitialErrors.at(i->first));
    }
    std::cout << "after names" << std::endl;

    if(fMinima.size() && fMaxima.size()){
        int i = 0;
        for(std::set<std::string>::iterator it = fParameterNames.begin();
            it != fParameterNames.end(); ++it){
            params.SetLimits(i++, fMinima[*it], fMaxima[*it]);
        }
    }
    std::cout << "after limits" << std::endl;

    if("Migrad" == fMethod)
        fMinimiser = new MnMigrad(fMinuitFCN, params);
    
    else if ("Minimize" == fMethod)
        fMinimiser = new MnMinimize(fMinuitFCN, params);
    
    else if ("Simplex" == fMethod)
        fMinimiser = new MnSimplex(fMinuitFCN, params);
    else 
        throw ValueError(Formatter() << "Minuit::Unrecognised method - "
                         << fMethod
                         );
}

void
Minuit::Fix(const std::string& name_){
    fFixedParameters.insert(name_);
}

void
Minuit::Release(const std::string& name_){
    fFixedParameters.erase(name_);
}

void
Minuit::SetMaxCalls(unsigned max_) {
    fMaxCalls = max_;
}

unsigned
Minuit::GetMaxCalls() const {
    return fMaxCalls;
}

const FitResult&
Minuit::Optimise(TestStatistic* testStat_){
    testStat_ -> RegisterFitComponents();

    Initialise(testStat_);

	std::cout << "TestStatistic: " << ToString(testStat_->GetParameterNames()) << std::endl;
    if(testStat_->GetParameterNames() != fParameterNames)
        throw LogicError(Formatter() << "Minuit config parameters don't match the test statistic!\n" 
                         << "TestStatistic: " << ToString(testStat_->GetParameterNames()) 
                         << "\nMinuit: " << ToString(fParameterNames)
                         );

    // fix the requested parameters
    // first work out which minuit index this refers to
    std::set<std::string>::iterator it = fFixedParameters.begin();
    for(; it != fFixedParameters.end(); ++it){
        size_t pos = std::distance(std::find(fParameterNames.begin(), fParameterNames.end(), *it), fParameterNames.begin());
        fMinimiser -> Fix(pos);
    }

    // defaults are same as ROOT defaults
    fnMin = fMinimiser -> operator()(fMaxCalls, fTolerance);

    fFitResult.SetBestFit(ContainerTools::CreateMap(fParameterNames, fMinimiser -> Params()));
    fFitResult.SetErrors(ContainerTools::CreateMap(fParameterNames, fMinimiser -> Errors()));
    fFitResult.SetValid(fnMin.IsValid());

    if(fMaximising)
        fFitResult.SetExtremeVal(-fnMin.Fval());
    else
        fFitResult.SetExtremeVal(fnMin.Fval());

    DenseMatrix covarianceMatrix = CalcCovarianceMatrix(fnMin);

    fFitResult.SetCovarianceMatrix(covarianceMatrix);

    fIsOptimised = true;

    return fFitResult;
}

void 
Minuit::SetTolerance(double tol_) {
  fTolerance = tol_;
}

double
Minuit::GetTolerance() const {
  return fTolerance;
}

FitResult
Minuit::GetFitResult() const{
    return fFitResult;
    
}

DenseMatrix
Minuit::CalcCovarianceMatrix(ROOT::Minuit2::FunctionMinimum functionMin) {
    ROOT::Minuit2::MnUserCovariance minCovariance = functionMin.UserCovariance();
    const std::vector<double>& minCovarianceVector = minCovariance.ROOT::Minuit2::MnUserCovariance::Data();
    int noVectorEntries = minCovarianceVector.size();
    double noRows = (sqrt((8*noVectorEntries) + 1) - 1) / 2;
    DenseMatrix minCovarianceMatrix (noRows, noRows);
    minCovarianceMatrix.SetSymmetricMatrix(minCovarianceVector);
    return minCovarianceMatrix;
}

std::vector<std::pair<double,double> > 
Minuit::GetContour(TestStatistic* testStat_, const std::string& parx,const std::string& pary, const double& upperContour){

    if(!fIsOptimised)
      throw LogicError("Minuit::GetContour : Tried to find contour without optermising function. See this->Optimise.");
    // Set fMinuitFCN upper contour
    fMinuitFCN.SetUp(upperContour);
    Initialise(testStat_);

    // ROOT::Minuit2::FunctionMinimum fnMin = fMinimiser -> operator()(fMaxCalls, fTolerance);
    ROOT::Minuit2::MnPrint printer;
    std::cout <<fnMin<<std::endl;
    ROOT::Minuit2::MnContours contours(fMinuitFCN,fnMin);
    
    std::cout << "here" << std::endl;
    std::vector<std::pair<double,double> > cont = contours(fMinimiser ->Index(parx.c_str()),fMinimiser ->Index(pary.c_str()), 20);
    std::cout << "there" << std::endl;
    // std::vector<std::pair<double,double> > cont = contours(0, 1, 20);
    ROOT::Minuit2::MnPlot plot;
    // cont.insert(cont.end(), cont4().begin(), cont4().end());
    plot(fnMin.UserState().Value(parx.c_str()), fnMin.UserState().Value(pary.c_str()), cont);
    return cont;


}

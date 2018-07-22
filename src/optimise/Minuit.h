/************************************************************************************/
/* Largely inspired by the similar class in RAT, written by P.G.Jones and M.Mottram */
/************************************************************************************/

#ifndef __OXSX_MINUIT__
#define __OXSX_MINUIT__
#include <Optimiser.h>
#include <string>
#include <vector>
#include <MinuitFCN.h>
#include <Minuit2/MnApplication.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnUserTransformation.h>
#include <Minuit2/MinimumState.h>
#include <Minuit2/MinimumSeed.h>
#include <FitResult.h>
#include <DenseMatrix.h>
#include <set>

class TestStatistic;

class Minuit : public Optimiser{
 public:
    Minuit() :  fMethod("Migrad"),
                fMinimiser(NULL), fMaxCalls(0), 
                fTolerance(0.1), fMaximising(false),
                fIsOptimised(false),
fnMin(ROOT::Minuit2::MinimumSeed(ROOT::Minuit2::MinimumState(2),ROOT::Minuit2::MnUserTransformation()),3)
  {
                
                // ROOT::Minuit2::MnUserTransformation mut;
                // ROOT::Minuit2::MinimumState ms(2);
                // ROOT::Minuit2::MinimumSeed mseed(ms,mut);
                // fnMin(ROOT::Minuit2::MinimumSeed(ROOT::Minuit2::MinimumState(2),ROOT::Minuit2::MnUserTransformation()),3);
                }
    ~Minuit();

    virtual const FitResult& Optimise(TestStatistic*);

    void Fix(const std::string& param_);
    void Release(const std::string& param_);

    std::vector< std::pair<double,double> > 
    GetContour(TestStatistic* testStat_,const std::string&,const std::string&, const double&);

    void SetMethod(const std::string&);
    std::string GetMethod() const;

    void SetInitialValues(const ParameterDict&);
    void SetInitialErrors(const ParameterDict&);

    void   SetUpperContourEdge(double);
    double GetUpperContourEdge() const;

    void SetMinima(const ParameterDict& minima_);
    ParameterDict GetMinima() const;

    void SetMaxima(const ParameterDict& maxima_);
    ParameterDict GetMaxima() const;
   
    void     SetMaxCalls(unsigned);
    unsigned GetMaxCalls() const;
    
    void   SetTolerance(double);
    double GetTolerance() const;

    void SetMaximising(bool b_) {fMaximising = b_;}
    bool GetMaximising() const  {return fMaximising;}

    FitResult GetFitResult() const;
    

 private:
    void Initialise(TestStatistic*);
    MinuitFCN   fMinuitFCN; // wrapper on evaluator so migrad can call it
    ParameterDict fInitialValues;
    ParameterDict fInitialErrors;

    ParameterDict fMinima;
    ParameterDict fMaxima;
    std::set<std::string>    fFixedParameters;

    unsigned fMaxCalls;
    double   fTolerance;
    bool fIsOptimised;

    std::string fMethod;
    ROOT::Minuit2::MnApplication* fMinimiser;

    std::set<std::string> fParameterNames; 
    // the order they are held in vectors for ROOT

    FitResult fFitResult;
    bool fMaximising;

    DenseMatrix CalcCovarianceMatrix(ROOT::Minuit2::FunctionMinimum);
    ROOT::Minuit2::FunctionMinimum fnMin;
};
#endif

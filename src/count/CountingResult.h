#ifndef __OXSX_COUNTING_RESULT__
#define __OXSX_COUNTING_RESULT__
#include <CutLog.h>
#include <Histogram.h>
#include <vector>
#include <string>

class CountingResult{
 public:
    CountingResult() : fObservedCounts(-1), fSignalEfficiency(-1) {}
    void AddBackground(double expectedRate_, const std::string& name_,
                       double error_, const CutLog&);
    void SetSignal(double sigEff_, const std::string& name_, const CutLog&);

    void SetDataLog(const CutLog&);
    void SetObservedCounts(int counts_);

    int    GetObservedCounts() const;
    double GetExpectedCounts() const;   
    double GetExpectedCountError() const;
    double GetSignalEfficiency() const;
    const std::vector<double>& GetExpectedRates() const;
    const std::vector<double>& GetExpectedRateErrors() const;

    void Print() const;
    void SaveAs(const std::string& filename_) const;
    std::string AsString() const;

    void SetSampleHist(const Histogram& hist_);
    const Histogram& GetSampleHist() const;

 private:
    int  fObservedCounts;
    double fSignalEfficiency;
    std::vector<CutLog>      fBackgroundLogs;
    std::vector<double>      fExpectedRates;
    std::vector<double>      fExpectedRateErrors;
    std::vector<std::string> fBackgroundNames;
    std::string              fSignalName;
    CutLog fDataLog;
    CutLog fSignalLog;
    Histogram fSampleHist;
};
#endif

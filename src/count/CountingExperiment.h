#ifndef __OXSX_COUNTING_EXPERIMENT__
#define __OXSX_COUNTING_EXPERIMENT__
#include <CutCollection.h>
#include <EventSystematicManager.h>
#include <CountingResult.h>
#include <string>

class DataSet;
class EventSystematic;
class CountingExperiment{
 public:
    void CountData(DataSet*);
    void CountBackgrounds();
    void CountSignal();

    void AddCut(const Cut&, const std::string& name_);
    void AddSystematic(EventSystematic*);
    void AddBackground(DataSet* mcData_, double rate_, const std::string& name_);
    void SetSignal(DataSet* mcData_, const std::string& name_);
    const CountingResult& GetCountingResult() const;

 private:
    CutCollection            fCuts;
    std::vector<DataSet*>    fBackgroundDataSets;
    std::vector<double>      fBackgroundNorms;
    DataSet*                 fTestDataSet;
    DataSet*                 fSignalDataSet;

    EventSystematicManager   fSystematics;    

    std::vector<std::string> fBackgroundNames;
    std::string              fSignalName;
    CountingResult           fResult;
};
#endif

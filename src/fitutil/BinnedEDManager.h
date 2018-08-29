/***************************************************************************************************/
/* A collection of BinnedPDFs with normalisations for calculating the probability of a given obs   */
/* given some rates and some pdfs                                                                  */
/***************************************************************************************************/
#ifndef __OXSX_BINNED_ED_MANAGER__
#define __OXSX_BINNED_ED_MANAGER__
#include <BinnedED.h>
#include <FitComponent.h>
#include <ParameterManager.h>
#include <vector>

class Event;
class SystematicManager;
class BinnedEDShrinker;
class BinnedEDManager : public FitComponent{
 public:
    BinnedEDManager() : fNPdfs(0), fName("norms") {}

	void   AddPdf(const BinnedED&);
    void   AddPdf(const std::string&, const BinnedED&);
    void   AddPdfs(const std::vector<BinnedED>&);

    double Probability(const Event&) const;
    double BinProbability(size_t) const;

    const std::vector<double>& GetNormalisations() const;
    void SetNormalisations(const std::vector<double>& normalisations_);

    void ApplySystematics(const SystematicManager& sysMan_);
    void ApplyShrink(const BinnedEDShrinker&);

    const BinnedED& GetOriginalPdf(size_t index_) const;
    unsigned GetNPdfs() const;
    size_t   GetNDims() const;

    // Make a fittable component - i.e. rescale the binned pdfs inside to fit
    void   SetParameter(const std::string& name_, double value);
    double GetParameter(const std::string& name_) const;

    void   SetParameters(const ParameterDict&);
    ParameterDict GetParameters() const;
    size_t GetParameterCount() const;

    std::set<std::string> GetParameterNames() const;
    void   RenameParameter(const std::string& old_, const std::string& new_);

    std::string GetName() const;
    void SetName(const std::string&);

	void ApplyOscillations();


 private:
    ParameterManager       fParameterManager;
    std::vector<BinnedED>  fOriginalPdfs;
    std::vector<BinnedED>  fWorkingPdfs;
    std::vector<double>    fNormalisations;
	std::vector<double>    fOscillationsP1;
	std::vector<double>    fOscillationsP2;
	std::vector<double>    fOscillationsP3;
    int                    fNPdfs;
    size_t fNDims;

    std::string fName; // component name

    void RegisterParameters();
	void RegisterParameters(const std::string&);

 protected:
	std::vector<std::string> fParameterTypes;

};
#endif

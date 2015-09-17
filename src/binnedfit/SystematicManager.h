/****************************************************************************************************/
/*  A collection of systematics, systematic manager controls the setting of the parameters          */
/* in each of them, parameters are ordered in the same way as the the systematics are added to the  */
/* collection                                                                                       */
/****************************************************************************************************/
#ifndef __SYSTEMATIC_MANAGER__
#define __SYSTEMATIC_MANAGER__
#include <vector>
#include <Systematic.h>

class SystematicManager{
 public:
    SystematicManager() {}
    ~SystematicManager();


    void Add(const Systematic&);
    const std::vector<Systematic>& GetSystematics() const;
    
    void SetParameters(const std::vector<double>& params_);
    const std::vector<double>& GetParameters() const;
    
 private:
    std::vector<size_t>     fParamCounts;
    std::vector<Systematic> fSystematics;
    std::vector<double>     fParams;
};
#endif


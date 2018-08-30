#ifndef __OXSX_PRIOR_MANAGER__
#define __OXSX_PRIOR_MANAGER__
#include <vector>
#include <map>
#include <string>
#include <stddef.h>
#include <Function.h>
#include <Prior.h>
#include <ParameterDict.h>

class PriorManager{
    public:
        PriorManager(){;}
        PriorManager& operator=(const PriorManager& other_);
        PriorManager(const PriorManager&);

        PriorManager(std::vector<Prior>& initList){
            priorList=initList;
        }

        void AddPrior(const Prior& prior);
        std::vector<Prior> GetPriors();
        int GetNPriors();
        double GetProbabilities( const ParameterDict& params_);
        double GetLogProbabilities( const ParameterDict& params_);

    private:
        std::vector<Prior> priorList;
};
#endif


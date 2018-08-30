#include <catch.hpp>
#include <iostream>
#include <PriorManager.h>
#include <Prior.h>
#include <ParameterDict.h>
#include <Function.h>

// Test function 
class Func: public Function
{
    public:
        Func(int num):Function(),  fNDims(num){;}
        double operator() (const std::vector<double>& params) const{
            double ans = params[0];
            for (int i = 1; i < params.size(); ++i) {
              ans -= params[i];
            }
            return ans;
        }
        void SetParameter(const std::string&, double){;}
        double GetParameter(const std::string&) const{;}
        void SetParameters(const ParameterDict&){;}
        ParameterDict GetParameters() const{;}
        size_t GetParameterCount() const{;}
        void RenameParameter(const std::string&, const std::string&){;}
        std::string GetName() const{;}
        void SetName(const std::string&){;}
        Function* Clone() const{
          return static_cast<Function*> (new Func(fNDims));
        }
        int GetNDims() const{;}
        std::set<std::string> GetParameterNames() const{;}


    private:
        int fNDims;
};

TEST_CASE("PRIORMANAGER TEST"){
    ParameterDict dict;
    dict["alpha"]  =10;
    dict["beta"]   =5;
    dict["gamma"]  =1;
    dict["delta"]  =1;
    dict["epsilon"]=1;
    Func* func = new Func(2);
    Func* func_2 = new Func(4);

    Prior* pri = new Prior();
    pri->SetFunction(func);
    pri->SetPrimary("alpha");
    pri->SetDependent("beta");
    Prior* pri_2 = new Prior();
    pri_2->SetFunction(func_2);
    pri_2->SetPrimary("beta");

    std::vector<std::string> deps;
    deps.push_back("gamma");
    deps.push_back("delta");
    deps.push_back("epsilon");
    pri_2->SetDependents(deps);

    std::vector<Prior> priors;
    priors.push_back(*pri);
    priors.push_back(*pri_2);
    PriorManager * man = new PriorManager(priors);
    
    PriorManager * man_2 = new PriorManager();
    man_2->AddPrior(*pri);
    man_2->AddPrior(*pri_2);

    SECTION("Prior :: Getters and Setter"){
        REQUIRE(man->GetNPriors() == 2);
        REQUIRE(man_2->GetNPriors() == 2);
    }
    SECTION("Prior :: Calling Probability"){
        REQUIRE(man->GetProbabilities(dict) == 7);
    }
    SECTION("Prior :: Checking error codes"){

        dict["delta"]=100;

        int flag = 0;
        try { 
          man->GetLogProbabilities(dict);
        }catch(const ParameterError& e_){
          flag=1;
        }

        REQUIRE( flag == 1 );
    }
}

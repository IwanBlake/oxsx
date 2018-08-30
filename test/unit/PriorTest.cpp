#include <catch.hpp>
#include <iostream>
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

TEST_CASE("PIROR TEST"){
    ParameterDict dict;
    dict["alpha"]  =10;
    dict["beta"]   =1;
    dict["gamma"]  =1;
    dict["delta"]  =1;
    dict["epsilon"]=1;
    Func* func = new Func(2);
    Func* func_2 = new Func(4);

    Prior* pir = new Prior();
   
    SECTION("Prior :: Getters and Setter"){
        pir->SetPrimary("alpha");
        REQUIRE(pir->GetPrimary() == "alpha");
        pir->SetDependent("beta");

        REQUIRE(pir->GetDependents().size() == 1);
        REQUIRE(pir->GetDependents()[0] == "beta");
        std::vector<std::string> deps;
        deps.push_back("beta");
        deps.push_back("gamma");
        deps.push_back("delta");
        deps.push_back("epsilon");
        pir->SetDependents(deps);
        REQUIRE(pir->GetDependents()[3] == "epsilon");
    }
    SECTION("Prior :: Calling Probability"){
        pir->SetFunction(func);
        pir->SetPrimary("alpha");
        // pir->SetPrimary("phi");
        pir->SetDependent("beta");

        REQUIRE(pir->Probability(dict) == 9);
        std::vector<std::string> deps;
        deps.push_back("beta");
        deps.push_back("gamma");
        deps.push_back("delta");
        deps.push_back("epsilon");
        // deps.push_back("phi");
        pir->SetDependents(deps);
        REQUIRE(pir->Probability(dict) == 6);
    }
    SECTION("Prior :: Checking error codes"){
        pir->SetFunction(func_2);
        pir->SetPrimary("phi");

        int flag = 0;
        try { 
          pir->Probability(dict);
        }catch(const ParameterError& e_){
          flag=1;
        }
 
        REQUIRE( flag == 1 );

        std::vector<std::string> deps;
        deps.push_back("beta");
        deps.push_back("gamma");
        deps.push_back("delta");
        deps.push_back("rho");
        pir->SetDependents(deps);

        pir->SetPrimary("alpha");
        try { 
          pir->Probability(dict);
        }catch(const ParameterError& e_){
          flag=2;
        }
 
        REQUIRE( flag == 2 );
    }
}

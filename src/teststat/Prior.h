#ifndef __OXSX_PIROR__
#define __OXSX_PIROR__
#include <vector>
#include <string>
#include <stddef.h>
#include <Function.h>
#include <Gaussian.h>
#include <ParameterDict.h>


class Prior{
    public:
        Prior(){;}
        ~Prior();
        Prior(const Prior&);
        Prior& operator=(const Prior&);
        Prior(const std::string& pri_, Function* func):
            fPrimary(pri_), function(func){;}
        Prior(const std::string& pri_,std::vector<std::string>& deps_ , Function* func):
            fPrimary(pri_), fDependence(deps_), function(func){;}
        Prior(const std::string& pri_,std::string& dep_ , Function* func):
            fPrimary(pri_), fDependence( std::vector<std::string>(1,dep_) ), function(func){;}


        void SetPrimary(const std::string&);
        std::string GetPrimary()const{return fPrimary;}
        void SetDependent(const std::string&);
        void SetDependents(std::vector<std::string>&);
        std::vector<std::string> GetDependents(){return fDependence;}
        void SetFunction(Function* func);
        Function* GetFunction();

        double Probability(const ParameterDict&);

    private:
        std::string fPrimary;
        std::vector<std::string> fDependence;
        Function* function;
};
#endif


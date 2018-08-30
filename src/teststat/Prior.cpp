#include <Prior.h>
#include <iostream>
#include <Exceptions.h>
#include <algorithm>
#include <ContainerTools.hpp>

using ContainerTools::ToString;
using ContainerTools::GetKeys;

Prior::Prior(const Prior& other_){
    //What do you need to copy over?
    fPrimary= other_.fPrimary;
    fDependence= other_.fDependence;
    function = dynamic_cast<Function*>(other_.function->Clone());
}

Prior& Prior::operator=(const Prior& other_)
{
    fPrimary= other_.fPrimary;
    fDependence= other_.fDependence;
    function = dynamic_cast<Function*>(other_.function->Clone());
    return *this;
}

Prior::~Prior(){
    delete function;
}

void
Prior::SetPrimary(const std::string& pri_){
  fPrimary = pri_;
}

void
Prior::SetDependent(const std::string& dep_){
  std::vector<std::string> dep;
  dep.push_back(dep_);
  fDependence = dep;
}

void
Prior::SetDependents(std::vector<std::string>& deps_){
  fDependence = deps_;
}

void
Prior::SetFunction(Function* func){
    function=func;
}

Function*
Prior::GetFunction(){
    return function;
}

double
Prior::Probability(const ParameterDict& param_){
    std::vector<double> parameters;
    if( param_.find(fPrimary) == param_.end())
        throw ParameterError(Formatter()<<"Prior:: The primary variable \'"<<fPrimary<<"\' isn't in the testStat. TestStat has: [ "<<ToString(GetKeys(param_))<<"]");
    parameters.push_back(param_.at(fPrimary));

    for (std::vector<std::string>::const_iterator para_ = fDependence.begin(); para_ != fDependence.end(); ++para_) {
        if( param_.find(*para_) == param_.end())
            throw ParameterError(Formatter()<<"Prior:: The dependence  \'"<<*para_<< "\' isn't in the testStat. TestStat has: [ "<<ToString(GetKeys(param_))<<"]");
        parameters.push_back(param_.at(*para_));
    }
    if(!function)
        throw LogicError(Formatter()<<"Prior:: function not set.");

    return function->operator()(parameters);
}

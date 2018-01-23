#include <Prior.h>
#include <math.h>
#include <PriorManager.h>
#include <iostream>
#include <Exceptions.h>

PriorManager::PriorManager(const PriorManager& other_ ){
  priorList =  other_.priorList;
}

PriorManager& PriorManager::operator=(const PriorManager& other_)
{
  priorList =  other_.priorList;
  return *this;
}

int
PriorManager::GetNPriors(){
  return priorList.size();
}
void
PriorManager::AddPrior(const Prior& prior){
  std::string pri = prior.GetPrimary();
  for (std::vector<Prior>::iterator i = priorList.begin(); i != priorList.end(); ++i) {
    if (i->GetPrimary()== pri)
      throw LogicError(Formatter()<<"PriorManager:: Two prior  priors with the same primary variable");
  }
  priorList.push_back(prior);
}

std::vector<Prior>
PriorManager::GetPriors(){
  return priorList;
}

double PriorManager::GetProbabilities( const ParameterDict& params_){
  double result=0;
  for (std::vector<Prior>::iterator iter = priorList.begin(); iter != priorList.end(); ++iter ) {
    result += (*iter).Probability(params_);
  }
  return result;

}


double PriorManager::GetLogProbabilities( const ParameterDict& params_){
  double result=0;
  for (std::vector<Prior>::iterator iter = priorList.begin(); iter != priorList.end(); ++iter ) {
    if((*iter).Probability(params_)<=0)
      throw ParameterError(Formatter()<<"PriorManager:: Attempting to log a negative value.");
    result += log((*iter).Probability(params_));
  }
  return result;
}

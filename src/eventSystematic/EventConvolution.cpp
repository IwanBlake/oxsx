#include <EventConvolution.h>
#include <PDF.h>
#include <JumpPDF.h>
#include <Exceptions.h>
#include <Rand.h>
#include <Event.h>
#include <iostream>

/////////////////////////////
// CONSTRUCTORS/DESTRUCTOR //
/////////////////////////////

EventConvolution::~EventConvolution(){
  delete fDist;
}

EventConvolution::EventConvolution(const EventConvolution& other_){
  if(other_.fDist)
    fDist = other_.fDist -> Clone();
  else
    fDist = NULL;
  fObservables = other_.fObservables;
}

EventConvolution
EventConvolution::operator=(const EventConvolution& other_){ 
  if(other_.fDist)
    fDist = other_.fDist -> Clone();
  else                                 
    fDist = NULL;
  fObservables = other_.fObservables;
  return *this;
}

/////////////
// Get/Set //
/////////////

void
EventConvolution::SetPDF(PDF* f_){
  if(!f_)
    fDist = NULL;

  else{
    if(f_->GetNDims() != 1)
      throw DimensionError("EventConvolution::SetFunction", 1, 
                           f_->GetNDims(), "Only implemented for 1D functions");
    fDist = static_cast<ConditionalPDF*>(new JumpPDF(f_));
  } 
}

void
EventConvolution::SetConditionalPDF(ConditionalPDF* c_){
  if(!c_)
    fDist = NULL;
  else
    fDist  = c_->Clone();
}

/////////////////////////////////////////////////////////////////////////////
// FIT COMPONENT INTERFACE : just forward the call to the underlying funcn //
/////////////////////////////////////////////////////////////////////////////

void
EventConvolution::MakeFittable(){
  if(!fDist)
    throw NULLPointerAccessError("EventConvolution::MakeFittable", 
                                 "Have you set the function?");
  fDist->MakeFittable();
}

std::vector<std::string> 
EventConvolution::GetParameterNames() const{
  if(!fDist)
    throw NULLPointerAccessError("EventConvolution::GetParameterNames", 
                                 "Have you set the sampling function?");
  return fDist->GetParameterNames();
}

std::vector<double> 
EventConvolution::GetParameters() const{
  if(!fDist)
    throw NULLPointerAccessError("EventConvolution::GetParameters", 
                                 "Have you set the sampling function?");
  return fDist->GetParameters();
}

size_t
EventConvolution::GetParameterCount() const{
  if(!fDist)
    throw NULLPointerAccessError("EventConvolution::GetParameterCount", 
                                 "Have you set the sampling function?");
  return fDist->GetParameterCount();
}

void
EventConvolution::SetParameters(const std::vector<double>& params_){
  if(!fDist)
    throw NULLPointerAccessError("EventConvolution::SetParameters", 
                                 "Have you set the sampling function?");
  return fDist->SetParameters(params_);
}

void 
EventConvolution::SetParameterNames(const std::vector<std::string>& names_){
  if(!fDist)
    throw NULLPointerAccessError("EventConvolution::SetParameters", 
                                 "Have you set the sampling function?");
  fDist->SetParameterNames(names_);
}

double 
EventConvolution::GetParameter(const std::string& name_) const{
    if(!fDist)
        throw NULLPointerAccessError("EventConvolution::SetParameters", 
                                 "Have you set the sampling function?");
    return fDist->GetParameter(name_);
}

void   
EventConvolution::SetParameter(const std::string& name_, double val_){
    if(!fDist)
        throw NULLPointerAccessError("EventConvolution::SetParameters", 
                                     "Have you set the sampling function?");
    fDist->SetParameter(name_, val_);
}

// Event Systematic Interface
Event
EventConvolution::operator()(const Event& event_){
  if(!fDist)
    throw NULLPointerAccessError("EventConvolution::operator()", 
                                 "Have you set the sampling function?");

  double newVal = fDist->Sample(event_.ToObsSet(fObservables)).at(0);

  std::vector<double> obs = event_.GetData();
  obs[fObservables.GetIndex(0)] = newVal;
  return Event(obs);
}


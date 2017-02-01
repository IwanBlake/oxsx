#include <AnalyticED.h>
#include <Exceptions.h>
#include <EventDistribution.h>
#include <PDF.h>

AnalyticED::AnalyticED(PDF* f_){
    fFunction = dynamic_cast<PDF*>(f_->Clone());
    fNorm     = 1;
}

AnalyticED::~AnalyticED(){
    delete fFunction;
}

AnalyticED::AnalyticED(const AnalyticED& other_){
    fNorm  = other_.fNorm;
    fObservables = other_.fObservables;
    fFunction = dynamic_cast<PDF*>(other_.fFunction->Clone());
}

EventDistribution*
AnalyticED::Clone() const{
    return static_cast<EventDistribution*>(new AnalyticED(*this));
}

double
AnalyticED::Probability(const std::vector<double>& vals_) const{    
    try{
        return fFunction->operator()(vals_)/fNorm;
    }
    catch(const DimensionError& e_){
        throw DimensionError(std::string("AnalyticED internal function ::") + e_.what());
                             
                             
    }
}

double
AnalyticED::Probability(const Event& event_) const{
    try{
        return Probability(event_.ToObsSet(fObservables));
    }

    catch(const RepresentationError& e_){
        throw RepresentationError("AnalyticED::Probability() failed with  "
                                  + std::string(e_.what()) + 
                                  " is the rep set correctly?");
    }
}

double 
AnalyticED::Integral() const{
    return fNorm;
}

void
AnalyticED::Normalise(){
    fNorm = 1;
}

void
AnalyticED::SetObservables(const ObsSet& rep_){
    fObservables = rep_;
}

ObsSet
AnalyticED::GetObservables() const {
    return fObservables;
}

unsigned
AnalyticED::GetNDims() const{
    return fFunction->GetNDims();
}

// Fitting this pdf to data means adjusting the underlying function
void
AnalyticED::MakeFittable(){
    fFunction -> MakeFittable();
}

std::vector<std::string> 
AnalyticED::GetParameterNames() const{
    std::vector<std::string> funcNames = fFunction->GetParameterNames();
    for(size_t i = 0; i < funcNames.size(); i++)
        funcNames[i] = "Analytic Dist: " + funcNames[i];
    return funcNames;
}

std::vector<double>
AnalyticED::GetParameters() const{
    return fFunction->GetParameters();
}

size_t 
AnalyticED::GetParameterCount() const{
    return fFunction->GetParameterCount();
}

void
AnalyticED::SetParameters(const std::vector<double>& params_){
    try{
        fFunction->SetParameters(params_);
    }
    catch(const ParameterCountError& e_){
        throw ParameterCountError(std::string("AnalyticED internal function : ") + e_.what());

    }
}

void
AnalyticED::SetParameterNames(const std::vector<std::string>& names_){
    try{
        fFunction->SetParameterNames(names_);
    }
    catch(const DimensionError& e_){
        throw DimensionError(std::string("AnalyticED internal function : ") + e_.what());
    }
}

double 
AnalyticED::GetParameter(const std::string& name_) const{
    try{
        return fFunction -> GetParameter(name_);
    }
    catch(const NotFoundError& e_){
        throw NotFoundError(std::string("AnalyticED internal function : ") + e_.what());
    }    
}

void   
AnalyticED::SetParameter(const std::string& name_, double val_){
    try{
        fFunction->SetParameter(name_, val_);
    }
    catch(const NotFoundError& e_){
        throw NotFoundError(std::string("AnalyticED internal function : ") + e_.what());
    }
    catch(const ValueError& e_){
        throw ValueError(std::string("AnalyticED internal function : ") + e_.what());
    }
}

#include <ContainerParameter.h>
#include <Exceptions.h>
#include <sstream>

template<typename Container> 
void  
ParameterManager::AddContainer(Container& cntr_, 
                               const std::string& sharedName_){
    std::stringstream ss;
    for(size_t i = 0; i < cntr_.size(); i++){
        ss << sharedName_ << " " << i;
        Add(new ContainerParameter<Container>(ss.str(), cntr_, i));
        ss.str("");
    }
}

template<typename Container> 
void  
ParameterManager::AddContainer(Container& cntr_, 
                               const std::vector<std::string>& names_){
    if(cntr_.size() != names_.size())
        throw DimensionError("ParameterManager::AddContainer", cntr_.size(),
                             names_.size(), "#names doesn't match #parameters");
    for(size_t i = 0; i < cntr_.size(); i++)
        Add(new ContainerParameter<Container>(names_.at(i), cntr_, i));
}

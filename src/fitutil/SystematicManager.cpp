#include <SystematicManager.h>
#include <Exceptions.h>
#include <Formatter.hpp>
#include <ContainerTools.hpp>


void
SystematicManager::Construct(){
    // Don't do anything if there are no systematics
    if(!fGroups.size())
        return;

    //Construct the response matrices.
    for (std::map<std::string,std::vector<Systematic*> >::const_iterator group=fGroups.begin(); group!= fGroups.end(); ++group ) {
        for(size_t i = 0; i < group->second.size(); i++)
            fGroups[group->first].at(i) -> Construct();
    }

    //This loop should construct all fGroups.
    for (std::map<std::string,std::vector<Systematic*> >::const_iterator group =fGroups.begin(); group != fGroups.end(); ++group ) {
        SparseMatrix resp = fGroups[group->first].at(0) -> GetResponse();
        for(size_t i = 1; i < group->second.size(); i++)
            resp *= fGroups[group->first].at(i) -> GetResponse();
        fTotalReponses[group->first]=resp;
    }
}

const SparseMatrix&
SystematicManager::GetTotalResponse(const std::string& groupName_) const{
    //should you construct?
    try{
        return fTotalReponses.at(groupName_);
    }
    catch(const std::out_of_range& e_){
        throw NotFoundError(Formatter()<<
                "SystematicManager :: name : "<< groupName_<<
                " not found in systematic map");
    }
}


void
SystematicManager::Add(Systematic* sys_, const std::string& groupName_){
    fGroups[groupName_].push_back(sys_);
    fNGroups = fGroups.size();
}

void
SystematicManager::UniqueSystematics(const std::vector<std::string>& syss_){
    std::vector<std::string> allname;
    for (int i = 0; i <syss_.size(); ++i) {
        if(fGroups.find( syss_.at(i) ) == fGroups.end())
            throw NotFoundError (
                    Formatter()<<"SystematicManager:: Systematic group "<<
                    syss_.at(i)<<
                    " not found in complete set of groups. " );

        std::vector<Systematic*> group = fGroups[syss_.at(i)];

        for (int j = 0; j < group.size(); ++j) {
            if(std::find(allname.begin(),allname.end(), group.at(j)->GetName() ) == allname.end())
                allname.push_back(group.at(j)->GetName());
            else
                throw NotFoundError(Formatter()<<
                        "SystematicManager :: Systematic: "<< 
                        group.at(i)->GetName()<<
                        " is in more than one group. BAD!!!");
        } 
    }
}

void
SystematicManager::AddDist(const BinnedED& pdf, const std::vector<std::string>& syss_){
    //Check whether all systematics in syss_ are unique.
    UniqueSystematics(syss_);
    fEDGroups[pdf.GetName()] = syss_;
}

void
SystematicManager::DistortEDs(std::vector<BinnedED>& fWorkingEDs_) const {
    for(size_t j = 0; j < fWorkingEDs_.size(); j++){
        const std::string name = fWorkingEDs_.at(j).GetName();

        //If default group exist then apply that set of systematics first regardless of name.
        if ( fGroups.find("default") != fGroups.end() )
            fWorkingEDs_[j].SetBinContents(GetTotalResponse("default").operator()(fWorkingEDs_.at(j).GetBinContents()));

        if(fEDGroups.find(name) == fEDGroups.end())
            continue;
            
        //Apply everything else.
        for (int i = 0; i < fEDGroups.at(name).size(); ++i) {
            //Check that the group has systematics init. 
            std::string groupName = fEDGroups.at(name).at(i);
            if (groupName == "default")
                continue;
            if (!fGroups.at( groupName ).size())
                throw LogicError(Formatter()<<"SystematicManager:: ED "<<
                        name
                        <<" has a systematic group of zero size acting on it");

            fWorkingEDs_[j].SetBinContents(GetTotalResponse(groupName).operator()(fWorkingEDs_.at(j).GetBinContents()));
        }
    }
}

//Getters and Setters

const size_t
SystematicManager::GetNSystematicsInGroup(const std::string& name_) const{
    try{        
        return fGroups.find(name_)->second.size();
    }
    catch(const std::out_of_range& e_){
        throw NotFoundError(Formatter()<<
                "SystematicManager :: name : "<< 
                name_ <<
                " not found in systematic map");
    }
}

const std::vector<std::string>
SystematicManager::GetGroupNames() const{
    std::vector<std::string> v;
    for(std::map<std::string,std::vector<Systematic*> >::const_iterator it = fGroups.begin(); it !=fGroups.end(); ++it) {
        v.push_back(it->first);
    }
    return v;
}

const std::vector<std::string>
SystematicManager::GetSystematicsNamesInGroup(const std::string& name) const{
    try{
        std::vector<Systematic*> sys = fGroups.at(name);
        std::vector<std::string> names;
        for (int i = 0; i < sys.size(); ++i) {
            names.push_back(sys.at(i)->GetName());    
        }
        return names;
    }
    catch(const std::out_of_range& e_){
        throw NotFoundError(Formatter()<<
                "SystematicManager :: name : "<< 
                name <<
                " not found in systematic map");
    }

}

const std::vector<Systematic*>&
SystematicManager::GetSystematicsInGroup(const std::string& name) const{
    try{
        return fGroups.at(name);
    }catch(const std::out_of_range& e_){
        throw NotFoundError(Formatter()<<
                "SystematicManager :: name : "<< 
                name <<
                " not found in systematic map");
    }
}

const std::vector<std::string>
SystematicManager::GetGroups(const std::string& name) const{
    std::vector<std::string> names;
    for (std::map<std::string,std::vector<Systematic*> >::const_iterator group =fGroups.begin(); group != fGroups.end(); ++group )
        names.push_back(group->first);

    return names;
}

const std::map<std::string, std::vector<Systematic*> >& 
SystematicManager::GetSystematicsGroup() const{
    return fGroups;
}

const size_t
SystematicManager::CountNSystematics() const{
    size_t NSystematics = 0;
    for (std::map<std::string,std::vector<Systematic*> >::const_iterator group =fGroups.begin(); group != fGroups.end(); ++group )
        NSystematics += group->second.size();
    return NSystematics;
}

const size_t
SystematicManager::GetNSystematics() const{
    return CountNSystematics();
}

const size_t
SystematicManager::GetNGroups() const{
    return fNGroups;
}

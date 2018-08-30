#include <FileManager.h>
#include <CutCollection.h>
#include <ROOTNtuple.h>
#include <Exceptions.h>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>     // std::distance
#include <glob.h>
#include <iostream>


FileManager::FileManager(const std::string& name_, const std::string& folder_, const std::string& treeName_){
    names.push_back(name_);
    folders.push_back(folder_);
    treeNames.push_back(treeName_);
}

FileManager::FileManager(const std::vector<std::string>& names_, const std::vector<std::string> & folders_, const std::vector<std::string>& treeNames_){
    //Check that vectors are all the same size. 
    if(names_.size()==folders_.size() && names_.size()==treeNames_.size()){
        names=names_;
        folders=folders_;
        treeNames=treeNames_;
    }else{
        throw RepresentationError(std::string("The vectors given to the FileManager are not the same size."));
    }
}


void
FileManager::AddFolder(const std::string& name_,const std::string& path_,const std::string& treeName_){
    names.push_back(name_);
    folders.push_back(path_);
    treeNames.push_back(treeName_);
}

// taken from http://stackoverflow.com/questions/8401777/simple-glob-in-c-on-unix-system.
inline std::vector<std::string> glob(const std::string& pat){
    using namespace std;
    glob_t glob_result;
    glob(pat.c_str(),GLOB_TILDE,NULL,&glob_result);
    vector<string> ret;
    for(unsigned int i=0;i<glob_result.gl_pathc;++i){
        ret.push_back(string(glob_result.gl_pathv[i]));
    }
    globfree(&glob_result);
    return ret;
}

void
FileManager::FillEDs(std::vector<BinnedED>& fPdfs_){

    for (size_t i = 0; i < fPdfs_.size(); ++i) {

        std::string name = fPdfs_.at(i).GetName();        
        std::cout << "filling pdf = "<< name << std::endl;
        try {
            std::vector<std::string>::iterator pos = std::find(names.begin(), names.end(), name);
            size_t kill = std::distance(names.begin(), pos);
            folders.at(kill);

            std::vector<std::string>::iterator iter = names.begin();

            while ((iter = std::find(iter, names.end(), name)) != names.end())
            {
                size_t index = std::distance(names.begin(), iter);
                // s1.find(s2) != std::string::npos
                std::vector<std::string> files;
                if(folders.at(index).find(std::string(".root")) != std::string::npos){
                    files.push_back(folders.at(index));
                }else{
                    files = glob(folders.at(index)+"*.root");
                }

                for (int file = 0; file < files.size(); ++file) {
                    ROOTNtuple temp(files.at(file),treeNames.at(i));
                    for (int j = 0; j < temp.GetNEntries(); ++j) {
                        if(cuts.PassesCuts(temp.GetEntry(j))){
                            fPdfs_.at(i).Fill(temp.GetEntry(j));
                        }
                    }//Entries loop
                }//File loop
                iter++;
            }
            }catch (const std::out_of_range& oor) {
                std::cerr << "Out of Range error: " << oor.what() << '\n';
                std::cerr << "ED name = \""<< name<< "\" not associated with a path in the FileManager "  << '\n';
            }


    }//loop over pdfs.
}

void
FileManager::AddCuts(const CutCollection& cuts_){
    cuts=cuts_;
}

CutCollection
FileManager::GetCuts() const{
    return cuts;
}


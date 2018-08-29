#include <BinnedEDManager.h>
#include <SystematicManager.h>
#include <BinnedEDShrinker.h>
#include <BinnedED.h>
#include <Exceptions.h>
#include <sstream>

unsigned
BinnedEDManager::GetNPdfs() const{
    return fOriginalPdfs.size();
}

size_t
BinnedEDManager::GetNDims() const{
    return fNDims;
}

double
BinnedEDManager::Probability(const Event& data_) const{
    double sum = 0;

    for(size_t i = 0; i < fWorkingPdfs.size(); i++){
        sum += fParameters.at(i) * fWorkingPdfs[i].Probability(data_);
    }

    return sum;
}

double
BinnedEDManager::BinProbability(size_t bin_) const{
    double sum = 0;
    try{
        for(size_t i = 0; i < fWorkingPdfs.size(); i++){
            sum += fParameters.at(i) * fWorkingPdfs.at(i).GetBinContent(bin_);

        }
    }
    catch(const std::out_of_range&){
        throw LogicError("BinnedEDManager:: Normalisation vector doesn't match pdf vector - are the normalisations set?");
    }
    return sum;
}

void
BinnedEDManager::SetNormalisations(const std::vector<double>& normalisations_){
    if (normalisations_.size() != fOriginalPdfs.size())
        throw LogicError("BinnedEDManager: number of norms doesn't match #pdfs");
    fParameters = normalisations_;
}

void
BinnedEDManager::ApplyOscillations(){
    // If there are no oscillation parameters defined don't do anything
    //  ( then, working pdfs = original pdfs from initialisation)

    // get parameters


    // do oscillations here... modify this code
    // for(size_t j = 0; j < fOriginalPdfs.size(); j++){
        // fWorkingPdfs[j] = fOriginalPdfs.at(j);
        // fWorkingPdfs[j].SetBinContents(sysMan_.GetTotalResponse().operator()(fOriginalPdfs.at(j).GetBinContents()));
    // }
}

void
BinnedEDManager::ApplySystematics(const SystematicManager& sysMan_){
    // If there are no systematics dont do anything
    //  ( working pdfs = original pdfs from initialisation)

    if(!sysMan_.GetSystematics().size())
        return;

    for(size_t j = 0; j < fOriginalPdfs.size(); j++){
        fWorkingPdfs[j] = fOriginalPdfs.at(j);
        fWorkingPdfs[j].SetBinContents(sysMan_.GetTotalResponse().operator()(fOriginalPdfs.at(j).GetBinContents()));
    }
}

const BinnedED&
BinnedEDManager::GetOriginalPdf(size_t index_) const{
    return fOriginalPdfs.at(index_);
}

// void
// BinnedEDManager::AddPdf(const BinnedED& pdf_){
	// AddPdf("norm", pdf_);
// }

void
BinnedEDManager::AddPdf(const std::string& type_, const BinnedED& pdf_){
	if (type_=="norm"){
		fParameters.push_back(0);
		fParameterTypes.push_back(type_);
		fOriginalPdfs.push_back(pdf_);
		fWorkingPdfs.push_back(pdf_);
		fNPdfs++;
		RegisterParameters(type_);
	}
	if (type_=="osc"){
		for (int i=0; i<3; i++){ //push back all vectors for the number of neutrino oscillation parameters
			fParameters.push_back(0);
			fParameterTypes.push_back(type_);
			fOriginalPdfs.push_back(pdf_);
			fWorkingPdfs.push_back(pdf_);
			fNPdfs++;
		}
		RegisterParameters(type_);
	}
}

// void
// BinnedEDManager::AddPdfs(const std::vector<BinnedED>& pdfs_){
    // for(size_t i = 0; i < pdfs_.size(); i++){
        // AddPdf(pdfs_.at(i));
    // }
    // RegisterParameters();
// }

const std::vector<double>&
BinnedEDManager::GetNormalisations() const{
    return fParameters;
}

void
BinnedEDManager::ApplyShrink(const BinnedEDShrinker& shrinker_){
    if (!shrinker_.GetBuffers().size())
        return;

    // only shrink if not already shrunk! FIXME: more obvious behaviour
    if (!fWorkingPdfs.size() || fWorkingPdfs.at(0).GetNBins() != fOriginalPdfs.at(0).GetNBins())
        return;

    for (size_t i = 0; i < fWorkingPdfs.size(); i++){
        fWorkingPdfs[i] = shrinker_.ShrinkDist(fWorkingPdfs.at(i));
        fWorkingPdfs[i].Normalise();
    }

}

////////////////////////////////
// Make this a fit component! //
////////////////////////////////

std::string
BinnedEDManager::GetName() const{
    return fName;
}

void
BinnedEDManager::SetName(const std::string& n_){
    fName = n_;
}

void
BinnedEDManager::RenameParameter(const std::string& old_, const std::string& new_){
    fParameterManager.RenameParameter(old_, new_);
}

void
BinnedEDManager::SetParameter(const std::string& name_, double value_){
    fParameterManager.SetParameter(name_, value_);
}

double
BinnedEDManager::GetParameter(const std::string& name_) const{
    return fParameterManager.GetParameter(name_);
}

void
BinnedEDManager::SetParameters(const ParameterDict& ps_){
    fParameterManager.SetParameters(ps_);
}

ParameterDict
BinnedEDManager::GetParameters() const{
    return fParameterManager.GetParameters();
}

size_t
BinnedEDManager::GetParameterCount() const{
    return fParameterManager.GetParameterCount();
}

std::set<std::string>
BinnedEDManager::GetParameterNames() const{
    return fParameterManager.GetParameterNames();
}

// void
// BinnedEDManager::RegisterParameters(){
	// RegisterParameters("norm");
// }

void
BinnedEDManager::RegisterParameters(const std::string& type_){
    fParameterManager.Clear();
	std::vector<std::string> parameterNames;
    for(size_t i = 0; i < fOriginalPdfs.size(); i++){
		if (fParameterTypes.at(i)=="norm"){
			parameterNames.push_back(fOriginalPdfs.at(i).GetName() + "_norm");
			std::cout << "here norm i: " << i << " fParameterTypes.at(i): " << fParameterTypes.at(i) << std::endl;
		}
		if (fParameterTypes.at(i)=="osc"){
			std::cout << "here osc i: " << i << " fParameterTypes.at(i): " << fParameterTypes.at(i) << std::endl;
			parameterNames.push_back(fOriginalPdfs.at(i).GetName() +"_delmsqr21");
			parameterNames.push_back(fOriginalPdfs.at(i).GetName() +"_sinsqrtheta12");
			parameterNames.push_back(fOriginalPdfs.at(i).GetName() +"_sinsqrtheta13");
			i=i+2;
		}
	}
	fParameterManager.AddContainer(fParameters, parameterNames);
}

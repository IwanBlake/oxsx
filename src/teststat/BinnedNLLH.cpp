#include <BinnedNLLH.h>
#include <math.h>
#include <DataSet.h>
#include <Exceptions.h>
#include <DistFiller.h>
#include <CutLog.h>
#include <Exceptions.h>
#include <Formatter.hpp>
#include <iostream>
#include <vector>

int stepcount = 0;

double 
BinnedNLLH::Evaluate(){
    if (stepcount % 10 == 0)
      std::cout<<stepcount<<std::endl; 
    stepcount += 1;
  
    if(!fDataSet && !fCalculatedDataDist) 
        throw LogicError("BinnedNNLH function called with no data set and no DataDist! set one of these first");
    
    if (!fCalculatedDataDist)
        BinData();
    
    if(!fAlreadyShrunk){
        fDataDist = fPdfShrinker.ShrinkDist(fDataDist);
        fAlreadyShrunk = true;
    }
    /*
    std::cout<<"\n 2) AlreadyShrunk"<<std::endl;
    for (size_t i = 0; i < fPdfManager.GetNormalisations().size(); i++){
      std::cout<<"Norm: "<<fPdfManager.GetNormalisations()[i]<<std::endl;;
    }
    for (size_t i = 0; i < 1; i++){
      std::cout<<"Integral: "<<fPdfManager.GetWorkingPdf(0).Integral()<<std::endl;;
    }
    */
    // Construct systematics 
    fSystematicManager.Construct(); 
    /*
    std::cout<<"\n 3) Constructed Sys"<<std::endl;
    for (size_t i = 0; i < fPdfManager.GetNormalisations().size(); i++){
      std::cout<<"Norm: "<<fPdfManager.GetNormalisations()[i]<<std::endl;;
    }
    for (size_t i = 0; i < 1; i++){
      std::cout<<"Integral: "<<fPdfManager.GetWorkingPdf(0).Integral()<<std::endl;;
    }
    */
    // Apply systematics
    fPdfManager.ApplySystematics(fSystematicManager);
    /*
    std::cout<<"\n 4) Applied Sys "<<std::endl;
    for (size_t i = 0; i < fPdfManager.GetNormalisations().size(); i++){
      std::cout<<"Norm: "<<fPdfManager.GetNormalisations()[i]<<std::endl;;
    }
    for (size_t i = 0; i < 1; i++){
      std::cout<<"Integral: "<<fPdfManager.GetWorkingPdf(0).Integral()<<std::endl;;
    }
    */
    
    //check num working pdfs == getNPdfs?
    //swap osc_loss vector with map to ensure getting correct loss for each pdf?

    //int numpdfs = fPdfManager.GetNPdfs();
    //double * osc_loss [numpdfs];
    std::vector<double> osc_loss;
    //std::vector<std::string> osc_loss_pdfnames;
    for (size_t i = 0 ; i < fPdfManager.GetNPdfs(); i++){
      bool foundoscgroup = false;
      std::string pdfname = fPdfManager.GetWorkingPdf(i).GetName();
      //std::vector<std::string> groupsined = fSystematicManager.GetGroupsfromED(pdfname);
      //for (size_t j = 0 ; j < groupsined.size(); j++){
      std::vector<std::string>::iterator it = std::find(fOscPdfs.begin(),fOscPdfs.end(),pdfname);
      if (it != fOscPdfs.end()){
	  //std::cout<<"ED: "<<pdfname<<" sys group: "<<groupsined[j]<<" CHANGE NORM!!"<<std::endl;
	  //std::cout<<"ED: "<<pdfname<<std::endl;
	foundoscgroup = true;
	//break;
      }//else
      //std::cout<<"ED: "<<pdfname<<" sys group: "<<groupsined[j]<<std::endl;
      //}
      if (foundoscgroup){
	//std::cout<<fPdfManager.GetWorkingPdf(i).GetName()<<" "<<fPdfManager.GetWorkingPdf(i).Integral()<<std::endl;
	osc_loss.push_back(fPdfManager.GetWorkingPdf(i).Integral());
	//osc_loss_pdfnames.push_back(pdfname + "_norm");
      }else//{
	osc_loss.push_back(1.);
      //osc_loss_pdfnames.push_back(pdfname);
      //}
      //std::cout<<" "<<std::endl;
    }
    //std::cout<<"\n"<<std::endl;
    
    //double osc_loss = fPdfManager.GetWorkingPdf(0).Integral();
    
    /*std::vector<double> oscnormalisations = fPdfManager.GetNormalisations();
    for (size_t i = 0; i < oscnormalisations.size(); i++)
      oscnormalisations[i] = oscnormalisations[i] * osc_loss;
    fPdfManager.SetNormalisations(oscnormalisations);
    
    std::cout<<"\n 5) changed norm "<<std::endl;
    for (size_t i = 0; i < fPdfManager.GetNormalisations().size(); i++){
      std::cout<<"Norm: "<<fPdfManager.GetNormalisations()[i]<<std::endl;;
    }
    for (size_t i = 0; i < 1; i++){
      std::cout<<"Integral: "<<fPdfManager.GetWorkingPdf(0).Integral()<<std::endl;;
    }
    */
    
    // Apply Shrinking
    fPdfManager.ApplyShrink(fPdfShrinker);
    
    /*std::cout<<"\n 6) AppliedShrink "<<std::endl;
    for (size_t i = 0; i < fPdfManager.GetNormalisations().size(); i++){
      std::cout<<"Norm: "<<fPdfManager.GetNormalisations()[i]<<std::endl;;
    }
    for (size_t i = 0; i < 1; i++){
      std::cout<<"Integral: "<<fPdfManager.GetWorkingPdf(0).Integral()<<std::endl;;
    }
    */
    //std::cout<<"\n"<<std::endl;
    // loop over bins and calculate the likelihood
    double nLogLH = 0;
    //std::cout<<"fdatadist bins: "<<fDataDist.GetNBins()<<std::endl;
    for(size_t i = 0; i < fDataDist.GetNBins(); i++){
      //std::cout<<fDataDist.GetBinContent(i)<<std::endl;
      //double prob = fPdfManager.BinProbability(i);
      //std::cout<<"oscloss"<<std::endl;
      double prob = fPdfManager.BinProbability(i, osc_loss);
        if(!prob)
            throw std::runtime_error("BinnedNLLH::Encountered zero probability bin!");
        nLogLH -= fDataDist.GetBinContent(i) *  log(prob);        
    }
    /*
    std::cout<<"\n 7) Calculated LH basic "<<std::endl;
    for (size_t i = 0; i < fPdfManager.GetNormalisations().size(); i++){
      std::cout<<"Norm: "<<fPdfManager.GetNormalisations()[i]<<std::endl;;
    }
    for (size_t i = 0; i < 1; i++){
      std::cout<<"Integral: "<<fPdfManager.GetWorkingPdf(0).Integral()<<std::endl;;
    }
    */
    // Extended LH correction
    const std::vector<double>& normalisations = fPdfManager.GetNormalisations();
    for(size_t i = 0; i < normalisations.size(); i++)
      nLogLH += osc_loss[i] * normalisations.at(i);
      //nLogLH += normalisations.at(i);
    
    /*std::cout<<"\n 8) Extended LH "<<std::endl;
    for (size_t i = 0; i < fPdfManager.GetNormalisations().size(); i++){
      std::cout<<"Norm: "<<fPdfManager.GetNormalisations()[i]<<std::endl;;
    }
    for (size_t i = 0; i < 1; i++){
      std::cout<<"Integral: "<<fPdfManager.GetWorkingPdf(0).Integral()<<std::endl;;
    }
    */
    
    //don't need to change this part, osc_loss parameters don't affect calculation?
    //get rid of finding functions and just insist constraints must be placed in the same order??
    // Constraints
    /*for(std::map<std::string, QuadraticConstraint>::iterator it = fConstraints.begin();
        it != fConstraints.end(); ++it){
        std::vector<std::string>::iterator pdfnameit = std::find(osc_loss_pdfnames.begin(), osc_loss_pdfnames.end(), it->first);
	if (pdfnameit != osc_loss_pdfnames.end()){
	  int pdfnameindex = std::distance(osc_loss_pdfnames.begin(), pdfnameit);
	  //std::cout<<"Constraint name: "<<it->first<<" pdf name + _norm: "<<osc_loss_pdfnames[pdfnameindex]<<std::endl;
	  
	  nLogLH += it->second.Evaluate(fComponentManager.GetParameter(it->first),osc_loss[pdfnameindex]);
	}
	else{
	  //std::cout<<"NOT OSCILLATED PDF NORM CONSTRAINT: "<<it->first<<std::endl;
	  nLogLH += it->second.Evaluate(fComponentManager.GetParameter(it->first));
	}
    }
    */
    
    // Constraints
    for(std::map<std::string, QuadraticConstraint>::iterator it = fConstraints.begin();
        it != fConstraints.end(); ++it)      
        nLogLH += it->second.Evaluate(fComponentManager.GetParameter(it->first));
    

    /*for (size_t i = 0; i < oscnormalisations.size(); i++)
    oscnormalisations[i] = (oscnormalisations[i]/osc_loss);
    fPdfManager.SetNormalisations(oscnormalisations);
    
    std::cout<<"\n 9) Changed norms back "<<std::endl;
    for (size_t i = 0; i < fPdfManager.GetNormalisations().size(); i++){
      std::cout<<"Norm: "<<fPdfManager.GetNormalisations()[i]<<std::endl;;
    }
    for (size_t i = 0; i < 1; i++){
      std::cout<<"Integral: "<<fPdfManager.GetWorkingPdf(0).Integral()<<std::endl;;
      }*/

    return nLogLH;
}

void
BinnedNLLH::BinData(){
    fDataDist =  BinnedED(fPdfManager.GetOriginalPdf(0)); // make a copy for same binning and data rep
    fDataDist.Empty();
    CutLog log(fCuts.GetCutNames());
    DistFiller::FillDist(fDataDist, *fDataSet, fCuts, log);
    fCalculatedDataDist = true;    
    fSignalCutLog = log;
}

void
BinnedNLLH::AddDist(const std::vector<BinnedED>& pdfs, const std::vector<std::vector<std::string> >& sys_){
    if (pdfs.size() != sys_.size())
       throw DimensionError(Formatter()<<"BinnedNLLH:: #sys_ != #group_");
    for (int i = 0; i < pdfs.size(); ++i)
        AddDist( pdfs.at(i), sys_.at(i) );
}

void
BinnedNLLH::AddDist(const std::vector<BinnedED>& pdfs, const std::vector<std::vector<std::string> >& sys_, const std::vector<bool> ifosc){
    if (pdfs.size() != sys_.size())
       throw DimensionError(Formatter()<<"BinnedNLLH:: #sys_ != #group_");
    for (int i = 0; i < pdfs.size(); ++i){
        AddDist( pdfs.at(i), sys_.at(i) );
	if (ifosc[i])
	    fOscPdfs.push_back(pdfs.at(i).GetName());
    }
}

void
BinnedNLLH::AddDist(const BinnedED& pdf_, const std::vector<std::string>& syss_){
    fPdfManager.AddPdf(pdf_);
    fSystematicManager.AddDist(pdf_,syss_);
}

void
BinnedNLLH::AddDist(const BinnedED& pdf_, const std::vector<std::string>& syss_, const bool ifosc){
    fPdfManager.AddPdf(pdf_);
    fSystematicManager.AddDist(pdf_,syss_);
    if (ifosc)
      fOscPdfs.push_back(pdf_.GetName());
}

void
BinnedNLLH::AddDist(const BinnedED& pdf_){
    fPdfManager.AddPdf(pdf_);
    fSystematicManager.AddDist(pdf_,"");
}

void
BinnedNLLH::AddDist(const BinnedED& pdf_, const bool ifosc){
    fPdfManager.AddPdf(pdf_);
    fSystematicManager.AddDist(pdf_,"");
    if (ifosc)
      fOscPdfs.push_back(pdf_.GetName());
}

void
BinnedNLLH::AddPdf(const BinnedED& pdf_){
    AddDist(pdf_);
}

void
BinnedNLLH::SetPdfManager(const BinnedEDManager& man_){
    fPdfManager = man_;
}

void
BinnedNLLH::SetSystematicManager(const SystematicManager& man_){
    fSystematicManager = man_;
}

void 
BinnedNLLH::AddSystematic(Systematic* sys_){
    fSystematicManager.Add(sys_);
}

void
BinnedNLLH::AddSystematic(Systematic* sys_, const std::string&  group_){
    fSystematicManager.Add(sys_, group_);
}

/*void
BinnedNLLH::AddSystematic(Systematic* sys_, const std::string&  group_, const bool ifosc){
    if (ifosc)
      fIfOscSystematics.push_back(group_);
    fSystematicManager.Add(sys_, group_);
    }*/

void
BinnedNLLH::SetDataSet(DataSet* dataSet_){
    fDataSet = dataSet_;
    fCalculatedDataDist = false;
}

DataSet*
BinnedNLLH::GetDataSet(){
    return fDataSet;
}

void
BinnedNLLH::SetDataDist(const BinnedED& binnedPdf_){
    fDataDist = binnedPdf_;
    fCalculatedDataDist = true;
}

BinnedED
BinnedNLLH::GetDataDist() const{
    return fDataDist;
}


void
BinnedNLLH::SetBuffer(size_t dim_, unsigned lower_, unsigned upper_){
    fPdfShrinker.SetBuffer(dim_, lower_, upper_);
}

std::pair<unsigned, unsigned>
BinnedNLLH::GetBuffer(size_t dim_) const{
    return fPdfShrinker.GetBuffer(dim_);
}

void
BinnedNLLH::SetBufferAsOverflow(bool b_){
    fPdfShrinker.SetUsingOverflows(b_);
}

bool
BinnedNLLH::GetBufferAsOverflow() const{
    return fPdfShrinker.GetUsingOverflows();
}

void
BinnedNLLH::AddSystematics(const std::vector<Systematic*> systematics_){
    for(size_t i = 0; i < systematics_.size(); i++)
        AddSystematic(systematics_.at(i));
}

void
BinnedNLLH::AddSystematics(const std::vector<Systematic*> sys_, const std::vector<std::string> & groups_){
    if (groups_.size() != sys_.size())
       throw DimensionError(Formatter()<<"BinnedNLLH:: #sys_ != #group_");
    for(size_t i = 0; i <sys_.size(); i++)
        AddSystematic(sys_.at(i), groups_.at(i));
}
/*
void
BinnedNLLH::AddSystematics(const std::vector<Systematic*> sys_, const std::vector<std::string> & groups_, const std::vector<bool> ifosc){
    if (groups_.size() != sys_.size())
       throw DimensionError(Formatter()<<"BinnedNLLH:: #sys_ != #group_");
    for(size_t i = 0; i <sys_.size(); i++){
        AddSystematic(sys_.at(i), groups_.at(i));
	if (ifosc[i])
	  fIfOscSystematics.push_back(groups_.at(i));
    }
}
*/
void
BinnedNLLH::SetNormalisations(const std::vector<double>& norms_){    
    fPdfManager.SetNormalisations(norms_);
}

std::vector<double>
BinnedNLLH::GetNormalisations() const{
    return fPdfManager.GetNormalisations();
}

void
BinnedNLLH::AddCut(const Cut& cut_){
    fCuts.AddCut(cut_);
}

void 
BinnedNLLH::SetCuts(const CutCollection& cuts_){
    fCuts = cuts_;
}

void 
BinnedNLLH::SetConstraint(const std::string& paramName_, double mean_, double sigma_){
    fConstraints[paramName_] = QuadraticConstraint(mean_, sigma_);
}


double
BinnedNLLH::GetSignalCutEfficiency() const{
    return fSignalCutEfficiency;
}

void
BinnedNLLH::SetSignalCutEfficiency(double eff_){
    fSignalCutEfficiency = eff_;
}

CutLog
BinnedNLLH::GetSignalCutLog() const{
    return fSignalCutLog;
}

void
BinnedNLLH::SetSignalCutLog(const CutLog& lg_){
    fSignalCutLog = lg_;
}

/////////////////////////////////////////////////////////
// Declare which objects should be adjusted by the fit //
/////////////////////////////////////////////////////////
void
BinnedNLLH::RegisterFitComponents(){
    fComponentManager.Clear();
    fComponentManager.AddComponent(&fPdfManager);
    
    //Because the limits are set by name they can be added in any order.
    const std::map<std::string, std::vector<Systematic*> > sys_ = fSystematicManager.GetSystematicsGroup();
    std::vector<std::string> alreadyAdded;
    for (std::map<std::string, std::vector<Systematic*> >::const_iterator group_ = sys_.begin(); group_ !=sys_.end(); ++group_) {
        for (int i = 0; i < group_->second.size(); ++i) {
            if( std::find( alreadyAdded.begin() , alreadyAdded.end() , group_->second.at(i)->GetName() ) == alreadyAdded.end() ){
                fComponentManager.AddComponent( group_->second.at(i) );
                alreadyAdded.push_back( group_->second.at(i)->GetName() );
            }
        }//End of group
    }//End of groups
}


void
BinnedNLLH::SetParameters(const ParameterDict& params_){
    try{
        fComponentManager.SetParameters(params_);
    }
    catch(const ParameterError& e_){
        throw ParameterError(std::string("BinnedNLLH::") + e_.what());
    }
}
                                             
ParameterDict
BinnedNLLH::GetParameters() const{
    return fComponentManager.GetParameters();
}

int
BinnedNLLH::GetParameterCount() const{
    return fComponentManager.GetTotalParameterCount();
}

std::set<std::string>
BinnedNLLH::GetParameterNames() const{
    return fComponentManager.GetParameterNames();
}

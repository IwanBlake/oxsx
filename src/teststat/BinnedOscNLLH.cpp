#include <BinnedOscNLLH.h>
#include <math.h>
#include <DataSet.h>
#include <Exceptions.h>
#include <DistFiller.h>
#include <CutLog.h>
#include <iostream>

double 
BinnedOscNLLH::Evaluate(){
    if(!fDataSet && !fCalculatedDataDist) 
        throw LogicError("BinnedOscNLLH function called with no data set and no DataDist! set one of these first");
    
    if (!fCalculatedDataDist)
        BinData();
    
    if(!fAlreadyShrunk){
        fDataDist = fPdfShrinker.ShrinkDist(fDataDist);
        fAlreadyShrunk = true;
    }
    
    // Apply oscillations here?
    //fPdfManager.ApplyOscillations();

    // Construct systematics
    fSystematicManager.Construct();

    // Apply systematics
    fPdfManager.ApplySystematics(fSystematicManager);

    // Apply Shrinking
    fPdfManager.ApplyShrink(fPdfShrinker);

    // loop over bins and calculate the likelihood
    double nLogLH = 0;
    for(size_t i = 0; i < fDataDist.GetNBins(); i++){
        double prob = fPdfManager.BinProbability(i);
        if(!prob)
            throw std::runtime_error("BinnedOscNLLH::Encountered zero probability bin!");
        nLogLH -= fDataDist.GetBinContent(i) *  log(prob);        
    }


    // Extended LH correction
    const std::vector<double>& normalisations = fPdfManager.GetNormalisations();
    for(size_t i = 0; i < normalisations.size(); i++)
        nLogLH += normalisations.at(i);
            
    // Constraints
    for(std::map<std::string, QuadraticConstraint>::iterator it = fConstraints.begin();
        it != fConstraints.end(); ++it)
        nLogLH += it->second.Evaluate(fComponentManager.GetParameter(it->first));

    return nLogLH;
}

double
BinnedOscNLLH::NuSurvProb(double nuE, double baseline, double delmsqr21, double sinsqrtheta12, double sinsqrtheta13) {
  double fSSqr2Theta12 = 0;//pow(sin(2.0 * TMath::ASin(sqrt(sinsqrtheta12))), 2.0);
  double fS4 = pow(sinsqrtheta13, 2.0);
  double fC4 = pow(1.0-sinsqrtheta13, 2.0);
  double scale = 1.267e3; // for nuE in [MeV] and baseline in [km]
  double sSqrDmBE = pow(sin(scale * delmsqr21 * baseline / nuE), 2.0);
  double fOscProb = (fC4 * (1.0 - fSSqr2Theta12 * sSqrDmBE) + fS4);
  return fOscProb;
}

void
BinnedOscNLLH::ntOscillate( const char* ntin, const char* ntout, double delmsqr21, double sinsqrtheta12, double sinsqrtheta13) {

  // TFile *fin = new TFile(ntin);
  // TTree *inT = (TTree*)fin->Get("nt");
  // int nentries = inT->GetEntries();
  // double survprob;
  // float ParKE;
  // float ReactorDistance;
  // inT->SetBranchAddress("ParKE", &ParKE);
  // inT->SetBranchAddress("ReactorDistance", &ReactorDistance);
  // TFile *fout = new TFile(ntout,"RECREATE");
  // TTree *outT = inT->CloneTree(0);

  // for (int i = 0; i < nentries ; i++){
    // inT->GetEntry(i);
    // survprob = NuSurvProb(ParKE, ReactorDistance, delmsqr21, sinsqrtheta12, sinsqrtheta13);
    // const double random = CLHEP::HepUniformRand();
    // std::cout<<"Rand:   ";
    // std::cout<<random<<std::endl;
    // std::cout<<"Prob:   ";
    // std::cout<<survprob<<std::endl;
    // if (survprob > random){
      // outT->Fill();
    // }
  // }
  // outT->Print();
  // outT->AutoSave();
  //delete fin;
  //delete fout;
}

void
BinnedOscNLLH::BinData(){
    fDataDist =  BinnedED(fPdfManager.GetOriginalPdf(0)); // make a copy for same binning and data rep
    fDataDist.Empty();
    CutLog log(fCuts.GetCutNames());
    DistFiller::FillDist(fDataDist, *fDataSet, fCuts, log);
    fCalculatedDataDist = true;    
    fSignalCutLog = log;
}

void
BinnedOscNLLH::SetPdfManager(const BinnedEDManager& man_){
    fPdfManager = man_;
}

void
BinnedOscNLLH::SetSystematicManager(const SystematicManager& man_){
    fSystematicManager = man_;
}

void
BinnedOscNLLH::AddPdf(const BinnedED& pdf_){
    fPdfManager.AddPdf("normalisation", pdf_);
}

void
BinnedOscNLLH::AddPdf(const std::string& name_, const BinnedED& pdf_){
    fPdfManager.AddPdf(pdf_);
}

void 
BinnedOscNLLH::AddSystematic(Systematic* sys_){
    fSystematicManager.Add(sys_);
}

void
BinnedOscNLLH::SetDataSet(DataSet* dataSet_){
    fDataSet = dataSet_;
    fCalculatedDataDist = false;
}

DataSet*
BinnedOscNLLH::GetDataSet(){
    return fDataSet;
}

void
BinnedOscNLLH::SetDataDist(const BinnedED& binnedPdf_){
    fDataDist = binnedPdf_;
    fCalculatedDataDist = true;
}

BinnedED
BinnedOscNLLH::GetDataDist() const{
    return fDataDist;
}


void
BinnedOscNLLH::SetBuffer(size_t dim_, unsigned lower_, unsigned upper_){
    fPdfShrinker.SetBuffer(dim_, lower_, upper_);
}

std::pair<unsigned, unsigned>
BinnedOscNLLH::GetBuffer(size_t dim_) const{
    return fPdfShrinker.GetBuffer(dim_);
}

void
BinnedOscNLLH::SetBufferAsOverflow(bool b_){
    fPdfShrinker.SetUsingOverflows(b_);
}

bool
BinnedOscNLLH::GetBufferAsOverflow() const{
    return fPdfShrinker.GetUsingOverflows();
}

void
BinnedOscNLLH::AddPdfs(const std::vector<BinnedED>& pdfs_){
  for(size_t i = 0; i < pdfs_.size(); i++)
    AddPdf(pdfs_.at(i));
}

void
BinnedOscNLLH::AddSystematics(const std::vector<Systematic*> systematics_){
  for(size_t i = 0; i < systematics_.size(); i++)
    AddSystematic(systematics_.at(i));
}

void
BinnedOscNLLH::SetNormalisations(const std::vector<double>& norms_){    
    fPdfManager.SetNormalisations(norms_);
}

std::vector<double>
BinnedOscNLLH::GetNormalisations() const{
    return fPdfManager.GetNormalisations();
}

void
BinnedOscNLLH::AddCut(const Cut& cut_){
    fCuts.AddCut(cut_);
}

void 
BinnedOscNLLH::SetCuts(const CutCollection& cuts_){
    fCuts = cuts_;
}

void 
BinnedOscNLLH::SetConstraint(const std::string& paramName_, double mean_, double sigma_){
    fConstraints[paramName_] = QuadraticConstraint(mean_, sigma_);
}


double
BinnedOscNLLH::GetSignalCutEfficiency() const{
    return fSignalCutEfficiency;
}

void
BinnedOscNLLH::SetSignalCutEfficiency(double eff_){
    fSignalCutEfficiency = eff_;
}

CutLog
BinnedOscNLLH::GetSignalCutLog() const{
    return fSignalCutLog;
}

void
BinnedOscNLLH::SetSignalCutLog(const CutLog& lg_){
    fSignalCutLog = lg_;
}

/////////////////////////////////////////////////////////
// Declare which objects should be adjusted by the fit //
/////////////////////////////////////////////////////////
void
BinnedOscNLLH::RegisterFitComponents(){
    fComponentManager.Clear();
    fComponentManager.AddComponent(&fPdfManager);
    for(size_t i = 0; i < fSystematicManager.GetSystematics().size(); i++)
        fComponentManager.AddComponent(fSystematicManager.GetSystematics().at(i));
}

void
BinnedOscNLLH::SetParameters(const ParameterDict& params_){
    try{
        fComponentManager.SetParameters(params_);
    }
    catch(const ParameterError& e_){
        throw ParameterError(std::string("BinnedOscNLLH::") + e_.what());
    }
}
                                             
                 
ParameterDict
BinnedOscNLLH::GetParameters() const{
    return fComponentManager.GetParameters();
}

int
BinnedOscNLLH::GetParameterCount() const{
    return fComponentManager.GetTotalParameterCount();
}

std::set<std::string>
BinnedOscNLLH::GetParameterNames() const{
    return fComponentManager.GetParameterNames();
}

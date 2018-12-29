#include <TFile.h>
#include <TMath.h>
#include <string>
#include <TNtuple.h>
#include <iostream>
#include <TObject.h>
#include <CLHEP/Random/Randomize.h>

double NuSurvProb(double nuE, double baseline, double delmsqr21, double sinsqrtheta12, double sinsqrtheta13)
{
  double fSSqr2Theta12 = pow(sin(2.0 * TMath::ASin(sqrt(sinsqrtheta12))), 2.0);
  double fS4 = pow(sinsqrtheta13, 2.0);
  double fC4 = pow(1.0-sinsqrtheta13, 2.0);
  double scale = 1.267e3; // for nuE in [MeV] and baseline in [km]
  double sSqrDmBE = pow(sin(scale * delmsqr21 * baseline / nuE), 2.0);
  double fOscProb = (fC4 * (1.0 - fSSqr2Theta12 * sSqrDmBE) + fS4);
  return fOscProb;
}

void ntOscillate( const char* ntin, const char* ntout, double delmsqr21, double sinsqrtheta12, double sinsqrtheta13) {

  TFile *fin = new TFile(ntin);
  TTree *inT = (TTree*)fin->Get("nt");
  int nentries = inT->GetEntries();
  double survprob;
  float ParKE;
  float ReactorDistance;
  inT->SetBranchAddress("ParKE", &ParKE);
  inT->SetBranchAddress("ReactorDistance", &ReactorDistance);
  
  std::stringstream fileout;
  //std::string namefile;
  //namefile.push_back(ntout);
  fileout<<ntout<<"ds21_"<<delmsqr21<<"_ss12_"<<sinsqrtheta12<<"_ss13_"<<sinsqrtheta13<<".root";
  TFile *fout = new TFile(fileout.str().c_str(),"RECREATE");
  TTree *outT = inT->CloneTree(0);
  
  for (int i = 0; i < nentries ; i++){
    inT->GetEntry(i);
    survprob = NuSurvProb(ParKE, ReactorDistance, delmsqr21, sinsqrtheta12, sinsqrtheta13);
    std::cout<<"Distance "<<ReactorDistance<<std::endl;
    const double random = CLHEP::HepUniformRand();
    if (survprob > random){
      outT->Fill();
    }
  }
  outT->Print();
  outT->AutoSave();
  delete fin;
  delete fout;
}

int main(int argc, char* argv[])
{
  if (argc != 6) {
    std::cout<<"5 arguments expected: \n 1: location of input unoscillated pruned ntuple \n 2: \
location/filename of output oscillated (UP TO .root!)\n 3: delmsqr21 \n 4: sinsqrtheta12 \n 5:sinsqrtheta13"<<std::endl;
  }else{
    const char* rootin = argv[1];
    const char* rootout = argv[2];

    double delmsqr21 = atof(argv[3]);
    double sinsqrtheta12 = atof(argv[4]);
    double sinsqrtheta13 = atof(argv[5]);
    ntOscillate(rootin, rootout, delmsqr21, sinsqrtheta12, sinsqrtheta13);
  }
}

/* Example shows how to create Prior objects which can be used to constrain
 * floated variable by applying a penalty term to the fitting function. The
 * magnitude of the penalty term is given by a class inheriting from the
 * Function class (interface).
 *
 * A Prior can depend solely on the primary variable as well as other terms
 * known the ComponentManager. Dependences can be added to the Prior with the
 * SetDependents member function.
 *
 * Once a Prior is set up it may be added to a fitting function by calling the
 * AddPrior() function on said fitting function.
 *
 */

#include <iostream>
#include <iomanip>
#include <Prior.h>
#include <Function.h>
#include <PriorManager.h>
#include <Gaussian.h>
#include <ParameterDict.h>
#include <math.h>
#include <TRandom3.h>
class Func: public Function
{
  public:
    Func(int num):Function(), fNdims(num){;}
    double operator() (const std::vector<double>& params) const{
      double exponent = 0;
      double mean;
      double stdDev =1;
      double nDevs;

      double diff = params[0];
      for (int i = 1; i < fNdims; ++i) {
        diff -= params[i];
      }
      nDevs = (diff)/stdDev;
      exponent += nDevs * nDevs;

      double norm = 1;
      norm *= sqrt(2*M_PI)*stdDev;

      return exp(- 0.5* exponent) / norm; 
    }
    void SetParameter(const std::string&, double){;}
    double GetParameter(const std::string&) const{;}
    void SetParameters(const ParameterDict&){;}
    ParameterDict GetParameters() const{;}
    size_t GetParameterCount() const{;}
    void RenameParameter(const std::string&, const std::string&){;}
    std::string GetName() const{;}
    void SetName(const std::string&){;}
    Function* Clone() const{
      return static_cast<Function*> (new Func(fNdims));
    }
    int GetNDims() const{;}
    std::set<std::string> GetParameterNames() const{;}

  private:
    int fNdims;
};

void GaussianPrior(ParameterDict dict);
void PriorWithDependence(ParameterDict dict);
void PriorWith2Dependence(ParameterDict dict);

int main(int argc, char *argv[])
{
  std::cout << "Starting prior example" << std::endl; 

  ParameterDict dict;
  dict["alpha"]  =10;
  dict["beta"]   =2;
  dict["gamma"]  =4;
  dict["delta"]  =5;
  dict["epsilon"]=9;

  GaussianPrior(dict);
  PriorWithDependence(dict);
  PriorWith2Dependence(dict);

  return 0;
}

void GaussianPrior(ParameterDict dict){
  // Setting a prior on a single floating parameter.
  Prior* pir = new Prior();
  // Set the function used to constrain the parameter. 
  // pir->SetFunction(new Gaussian(4,0.3));
  pir->SetFunction(new Gaussian(4,3));
  // Set the primary variable of by name. 
  pir->SetPrimary("alpha");

  // Then you're done.
  // To add to a fitting function e.g. BinnedNLLH. Use
  // the AddPrior function. e.g. lh.AddPrior(&pir); 

  {
    // This is a simple optimisation loop. Just to give some kind of output
    // showing that the parameter "alpha" changes according to the output of
    // the probability method which uses the underline function to put a numerical value on the difference between the current values and the required value.   

    double perviousProb = 0;
    bool Qup =false;
    TRandom3 rand;

    std::cout << "\nGaussian is centred about 4" << std::endl; 
    std::cout<<"Starting dummy optimisation"<<std::endl;
    for (int i = 0; i < 15; ++i) {
      if(i==0){
        double randd= rand.Integer(1);
        dict["beta"] = dict["beta"]+pow(1,randd);
        if(randd<0)
          Qup = false;
        else
          Qup = true;
      }
      double prob =pir->Probability(dict); 
      std::cout<<"alpha = "<<std::setprecision(3)<< dict["alpha"]<< std::fixed 
        << std::setw(8)<<std::setprecision(6)
        << " prob = "<< prob << std::setprecision(3)
        <<"\t => Term added to a BinnedNLLH = "<<-log(prob)<<std::endl;
      if(prob < perviousProb){
        if(Qup == true){
          dict["alpha"]  =dict["alpha"]-1;
          Qup= false;
        }else{
          dict["alpha"]  =dict["alpha"]+1;
          Qup= true;
        }

      }else{
        if(Qup == true){
          dict["alpha"]  =dict["alpha"]+1;
        }else{
          dict["alpha"]  =dict["alpha"]-1;
        }
      }
      perviousProb = prob;
    }
  }
}

void PriorWithDependence(ParameterDict dict){
  // For priors that depend on a single other floated parameters.

  // The function relating the primary variable value to its dependences can be
  // any user defined function that inherits from Function. 
  //
  // The vector that the function acts on must be of the form {fPrimary,dep} 
  //
  // (fyi Func is a gaussian constraint centred on the dependent variable.)
  Func* func= new Func(2);
  Prior* pir = new Prior();
  pir->SetFunction(func);
  // Set the primary variable of by name. 
  pir->SetPrimary("alpha");
  // Set the dependent variable by passing the parameter name. It
  // is important that the order matches that of the dependants in the function
  // definition.

  pir->SetDependent("beta");

  //Now everything is setup
  //To add to a fitting function e.g. BinnedNLLH. Use
  //the AddPrior function. e.g. lh.AddPrior(&pir); 

  {
    // This is a simple optimisation loop. Just to give some kind of output
    // showing that changing the ParameterDict passed to the probability
    // method returns 

    double perviousProb = 0;
    bool Qup =false;
    TRandom3 rand;

    std::cout<<"\nStarting single dependent dummy optimisation"<<std::endl;
    for (int i = 0; i < 15; ++i) {
      if(i==0){
        double randd= rand.Integer(1);
        dict["beta"] = dict["beta"]+pow(1,randd);
        if(randd<0)
          Qup = false;
        else
          Qup = true;
      }
      double prob =pir->Probability(dict); 
      std::cout<<"alpha = "<<std::setprecision(3)<< dict["alpha"]<< std::fixed 
        << std::setw(10)<< " beta = "<< dict["beta"] <<std::setprecision(6)
        << " prob = "<< prob << std::setprecision(3)
        <<"\t => Term added to a BinnedNLLH = "<< -log(prob)<<std::endl;
      if(prob < perviousProb){
        if(Qup == true){
          dict["alpha"]  =dict["alpha"]-1;
          Qup= false;
        }else{
          dict["alpha"]  =dict["alpha"]+1;
          Qup= true;
        }

      }else{
        if(Qup == true){
          dict["alpha"]  =dict["alpha"]+1;
        }else{
          dict["alpha"]  =dict["alpha"]-1;
        }
      }
      perviousProb = prob;
    }
  }
}

void PriorWith2Dependence(ParameterDict dict){
  // The function relating the primary variable value to its dependences can be
  // any user defined function that inherits from Function. 
  //
  // The vector that the function acts on must be of the form {fPrimary,dep_1,dep_2,...} 
  //
  // (FYI Func is a gaussian constraint centred on the sum dependent variables. )
  Func* func= new Func(3);
  // Set the dependent variable by passing a vector of parameter names. It
  // is important that the order matches that of the dependants in the function
  // definition.
  std::vector<std::string> deps;
  deps.push_back("gamma");
  deps.push_back("delta");

  // beta being the primary variable.
  Prior* pir2 = new Prior("beta",deps,func);

  {
    double perviousProb = 0;
    bool Qup =false;
    TRandom3 rand;
    std::cout<<"\nStarting double dependences dummy optimisation"<<std::endl;
    for (int i = 0; i < 15; ++i) {
      if(i==0){
        double randd= rand.Integer(1);
        dict["beta"] = dict["beta"]+pow(1,randd);
        if(randd<0)
          Qup = false;
        else
          Qup = true;
      }
      double prob =pir2->Probability(dict); 
      std::cout<<"beta = "<<std::setprecision(3)<< dict["beta"]<< std::fixed 
        << std::setw(10)<< " gamma + delta = "<< dict["gamma"] + dict["delta"] << std::setprecision(6)
        << " prob = "<< prob << std::setprecision(3)
        <<"\t => Term added to a BinnedNLLH = "<< -log(prob)<<std::endl;
      if(prob < perviousProb){
        if(Qup == true){
          dict["beta"]  =dict["beta"]-1;
          Qup= false;
        }else{
          dict["beta"]  =dict["beta"]+1;
          Qup= true;
        }

      }else{
        if(Qup == true){
          dict["beta"]  =dict["beta"]+1;
        }else{
          dict["beta"]  =dict["beta"]-1;
        }
      }
      perviousProb = prob;
    }
  }
}

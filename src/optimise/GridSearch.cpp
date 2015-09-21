#include <GridSearch.h>
#include <TestStatistic.h>

void 
GridSearch::SetMinima(const std::vector<double>& minima_){
    fMinima = minima_;
}

void 
GridSearch::SetMaxima(const std::vector<double>& maxima_){
    fMaxima = maxima_;
}

void
GridSearch::SetStepSizes(const std::vector<double>& steps_){
    fStepSizes = steps_;
}

std::vector<double> 
GridSearch::GetMaxima() const{
    return fMaxima;
}

std::vector<double> 
GridSearch::GetMinima() const{
    return fMinima;
}

std::vector<double>
GridSearch::GetStepSizes() const{
    return fStepSizes;
}

void 
GridSearch::Optimise(){
    // list of rates followed by list of systematics
    fBestFit.resize(pTestStatistic -> GetNParams());
    fMinVal = 0;

    // start at min value
    fParams = fMinima;
    while(Increment(0)){
        // calculate the new value
        // if bigger, grab this as new best fit
        pTestStatistic->SetParams(fParams);

        double currentVal = pTestStatistic->Evaluate();

        if (currentVal < fMinVal || fMinVal == 0){
                fMinVal = currentVal;
                fBestFit = fParams;
        }
    } 
}

bool 
GridSearch::Increment(size_t index_){
    fParams[index_] += fStepSizes.at(index_);

    // wrap around past the maximum
    if (fParams[index_] > fMaxima.at(index_)){
        fParams[index_] = fMinima.at(index_);

        // if its the last index no rippling to do
        if (index_ == fStepSizes.size() - 1)
            return false;

        // ripple up to next index
        if (Increment(index_ + 1)) // propagate the false down
            return true;
        else
            return false;

    }
    return true;
}

std::vector<double> 
GridSearch::GetBestFit() const{
    return fBestFit;
}

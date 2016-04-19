#include <PdfAxis.h>
#include <algorithm>
#include <Exceptions.h>
#include <iostream>


PdfAxis::PdfAxis(const std::string& name_, double min_, double max_, size_t nBins_,
                 const std::string& latexName_){
    fName = name_;
    fLatexName = latexName_;
    if(fLatexName == "")
        fLatexName = name_;

    fMin = min_;
    fMax = max_;

    if (fMin >= fMax || !nBins_)
        throw LogicError("PdfAxis: Invalid bin specification: min > max or nbins = 0");

    fNBins = nBins_;
    double binWidth =  double(fMax - fMin) /fNBins;

    fBinLowEdges.resize(fNBins, 0);    
    fBinHighEdges.resize(fNBins, 0);
    fBinCentres.resize(fNBins, 0);
	fBinWidths.resize(fNBins, binWidth);

    for(size_t i = 0; i < fNBins; i++){
        fBinLowEdges[i]  = fMin  + i          * binWidth;
        fBinCentres[i]   = fMin  + (i + 0.5)  * binWidth;
        fBinHighEdges[i] = fMin  + (i+1)      * binWidth;
    }
}

PdfAxis::PdfAxis(const std::string& name_, const std::vector<double>& lowEdges_,
                 const std::vector<double>& highEdges_, const std::string& latexName_){

    if (highEdges_.size() != lowEdges_.size() || !lowEdges_.size())
        throw LogicError("PdfAxis: Invalid bins. Must have the same number of low and high edges!");

    fName = name_;
    fBinLowEdges  = lowEdges_;
    fBinHighEdges = highEdges_;

    for(size_t i = 0; i < fBinLowEdges.size(); i++){
        if(i && fBinLowEdges.at(i-1) > fBinLowEdges.at(i))
            throw LogicError("PdfAxis: Low edges not ordered!");
        if(fBinLowEdges.at(i) > fBinHighEdges.at(i))
            throw LogicError("PdfAxis: Bin Low Edge is bigger than equivilent high edge!");
    }

    fNBins = fBinLowEdges.size();

    fBinCentres.resize(fNBins, 0);
	fBinWidths.resize(fNBins, 0);

    for(size_t i = 0; i < fBinLowEdges.size(); i++){
	    fBinCentres[i] = 0.5 * (fBinLowEdges.at(i) + fBinHighEdges.at(i));
		fBinWidths[i]  = fBinHighEdges.at(i) - fBinLowEdges.at(i);
	}

    fLatexName = latexName_;
    fMin = fBinLowEdges.at(0);
    fMax = fBinHighEdges.back();

    if("" == fLatexName)
        fLatexName = name_;
}


size_t 
PdfAxis::FindBin(double value_) const{
  size_t insertIndex = std::lower_bound(fBinHighEdges.begin(), fBinHighEdges.end(), value_) - fBinHighEdges.begin();
    if (insertIndex == fNBins)
        return insertIndex-1;
    return insertIndex;
}

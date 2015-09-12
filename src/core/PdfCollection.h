/***************************************************************************************************/
/* A collection of PDFs with normalisations for calculating the probability of a given observation */
/* given some rates and some pdfs                                                                  */
/***************************************************************************************************/

#ifndef __PDF_COLLECTION__
#define __PDF_COLLECTION__
#include <vector>
#include <BinnedPdf.h>

class DataHandler;
class Systematic;

class PdfCollection{
 public:
    PdfCollection() {}
    ~PdfCollection();
    size_t GetNDims() const;
    double Probability(const DataHandler&) const;
    void   ApplySystematics(const std::vector<Systematic>&);

    std::vector<double> GetNormalisations() const;
    void SetNormalisations(const std::vector<double>& normalisations_);

 private:
    std::vector<BinnedPdf> fOriginalPdfs;
    std::vector<BinnedPdf> fWorkingPdfs;
    std::vector<double>    fNormalisations;
    size_t fNDims;
};
#endif

/***********************************************************************************************/
/* A quadratic constraint on a fit parameter, for log likelihood and Chi-Square tests this is  */
/* equivlent to a gaussian contraint.                                                          */
/***********************************************************************************************/
#ifndef __OXSX_QUADRATIC_CONSTRAINT__
#define __OXSX_QUADRATIC_CONSTRAINT__

class QuadraticConstraint{
 public:
 QuadraticConstraint(){}
 QuadraticConstraint(double mean_, double width_) : 
    fMean(mean_), fWidth(width_) {}

    double Evaluate(double val_) const {
        return (val_ - fMean) * (val_ - fMean) / (2 * fWidth * fWidth);
    }
    
    /*double Evaluate(double val_, double osc_loss) const {
      double frac = fWidth/fMean;
      
      return (val_ - fMean) * (val_ - fMean) / (2 * fWidth * fWidth);
      }*/
    
 private:
    double fMean;
    double fWidth;
};
#endif

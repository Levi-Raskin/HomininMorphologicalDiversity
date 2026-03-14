#ifndef HomininDisparityAnalysis_hpp
#define HomininDisparityAnalysis_hpp

#include <string>
#include <vector>

class HomininDisparityAnalysis{
    public:
                                    HomininDisparityAnalysis(std::string treeIn, std::string treeOut, std::string shapeOut);
        void                        run(int numAdditionalTaxa, //number of additional taxa to simulate
                                        int nreps, // number of repetitions per simulated condition
                                        double lddmmSigma, // LDDMM inter-landmark correlation
                                        double lddmmAlpha, // LDDMM landmark rate
                                        double betaDistAlpha, // beta dist alpha parameter for tree simulations
                                        double betaDistBeta, // beta dist beta parameter for tree simulations
                                        int numLandmarks, // vector of the number of landmarks tested
                                        int nthreads);
    private:
        std::vector<std::string>    backboneTreeNewicks;
        std::string                 shapeOutfile;
        std::string                 treeOutfile;
};

#endif

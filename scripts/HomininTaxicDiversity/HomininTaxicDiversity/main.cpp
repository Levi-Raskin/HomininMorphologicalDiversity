#include <filesystem>
#include <iostream>
#include <vector>

#include "HomininDisparityAnalysis.hpp"

int main(int argc, const char* argv[]) {
    // Simulation parameters:
    //constant parameters
    int numReps = 100;
    double lddmmSigma = 1.0;
    double betaDistAlpha = 1.0;
    double betaDistBeta = 1.0;
    int threads = 10;
    
    //changing parameters
    std::vector<double> alphas = {0.1, 0.2, 0.3, 0.4};
    std::vector<int> numAdditionalTaxa = {3, 5, 8, 11, 15, 20, 25, 50, 100, 144}; //abstract explicitly mentions 5, 11, 25,50, and 144
//    std::vector<int> numLandmarks = {10, 25, 50};
    std::vector<int> numLandmarks = {25, 50};
    
    std::string base ="/Users/levir/Documents/GitHub/HomininTaxicDiversity/results/";
    std::string treeIn = "/Users/levir/Documents/GitHub/HomininTaxicDiversity/data/sampledTrees.tsv";

    int totalIter = numLandmarks.size() * alphas.size() * numAdditionalTaxa.size();
    int iter = 0;
    auto loopStart = std::chrono::steady_clock::now();

    for(int lm : numLandmarks){
        for(double a : alphas){
            for(int nt : numAdditionalTaxa){
                std::cout << "Current simulation:  " << "\n";
                std::cout << "Num landmarks:       " << lm << "\n";
                std::cout << "LDDMM alpha:         " << a << "\n";
                std::cout << "Num additional taxa: " << nt << "\n";
                
                // ETA code
                if(iter > 0){
                    auto now = std::chrono::steady_clock::now();
                    double elapsed = std::chrono::duration<double>(now - loopStart).count();
                    double avgPerIter = elapsed / iter;
                    double eta = avgPerIter * (totalIter - iter);
                    int etaMin = (int)(eta / 60);
                    int etaSec = (int)(eta) % 60;
                    std::cout << "[" << iter << "/" << totalIter << "] ETA: " << etaMin << "m " << etaSec << "s" << " | Average per iteration: " << avgPerIter  <<  "\n";
                }
                
                int numAddTaxa = nt;
                double lddmmAlpha = a;
                int numLM = lm;
                
                
                //helperlambda for tree out/shape out
                auto fmt_path = [](int lm, int nt, double a) {
                    char buf[256];
                    std::snprintf(buf, sizeof(buf), "numLandmarks%d/numAdditionalTaxa%d/alpha%.2f/", lm, nt, a);
                    return std::string(buf);
                };
                
                std::string treeOut =   base +
                                        "simulatedTrees/" +
                                        fmt_path(lm , nt , a);
                std::string shapeOut =  base +
                                        "simulatedShapes/"+
                                        fmt_path(lm , nt , a);
                                        
                std::error_code ec;
                std::filesystem::create_directories(treeOut, ec);
                if(ec) std::cerr << "Failed to create treeOut: " << ec.message() << "\n";
                std::filesystem::create_directories(shapeOut, ec);
                if(ec) std::cerr << "Failed to create shapeOut: " << ec.message() << "\n";
                
                HomininDisparityAnalysis disp(treeIn, treeOut, shapeOut);
                disp.run(numAddTaxa, numReps, lddmmSigma, lddmmAlpha, betaDistAlpha, betaDistBeta, numLM, threads);
                iter++;
            }
        }
    }
    
    return 0;
}

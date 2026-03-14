#include "Eigen/Core"
#include "HomininDisparityAnalysis.hpp"
#include "ReadTSV.hpp"
#include "SimShapeLDDMMv2.hpp"
#include "Tree.hpp"
#include "Utility.hpp"
#include "WriteTSV.hpp"

#include <iostream>

HomininDisparityAnalysis::HomininDisparityAnalysis(std::string treeIn, std::string treeOut, std::string shapeOut) : treeOutfile(treeOut), shapeOutfile(shapeOut){

    ReadTSV r = ReadTSV(treeIn, true, true, true);

    for(int i = 0; i < r.getReadStringData().size(); i++)
        backboneTreeNewicks.push_back(r.getReadStringData()[i][0]);
}

void HomininDisparityAnalysis::run( int numAdditionalTaxa, //number of additional taxa to simulate
                                    int nreps, // number of repetitions per simulated condition
                                    double lddmmSigma, // LDDMM inter-landmark correlation
                                    double lddmmAlpha, // LDDMM landmark rate
                                    double betaDistAlpha, // beta dist alpha parameter for tree simulations
                                    double betaDistBeta, // beta dist beta parameter for tree simulations
                                    int numLandmarks, // vector of the number of landmarks tested
                                    int nthreads){

        Eigen::MatrixXd rootShape2D = Utility::Shapes::generateUnitCirclePoints(numLandmarks);
        Eigen::MatrixXd rootShape3D = Utility::Shapes::generateUnitSpherePoints(numLandmarks);
        
        int ntrees = backboneTreeNewicks.size();
        //tree file helper
            auto fmt_path = [](int t, double alpha, double beta, int nt) {
                char buf[512];
                std::snprintf(buf, sizeof(buf), "TreeIndex%d_TreeAlpha%.2f_TreeBeta%.2f_NumAdditionalTaxa%d.tsv", t, alpha, beta, nt);
                return std::string(buf);
            };
        //shape file helper
        auto fmt_shape_path = [](const char* dim, int t, double alpha, double beta, int nt, int lm, double rate, double sig, int rep) {
            char buf[512];
            std::snprintf(buf, sizeof(buf), "%s_TreeIndex%d_TreeAlpha%.2f_TreeBeta%.2f_NumAdditionalTaxa%d_numLandmarks%d_lddmmAlpha%.2f_Sigma%.2f_rep%d",
                          dim, t, alpha, beta, nt, lm, rate, sig, rep);
            return std::string(buf);
        };
        
        #pragma omp parallel for num_threads(nthreads)
        for(int t = 0; t < ntrees; t++){
            //create a new tree object from the newick string
            Tree fixedTree  = Tree(backboneTreeNewicks[t]);
            
            //drop all non-pan great ape tips from these trees
            fixedTree.dropTip("Pongo_pygmaeus");
            fixedTree.dropTip("Papio");
            fixedTree.dropTip("Hylobates");
            fixedTree.dropTip("Colobus_guereza");
            fixedTree.dropTip("Gorilla_gorilla");
            fixedTree.forceBinary();
            fixedTree.initializeDownPassSequence();
            
            WriteTSV w = WriteTSV();
            std::string treeFile =treeOutfile + fmt_path(t, betaDistAlpha, betaDistBeta, numAdditionalTaxa);
            w.addFilepath( treeFile , true);
            
            SimShapeLDDMMv2 twoD = SimShapeLDDMMv2(rootShape2D);
            SimShapeLDDMMv2 threeD = SimShapeLDDMMv2(rootShape3D);
            
            for(int i = 0; i < nreps; i++){
                Tree tree = fixedTree;
                tree.addNTipsBetaDistTime(  numAdditionalTaxa,
                                            betaDistAlpha,
                                            betaDistBeta);
                w.appendDataTSV({tree.getNewickString()});
                
                //2D simulations
                twoD.runSimulation(&tree, lddmmAlpha, lddmmSigma);
                twoD.writeShapes(shapeOutfile + fmt_shape_path("2D", t, betaDistAlpha, betaDistBeta, numAdditionalTaxa, numLandmarks, lddmmAlpha, lddmmSigma, i), true);
                
                //3D simulations
                threeD.runSimulation(&tree, lddmmAlpha, lddmmSigma);
                threeD.writeShapes(shapeOutfile + fmt_shape_path("3D", t, betaDistAlpha, betaDistBeta, numAdditionalTaxa, numLandmarks, lddmmAlpha, lddmmSigma, i), true);
                }
                
            //compressing tree output
            std::string cmd = "gzip -f \"" + treeFile + "\"";
            int result = std::system(cmd.c_str());
            if (result != 0)
                std::cerr << "Warning: gzip compression failed for " << treeFile << std::endl;
        }
}

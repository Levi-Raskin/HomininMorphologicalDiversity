#include "Eigen/Dense"
#include "Msg.hpp"
#include "Node.hpp"
#include "SDE.hpp"
#include "SimShapeLDDMMv2.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "Utility.hpp"
#include "ThreadPool.hpp"
#include "TicToc.hpp"
#include "Tree.hpp"
#include "WriteTSV.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include </usr/local/include/omp.h>

SimShapeLDDMMv2::SimShapeLDDMMv2(const Eigen::MatrixXd& lm) : rootShape(lm), alpha(0.01), sigma(0.02), d((int)lm.cols()), n((int)lm.rows()){
    //instantiating kernal parmas alpha and sigma with reasonable values
    //j is the number of background grid items
    //assuming lm is nxd where d is the dimensionality of the landmarks and n is the # of landmarks
    //xmin/max and ymin/max are grid bounds
    efficientDiffMat = Eigen::MatrixXd::Zero(d * n, d * n);
    scratchShapeMatrix.resize(n, d);
    
    trueRoot.resize(n, d);
    currShape.resize(n, d);
    w.resize(n*d);
    incrementVec.resize(n*d);
    scratch.resize(n, d);
    scratch2.resize(n, d);
}

void SimShapeLDDMMv2::calculateDiffusionMatrixEfficient(const Eigen::MatrixXd& shape){
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {   // j starts at i, not 0
            double kernel_val = computeKernel(shape, i, j);
            for (int k = 0; k < d; k++) {
                efficientDiffMat(d*i + k, d*j + k) = kernel_val;
                efficientDiffMat(d*j + k, d*i + k) = kernel_val;
            }
        }
    }
}

double SimShapeLDDMMv2::computeKernel(const Eigen::MatrixXd& shape, int i, int j){
    double r_val = (1e-7 + (shape.row(i) - shape.row(j)).squaredNorm()) / sigma;
    return alpha * 2.0 * (3.0 + 3.0 * r_val + r_val * r_val) * std::exp(-r_val);
}

bool SimShapeLDDMMv2::doLinesIntersect(const Eigen::Vector2d &p1, const Eigen::Vector2d &p2,
                      const Eigen::Vector2d &q1, const Eigen::Vector2d &q2) {
    auto cross2d = [](const Eigen::Vector2d& a, const Eigen::Vector2d& b) {
        return a.x() * b.y() - a.y() * b.x();
    };

    Eigen::Vector2d r = p2 - p1;
    Eigen::Vector2d s = q2 - q1;
    Eigen::Vector2d qp = q1 - p1;

    double rxs    = cross2d(r, s);
    double qpxr   = cross2d(qp, r);
    double t_num  = cross2d(qp, s);

    if (std::abs(rxs) < 1e-10) {
        // Parallel: only intersect if collinear and overlapping
        if (std::abs(qpxr) > 1e-10) return false; // parallel but not collinear
        double rr = r.squaredNorm();
        double t0 = qp.dot(r) / rr;
        double t1 = t0 + s.dot(r) / rr;
        double lo = std::min(t0, t1);
        double hi = std::max(t0, t1);
        return hi >= 0.0 && lo <= 1.0;
    }

    double t = t_num  / rxs;
    double u = qpxr   / rxs;
    return t >= 0.0 && t <= 1.0 && u >= 0.0 && u <= 1.0;
}

bool SimShapeLDDMMv2::doSegmentsIntersect3D(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2,
                           const Eigen::Vector3d& q1, const Eigen::Vector3d& q2) {
    // Direction vectors
    Eigen::Vector3d u  = p2 - p1;
    Eigen::Vector3d v  = q2 - q1;
    Eigen::Vector3d w0 = p1 - q1;

    double a = u.dot(u);
    double b = u.dot(v);
    double c = v.dot(v);
    double d = u.dot(w0);
    double e = v.dot(w0);
    double D = a * c - b * b;

    double sc, tc;

    if (D < 1e-10) {
        sc = 0.0;
        tc = (b > c ? d / b : e / c);
    } else {
        sc = (b * e - c * d) / D;
        tc = (a * e - b * d) / D;
        // Reject immediately if parameters are out of segment bounds
        if (sc < 0.0 || sc > 1.0 || tc < 0.0 || tc > 1.0) return false;
    }

    Eigen::Vector3d closest = (p1 + sc * u) - (q1 + tc * v);
    return closest.squaredNorm() < 1e-12; // squaredNorm avoids sqrt
}

bool SimShapeLDDMMv2::checkLandmarksCrossed(const Eigen::MatrixXd& s) {
    int n = (int)s.rows();
    if (n < 4) return false;

    for (int i = 0; i < n; i++) {
        int i_next = (i + 1) % n;
        for (int j = i + 2; j < n; j++) {
            if (i == 0 && j == n - 1) continue;

            int j_next = (j + 1) % n;

            if (d == 2) {
                if (doLinesIntersect(s.row(i), s.row(i_next), s.row(j), s.row(j_next)))
                    return true;
            } else {
                if (doSegmentsIntersect3D(s.row(i), s.row(i_next), s.row(j), s.row(j_next)))
                    return true;
            }
        }
    }
    return false;
}


void SimShapeLDDMMv2::runSimulation(Tree* t, double a, double s){
    tree = t;
    dpseq = tree->getDownPassSequence();
    
    int maxNodeIdx = 0;
    for(Node* n : dpseq)
        if(n->getIndex() > maxNodeIdx)
            maxNodeIdx = n->getIndex();
    nodeShapes.clear();
    nodeShapes.resize(maxNodeIdx + 1);
    alpha = a;
    sigma = s;
    
    int gotocounter = 0;
    
    simissue:
    gotocounter++;
    if(gotocounter > 100)
        Msg::error("you're really struggling on this tree");
    
    for (auto i = dpseq.rbegin(); i != dpseq.rend(); i++){
        Node* n = *i;
        if(n == tree->getRoot())
            nodeShapes[n->getIndex()] = rootShape;
        else{
            Node* nAnc = n->getAncestor();
            scratchShapeMatrix = nodeShapes[nAnc->getIndex()];
            
            double bl = tree->getBranchLength(n, nAnc);
            if(bl == 0){
                nodeShapes[n->getIndex()] = scratchShapeMatrix;
            }else{
                try{
                    nodeShapes[n->getIndex()] = simulateShape(scratchShapeMatrix, bl);
                }catch (const std::runtime_error& e) {
                    std::cout << "Caught runtime_error: " << e.what() << std::endl;
                    goto simissue;
                }
            }
            
        }
    }
}

Eigen::MatrixXd SimShapeLDDMMv2::simulateShape(const Eigen::MatrixXd& s, double t){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    
    bool whileLoopExcept = false;
    int bigWhileCnt = 0;
    int numSteps = (int)(1000 * t);
    double timeStep = t / numSteps;
    
    trueRoot = s;
    
    do{
        currShape = trueRoot;
        whileLoopExcept = false;
        calculateDiffusionMatrixEfficient(currShape);
        
        for(int i = 0 ; i<numSteps; i++){
            int whileCnt = 0;
            do{
                for(int j = 0; j < d*n; j++)
                    w(j) = Probability::Normal::rv(&rng, 0, timeStep);
                    
                incrementVec.noalias() = efficientDiffMat * w;
                
                Eigen::Map<Eigen::MatrixXd> scratch_map(incrementVec.data(), n, d);
                scratch2 = currShape + scratch_map;
                
                whileCnt++;
                if(whileCnt > 500){
                    whileLoopExcept = true;
                    break;
                }
            }while(checkLandmarksCrossed(scratch2) == true);
            
            if(whileLoopExcept == true)
                break;
            currShape = scratch2;
            calculateDiffusionMatrixEfficient(currShape);
        }
        
        bigWhileCnt++;
        if(bigWhileCnt > 10)
            throw std::runtime_error("big while loop running forever");
        if(bigWhileCnt > 2)
            Msg::warning("struggling to simulate lddmm shapes without crossing: in big while loop");
    }while(whileLoopExcept == true);

    return currShape;
}

void SimShapeLDDMMv2::writeShapes(std::string fp, bool compress){
    std::vector<std::string> rn;
    
    int numLandmarks = (int)rootShape.rows();
    int numTips = tree->getNumTaxa();
    Eigen::MatrixXd plot(numLandmarks * numTips, d + 1);

    int rowIdx = 0;

    for (Node* n : dpseq) {
        if (n->getIsTip()) {
            scratchShapeMatrix = nodeShapes[n->getIndex()];
            for (int i = 0; i < numLandmarks; i++) {
                plot(rowIdx, 0) = i + 1;           // landmark index (1-based)
                for(int j = 0; j < d; j++){
                    plot(rowIdx, j+1) = scratchShapeMatrix(i, j);     // x
                }
                rowIdx++;
                rn.push_back(n->getName());
            }
        }
    }
    
    std::string out = fp + "nodeShapes.tsv";
    WriteTSV w = WriteTSV();
    w.addFilepath(out, true);
    w.addRownamesTSV(rn);
    w.appendDataTSV(plot);
    w.closeTSV();
    
    if(compress == true){
        std::string cmd = "gzip -f \"" + out + "\"";
        int result = std::system(cmd.c_str());
        if (result != 0)
            std::cerr << "Warning: gzip compression failed for " << out << std::endl;
    }
}

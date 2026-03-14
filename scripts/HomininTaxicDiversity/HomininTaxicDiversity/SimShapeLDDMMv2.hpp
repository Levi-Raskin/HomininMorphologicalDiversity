#ifndef SimShapeLDDMMv2_hpp
#define SimShapeLDDMMv2_hpp

#include "Eigen/Dense"
#include <map>
#include <string>

class Node;
class Tree;

class SimShapeLDDMMv2{
    public:
                                            SimShapeLDDMMv2(const Eigen::MatrixXd& lm);
        void                                runSimulation(Tree* t, double a, double s);
        void                                writeShapes(std::string fp, bool compress);
    private:
        void                                calculateDiffusionMatrixEfficient(const Eigen::MatrixXd& shape);
        bool                                checkLandmarksCrossed(const Eigen::MatrixXd& s);
        double                              computeKernel(const Eigen::MatrixXd& shape, int i, int j);
        bool                                doLinesIntersect(       const Eigen::Vector2d &p1,
                                                                    const Eigen::Vector2d &p2,
                                                                    const Eigen::Vector2d &q1,
                                                                    const Eigen::Vector2d &q2);
        bool                                doSegmentsIntersect3D(  const Eigen::Vector3d& p1,
                                                                    const Eigen::Vector3d& p2,
                                                                    const Eigen::Vector3d& q1,
                                                                    const Eigen::Vector3d& q2);
        Eigen::MatrixXd                     simulateShape(const Eigen::MatrixXd& s, double t);
        std::vector<Eigen::MatrixXd>        nodeShapes; //indexed by node idx
        std::vector<Node*>                  dpseq;
        Eigen::MatrixXd                     efficientDiffMat;
        Eigen::MatrixXd                     rootShape;
        Eigen::MatrixXd                     scratchShapeMatrix;
        //simulation scratch objs
        Eigen::MatrixXd                     trueRoot;
        Eigen::MatrixXd                     currShape;
        Eigen::VectorXd                     w;
        Eigen::VectorXd                     incrementVec;
        Eigen::MatrixXd                     scratch;
        Eigen::MatrixXd                     scratch2;
        Tree*                               tree;
        double                              alpha;
        double                              sigma;
        int                                 d;
        int                                 n;
};

#endif

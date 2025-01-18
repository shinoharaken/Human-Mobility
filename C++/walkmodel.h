#ifndef WALKMODEL_H
#define WALKMODEL_H

#define MAX_L_RATE 0.7
#define LAMBDA 1E20


#define isFix 1 //基点固定か否か
#define GAMMA 0.0
#define ALPHA 1.0


#include "utility.h"
using namespace std;


class WalkModel
{
public:
    WalkModel();

    WalkModel(int dimension=2,double range=1.0,double gamma=GAMMA);

    void initialize();
    void setRandomPosition();
    std::vector<double> calcNextPosition();

    std::vector<double> X();
    void setX(std::vector<double> x);
    double l();

    std::vector<double> targetPoint();
    std::vector<double> basePoint();


    double distanceFromBasepoint();
    double distanceFromTarget();




private:

    Utility *util;

    double m_lmax;

    int m_dimension;
    double m_range;
    double  m_l;

    double m_gamma;


    std::vector<double> m_deltaX;
    std::vector<double>  m_X;
    std::vector<double>  m_basepoint;
    std::vector<double> m_target;



    std::vector<double> calcDeltaX(std::vector<double> rr);
    std::vector< double> samplingInverseWalk();
    std::vector< double> generateTarget();
    std::vector<double> dataSamplingExponential();



};

#endif // WALKMODEL_H

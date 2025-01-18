#include "walkmodel.h"

WalkModel::WalkModel()
{
    util=new Utility;

    m_dimension=2;
    m_range=1.0;
    m_gamma=GAMMA;

    initialize();
}

WalkModel::WalkModel(int dimension,double range,double gamma)
{
    util=new Utility;

    m_dimension=dimension;
    m_range=range;
    m_gamma=gamma;


    initialize();
}

void WalkModel::initialize()
{

    m_lmax=MAX_L_RATE*m_range;
    m_X.clear();

    m_deltaX.clear();
    m_l=0.0;
    m_X.resize(m_dimension);
    m_deltaX.resize(m_dimension,1.0);
    m_basepoint.resize(m_dimension);
    m_target.resize(m_dimension);

    for(int d=0;d<m_dimension;d++)
    {

        m_basepoint[d]=0.0;
    }


    setRandomPosition();


}

void WalkModel::setX(std::vector<double> x)
{
    m_X=x;
}




std::vector<double>  WalkModel::X()
{
    return m_X;
}

std::vector<double>  WalkModel::targetPoint()
{
    return m_target;
}

std::vector<double>  WalkModel::basePoint()
{
    return m_basepoint;
}

double  WalkModel::l()
{
    return m_l;
}

double WalkModel::distanceFromTarget()
{
    return util->calcVectorNorm(m_target);
}

double WalkModel::distanceFromBasepoint()
{
    return util->calcVectorNorm(m_basepoint,this->X());
}

void WalkModel::setRandomPosition()
{


    double range=m_range;
    double r=range;
    while(r>=range)
    {
        for(int d=0;d<m_dimension;d++)
        {
            m_X[d]=-range+2.0*range*util->drnd();
        }
        r=util->calcVectorNorm(m_X);
    }


    for(int d=0;d<m_dimension;d++)
    {
        m_X[d]+=m_basepoint[d];
    }
}
std::vector< double> WalkModel::calcNextPosition()
{

    std::vector<double> newX(m_dimension,0.0);
    std::vector<double>  deltaxx(m_dimension,0.0);


    deltaxx=this->samplingInverseWalk();


    double distance=util->calcVectorNorm(deltaxx);
    if(distance>m_lmax)
    {
        double kk=m_lmax/distance;
        for(int d=0;d<m_dimension;d++)
            deltaxx[d]=kk*deltaxx[d];
    }
    else if(isnan(distance) || isinf(distance))
    {

        std::vector<double> deltaxx1(m_dimension,0.0);
        return deltaxx1;
    }


    for(int d=0;d<m_dimension;d++)
    {

        m_deltaX[d]=deltaxx[d];
        newX[d]=m_X[d]+m_deltaX[d];
    }


    for(int d=0;d<m_dimension;d++)
    {
        m_deltaX[d]=newX[d]-m_X[d];
    }
    m_l=util->calcVectorNorm(m_deltaX);

    return newX;
}

std::vector< double> WalkModel::samplingInverseWalk()
{

    std::vector<double>  deltaxx(m_dimension,0.0);
    m_target=generateTarget();


    std::vector< double> rr(m_dimension,0.0);
    for(int d=0;d<m_dimension;d++)
    {
        rr[d]=(m_target[d]-m_X[d]);
    }
    double al=ALPHA;


    for(int d=0;d<m_dimension;d++)
    {
        rr[d]=al*rr[d];
    }

    std::vector< std::vector< double> > aa;
    std::vector< std::vector< double> > aa1;
    std::vector< std::vector< double> > AA;
    std::vector< std::vector< double> > invAA;


    std::vector< std::vector< double> > baseVec;
    baseVec.resize(m_dimension);

    for(int d=0;d<m_dimension;d++)
    {
        std::vector< double> tmp(m_dimension,0.0);

        tmp=util->pointOnHypersurface(1.0,m_dimension);

        baseVec[d]=tmp;
    }
    baseVec=util->baseOrthogonalization(baseVec);//新たな直交基底作成

    AA.resize(m_dimension);
    for(int i=0;i<m_dimension;i++)
    {
        std::vector< double> tmp(m_dimension,0.0);
        AA[i]=tmp;
    }


    for(int i=0;i<m_dimension;i++)
    {
        for(int j=0;j<m_dimension;j++)
        {
            AA[i][j]=baseVec[j][i];
        }
    }

    double detAA;
    invAA=util->calcInverseMatrix(AA,1.e-6,&detAA);

    aa.resize(m_dimension);
    for(int d=0;d<m_dimension;d++)
    {
        aa[d].push_back(rr[d]);
    }



    aa1=util->MatrixMultiplication(invAA,aa);
    for(int d=0;d<m_dimension;d++)
    {
        rr[d]=aa1[d][0];
    }


    deltaxx=this->calcDeltaX(rr);


    for(int d=0;d<m_dimension;d++)
    {
        aa1[d][0]=deltaxx[d];
    }

    aa=util->MatrixMultiplication(AA,aa1);
    for(int d=0;d<m_dimension;d++)
    {
        deltaxx[d]=aa[d][0];
    }





    return deltaxx;
}
std::vector<double> WalkModel::calcDeltaX(std::vector<double> rr)
{

    std::vector<double> deltaxx(m_dimension,0.0);
    std::vector<double> eta(m_dimension,1.0/(double)m_dimension);
    std::vector<double> bb(m_dimension,1.0/(double)m_dimension);

    double r=util->calcVectorNorm(rr);


    if(r==0.0)
    {
        return deltaxx;
    }


    for(int d=0;d<m_dimension;d++)
    {
        bb[d]=util->drnd();
    }
    bb=util->normalize(bb);




    for(int d=0;d<m_dimension;d++)
    {
        eta[d]=pow(bb[d],1.0-fabs(m_gamma))*pow(rr[d]*rr[d]/(r*r),m_gamma);

    }




    for(int d=0;d<m_dimension;d++)
    {

        if(rr[d]!=0.0)
        {
            deltaxx[d]=eta[d]*r*r/rr[d];
        }
        else
            deltaxx[d]=0.0;
    }

    return deltaxx;
}



std::vector< double> WalkModel::generateTarget()
{
    std::vector<double> target;

    target=this->dataSamplingExponential();

    double distance=util->calcVectorNorm(target);
    while(distance>m_range)
    {
        target=this->dataSamplingExponential();
        distance=util->calcVectorNorm(target);
    }

    if(!isFix)
    {
        m_basepoint=m_X;
    }

    for(int d=0;d<m_dimension;d++)
    {
        target[d]+=m_basepoint[d];
    }

    return target;
}


vector<double> WalkModel::dataSamplingExponential()
{
    int dim=m_dimension;
    vector<double> data(dim,0.0);


    data=util->samplingExponential(LAMBDA,m_dimension);


    if(std::isinf(data[0]) || std::isinf(data[1]))
    {
        cout<<"dat error!"<<endl;
        cout<<"data="<<data[0]<<endl;

        vector<double> data1(dim,0.0);
        return data1;

    }

    return data;
}



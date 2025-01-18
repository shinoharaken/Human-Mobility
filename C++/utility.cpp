#include "utility.h"
using namespace std;

std::random_device seed_gen;
std::mt19937 engine(seed_gen());
std::uniform_real_distribution<> uniformdist(0.0, 1.0);
std::normal_distribution<> normaldist(0.0, 1.0);
Utility::Utility()
{
}

Utility::~Utility()
{
}


double Utility::drnd()
{
    return uniformdist(engine);
}
double Utility::innerProduct(std::vector< double> x,std::vector< double> y)
{
    double rt=0.0;
    if(x.size()!=y.size())
    {
        return rt;
    }

    for(int i=0;i<x.size();i++)
    {
        rt+=x[i]*y[i];
    }

    return rt;
}


double Utility::samplingExponential1D( double lambda )
{
    double y=this->drnd();
    double x=-log(1.0-y)/lambda;
    //    double x=-log(1.0-y)*lambda;

    return x;
}

std::vector< double> Utility::pointOnHypersurface(double r,int dimension)
{
    std::vector< double> data(dimension,0.0);


    std::vector< double> rad(dimension,0.0);
    for(int d=0;d<rad.size();d++)
    {
        rad[d]=this->drnd()*M_PI;
        if(d==rad.size()-2)
            rad[d]=2.0*this->drnd()*M_PI;
        if(d==rad.size()-1)
            rad[d]=0.0;
    }
    for(int d1=0;d1<dimension;d1++)
    {
        data[d1]=r*cos(rad[d1]);
        for(int d2=0;d2<d1;d2++)
        {
            data[d1]=data[d1]*sin(rad[d2]);
        }
    }

    return data;

}
std::vector< double> Utility::samplingExponential(double lambda,int dimension)
{
    std::vector< double> deltaxx(dimension,0.0);


    double l=this->samplingExponential1D(lambda);

    deltaxx=this->pointOnHypersurface(l,dimension);

    if(dimension==1)
    {
        deltaxx[0]=fabs(deltaxx[0]);
    }
    return deltaxx;
}

double Utility::calcVectorNorm(std::vector<double> data,double p)
{
    double sum=0.0;
    double rt=0.0;

    for(int i=0;i<data.size();i++)
    {
        sum+=std::pow(data.at(i),p);
    }

    rt= std::pow(sum,1.0/p);

    return rt;
}
double Utility::calcVectorNorm(std::vector<double> data0,std::vector<double> data1,double p)
{
    std::vector<double> data(data0.size(),0.0);
    for(int i=0;i<data.size();i++)
    {
        data[i]=data1[i]-data0[i];
    }

    return this->calcVectorNorm(data,p);
}
std::vector< double> Utility::vectorNormalization(std::vector< double> data)
{
    std::vector< double> rt(data.size(),0.0);
    double nm=this->calcVectorNorm(data);
    if(nm<=0.0)
        return rt;

    for(int i=0;i<data.size();i++)
    {
        rt[i]=data[i]/nm;

    }

    return rt;
}

std::vector< double> Utility::normalize(std::vector< double> data)
{
    std::vector< double> rt;
    rt.clear();
    double sum1=0.0;
    for(int i=0;i<data.size();i++)
    {
        sum1+=data[i];
    }

    for(int i=0;i<data.size();i++)
    {
        if(sum1!=0.0)
        {
            rt.push_back(data.at(i)/sum1);
        }
        else
        {
            rt.push_back(1.0/(double)data.size());

        }
    }

    return rt;
}


std::vector< std::vector<double> > Utility::calcInverseMatrix(std::vector< std::vector<double> > input, double eps, double *det)
{
    std::vector< std::vector<double> > output;
    output.clear();
    output.resize(input.size());

    if(input.size()==1)
    {
        output[0].resize(1);
        *det=fabs(input[0][0]);
        output[0][0]=1.0/input[0][0];
        return output;
    }
    int m=input.size();
    int l=input[0].size();
    int n=m*l;
    double *mat = (double *)malloc(n * sizeof(double));

    int ct=0;
    for(int i=0;i<input.size();i++)
    {
        output[i].resize(input[i].size());
        for(int j=0;j<input[i].size();j++)
        {
            mat[ct]=input[i][j];
            ct++;
        }

    }
    *det=this->minver(mat,l,m,eps);
    ct=0;
    for(int i=0;i<output.size();i++)
    {
        for(int j=0;j<output[i].size();j++)
        {
            output[i][j]=mat[ct];
            ct++;
        }

    }

    free(mat);

    return output;
}

void Utility::save2(std::vector< std::vector<double> > ary, std::string filename)
{
    std::string sep=",";
    std::ofstream file(filename);

    std::stringstream ss;

    for(int i=0;i<ary.size();i++)
    {
        for(int j=0;j<ary[i].size();j++)
        {
            if(j>0)
                ss<<sep<<ary[i][j];
            else
                ss<<ary[i][j];
        }
        ss<<std::endl;

    }

    file << ss.str();

    file.close();
}

void Utility::save2(std::vector< std::vector<int> > ary, std::string filename)
{
    std::string sep=",";
    std::ofstream file(filename);

    std::stringstream ss;

    for(int i=0;i<ary.size();i++)
    {
        for(int j=0;j<ary[i].size();j++)
        {
            if(j>0)
                ss<<sep<<ary[i][j];
            else
                ss<<ary[i][j];
        }
        ss<<std::endl;

    }

    file << ss.str();

    file.close();
}

void Utility::save1(std::vector<double> ary, std::string filename)
{

    std::ofstream file(filename);

    std::stringstream ss;

    for(int i=0;i<ary.size();i++)
    {
        ss<<ary[i];
        ss<<std::endl;

    }

    file << ss.str();

    file.close();
}

void Utility::save1(std::vector<int> ary, std::string filename)
{

    std::ofstream file(filename);

    std::stringstream ss;

    for(int i=0;i<ary.size();i++)
    {
        ss<<ary[i];
        ss<<std::endl;

    }

    file << ss.str();

    file.close();
}

std::vector< std::vector<double> > Utility::MatrixMultiplication(std::vector< std::vector<double> > xx,std::vector< std::vector<double> >yy)
{
    std::vector< std::vector<double> > output;
    int nn=yy.size();


    if(xx[0].size()!=yy.size())
    {
        std::cout<<"MatrixMultiplication size error!! "<<xx[0].size()<<" "<<yy.size()<<std::endl;
        return output;
    }

    output.resize(xx.size());
    for(int i=0;i<xx.size();i++)
    {
        output[i].resize(yy[0].size());
    }


    for(int i=0;i<output.size();i++)
    {
        for(int j=0;j<output[i].size();j++)
        {
            double sum=0.0;
            for(int k=0;k<nn;k++)
            {
                sum+=xx[i][k]*yy[k][j];
            }

            output[i][j]=sum;

        }
    }
    return output;
}

double Utility::minver(double a[], int l, int m, double eps)
{
    int i, iw, j, k, *p, r, s, t, u, v, *work;
    double api, pivot, *q, *q1, w, w1, wmax;

    if(m < 2 || l < m || eps <= 0.)
    {
        fprintf(stderr, "Error : Illegal parameter  in minver()\n");
        return 0.;
    }
    work = (int *)malloc(m * sizeof(int));
    if(work == NULL)
    {
        fprintf(stderr, "Error : Out of memory  in minver()\n");
        return 0.;
    }
    w1 = 1.;
    for(i = 0, p = work; i < m; i++)	*p++ = i;
    for(k = 0, u = 0; k < m; k++, u += l)
    {
        wmax = 0.;
        for(i = k; i < m; i++)
        {
            w = fabs(*(a + i * l + k));
            if(w > wmax)
            {
                wmax = w;
                r = i;
            }
        }
        api = fabs(pivot = *(a + r * l + k));
        if(api < eps)
        {
            //            fprintf(stderr, "Error : api < eps  in minver()\n");
            free((char *)work);
            return w1;
        }
        w1 *= pivot;
        v = r * l;
        if(r != k)
        {
            w1 = -w1;
            this->swapi(work + k, work + r);
            for(j = 0, q = a + u, q1 = a + v; j < m; j++)	this->swapd(q++, q1++);
        }
        for(i = 0, q = a + u; i < m; i++)	*q++ /= pivot;
        for(i = 0, v = 0; i < m; i++, v += l)
        {
            if(i != k)
            {
                s = v + k;
                w = *(a + s);
                if(w != 0.)
                {
                    for(j = 0, q = a + u, q1 = a + v; j < m; j++, q++, q1++)
                        if (j != k)	*q1 -= w * *q;
                    *(a + s) = - w / pivot;
                }
            }
        }
        *(a + u + k) = 1. / pivot;
    }
    for(i = 0; i < m; i++)
    {
        for(;;)
        {
            k = *(work + i);
            if(k == i)	break;
            this->swapi(work + k, work + i);
            for(j = 0, u = 0; j < m; j++, u += l)	this->swapd(a + u + i, a + u + k);
        }
    }
    free((char *)work);
    return w1;
}

void Utility::swapi(int *a, int *b)
{
    int w;

    w = *a;
    *a = *b;
    *b = w;
    return;
}

void Utility::swapd(double *a, double *b)
{
    double w;

    w = *a;
    *a = *b;
    *b = w;
    return;
}

std::vector< std::vector< double> > Utility::baseOrthogonalization(std::vector< std::vector< double> > baseVec)
{
    std::vector< std::vector< double> > uu;
    int dim=baseVec.size();
    uu.resize(dim);
    for(int k=0;k<dim;k++)
    {
        std::vector< double> tmp(dim,0.0);
        uu[k]=tmp;
    }


    uu[0]=this->vectorNormalization(baseVec[0]);
    for(int k=1;k<dim;k++)
    {
        std::vector< double> vv(dim,0.0);
        vv=baseVec[k];
        for(int i=0;i<k;i++)
        {
            double ip=this->innerProduct(baseVec[k],uu[i]);
            for(int d1=0;d1<dim;d1++)
            {
                vv[d1]-=ip*uu[i][d1];
            }
        }

        uu[k]=this->vectorNormalization(vv);

    }
    return uu;
}

vector< vector<double> > Utility::choleskyDecomp(vector< vector<double> > m)
{
    vector< vector<double>> l;
    l.clear();

    if(m.size() <= 0 || m.size() != m[0].size()) {
        return l;
    }

    int n=m.size();
    int i, j, k;

    for(i = 1; i < n; i++) {
        for(j = 0; j < i; j++) {
            if(m[i][j] != m[j][i]) {
                //                qDebug()<<"choleskyDecomp error";
                m[i][j] = m[j][i];
                //                return l;
            }
        }
    }

    double ld, lld;
    vector<double> d;
    d.resize(n);

    l.resize(n);

    for(i = 0; i < n; i++) {
        l[i].resize(i+1);
    }

    l[0][0] = 1;
    d[0] = m[0][0];

    for(i = 0; i < n; i++) {
        for(j = 0; j < i; j++) {
            lld = m[i][j];
            for(int k = 0; k < j; k++) {
                lld -= l[i][k] * l[j][k] * d[k];
            }
            l[i][j] = lld / d[j];
        }

        ld = m[i][i];
        for(k = 0; k < i; k++) {
            ld -= l[i][k] * l[j][k] * d[k];
        }
        l[i][j] = 1;
        d[i] = ld;
    }

    for(j = 0; j < n; j++) {
        d[j] = sqrt(d[j]);
    }

    for(i = 0; i < n; i++) {
        for(j = 0; j <= i; j++) {
            l[i][j] *= d[j];
        }
    }
    return l;
}



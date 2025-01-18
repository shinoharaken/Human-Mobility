#ifndef Utility_H
#define Utility_H


#include <string>
#include <sstream>
#include <time.h>
#include <iostream>
#include <fstream>
#include <random>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <dirent.h>


using namespace std;

class Utility

{
public:

    Utility();
    ~Utility();


    std::vector< double> pointOnHypersurface(double r,int dimension);



    void save1(std::vector<double> ary, std::string filename);
    void save1(std::vector<int> ary, std::string filename);
    void save2(std::vector< std::vector<double> > ary, std::string filename);
    void save2(std::vector< std::vector<int> > ary, std::string filename);
    double drnd();
    std::vector< std::vector< double> > baseOrthogonalization(std::vector< std::vector< double> >);

    std::vector< std::vector<double> > MatrixMultiplication(std::vector< std::vector<double> > xx,std::vector< std::vector<double> >yy);
    std::vector< std::vector<double> > calcInverseMatrix(std::vector< std::vector<double> > input, double eps, double *det);
    std::vector< double> normalize(std::vector< double> data);

    double calcVectorNorm(std::vector<double> data,double p=2.0);
    double calcVectorNorm(std::vector<double> data0,std::vector<double> data1,double p=2.0);
    double samplingExponential1D( double lambda );
    std::vector< double> samplingExponential(double lambda,int dimension);
    std::vector< double> vectorNormalization(std::vector< double> data);
    double innerProduct(std::vector< double> x,std::vector< double> y);


private:
    double minver(double a[], int l, int m, double eps);
    void swapi(int *a, int *b);
    void swapd(double *a, double *b);

    vector< vector<double> > choleskyDecomp(vector< vector<double> > m);
};

#endif // Utility_H

//#include <QCoreApplication>


#include "walkmodel.h"


#define RANGE 1.0

#define TMAX 100000


#define DIMENSION 2


#define TRYMAX 1000000

using namespace std;
Utility *util;


int main(int argc, char *argv[])
{

    util=new Utility;
    for(int tr=0;tr<TRYMAX;tr++)
    {


        WalkModel *agent=new WalkModel(DIMENSION,RANGE,GAMMA);
        agent->initialize();
        for(int t=0;t<TMAX;t++)
        {

            vector<double> cpos=agent->calcNextPosition();
            if(util->calcVectorNorm(cpos)<=RANGE)
            {
                agent->setX(cpos);
            }
            else
            {
                break;
            }

        }//t

    }//tr


    return 1;


}







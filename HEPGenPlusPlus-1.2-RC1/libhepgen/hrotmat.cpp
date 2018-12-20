#include "hrotmat.h"


void HRotMat::resetMatrix()
{
    for (unsigned int i =0; i < 3; i++)
        for (unsigned int j =0; j < 3; j++)
        {
            if (i==j)
                components[i][j] = 1;
            else
                components[i][j] = 0;
        }

}



void HRotMat::setFromThetaPhiVector(double _theta, double _phi, const HVector3& _axis)
{
    //make a copy, we want to change it up!
    HVector3 axis = _axis;
    //normalize the axis
    axis.normalize(1.0);

    HVector3 helperx = HVector3(cos(_theta),sin(_theta)*cos(_phi),sin(_theta)*sin(_phi));
    HVector3 helperz = helperx.crossProduct(axis);
    HVector3 helpery = helperz.crossProduct(helperx);


    //normalize them
    helperx.normalize(1.0);
    helpery.normalize(1.0);
    helperz.normalize(1.0);


    //components[COLUMN][ROW]
    for (unsigned int i = 0; i < 3; i++)
    {
        components[0][i] = helperx.at(i);
        components[1][i] = helpery.at(i);
        components[2][i] = helperz.at(i);
    }
}


void HRotMat::printDebug()
{
    printf("-------------------------------------------------------\n");
    for (int i =0; i<3; i++)
        printf("%f \t %f \t %f \n",components[i][0],components[i][1],components[i][2]);
    printf("-------------------------------------------------------\n");
}





HVector3 HRotMat::rotateVector(HVector3 _in)
{
    double ret[3];
    double inc[3];
    inc[0]=_in.X();
    inc[1]=_in.Y();
    inc[2]=_in.Z();
    memset(&ret[0],0,sizeof(double)*3);
    //keep columns fixed, iterate overs rows "row times column"
    for (unsigned int i=0; i<3; i++)
    {
        for (unsigned int j=0; j<3; j++)
            ret[i] += inc[j]*components[j][i];
    }
    return HVector3(ret[0],ret[1],ret[2]);
}

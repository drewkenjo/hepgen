
#include <iostream>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <fstream>

#include <string>
#include <sys/stat.h>
#include <unistd.h>

int main (int argc, char** argv){
    if (argc < 2){
        printf("Use me like: %s outFile \n",argv[0]);
        return -1;
    }
    
    for (double qsq = 2.0; qsq < 16; qsq++){
        std::cout << qsq << std::endl;
        for (double xbjC = 0.005; xbjC < 0.35; xbjC+=0.01){
//             for(double tprimCount = -0.1; tprimCount > -0.75; tprimCount-= 0.05){
                double tprimCount = -0.30;
                char buffer[200];
                sprintf(buffer,"./libGKPi0/gkCrossSection %f %f %f | grep Final | cut -d':' -f2 >> %s",qsq,xbjC,tprimCount,argv[1]);
                 std::cout << "running: " << buffer << std::endl;
                system(buffer);
            }
        }
    
    
}

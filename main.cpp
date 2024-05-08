#include "pointGroup.h"
#include <iostream>


// (1) Based on algorithms:
//   Marcus Johansson and Valera Veryazov,
//   J. Cheminform. 2017,9,8. (doi:10.1186/s13321-017-0193-3)
// (2) compile: qmake && make
// (3) more information please see: https://symotter.org/gallery

int main(int argc, char *argv[])
{
    std::string fileName="/home/zhangfq/project/new/symmetrizer/example/benzene.xyz";
    PointGroup p(0.01);
    p.loadFile(fileName);
    return 1;
}

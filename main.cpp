#include <iostream>
#include <string>
#include <vector>
#include <eigen3/Eigen/Eigen>
#include "EllipSFA.hpp"

using cMat = std::vector<std::vector<std::complex<double>>>;

int main()
{
    cMat result(2, cVec(2));
    result[0][0] = 1;
    result[1][0] = 3;
    result[0][1] = 2;
    result[1][1] = 3;

    EllipticSFA test = EllipticSFA("test");
    return 0;
}
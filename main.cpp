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

    EllipticSFA test = EllipticSFA();
    std::string fileName = "testfil.txt";

    //test.loadInputFile();

    double omega = 0.057;
    double Ip = 0.579;
    double Up = 0.55;
    double cep = M_PI/2.0;
    double epsilon = M_PI/4.0;
    int N = 4;

    EllipticSFA sfa = EllipticSFA(Ip, Up, N, cep, omega, epsilon);
    dVec p = {0.3, 0.2, 0.0};
    cVec ts = sfa.getSaddleTimes(p);
    for(dcmplx t : ts)
    {
        std::cout << t << "\n";
    }
    return 0;
}
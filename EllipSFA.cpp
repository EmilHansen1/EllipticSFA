#include <iostream>
#include <eigen3/Eigen/Eigen>
#include "EllipSFA.hpp"

EllipticSFA::EllipticSFA(){}

EllipticSFA::EllipticSFA(double Ip, double Up, int N, double phi, double omega, double epsilon)
{

}


cVec EllipticSFA::eField(dcmplx t)
{
    double factor = 2.0 * Up * omega / std::sqrt(1.0 + std::pow(epsilon, 2));
    cVec result = {factor * std::sin(omega  * t), -factor * epsilon * std::cos(omega * t)};
    return result;
}

cVec EllipticSFA::aField(dcmplx t)
{
    double factor = 2.0 * Up / std::sqrt(1.0 + std::pow(epsilon, 2));
    cVec result = {factor * std::cos(omega * t), factor * epsilon * std::sin(omega * t)};
    return result;
}


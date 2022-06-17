/**
 * SFA direct ionisation solver for elliptically polarised light
 * @author Emil Hansen
 * @date 2022
*/

#ifndef ELLIPTIC_SFA_H
#define ELLIPTIC_SFA_H

#include <complex>

/* Shorthands */
using dcmplx = std::complex<double>;
using dVec = std::vector<double>;
using cVec = std::vector<std::complex<double>>;

/**
 * Where the magic happens
 */
class EllipticSFA
{
public:

    /* Field variables */
    double Ip, Up, phi, omega, epsilon;
    int N;

    /**
     * Empty CTOR
     */
    EllipticSFA();

    /**
     * Ze good CTOR
     * @param Ip
     * @param Up
     * @param N
     * @param phi
     * @param omega
     * @param epsilon
     */
    EllipticSFA(double Ip, double Up, int N, double phi, double omega, double epsilon);

    dcmplx envelope(dcmplx t)
    {
        return std::pow(std::sin(omega * t / (2.0 * N)), 2);
    }

    /**
     * The electric field
     * @param t complex time
     * @return the value of the electric field at complex time t
     */
    cVec eField(dcmplx t);

    /**
     * The vector potential
     * @param t complex time
     * @return the value of the vector potential at complex time t
     */
    cVec aField(dcmplx t);


private:

};



#endif


/**
 * SFA direct ionisation solver for elliptically polarised light
 * @author Emil Hansen
 * @date 2022
*/

#ifndef ELLIPTIC_SFA_H
#define ELLIPTIC_SFA_H

#include <complex>
using dcmplx = std::complex<double>;

/**
 * Where the magic happens
 */
class EllipticSFA
{
public:

    /**
     * Empty CTOR
     */
    EllipticSFA();


    EllipticSFA(double Ip, double Up, int N, double phi, double omega, std::string target, );


private:

};



#endif


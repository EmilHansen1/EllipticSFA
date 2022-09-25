/**
 * SFA direct ionisation solver for elliptically polarised light
 * @authors Emil Hansen & Mads Carlsen
 * @date 2022
*/

#ifndef ELLIPTIC_SFA_H
#define ELLIPTIC_SFA_H

#include <complex>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <libconfig.h++>

/* Shorthands and constants */
using dcmplx = std::complex<double>;
using dVec = std::vector<double>;
using cVec = std::vector<std::complex<double>>;
using dMat = std::vector<std::vector<double>>;
using cMat = std::vector<std::vector<std::complex<double>>>;
constexpr dcmplx I(0., 1.);

/**
 * Where the magic happens
 */
class EllipticSFA
{
public:

    /* Field variables */
    double Ip, Up, cep, omega, epsilon;
    int N;
    std::string target;

    /**
     * Empty CTOR
     */
    EllipticSFA();

    /**
     * Ze good CTOR
     * @param Ip
     * @param Up
     * @param N
     * @param cep carrier envelope phase of the field
     * @param omega frequency of field
     * @param epsilon elliptic parameter
     */
    EllipticSFA(double Ip, double Up, int N, double cep, double omega, double epsilon);

    /**
     * Ze real good CTOR
     * @param inputFileName file name of the input file
     */
    EllipticSFA(std::string inputFileName);

    /**
     * The electric field, minus the derivative of aField
     * @param t complex time
     * @return the value of the electric field at complex time t
     */
    cVec eField(dcmplx t);

    /**
     * The vector potential for elliptical field
     * @param t complex time
     * @return the value of the vector potential at complex time t
     */
    cVec aField(dcmplx t);

    /**
     *
     * @param t complex time
     * @param p momentum vector (px, py, pz)
     * @return the value of the action S at a given time and momentum
     */
    dcmplx action(dcmplx t, dVec pVec);

    /**
     * Derivative of the action with respect to time
     * @param t complex time
     * @param pVec momentum vector (px, py, pz)
     */
    dcmplx actionDt(dcmplx t, dVec pVec);

    /**
     * Double derivative of the ation with respect to time
     * @param t complex time
     * @param pVec momentum vector (px, py, pz)
     */
    dcmplx actionDDt(dcmplx t, dVec pVec);

    /**
     * Function for finding the saddle point times for a given momentum
     * @param pVec momentum vector (px, py, pz)
     * @return vector of complex<double> saddle point times
     */
    cVec getSaddleTimes(dVec pVec);

    /**
     * The transition amplitude for a given final momentum pVec
     * @param pVec momentum vector (px, py, pz)
     */
    dcmplx transAmp(dVec pVec);

    /**
     * The initial state matrix element found in the SFA transition amplitude
     * @param pVec momentum vector (px, py, pz)
     * @param ts complex saddle point time
     * @return
     */
    dcmplx getMatrixElement(dVec pVec, dcmplx ts);

    /**
     * Transition amplitude calculated for a XY momentum grid spanned by pxList and pyList
     * @param pxList
     * @param pyList
     */
    dMat transAmpXY(dVec pxList, dVec pyList, double pz);

    /**
     * Function that saves a given matrix to a file
     * @param fileName name of data file
     * @param mat matrix to save
     */
    void saveMatrixToFile(std::string fileName, dMat& mat);

    /**
     * Template function used to load the variables using libconfig
     * @tparam T type of the variable to be loaded
     * @param paramName name of the parameter in the input file
     * @param cfg input file object from libconfig
     * @param defaultValue default value of the parameter, if not found in the input file
     * @return value of the parameter specified in the input file
     */
    template<typename T> T loadParam(std::string paramName, T defaultValue, libconfig::Config& cfg);

    /**
     * Loads all the class variables and settings from the specified input file
     * @param fileName file name of input file
     */
    void loadInputFile(std::string fileName);

private:

};



#endif


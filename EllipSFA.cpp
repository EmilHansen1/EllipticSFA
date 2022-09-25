#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Eigen>
#include "EllipSFA.hpp"
#include <eigen3/unsupported/Eigen/Polynomials>

EllipticSFA::EllipticSFA(){}

EllipticSFA::EllipticSFA(double Ip, double Up, int N, double phi, double omega, double epsilon)
{
    this->Ip = Ip;
    this->Up = Up;
    this->N = N;
    this->cep = phi;
    this->omega = omega;
    this->epsilon = epsilon;
}

EllipticSFA::EllipticSFA(std::string inputFileName)
{
    loadInputFile(inputFileName);
}

cVec EllipticSFA::eField(dcmplx t)
{
    dcmplx factor = -1. * std::sqrt(2.*Up) * std::sin(omega*t/(2.*N)) * omega / (double)N;
    dcmplx cosCep = std::cos(omega*t+cep);
    dcmplx sinCep = std::cos(omega*t+cep);
    dcmplx cosN = std::cos(omega*t/(2.*N));
    dcmplx sinN = std::sin(omega*t/(2.*N));
    cVec result = {factor * std::cos(epsilon/2.) * (cosCep * cosN - (double)N * sinN * sinCep),
                   factor * std::sin(epsilon/2.) * ((double)N * sinN * cosCep + sinCep * cosN)};
    return result;
}

cVec EllipticSFA::aField(dcmplx t)
{
    dcmplx factor = std::sqrt(2.*Up) * std::pow(std::sin(omega*t / (2.*N)), 2);
    cVec result = {factor * std::cos(omega * t + cep) * std::cos(epsilon/2.),
                   factor * std::sin(omega * t + cep) * std::sin(epsilon/2.)};
    return result;
}

dcmplx EllipticSFA::action(dcmplx t, dVec pVec)
{
    double px = pVec[0]; double py = pVec[1]; double pz = pVec[2];
    double cp = std::cos(cep);
    double sp = std::sin(cep);
    dcmplx co = std::cos(omega * t);
    dcmplx so = std::sin(omega * t);
    double ce = std::cos(epsilon);
    double ce2 = std::cos(epsilon/2.);
    double se2 = std::sin(epsilon/2.);
    dcmplx con = std::cos(omega * t / (double)N);
    dcmplx son = std::sin(omega * t / (double)N);
    double A0 = std::sqrt(2.*Up);

    double frontFac = 1. / (64. * std::pow(N,4) * omega - 80. * std::pow(N,2)*omega + 16. * omega);
    dcmplx term1 = 16. * (N-0.5) * std::pow(A0,2) * std::pow(N,2) * (cp * sp * std::pow(co,2) + so *
            (std::pow(cp,2) - 0.5) * co - sp * cp/2) * (N+0.5) * ce * std::pow(con,2);

    dcmplx term2 = 4. * A0 * N * con * (((-4. * std::pow(cp,2) * ce + 2.*ce) * std::pow(co,2) + 4.*cp*sp*ce*co*so +
            2.*std::pow(cp,2)*ce + std::pow(N,2) - ce - 1.) * (N-0.5) * A0 * (N+0.5) * son - 8.*(2.*(N-0.5) * px *
                    (so*cp + co*sp) * (N+0.5) * ce2 - 2.*(N-0.5) * (co*cp - so*sp) * py * (N+0.5) * se2 + A0*(N-1.) *
                    (cp*sp * std::pow(co,2) + so*(std::pow(cp,2) - 0.5)*co - sp*cp/2.) * (N+1.) * ce) * (double)N);

    dcmplx term3 = -16.*A0*N*son * (-4.*(N-0.5) * (co*cp - so*sp) * px * (N+0.5)*ce2 - 4.*(N-0.5)*py*
            (so*cp + co*sp) * (N+0.5)*se2 + A0*(N-1.)*((-std::pow(cp,2)*ce + ce/2) * std::pow(co,2) + cp*sp*ce*co*so +
            std::pow(cp,2) * ce/2. + std::pow(N,2) - ce/4. - 1./4.)*(N+1.));

    dcmplx term4 = 64. * (N-0.5) * (sp * (std::pow(N,2)-1.)*co + cp*(std::pow(N,2)-1.)*so + sp) * A0 * px * (N+0.5)*ce2;
    dcmplx term5 = -64.*(N-0.5)*A0*(cp*(std::pow(N,2)-1.)*co + (1.-std::pow(N,2))*sp*so + cp) * py * (N+0.5) * se2;
    dcmplx term6 = 16. * (N-0.5) * std::pow(A0,2) *(std::pow(N,2) - 3./2.) * sp * cp * (N+0.5) * ce * std::pow(co,2);
    dcmplx term7 = 16. * (N-0.5) * std::pow(A0,2) * (std::pow(N,2)-3./2.) * so * (std::pow(cp,2) - 0.5) * (N+0.5) * ce * co;
    dcmplx term8 = -8.*(std::pow(A0,2) * ce * sp * (std::pow(N,2)-3./4.) * cp - (3.*(N-0.5)*(std::pow(A0,2) +
            16.*std::pow(px,2)/3. + 16.*std::pow(py,2)/3. + 16.*std::pow(pz,2)/3.)*omega*(N+0.5)*t)/2.) * (N-1.)*(N+1.);
    return Ip * t + 0.5 * frontFac * (term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8);
}

dcmplx EllipticSFA::actionDt(dcmplx t, dVec pVec)
{
    double px = pVec[0]; double py = pVec[1]; double pz = pVec[2];
    dcmplx sinFac = std::sqrt(2.*Up) * std::pow(std::sin(omega*t/(2.*N)),2);
    return Ip + 0.5 * (std::pow(pz,2) + std::pow(px + sinFac * std::cos(omega*t+cep) * std::cos(epsilon/2.),2) +
        std::pow(py + sinFac * std::sin(omega*t+cep) * std::sin(epsilon/2.),2));
}

dcmplx EllipticSFA::actionDDt(dcmplx t, dVec pVec)
{
    double px = pVec[0]; double py = pVec[1];
    dcmplx sinN = std::sin(omega*t/(2.*N));
    dcmplx cosN = std::cos(omega*t/(2.*N));
    dcmplx sinN2 = std::pow(sinN,2);
    dcmplx cosCep = std::cos(omega*t + cep);
    dcmplx sinCep = std::sin(omega*t + cep);
    double cosE = std::cos(epsilon/2.);
    double sinE = std::sin(epsilon/2.);
    double A0 = std::sqrt(2. * Up);
    dcmplx factor1 = 2.*(px + A0*sinN2 * cosCep * cosE);
    dcmplx factor2 = A0 * sinN * cosCep * cosE * omega * cosN/(double)N - A0 * sinN2 * omega * sinCep * cosE;
    dcmplx factor3 = 2.*(py + A0 * sinN2 * sinCep * sinE);
    dcmplx factor4 = A0 * sinN * sinCep * sinE * omega * cosN/(double)N + A0 * sinN2 * omega * cosCep * sinE;
    return factor1 * factor2 + factor3 * factor4;
}

cVec EllipticSFA::getSaddleTimes(dVec pVec)
{
    // Abbrevations for (slightly) more readability
    dcmplx ec = std::exp(I*cep);
    dcmplx ec2 = ec*ec;
    double px = pVec[0], py = pVec[1], pz = pVec[2];
    double A = std::sqrt(2. * Up);
    double C = std::cos(epsilon/2.0);
    double S = std::sin(epsilon/2.0);
    double C2 = C*C;
    double S2 = S*S;

    // Create the polynomial, of order 4*N + 4, and fill in the coefficients
    Eigen::VectorXcd coefficients(4*N + 5);
    coefficients.setZero(); // SKAL DEN VÆRE HER?

    coefficients(0) = (C2 - S2)*Up*ec2/64.0;
    coefficients(1) = (S2 - C2)*Up*ec2/16.0;
    coefficients(2) = 3.0*(C2 - S2)*Up*ec2/32.0;
    coefficients(3) = (S2 - C2)*Up*ec2/16.0;
    coefficients(4) = (C2 - S2)*Up*ec2/64.0;
    coefficients(N + 1) += (-C*px + I*S*py)*A*ec/8.0;
    coefficients(N + 2) += (C*px - I*S*py)*A*ec/4.0;
    coefficients(N + 3) += (-C*px + I*S*py)*A*ec/8.0;
    coefficients(2*N + 0) += (C2 + S2)*Up/32.0;
    coefficients(2*N + 1) += -(C2 + S2)*Up/8.0;
    coefficients(2*N + 2) += (px*px + py*py + pz*pz)/2.0 + Ip + 3.0*(C2 + S2)*Up/16.0;
    coefficients(2*N + 3) += -(C2 + S2)*Up/8.0;
    coefficients(2*N + 4) += (C2 + S2)*Up/32.0;
    coefficients(3*N + 1) += -(C*px + I*S*py)*A/(8.0*ec);
    coefficients(3*N + 2) += (C*px + I*S*py)*A/(4.0*ec);
    coefficients(3*N + 3) += -(C*px + I*S*py)*A/(8.0*ec);
    coefficients(4*N + 0) += (C2 - S2)*Up/(64.0*ec2);
    coefficients(4*N + 1) += (S2 - C2)*Up/(16.0*ec2);
    coefficients(4*N + 2) += 3.0*(C2 - S2)*Up/(32.0*ec2);
    coefficients(4*N + 3) += (S2 - C2)*Up/(16.0*ec2);
    coefficients(4*N + 4) += (C2 - S2)*Up/(64.0*ec2);

    // Reverse the coefficients and solve the polynomial
    Eigen::PolynomialSolver<dcmplx, Eigen::Dynamic> solver;
    solver.compute(coefficients);
    const Eigen::PolynomialSolver<dcmplx, Eigen::Dynamic>::RootsType &roots = solver.roots();

    // Find corresponding times
    cVec ts;
    for(dcmplx r : roots)
    {
        dcmplx t = I*((double) N)*std::log(r)/omega;
        if(t.imag() < 0) continue;
        if(t.real() < 0) t += 2.0*N*M_PI/omega;
        ts.push_back(t);
    }

    // Finally sort according to real part
    std::sort(ts.begin(), ts.end(), [&](dcmplx a, dcmplx b) -> bool { return a.real() < b.real(); });
    return ts;
}

dcmplx EllipticSFA::transAmp(dVec pVec)
{
    cVec saddleTimes = getSaddleTimes(pVec);
    dcmplx matrixElement = 0.;
    dcmplx result = 0.;
    for (dcmplx ts : saddleTimes){
        matrixElement = getMatrixElement(pVec, ts);
        result += std::sqrt(2. * M_PI * I / actionDDt(ts, pVec)) * std::exp(I * action(ts, pVec)) * matrixElement;
    }
    return result;
}

dMat EllipticSFA::transAmpXY(dVec pxList, dVec pyList, double pz)
{
    int sizeX = pxList.size(), sizeY = pyList.size();
    dVec pVec;
    // First initialize the result matrix with the correct size
    dMat result(sizeX, dVec(sizeY));

    // Then fill it with the transition amplitudes
    for (int i = 0; i < sizeX; ++i)
    {
        for (int j = 0; j < sizeY; ++j)
        {
            pVec = {pxList[i], pyList[j], pz};
            result[i][j] = (double) std::pow(std::abs(transAmp(pVec)), 2);
        }
    }
    return result; 
}

dcmplx EllipticSFA::getMatrixElement(dVec pVec, dcmplx ts)
{
    return 1.0;
}

void EllipticSFA::saveMatrixToFile(std::string fileName, dMat& mat)
{
    std::ofstream ost{fileName};
    for (dVec row : mat)
    {
        for (int i = 0; i < row.size(); ++i)
        {
            if(i != 0)
            {
                ost << " ";
            }
            ost << row[i];
        }
        ost << "\n";
    }
    ost.close();
}

template<typename T>
T EllipticSFA::loadParam(std::string paramName, T defaultValue, libconfig::Config& cfg)
{
    T param;
    try
    {
        param = cfg.lookup(paramName.c_str());
        std::cout << paramName + " : " << param << "\n";
    }
    catch(const libconfig::SettingNotFoundException &nfex)
    {
        param = defaultValue;
        std::cout << "No setting for " + paramName + " found in configuration file. Using standard value of " << param << "." << std::endl;
    }
    return param;
}

void EllipticSFA::loadInputFile(const std::string fileName)
{
    double lambda, I0;
    libconfig::Config cfg;

    // Open the file and make the Config object
    try{
        cfg.readFile(fileName.c_str()); //(char*) &fileName
    }
    catch(const libconfig::FileIOException &fioex)
    {
        std::cerr << "I/O error while reading file." << std::endl;
    }

    // Load the parameters or set their default value
    lambda = loadParam("lambda", 800., cfg);
    I0 = loadParam("I0", 1.e14, cfg);
    cep = loadParam("CEP", 0., cfg);
    epsilon = loadParam("epsilon", 0., cfg);
    N = loadParam("N", 2, cfg);
    Ip = loadParam("Ip", 0.5, cfg);
    target = loadParam("target", "default", cfg);

    // Calculate omega and Up in a.u.  (We SHOULD add more significant digits in this conversion)
    omega = 2.*M_PI*137.036 / (lambda * 1.e-9 / (5.29177e-11));
    double E_max = std::sqrt(I0 / (3.50945e16));  // Max electric field amplitude in a.u.
    Up = std::pow(E_max,2) / (4*std::pow(omega,2));
}
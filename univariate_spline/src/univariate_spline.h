//
// Created by Aleksandr on 04-Sep-18.
//

#ifndef UNIVARIATE_SPLINE_H_
#define UNIVARIATE_SPLINE_H_

#include <vector>

class UnivariateSpline {

public:

    UnivariateSpline(const unsigned int& degree, const unsigned int& knotsCount);

    void initializeUniformKnots(const std::vector<double> &x);
    void computingCoefficients(const std::vector<double> &x, const std::vector<double> &y,
                               const std::vector<double> &w);

    unsigned int getKnotsCount() const;
    unsigned int getDegree() const;

    double getValue(const double& x) const;
    double getAlpha(const double &x, const unsigned int &knotId, const unsigned int &deg) const;
    std::vector<double> getKnots() const;


private:

    std::vector<double> knots;
    std::vector<double> coefficients;

    unsigned int degree;
    unsigned int knotsCount;
    unsigned int internalKnotsCount;

    void initializeBSplines(const double& x, std::vector<double>& bSplines);
    void computingMatrixA(std::vector<std::vector<double>>& A, const unsigned int& sizeMatrix,
                                  const std::vector<double>& x, const std::vector<double>& y,
                                  const std::vector<double>& w);

    int getLeftKnotIndex(const double &x, const int &minId) const;
    double deBoorAlgorithm(const double& x) const;
};


#endif //UNIVARIATE_SPLINE_H_

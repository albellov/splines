//
// Created by Aleksandr on 04-Sep-18.
//

#ifndef SPLINE_2D_SPLINE_H
#define SPLINE_2D_SPLINE_H

#include <vector>

class Spline {

public:

    Spline(const unsigned int& degree, const unsigned int& knotsCount);

    void initializeUniformKnots(const std::vector<double> &x);
    void computingCoefficients(const std::vector<double> &x, const std::vector<double> &y,
                               const std::vector<double> &w);

    unsigned int getKnotsCount() const;
    unsigned int getDegree() const;

    double getValue(const double& x) const;
    double getAlpha(const double &x, const unsigned int &knotId, const unsigned int &deg) const;


private:

    std::vector<double> knots;
    std::vector<double> coefficients;

    unsigned int degree;
    unsigned int knotsCount;
    unsigned int internalKnotsCount;

    double deBoorAlgorithm(const double& x) const;
    void initializeBSplines(const double& x, std::vector<double>& bSplines);
    void computingMatrixA(std::vector<std::vector<double>>& A, const unsigned int& sizeMatrix,
                                  const std::vector<double>& x, const std::vector<double>& y,
                                  const std::vector<double>& w);

    int getLeftKnotIndex(const double &x, const int &minId) const;
};


#endif //SPLINE_2D_SPLINE_H

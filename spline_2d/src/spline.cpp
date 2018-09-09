//
// Created by Aleksandr on 04-Sep-18.
//

#include <iostream>
#include <cmath>

#include "spline.h"
#include "../tools/tools.h"


// public:

Spline::Spline(const unsigned int& degree, const unsigned int& knotsCount){
    this->degree = degree;
    this->knotsCount = knotsCount;
    this->internalKnotsCount = knotsCount - 2;
}


void Spline::initializeUniformKnots(const std::vector<double> &x) {
    const auto points_count = static_cast<const unsigned int>(x.size());
    unsigned int unique_points = 0;
    const unsigned int oldKnotsCount = knotsCount;

    knots.resize(oldKnotsCount + 2 * degree);

    knotsCount = static_cast<unsigned int>(knots.size());

    for(int i = 0; i < degree + 1; i++){
        knots[i] = x[0];
        knots[knotsCount - 1 - i] = x[points_count-1];
    }

    for(int i = 1; i < points_count - 1; i++){
        if (x[i] != x[i-1]){
            unique_points ++;
        }
    }

    // Number of knots should be less than n - k for n points with unique x.
    if(unique_points == 0 || unique_points - degree < oldKnotsCount){
        return;
    }

    // Put the knots between the desired points.
    double points_for_knot = (double)unique_points / (oldKnotsCount - 1);

    unsigned int counterDiffPoints = 0;
    unsigned int ind = 1;
    bool isCoincidingPoints = x[0] == x[1];

    for (int knotId = 1; knotId < oldKnotsCount - 1; knotId++){
        for (; counterDiffPoints < points_for_knot * knotId || isCoincidingPoints; ind++){
            isCoincidingPoints = x[ind] == x[ind-1];

            if (x[ind] != x[ind-1]){
                counterDiffPoints++;
            }
        }
        knots[knotId + degree] = 0.5 * (x[ind] + x[ind-1]);
    }
}


void Spline::computingCoefficients(const std::vector<double> &x, const std::vector<double> &y,
                                   const std::vector<double> &w) {

    const unsigned int sizeMatrix = internalKnotsCount + degree + 1;

    std::vector<std::vector<double>> A(sizeMatrix, std::vector<double>(sizeMatrix, 0.0));
    std::vector<std::vector<double>> L(sizeMatrix, std::vector<double>(sizeMatrix, 0.0));

    std::vector<double> bSplines(degree + 1, 0.0);
    coefficients.resize(sizeMatrix);

    computingMatrixA(A, sizeMatrix, x, y, w);
    decompositionOfCholesky(A, L);

    // Solve L * b = r, where L_T * c = b and L is the lower-triangular matrix.
    // Vector "b" and vector "r" replaces vector coefficients.
    for (int i = 0; i < sizeMatrix; ++i) {
        for (int j = 0; j < i; ++j) {
            coefficients[i] -= L[i][j] * coefficients[j];
        }
        coefficients[i] /= L[i][i];
    }

    // Solve L_T * c = b, where L_T is the upper-triangular matrix.
    for (int i = sizeMatrix - 1; i >= 0; i--) {
        for (int j = sizeMatrix - 1; j > i; j--) {
            coefficients[i] -= L[j][i] * coefficients[j];
        }
        coefficients[i] /= L[i][i];
    }
}


unsigned int Spline::getKnotsCount() const{
    return internalKnotsCount + 2;
}


unsigned int Spline::getDegree() const{
    return degree;
}


std::vector<double> Spline::getKnots() const{

    std::vector<double> differentKnots(internalKnotsCount + 2);

    for (int i = 0; i < internalKnotsCount + 2; i++){
        differentKnots[i] = knots[i + degree];
    }

    return differentKnots;
}


// private:

void Spline::computingMatrixA(std::vector<std::vector<double>>& A, const unsigned int& sizeMatrix,
                              const std::vector<double>& x, const std::vector<double>& y,
                              const std::vector<double>& w){

    std::vector<double> bSplines(degree + 1, 0.0);
    coefficients.resize(sizeMatrix);

    int leftBoundId = 0;
    double wPow2;

    for (unsigned int i = 0; i < x.size(); ++i) {
        leftBoundId = getLeftKnotIndex(x[i], leftBoundId);

        initializeBSplines(x[i], bSplines);

        for (unsigned int j = 0; j < degree + 1; j++){
            wPow2 = w[i]*w[i];
            for (unsigned int k = 0; k < j + 1; k++){
                A[j + leftBoundId - degree][k + leftBoundId - degree] += wPow2 * bSplines[j] * bSplines[k];
            }
            coefficients[j + leftBoundId - degree] += wPow2 * y[i] * bSplines[j];
        }
    }
    reflectionOfMatrix(A, sizeMatrix);
}


void Spline::initializeBSplines(const double& x, std::vector<double>& bSplines){

    const int leftBoardId = getLeftKnotIndex(x, 0);

    bSplines[degree] = 1;

    unsigned int v;
    double alphaN, alphaNN;

    // TODO: Optimization
    // Replacement of recurrent calculation of B-splines.
    for (unsigned int currentDegree = 1; currentDegree < degree + 1; currentDegree++){
        v = leftBoardId - currentDegree + 1;

        alphaNN = getAlpha(x, v, currentDegree);
        bSplines[degree - currentDegree] = (1-alphaNN) * bSplines[degree - currentDegree + 1];
        for (int j = degree - currentDegree + 1; j < degree; j++){
            alphaN = alphaNN;
            v++;
            alphaNN = getAlpha(x, v, currentDegree);
            bSplines[j] = alphaN * bSplines[j] + (1 - alphaNN) * bSplines[j+1];
        }
        bSplines[degree] *= alphaNN;
    }
}


// For B-Spline.
double Spline::getAlpha(const double &x, const unsigned int &knotId, const unsigned int &deg) const{
    return (x - knots[knotId]) / (knots[knotId + deg] - knots[knotId]);
}


// Get S(x).
double Spline::getValue(const double &x) const {
    if (x < knots[0] || x > knots[knotsCount-1]){
        return 0;
    }
    return deBoorAlgorithm(x);
}


double Spline::deBoorAlgorithm(const double& x) const{

    int leftBound = getLeftKnotIndex(x, 0);
    if (leftBound < 0){
        return 0;
    }

    std::vector<double> values(degree + 1, 0.0);
    double alpha;

    for (int i = 0; i < degree + 1; i++){
        values[i] = coefficients[i + leftBound - degree];
    }
    for (int i = 1; i < degree + 1; i++){
        for (auto j = static_cast<unsigned int>(leftBound); j > leftBound - degree + i - 1; j--){
            alpha = getAlpha(x, j, degree - i + 1);
            values[j - leftBound + degree] = alpha * values[j - leftBound + degree] +
                                             (1 - alpha) * values[j - leftBound + degree - 1];
        }
    }
    return values[degree];
}


int Spline::getLeftKnotIndex(const double &x, const int &minId) const{
    if(x < knots[0] or x > knots[knotsCount-1]){
        return -1;
    }

    if (x == knots[knotsCount - 1]){
        return knotsCount - degree - 2;
    }

    for(int ind = minId; ind < knotsCount + degree - 2;){
        if(knots[ind] <= x && knots[ind + 1] > x){
            return ind;
        }
        else {
            ind++;
        }
    }
    return minId;
}

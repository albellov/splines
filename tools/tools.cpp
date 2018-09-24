//
// Created by Aleksandr on 06-Sep-18.
//
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

#include "tools.h"


template <class Container>
void split(const std::string& str, Container& cont, char delim = ' ');

// A is a decomposition of the form A = L*L_T, where L is the lower-triangular matrix.
void decompositionOfCholesky(const std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& L){
    for (int i = 0; i < L.size(); ++i) {
        for (int j = 0; j <= i; ++j) {
            if (i == j) {
                L[i][j] = A[i][j];
                for (int k = 0; k < j; ++k)
                    L[i][j] -= L[i][k] * L[i][k];
                L[i][j] = std::sqrt(L[i][j]);
            }

            else {
                L[i][j] = A[i][j];
                for (int k = 0; k < j; ++k)
                    L[i][j] -= L[i][k] * L[j][k];
                L[i][j] /= L[j][j];
            }
        }
    }
}


void reflectionOfMatrix(std::vector<std::vector<double>>& matrix, const unsigned int& size){
    for (int i = 0; i < size; i++){
        for (int j = 0; j < i; j++){
            matrix[j][i] = matrix[i][j];
        }
    }
}


void writingResult2d(const double &leftBound, const double &rightBound, const unsigned int &stepCount,
                     const UnivariateSpline &spl, const std::string &out_name,
                     const std::vector<std::string> &vals_name) {

    std::ofstream o_file(out_name);
    o_file << vals_name[0] << "," << vals_name[1] + "_spl" << std::endl;

    const double step = (rightBound - leftBound) / stepCount;

    double currentPoint;
    for (unsigned int i = 0; i < stepCount + 1; ++i) {

        currentPoint = leftBound + step*i;
        o_file << currentPoint << "," << spl.getValue(currentPoint) << std::endl;
    }
    o_file.close();


    o_file.open(out_name + "_knots");
    o_file << vals_name[0] << "," << vals_name[1] + "_knot" << std::endl;

    std::vector<double> knots = spl.getKnots();
    for (unsigned int i = 0; i < spl.getKnotsCount(); ++i) {

        currentPoint = knots[i];
        o_file << currentPoint << "," << spl.getValue(currentPoint) << std::endl;
    }
    o_file.close();
}

void writingResult3d(const double &leftBound, const double &rightBound, const unsigned int &stepCount,
                     const BivariateSpline &spl, const std::string &out_name, const std::vector<std::string> &vals_name) {

    std::ofstream o_file(out_name);
    o_file << vals_name[0] << "," << vals_name[1] + "_spl" << std::endl;

    const double step = (rightBound - leftBound) / stepCount;

    double currentPoint;
    for (unsigned int i = 0; i < stepCount + 1; ++i) {

        currentPoint = leftBound + step*i;
        o_file << currentPoint << "," << spl.getValue(currentPoint) << std::endl;
    }
    o_file.close();


    o_file.open(out_name + "_knots");
    o_file << vals_name[0] << "," << vals_name[1] + "_knot" << std::endl;

    std::vector<double> knots = spl.getKnots();
    for (unsigned int i = 0; i < spl.getKnotsCount(); ++i) {

        currentPoint = knots[i];
        o_file << currentPoint << "," << spl.getValue(currentPoint) << std::endl;
    }
    o_file.close();
}


void readingData2d(std::vector<double> &x, std::vector<double> &y, std::vector<double> &w,
                   std::vector<std::string> &vals_name, const std::string &inputFilename) {

    std::ifstream i_file(inputFilename);

    std::vector<std::string> strValues;
    std::string line;
    double wValue;

    if (!i_file.is_open())
    {
        std::cout << "In file doesn't exist!" << std::endl;
        exit(1);
    }

    std::getline(i_file, line);
    split(line, vals_name, ',');

    while (true)
    {
        std::getline(i_file, line);
        if (i_file.eof()) {
            break;
        }

        split(line, strValues, ',');

        // TODO: Add Weights.
        wValue = 1.0;

        x.push_back(stod(strValues[0], nullptr));
        y.push_back(stod(strValues[1], nullptr));
        w.push_back(wValue);
    }
    i_file.close();
}

void readingData3d(std::vector<double> &x, std::vector<double> &y, std::vector<double> &z, std::vector<double> &w,
                   std::vector<std::string> &vals_name, const std::string &inputFilename) {

    std::ifstream i_file(inputFilename);

    std::vector<std::string> strValues;
    std::string line;
    double wValue;

    if (!i_file.is_open())
    {
        std::cout << "File "<< inputFilename << " doesn't exist!" << std::endl;
        exit(1);
    }

    std::getline(i_file, line);
    split(line, vals_name, ',');

    while (true)
    {
        std::getline(i_file, line);
        if (i_file.eof()) {
            break;
        }

        split(line, strValues, ',');

        // TODO: Add Weights.
        wValue = 1.0;

        x.push_back(stod(strValues[0], nullptr));
        y.push_back(stod(strValues[1], nullptr));
        z.push_back(stod(strValues[2], nullptr));
        w.push_back(wValue);
    }
    i_file.close();
}

template <class Container>
void split(const std::string& str, Container& cont, char delim)
{
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delim)) {
        cont.push_back(token);
    }
}
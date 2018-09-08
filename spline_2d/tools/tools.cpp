//
// Created by Aleksandr on 06-Sep-18.
//
#include <cmath>
#include <fstream>
#include <iostream>

#include "tools.h"
#include "../spline.h"


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


void writing_result(const double& leftBound, const double& rightBound, const unsigned int& stepCount,
                    const Spline& spl, const std::string& out_name, const std::vector<std::string>& vals_name) {

    std::ofstream o_file(out_name);
    o_file << vals_name[0] << "," << vals_name[1] + "_spl" << std::endl;

    const double step = (rightBound - leftBound) / stepCount;

    double currentPoint;
    for (unsigned int i = 0; i < stepCount + 1; ++i) {

        currentPoint = leftBound + step*i;
        o_file << currentPoint << "," << spl.getValue(currentPoint) << std::endl;
    }
    o_file.close();
}


void reading_data(std::vector<double>& x, std::vector<double>& y, std::vector<double>& w,
                  std::vector<std::string>& vals_name, const std::string& inputFilename) {

    std::ifstream i_file(inputFilename);
    std::string line;

    if (!i_file.is_open())
    {
        std::cout << "In file doesn't exist!" << std::endl;
        exit(1);
    }

    std::getline(i_file, line);
    vals_name = { line.substr(0,line.find(',')), line.substr(line.find(',') + 1,line.size() - 1) };

    while (true)
    {
        std::getline(i_file, line);
        if (i_file.eof()) {
            break;
        }

        //TODO: Refactoring
        x.push_back(strtod(line.substr(0, line.find(',')).c_str(), nullptr));
        y.push_back(strtod(line.substr(line.find(',') + 1, line.size() - 1).c_str(), nullptr));
        w.push_back(1.0);
    }
    i_file.close();
}
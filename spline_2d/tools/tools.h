//
// Created by Aleksandr on 06-Sep-18.
//

#ifndef SPLINE_2D_TOOLS_H
#define SPLINE_2D_TOOLS_H

#include "../src/spline.h"



void reading_data(std::vector<double>& x, std::vector<double>& y, std::vector<double>& w,
                  std::vector<std::string>& vals_name, const std::string& inputFilename);
void writing_result(const double& leftBound, const double& rightBound, const unsigned int& stepCount,
                    const Spline& spl, const std::string& out_name, const std::vector<std::string>& vals_name);

void decompositionOfCholesky(const std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& L);
void reflectionOfMatrix(std::vector<std::vector<double>>& matrix, const unsigned int& size);

#endif //SPLINE_2D_TOOLS_H

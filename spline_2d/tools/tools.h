//
// Created by Aleksandr on 06-Sep-18.
//

#ifndef SPLINE_2D_TOOLS_H
#define SPLINE_2D_TOOLS_H

#include "../spline.h"

void get_coefficients(const std::vector<double> &x, const std::vector<double> &y, const Spline &spl,
                      std::vector<double> &c);

#endif //SPLINE_2D_TOOLS_H

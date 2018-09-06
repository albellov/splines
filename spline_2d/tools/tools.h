//
// Created by Aleksandr on 06-Sep-18.
//

#ifndef SPLINE_2D_TOOLS_H
#define SPLINE_2D_TOOLS_H

#include "../spline.h"


void getCoefficients(const std::vector<double> &x, const std::vector<double> &y, const Spline &spl,
                     std::vector<double> &c);
void reading_data(std::vector<double>&, std::vector<double>&, std::vector<std::string>&, const std::string&);
void writing_base(const std::vector<double>&, const Spline&, const std::string&, const std::string&, const unsigned int&);
void writing_result(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&,
                    const Spline&, const std::string&, const std::vector<std::string>&, const unsigned int&);


#endif //SPLINE_2D_TOOLS_H

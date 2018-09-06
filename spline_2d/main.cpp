//
// Created by Aleksandr on 04-Sep-18.
//

#include <cmath>
#include <string>
#include <iostream>
#include <vector>

#include "spline.h"
#include "tools/tools.h"


int main(int agrc, char* argv[])
{
    if (agrc < 5){
        std::cout << "Args: <in filename> <out filename> <spline knots count> <spline degree>" << std::endl;
        return 1;
    }

    const std::string in_name = argv[1];
    const std::string out_name = argv[2];
    const unsigned int knots_count = strtoul(argv[3], nullptr, 10);
    const unsigned int spline_degree = strtoul(argv[4], nullptr, 10);

    const unsigned int basis_count = knots_count + spline_degree - 1;
    const unsigned int basis_steps_count = 1000;

    std::vector<std::string> vals_name;
    std::vector<double> x, y;
    std::vector<double> c(basis_count, 0.0);

    reading_data(x, y, vals_name, in_name);

    const Spline spl(spline_degree, knots_count, x);

    getCoefficients(x, y, spl, c);

    // Constants found, spline ready.
    //writing_base(x, spl, out_name, vals_name[0], basis_steps_count);
    //writing_result(x, y, c, spl, out_name, vals_name, basis_steps_count);

    return 0;
}
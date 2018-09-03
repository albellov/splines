//
// Created by Aleksandr on 03-Sep-18.
//

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <random>


int main(int agrc, char* argv[])
{
    std::string f_type;
    double begin_val;
    double end_val;
    int step_count;
    std::string out_name;

    if (agrc > 5)
    {
        f_type = argv[1];
        begin_val = strtof(argv[2], nullptr);
        end_val = strtof(argv[3], nullptr);
        step_count = strtol(argv[4], nullptr, 10);
        out_name = argv[5];
        if (f_type != "linear"&&f_type != "exp"&&f_type != "sin"&&f_type != "cos")
        {
            std::cout << "Unknown function type: <" << f_type << ">" << std::endl
                      << "Applied vals - 'linear' 'exp', 'sin', 'cos'" << std::endl;
            return -1;
        }
        if (end_val - begin_val < 0)
        {
            std::cout << "Error: begin val > end val" << std::endl;
            return -1;
        }
    }
    else
    {
        std::cout << "Args: <func type*> <begin val> <end val> <step count> <out filename>" << std::endl
                  << "*: Applied vals - 'linear' 'exp', 'sin', 'cos'" << std::endl;
        return -1;
    }

    std::ofstream o_f(out_name);
    std::stringstream ss;
    std::random_device rd;
    std::default_random_engine mt(rd());
    std::uniform_real_distribution<double> noise(-0.5, 0.5);

    ss << "X, Y" << std::endl;
    o_f << ss.rdbuf();
    ss.clear(); ss.str("");

    double h = (end_val - begin_val) / (double)step_count;
    double x = begin_val, y = 0.0;

    for (int i = 0; i < step_count; ++i)
    {
        if (f_type == "linear")
            y = x;
        if (f_type == "exp")
            y = exp(x);
        if (f_type == "sin")
            y = sin(x);
        if (f_type == "cos")
            y = cos(x);

        y += noise(mt);

        ss << std::scientific << x << ", " << std::scientific << y << std::endl;
        o_f << ss.rdbuf();
        ss.clear(); ss.str("");
        x += h;
    }

    return 1;
}

//
// Created by Aleksandr on 03-Sep-18.
//

#include <fstream>
#include <string>
#include <iostream>
#include <random>


int main(int argc, char* argv[]){

    if (argc < 6){
        std::cout << "Args: <func type*> <begin val> <end val> <step count> <out filename>" << std::endl
                  << "*: Applied vals - 'linear' 'exp', 'sin', 'cos'" << std::endl;
        return 1;
    }

    const std::string f_type = argv[1];
    const double begin_val = strtod(argv[2], nullptr);
    const double end_val = strtod(argv[3], nullptr);
    const int step_count = strtol(argv[4], nullptr, 10);
    const std::string out_name = argv[5];

    if (f_type != "linear"&&f_type != "exp"&&f_type != "sin"&&f_type != "cos"&&f_type != "spec")
    {
        std::cout << "Unknown function type: <" << f_type << ">" << std::endl
                  << "Applied vals - 'linear' 'exp', 'sin', 'cos'" << std::endl;
        return 1;
    }
    if (end_val - begin_val < 0)
    {
        std::cout << "Error: begin_val > end_val" << std::endl;
        return 1;
    }

    const float sigma = 0.2;

    std::ofstream o_f(out_name);
    std::random_device rd;
    std::default_random_engine mt(rd());
    std::uniform_real_distribution<double> noise(0, 1);

    o_f << "X, Y" << std::endl;

    double h = (end_val - begin_val) / (double)step_count;
    double x = begin_val, y = 0.0;

    for (int i = 0; i < step_count; ++i)
    {

        x = i + 1./(i+1) * noise(mt);
        if (f_type == "linear")
            y = x;
        if (f_type == "exp")
            y = std::exp(x);
        if (f_type == "sin")
            y = std::sin(x);
        if (f_type == "cos")
            y = std::cos(x);
        if (f_type == "spec")
            y = std::sin(0.2*x);

        y += noise(mt)*15/(10+x);

        o_f << std::scientific << x << ", " << std::scientific << y << std::endl;
        x += h;
    }

    return 0;
}
//
// Created by Aleksandr on 03-Sep-18.
//

#include <fstream>
#include <string>
#include <iostream>
#include <random>
#include <ctime>


double getY(const double&, const std::string&);
void printAllowFunctionTypes(const std::vector<std::string>& allowFunctionTypes);
void checkFunctionTypes(const std::vector<std::string>& allowFunctionTypes, const std::string& functionType);

int main(int argc, char* argv[]){

    const std::vector<std::string> allowFunctionTypes = {"linear", "exp", "sin", "cos", "spec"};

    if (argc < 6){
        std::cout << "Args: <func type*> <begin val> <end val> <step count> <out filename>" << std::endl;
        printAllowFunctionTypes(allowFunctionTypes);
        return 1;
    }

    const std::string functionType = argv[1];
    const double beginValue = strtod(argv[2], nullptr);
    const double endValue = strtod(argv[3], nullptr);
    const int stepCount = strtol(argv[4], nullptr, 10);
    const std::string outputFilename = argv[5];

    checkFunctionTypes(allowFunctionTypes, functionType);

    if (endValue - beginValue < 0)
    {
        std::cout << "Error: begin value > end value" << std::endl;
        return 1;
    }

    std::ofstream o_f(outputFilename);

    o_f << "X,Y" << std::endl;


    std::mt19937_64 mt(static_cast<unsigned long long int>(time(nullptr)));
    std::uniform_real_distribution<> noise(-1, 1);

    double h = (endValue - beginValue) / (double)stepCount;
    double x = beginValue, y;

    // const double sigma = 0.2;

    for (int i = 0; i < stepCount; ++i)
    {
        y = getY(x, functionType) + noise(mt)*15/(10+x);
        o_f << std::scientific << x << "," << std::scientific << y << std::endl;
        x += h;
    }

    return 0;
}


double getY(const double& x, const std::string& functionType){

    double y = 0;

    if (functionType == "linear")
        y = x;
    else if (functionType == "exp")
        y = std::exp(x);
    else if (functionType == "sin")
        y = std::sin(x);
    else if (functionType == "cos")
        y = std::cos(x);
    else if (functionType == "spec")
        y = std::sin(0.2 * x);

    return y;
}

void checkFunctionTypes(const std::vector<std::string>& allowFunctionTypes, const std::string& functionType){
    bool isAllowFunctionType = false;

    for (const std::string &allowFunctionType : allowFunctionTypes) {
        if (functionType == allowFunctionType){
            isAllowFunctionType = true;
            break;
        }
    }
    if (!isAllowFunctionType){
        std::cout << "Unknown function type: <" << functionType << ">" << std::endl;
        printAllowFunctionTypes(allowFunctionTypes);
        exit(1);
    }
}


void printAllowFunctionTypes(const std::vector<std::string>& allowFunctionTypes){
    std::cout << "Applied vals - ";
    for (const std::string &allowFunctionType : allowFunctionTypes){
        std::cout << "'" << allowFunctionType << "', ";
    }
    std::cout << std::endl;
}
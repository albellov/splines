//
// Created by Aleksandr on 03-Sep-18.
//

#include <fstream>
#include <string>
#include <iostream>
#include <random>
#include <ctime>


void dataGeneration2d(const std::vector<std::string>& allowFunctionTypes, char **argv);
void dataGeneration3d(const std::vector<std::string>& allowFunctionTypes, char **argv);

double getY(const double& x, const std::string& functionType);
double getZ(const double& x, const double& y, const std::string& functionType);
void printAllowFunctionTypes(const std::vector<std::string>& allowFunctionTypes, const int& dim);
void checkFunctionTypes(const std::vector<std::string>& allowFunctionTypes,
                        const std::string& functionType, const int& dim);


int main(int argc, char* argv[]){

    const std::vector<std::vector<std::string>> allowFunctionTypes = { {"linear", "exp", "sin", "cos", "spec"},
                                                                       {"spec"}
                                                                       };

    if (argc != 7 && argc != 9){
        std::cout << "Args: <dim> <func type> <begin val x> <end val x>\n" <<
                     "      [<begin val y> <end val y>]** <step count> <out filename>\n" <<
                     "      [for dim = 3]" << std::endl;
        printAllowFunctionTypes(allowFunctionTypes[0], 2);
        printAllowFunctionTypes(allowFunctionTypes[1], 3);
        return 1;
    }

    //std::vector<std::string> allowFunctionTypes;
    const int dim = strtol(argv[1], nullptr, 10);

    if (dim == 2){
        dataGeneration2d(allowFunctionTypes[0], argv);
    }
    else if (dim == 3){
        dataGeneration3d(allowFunctionTypes[1], argv);
    }
    else{
        std::cout << "Error: dim may be 2 or 3" << std::endl;
        return 1;
    }

    return 0;
}


void dataGeneration2d(const std::vector<std::string>& allowFunctionTypes, char **argv){

    const int dim = 2;

    const std::string functionType = argv[2];
    const double beginX = strtod(argv[3], nullptr);
    const double endX = strtod(argv[4], nullptr);
    const int stepCount = strtol(argv[5], nullptr, 10);
    const std::string outputFilename = argv[6];

    checkFunctionTypes(allowFunctionTypes, functionType, dim);

    if (endX - beginX < 0)
    {
        std::cout << "Error: begin value > end value" << std::endl;
        exit(1);
    }

    std::ofstream o_f(outputFilename);

    o_f << "X,Y" << std::endl;

    std::mt19937_64 mt(static_cast<unsigned long long int>(time(nullptr)));
    std::uniform_real_distribution<> noise(-1, 1);

    double h = (endX - beginX) / (double)stepCount;
    double x = beginX, y;

    for (int i = 0; i < stepCount; ++i)
    {
        y = getY(x, functionType) + noise(mt)*15/(10+x);
        o_f << std::scientific << x << "," << std::scientific << y << std::endl;
        x += h;
    }
}

void dataGeneration3d(const std::vector<std::string>& allowFunctionTypes, char **argv){

    const int dim = 3;

    const std::string functionType = argv[2];
    const double beginX = strtod(argv[3], nullptr);
    const double endX = strtod(argv[4], nullptr);
    const double beginY = strtod(argv[5], nullptr);
    const double endY = strtod(argv[6], nullptr);
    const int stepCount = strtol(argv[7], nullptr, 10);
    const std::string outputFilename = argv[8];

    checkFunctionTypes(allowFunctionTypes, functionType, dim);

    if (endX - beginX < 0 || endY - beginY < 0)
    {
        std::cout << "Error: begin value > end value" << std::endl;
        exit(1);
    }

    std::ofstream o_f(outputFilename);

    o_f << "X,Y,Z" << std::endl;

    std::mt19937_64 mt(static_cast<unsigned long long int>(time(nullptr)));
    std::uniform_real_distribution<> noise(-1, 1);

    // TODO: Add different grids.
    const int stepCountX = stepCount;
    const int stepCountY = stepCount;

    double hX = (endX - beginX) / (double)stepCountX;
    double hY = (endY - beginY) / (double)stepCountY;

    double y = beginY, z;

    for (int j = 0; j < stepCountY; ++j)
    {
        double x = beginX;
        for (int i = 0; i < stepCountX; ++i){
            z = getZ(x, y, functionType) + noise(mt)*15/(10+x);
            o_f << std::scientific << x << "," << std::scientific << y << "," << std::scientific << z << std::endl;
            x += hX;
        }
        y += hY;
    }
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

double getZ(const double& x, const double& y, const std::string& functionType){

    double z = 0;

    if (functionType == "linear")
        z = x;
    else if (functionType == "exp")
        z = std::exp(x);
    else if (functionType == "sin")
        z = std::sin(x);
    else if (functionType == "cos")
        z = std::cos(x);
    else if (functionType == "spec")
        z = std::sin(0.2 * x);

    return z;
}

void checkFunctionTypes(const std::vector<std::string>& allowFunctionTypes,
                        const std::string& functionType, const int& dim){

    for (const std::string &allowFunctionType : allowFunctionTypes) {
        if (functionType == allowFunctionType) {
            return;
        }
    }

    std::cout << "Unknown function type: <" << functionType << ">" << std::endl;
    printAllowFunctionTypes(allowFunctionTypes, dim);
    exit(1);
}


void printAllowFunctionTypes(const std::vector<std::string>& allowFunctionTypes, const int& dim){
    std::cout << "Allow function types for dim = " << dim << ": ";
    for (const std::string &allowFunctionType : allowFunctionTypes){
        std::cout << "'" << allowFunctionType << "', ";
    }
    std::cout << std::endl;
}
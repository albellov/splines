//
// Created by Aleksandr on 04-Sep-18.
//

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <vector>

#include "spline.h"


void making_monospace_knots(std::vector<double>&, const std::vector<double>&, const unsigned int&);
void reading_data(std::vector<double>&, std::vector<double>&, std::vector<std::string>&, const std::string&);
void writing_base(const std::vector<double>&, const Spline&, const std::string&, const std::string&, const unsigned int&);
void method_Cholesky(const std::vector<double> &, const std::vector<double> &, const Spline &, std::vector<double> &);
void writing_result(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&,
                    const Spline&, const std::string&, const std::vector<std::string>&, const unsigned int&);
double get_basis_step(const std::vector<double>&, const unsigned int&);


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
    std::vector<double> knots;
    std::vector<double> c(basis_count, 0.0);

    reading_data(x, y, vals_name, in_name);

    making_monospace_knots(knots, x, knots_count);
    const Spline spl(spline_degree, knots);

    method_Cholesky(x, y, spl, c);

    // Constants found, spline ready.
    writing_base(x, spl, out_name, vals_name[0], basis_steps_count);
    writing_result(x, y, c, spl, out_name, vals_name, basis_steps_count);

    return 0;
}


/*
	Solve A*c = r, where

               | (N_o, N_o)			...			(N_o, N_(g+k))	|
	A = E_T*E =|	...				...			...				|
			   | (N_(g+k), N_o)		...		(N_(g+k), N_(g+k))	|


				|	(N_o, y)	|
	r = E_T*y = |	  ...		|
				| ((N_(g+k), y)	|

	k = spline_degree,
	g = knots_count
*/
void method_Cholesky(const std::vector<double> &x, const std::vector<double> &y, const Spline &spl,
                     std::vector<double> &c) {

    const unsigned int knots_count = spl.get_knots_count() - 1;
    const unsigned int spline_degree = spl.get_degree();
    const unsigned int basis_count = knots_count + spline_degree;// -1;

    std::vector<std::vector<double>> A(basis_count, std::vector<double>(basis_count, 0.0));
    std::vector<double> r(basis_count, 0.0);

    for (unsigned int i = 0; i < x.size(); ++i) {
        for (unsigned int j = 1; j < basis_count + 1; ++j)
        {
            r[j - 1] += y[i] * spl.get_basis_val(x[i], j, spline_degree);
            for (unsigned int k = 1; k < knots_count + spline_degree; ++k) {
                A[j - 1][k - 1] += spl.get_basis_val(x[i], j, spline_degree) * spl.get_basis_val(x[i], k, spline_degree);
            }
        }
    }


    std::vector<std::vector<double>> L(basis_count, std::vector <double>(basis_count, 0.0));
    for (int i = 0; i<basis_count; ++i)
        for (int j = 0; j <= i; ++j)
        {
            if (i == j)
            {
                L[i][j] = A[i][j];
                for (int k = 0; k<j; ++k)
                    L[i][j] -= L[i][k] * L[i][k];
                L[i][j] = std::sqrt(L[i][j]);
            }
            else
            {
                L[i][j] = A[i][j];
                for (int k = 0; k<j; ++k)
                    L[i][j] -= L[i][k] * L[j][k];
                L[i][j] /= L[j][j];
            }
        }

    std::vector <double> b(basis_count, 0.0);
    for (int i = 0; i < basis_count; ++i) {

        double sum_buf = r[i];

        for (int j = 0; j < i; ++j) {
            sum_buf -= L[i][j] * b[j];
        }

        b[i] = sum_buf / L[i][i];
    }

    std::vector<std::vector <double>> L_s(basis_count, std::vector <double>(basis_count, 0.0));
    for (int i = 0; i < basis_count; ++i) {
        for (int j = 0; j < basis_count; ++j) {
            L_s[i][j] = L[j][i];
        }
    }

    for (int i = 0; i < basis_count; ++i) {
        double sum_buf = b[basis_count - 1 - i];
        for (int j = 0; j < i; ++j) {
            sum_buf -= L_s[basis_count - 1 - i][basis_count - 1 - j] * c[basis_count - 1 - j];
        }
        c[basis_count - 1 - i] = sum_buf / L_s[basis_count - 1 - i][basis_count - 1 - i];
    }
}


void writing_result(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& c,
                    const Spline& spl, const std::string& out_name, const std::vector<std::string>& vals_name,
                    const unsigned int& basis_steps_count) {

    const unsigned int knots_count = spl.get_knots_count() - 1;
    const unsigned int spline_degree = spl.get_degree();
    const double h = get_basis_step(x, basis_steps_count);
    const double x_min = x[0];

    std::ofstream o_file(out_name);
    o_file << vals_name[0] << "," << vals_name[1] + "_orig" << "," << vals_name[1] + "_spl" << std::endl;
    double spl_val;

    for (unsigned int i = 0; i < x.size(); ++i) {
        o_file << x[i] << "," << y[i];
        spl_val = 0.0;
        for (unsigned int j = 1; j < knots_count + spline_degree; ++j) {
            spl_val += c[j - 1] * spl.get_basis_val(x_min + i*h, j, spline_degree);
        }
        o_file << "," << spl_val << std::endl;
    }
    o_file.close();
}


void writing_base(const std::vector<double>& x, const Spline& spl, const std::string& out_name,
                  const std::string& first_val, const unsigned int& basis_steps_count) {

    const unsigned int knots_count = spl.get_knots_count() - 1;
    const unsigned int spline_degree = spl.get_degree();
    const double x_min = x[0];

    std::ofstream o_file(out_name.substr(0, out_name.find('.')) + "_basis" + out_name.substr(out_name.find('.'), out_name.size() - 1));

    o_file << first_val;
    for (int i = 1; i < knots_count + spline_degree; ++i) {
        o_file << ", B_" << i;
    }
    o_file << std::endl;

    // Writing basis for vis.
    const double h = get_basis_step(x, basis_steps_count);
    for (int j = 0; j < basis_steps_count; ++j)
    {
        o_file << x_min + j*h;
        for (unsigned int i = 1; i < knots_count + spline_degree; ++i) {
            o_file << ", " << spl.get_basis_val(x_min + j*h, i, spline_degree);
        }
        o_file << std::endl;
    }
    o_file.close();
}


void making_monospace_knots(std::vector<double>& knots, const std::vector<double>& x, const unsigned int& knots_count){

    const double x_min = x[0];
    const double x_max = x[x.size() - 1];

    const double knots_step = (x_max - x_min) / (double)knots_count;

    for (int i = 0; i < knots_count; ++i) {
        knots.push_back(x_min + i*knots_step);
    }
    knots.push_back(x_max + 0.001);
}

void reading_data(std::vector<double>& x, std::vector<double>& y, std::vector<std::string>& vals_name,
                  const std::string& in_name) {

    std::ifstream i_file(in_name);
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
    }
    i_file.close();
}


double get_basis_step(const std::vector<double>& x, const unsigned int& basis_steps_count) {
    return (x[x.size() - 1] - x[0]) / (double)basis_steps_count;
}

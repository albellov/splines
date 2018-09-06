//
// Created by Aleksandr on 06-Sep-18.
//
#include <cmath>
#include <fstream>
#include <iostream>

#include "tools.h"

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

void decompositionOfCholesky(std::vector<std::vector<double>>& L, const std::vector<std::vector<double>>& A){
    for (int i = 0; i < L.size(); ++i) {
        for (int j = 0; j <= i; ++j) {
            if (i == j) {
                L[i][j] = A[i][j];
                for (int k = 0; k < j; ++k)
                    L[i][j] -= L[i][k] * L[i][k];
                L[i][j] = std::sqrt(L[i][j]);
            } else {
                L[i][j] = A[i][j];
                for (int k = 0; k < j; ++k)
                    L[i][j] -= L[i][k] * L[j][k];
                L[i][j] /= L[j][j];
            }
        }
    }
}


void getCoefficients(const std::vector<double> &x, const std::vector<double> &y, const Spline &spl,
                     std::vector<double> &c) {

    double sum_buf = 0;
    const unsigned int knots_count = spl.getKnotsCount() - 1;
    const unsigned int spline_degree = spl.getDegree();
    const unsigned int basis_count = knots_count + spline_degree -1;

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
    decompositionOfCholesky(L, A);

    std::vector <double> b(basis_count, 0.0);
    for (int i = 0; i < basis_count; ++i) {

        sum_buf = r[i];

        for (int j = 0; j < i; ++j) {
            sum_buf -= L[i][j] * b[j];
        }

        b[i] = sum_buf / L[i][i];
    }

    for (int i = 0; i < basis_count; ++i) {
        sum_buf = b[basis_count - 1 - i];
        for (int j = 0; j < i; ++j) {
            sum_buf -= L[basis_count - 1 - j][basis_count - 1 - i] * c[basis_count - 1 - j];
        }
        c[basis_count - 1 - i] = sum_buf / L[basis_count - 1 - i][basis_count - 1 - i];
    }
}

double get_basis_step(const std::vector<double>& x, const unsigned int& basis_steps_count) {
    return (x[x.size() - 1] - x[0]) / (double)basis_steps_count;
}

void writing_result(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& c,
                    const Spline& spl, const std::string& out_name, const std::vector<std::string>& vals_name,
                    const unsigned int& basis_steps_count) {

    const unsigned int knots_count = spl.getKnotsCount() - 1;
    const unsigned int spline_degree = spl.getDegree();
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

    const unsigned int knots_count = spl.getKnotsCount() - 1;
    const unsigned int spline_degree = spl.getDegree();
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
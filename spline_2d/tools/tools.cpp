//
// Created by Aleksandr on 06-Sep-18.
//
#include "tools.h"
#include <cmath>


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

void get_coefficients(const std::vector<double> &x, const std::vector<double> &y, const Spline &spl,
                      std::vector<double> &c) {

    const unsigned int knots_count = spl.get_knots_count() - 1;
    const unsigned int spline_degree = spl.get_degree();
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
//
// Created by Aleksandr on 04-Sep-18.
//

#ifndef SPLINE_2D_SPLINE_H
#define SPLINE_2D_SPLINE_H

#include <vector>

class Spline
{
private:
    std::vector <double> knots;
    unsigned int degree;
    unsigned int knots_count;

    double alpha(const double& x, const unsigned int& knot_id, const unsigned int& deg) const;

public:
    Spline(const unsigned int& sp_deg, const std::vector<double>& knots_vals);

    // Recurrent process of getting basis value.
    double get_basis_val(const double& x, const unsigned int& knot_id, const unsigned int& deg) const;

    unsigned int get_knots_count() const;
    unsigned int get_degree() const;
};



#endif //SPLINE_2D_SPLINE_H

//
// Created by Aleksandr on 04-Sep-18.
//

#include "spline.h"


double Spline::alpha(const double& x, const unsigned int& knot_id, const unsigned int& deg) const{
    if (knots[knot_id] == knots[knot_id + deg])
        return 0;
    else
        return (x - knots[knot_id]) / (knots[knot_id + deg] - knots[knot_id]);
}


Spline::Spline(const unsigned int& sp_deg, const std::vector<double>& knots_vals){
    degree = sp_deg;
    knots_count = static_cast<unsigned int>(knots_vals.size());

    // Knots duplicates spline_degree times in the first and the last positions.
    knots.resize(knots_count + 2 * degree);

    for (int i = 0; i < degree; ++i)
    {
        knots[i] = knots_vals[0];
        knots[knots_count + degree + i] = knots_vals[knots_count - 1];
    }
    for (int i = 0; i < knots_count; ++i)
        knots[degree + i] = knots_vals[i];
}


double Spline::get_basis_val(const double& x, const unsigned int& knot_id, const unsigned int& deg) const{
    if ((x < knots[knot_id]) || (x >= knots[knot_id + deg + 1]))
        return 0;
    if (deg == 1)
    {
        if ((x >= knots[knot_id]) && (x < knots[knot_id + 1]))
            return 1;
        else
            return 0;
    }
    double alpha_n = alpha(x, knot_id, deg - 1);
    double alpha_nn = (1 - alpha(x, knot_id + 1, deg - 1));
    double bas_val = 0;

    if (alpha_n > 0.0) {
        bas_val = alpha_n*get_basis_val(x, knot_id, deg - 1);
    }
    if (alpha_nn > 0.0) {
        bas_val += alpha_nn*get_basis_val(x, knot_id + 1, deg - 1);
    }

    return bas_val;
}

unsigned int Spline::get_knots_count() const{
    return knots_count;
}

unsigned int Spline::get_degree() const {
    return degree;
}
//
// Created by Aleksandr on 04-Sep-18.
//

#include <iostream>
#include "spline.h"


double Spline::getAlpha(const double &x, const unsigned int &knot_id, const unsigned int &deg) const{
    if (knots[knot_id] == knots[knot_id + deg])
        return 0;
    else
        return (x - knots[knot_id]) / (knots[knot_id + deg] - knots[knot_id]);
}


void Spline::initializationOfKnots(const std::vector<double> &x) {
    const auto points_count = static_cast<const unsigned int>(x.size());
    unsigned int unique_points = 0;
    unsigned int counter;

    knots.resize(knots_count);
    knots[0] = x[0];
    knots[knots_count-1] = x[points_count-1];

    for(int i = 1; i < points_count; i++){
        if(x[i] != x[i-1]){
            unique_points ++;
        }
    }

    // Number of knots should be less than n - k for n points with unique x.
    if(unique_points == 0 || unique_points - degree < knots_count){
        return;
    }

    unsigned int points_for_knot = unique_points / (knots_count - 1);

    // TODO: Declaration
    for(int i = 1, k = 1, j = 0; i < knots_count - 1; i++){
        for(; j < points_for_knot * i || x[k] == x[k-1]; k++){
            if(x[k]!=x[k-1]){
                j++;
            }
        }
        knots[i] = 0.5 * (x[k] + x[k-1]);
    }
}


Spline::Spline(const unsigned int& degree, const unsigned int& knots_count, const std::vector<double>& x){
    this->degree = degree;
    this->knots_count = knots_count;
    initializationOfKnots(x);
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
    double alpha_n = getAlpha(x, knot_id, deg - 1);
    double alpha_nn = (1 - getAlpha(x, knot_id + 1, deg - 1));
    double bas_val = 0;

    if (alpha_n > 0.0) {
        bas_val = alpha_n*get_basis_val(x, knot_id, deg - 1);
    }
    if (alpha_nn > 0.0) {
        bas_val += alpha_nn*get_basis_val(x, knot_id + 1, deg - 1);
    }

    return bas_val;
}

unsigned int Spline::getKnotsCount() const{
    return knots_count;
}

unsigned int Spline::getDegree() const {
    return degree;
}
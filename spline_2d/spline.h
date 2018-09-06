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

    double getAlpha(const double &x, const unsigned int &knot_id, const unsigned int &deg) const;
    void initializationOfKnots(const std::vector<double> &x);

public:
    Spline(const unsigned int& degree, const unsigned int& knots_count, const std::vector<double>& x);

    // Recurrent process of getting basis value.
    double get_basis_val(const double& x, const unsigned int& knot_id, const unsigned int& deg) const;

    unsigned int getKnotsCount() const;
    unsigned int getDegree() const;
};


#endif //SPLINE_2D_SPLINE_H

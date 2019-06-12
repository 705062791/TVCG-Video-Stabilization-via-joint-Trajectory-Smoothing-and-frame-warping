#include "Quad.h"
#include <cmath>

Bundled::Quad::Quad(const v2d & V00, const v2d & V01, const v2d & V10, const v2d & V11)
{
    V00_[0] = V00[0];
    V00_[1] = V00[1];

    V01_[0] = V01[0];
    V01_[1] = V01[1];

    V10_[0] = V10[0];
    V10_[1] = V10[1];

    V11_[0] = V11[0];
    V11_[1] = V11[1];
}

bool Bundled::Quad::GetBilinearCoordinates(const v2d & pt, double coefficients[4])
{
    double a_x = V00_[0] - V01_[0] - V10_[0] + V11_[0];
    double b_x = -V00_[0] + V01_[0];
    double c_x = -V00_[0] + V10_[0];
    double d_x = V00_[0] - pt[0];

    double a_y = V00_[1]  - V01_[1] - V10_[1] + V11_[1];
    double b_y = -V00_[1] + V01_[1];
    double c_y = -V00_[1] + V10_[1];
    double d_y = V00_[1]  - pt[1];

    double bigA = -a_y*b_x + b_y*a_x;
    double bigB = -a_y*d_x - c_y*b_x + d_y*a_x +b_y*c_x;
    double bigC = -c_y*d_x + d_y*c_x;

    double tmp1 = -1;
    double tmp2 = -1;
    double tmp3 = -1;
    double tmp4 = -1;

    double k1 = 0,k2 = 0;

    if(bigB*bigB - 4*bigA*bigC >= 0.0)
    {
        if ( fabs(bigA) >= 0.000001)
        {
            tmp1 = ( -bigB + sqrt(bigB*bigB - 4*bigA*bigC) ) / ( 2*bigA );
            tmp2 = ( -bigB - sqrt(bigB*bigB - 4*bigA*bigC) ) / ( 2*bigA );
        }
        else
        {
            tmp1 = -bigC/bigB;
        }

        if ( tmp1 >= -0.999999 && tmp1 <= 1.000001)
        {
            tmp3 = -(b_y*tmp1 + d_y) / (a_y*tmp1 + c_y);
            tmp4 = -(b_x*tmp1 + d_x) / (a_x*tmp1 + c_x);
            if(tmp3 >= -0.999999 && tmp3 <= 1.000001)
            {
                k1 = tmp1;
                k2 = tmp3;
            }
            else if (tmp4 >= -0.999999 && tmp4 <= 1.000001)
            {
                k1 = tmp1;
                k2 = tmp4;
            }
        }

        if ( tmp2 >= -0.999999 && tmp2 <= 1.000001)
        {
            if ( tmp3 >= -0.999999 && tmp3 <= 1.000001)
            {
                k1 = tmp2;
                k2 = tmp3;
            }
            else if ( tmp4 >= -0.999999 && tmp4 <= 1.000001)
            {
                k1 = tmp2;
                k2 = tmp4;
            }
        }
    }


    if ( k1>=-0.999999 && k1<=1.000001 && k2>=-0.999999 && k2<=1.000001)
    {
        coefficients[0] = (1.0-k1)*(1.0-k2);
        coefficients[1] = k1*(1.0-k2);
        coefficients[2] = (1.0-k1)*k2;
        coefficients[3] = k1*k2;

        return true;
    }
    else
    {
        return false;
    }
}

Bundled::Quad & Bundled::Quad::operator=(const Bundled::Quad & qd)
{
    this->V00_ = qd.V00_;
    this->V01_ = qd.V01_;
    this->V10_ = qd.V10_;
    this->V11_ = qd.V11_;

    return (*this);
}



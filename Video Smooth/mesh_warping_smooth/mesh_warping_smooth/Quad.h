#ifndef QUAD_H_INCLUDED
#define QUAD_H_INCLUDED

#include "common.h"

namespace Bundled
{
    class Quad
    {
    public:
        Quad(const v2d & V00, const v2d & V01, const v2d & V10, const v2d & V11);

        Quad & operator=(const Quad & qd);

        bool GetBilinearCoordinates(const v2d & pt, double coefficients[4]);

    //private:
        v2d V00_, V01_, V10_, V11_;
    };
}

#endif // QUAD_H_INCLUDED

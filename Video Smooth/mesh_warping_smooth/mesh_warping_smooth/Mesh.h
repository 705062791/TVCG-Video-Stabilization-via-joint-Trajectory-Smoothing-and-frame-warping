#ifndef MESH_H_INCLUDED
#define MESH_H_INCLUDED

#include "common.h"
#include "Quad.h"

namespace Bundled
{
    class Mesh
    {
    public:
        Mesh(int rows, int cols, double qW, double qH);

        int             GetMeshWidth    () {return meshWidth;}
        int             GetMeshHeight   () {return meshHeight;}

        v2d             GetVertex       (int i, int j);
        void            SetVertex       (int i, int j, const v2d & pt);

        Bundled::Quad   getQuad         (int i, int j);

    private:
        int             imgRows, imgCols;
        int             meshWidth, meshHeight;
        double          quadWidth, quadHeight;

        TensorDou       xMat, yMat;
    };
}

#endif // MESH_H_INCLUDED

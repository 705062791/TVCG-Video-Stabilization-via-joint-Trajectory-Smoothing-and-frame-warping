#include <vector>
#include "Mesh.h"
#include "common.h"

Bundled::Mesh::Mesh(int rows, int cols, double qW, double qH)
{
    imgRows     =   rows;
    imgCols     =   cols;
    quadWidth   =   qW;
    quadHeight  =   qH;

    double x    = 0, y = 0; // 0
    int i       = 0, j = 0;

    int ww  = (int)(imgCols / quadWidth) + 10;
    int hh  = (int)(imgRows / quadHeight) + 10;

    double * xset   = new double[ww];
    double * yset   = new double[hh];

    while(imgCols - 1 - x > 0.5 * quadWidth)
    {
        xset[i] = x;
        x       = x + quadWidth;
        i       = i + 1;
    }
    xset[i] = imgCols - 1; // -1

    while(imgRows - 1 - y > 0.5 * quadHeight)
    {
        yset[j] = y;
        y       = y + quadHeight;
        j       = j + 1;
    }
    yset[j] = imgRows - 1; // -1

    meshWidth   = i+1;                          // the number of vertex of a horizontal line
    meshHeight  = j+1;


    xMat.set(meshWidth, meshHeight);
    yMat.set(meshWidth, meshHeight);

    for(i = 0;i < meshHeight; i++)
    {
        for(j = 0;j < meshWidth; j++)
        {
            xMat.at(i,j) = xset[j];
            yMat.at(i,j) = yset[i];
        }
    }

    delete [] xset;
    delete [] yset;
}

v2d Bundled::Mesh::GetVertex(int i, int j)
{
    double x = xMat.at(i,j);
    double y = yMat.at(i,j);
    return _v2d_(x,y);
}

void Bundled::Mesh::SetVertex(int i, int j, const v2d & pt)
{
    xMat.at(i,j) = pt[0];
    yMat.at(i,j) = pt[1];
}

Bundled::Quad Bundled::Mesh::getQuad(int i, int j)
{
    v2d V00,V01,V10,V11;
    V00 = GetVertex(i,j);
    V01 = GetVertex(i,j+1);
    V10 = GetVertex(i+1,j);
    V11 = GetVertex(i+1,j+1);

    return Bundled::Quad(V00,V01,V10,V11);
}

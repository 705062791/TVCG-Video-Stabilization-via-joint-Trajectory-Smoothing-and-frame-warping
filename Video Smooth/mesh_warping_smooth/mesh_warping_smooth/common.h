#pragma once
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <opencv2\opencv.hpp>
#include <string>
#include <libiv.h>
using namespace std;
using namespace cv;


typedef struct point_path
{
	vector<v2d> PointPath;
	int start_frame;
	int end_frame;
};


//typedef Eigen::Triplet<double> T_d;
//typedef vector<T_d> Vec_Triplet;
typedef vector<point_path> set_of_path;




namespace Bundled
{

	typedef LibIV::Memory::Array::FastArray1D<int>              ArrayInt;
	typedef LibIV::Memory::Array::FastArray1D<double>           ArrayDou;
	typedef LibIV::Memory::Array::FastArray1D<v2d>              ArrayV2d;
	typedef LibIV::Memory::Array::FastArray1D<v3d>              ArrayV3d;

	typedef LibIV::Memory::Array::FastArray2D<int>              TensorInt;
	typedef LibIV::Memory::Array::FastArray2D<double>           TensorDou;
	typedef LibIV::Memory::Array::FastArray2D<v3d>              TensorV3d;

	typedef LibIV::Memory::Array::FastArray2D<cv::Mat>          TensorMatrix3d;
	typedef LibIV::Memory::Array::FastArray1D<TensorMatrix3d>   CubeMatrix3d;

	typedef LibIV::Memory::Array::FastArray2D<TensorDou>        FourDim;
}
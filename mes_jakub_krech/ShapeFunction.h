#pragma once

#ifndef EIGEN
#define EIGEN
#include <Eigen>
#endif

class ShapeFunction
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	Eigen::Vector4d N;

	ShapeFunction(double, double, double, double);
	~ShapeFunction();
};


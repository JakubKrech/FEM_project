#pragma once

#ifndef EIGEN
#define EIGEN
#include <Eigen>
#endif

class InterpolationCoordinates
{
public:
	double x, y;

	InterpolationCoordinates(double, double);
	~InterpolationCoordinates();
};


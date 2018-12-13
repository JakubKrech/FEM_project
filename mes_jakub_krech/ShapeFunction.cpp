#include "ShapeFunction.h"


ShapeFunction::ShapeFunction(double n1, double n2, double n3, double n4)
{
	N(0) = n1;
	N(1) = n2;
	N(2) = n3;
	N(3) = n4;
}


ShapeFunction::~ShapeFunction()
{
}

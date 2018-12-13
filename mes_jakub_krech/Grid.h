#pragma once
#include "Element.h"
#include <vector>


#ifndef EIGEN
#define EIGEN
#include <Eigen>
#endif

#define matrix_size 4

class Grid
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	std::vector<Node*> nodes;
	std::vector<Element*> elements;

	std::vector<ShapeFunction*> ShapeFunc;
	Eigen::Matrix<double, matrix_size, matrix_size> ShFnDerivativesToKsi;
	Eigen::Matrix<double, matrix_size, matrix_size> ShFnDerivativesToEta;

	static Grid* createGridFromDataFile(std::string fileName);

	Grid(double H, double L, int nH, int nL);
	~Grid();

	std::vector<ShapeFunction*> calculateShapeFunctions(Node*, Node*, Node*, Node*);
	Eigen::Matrix<double, matrix_size, matrix_size> calculateShapeFunctionDerivativesToKsi(Node*, Node*, Node*, Node*);
	Eigen::Matrix<double, matrix_size, matrix_size> calculateShapeFunctionDerivativesToEta(Node*, Node*, Node*, Node*);

	void print_nodes();
	void print_elements();
	void print_element_by_index(int);
	void print_shape_function_derivatives();
};
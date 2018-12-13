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

	Eigen::Matrix<double, matrix_size*matrix_size, matrix_size*matrix_size> Global_Matrix_H;
	Eigen::Matrix<double, matrix_size*matrix_size, matrix_size*matrix_size> Global_Matrix_C;
	Eigen::Matrix<double, 16, 1> Global_Matrix_P;

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

	void calculate_global_matrix_H();
	void print_global_matrix_H();

	void calculate_global_matrix_C();
	void print_global_matrix_C();

	void calculate_global_matrix_P();
	void print_global_matrix_P();

	void calculate_next_iterations(int);
};
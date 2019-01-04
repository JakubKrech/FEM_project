#pragma once
#include "Node.h"
#include "ShapeFunction.h"
#include "InterpolationCoordinates.h"
#include <vector>

#ifndef EIGEN
#define EIGEN
#include <Eigen>
#endif

#define matrix_size 4

class Element
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	static int global_id;
	int id;

	std::vector<Node*> Nodes;
	std::vector<ShapeFunction*> ShapeFunctions;
	//std::vector<InterpolationCoordinates*> InterpolCoord;

	Eigen::Matrix<double, matrix_size, matrix_size> ShapeFunctionsDerivativesToKsi;
	Eigen::Matrix<double, matrix_size, matrix_size> ShapeFunctionsDerivativesToEta;

	Eigen::Matrix<double, matrix_size, matrix_size> JacobianMatrix;
	Eigen::Vector4d JacobianDet;
	Eigen::Matrix<double, matrix_size, matrix_size> InverseJacobianMatrix;

	Eigen::Matrix<double, matrix_size, matrix_size> dN_dx;
	Eigen::Matrix<double, matrix_size, matrix_size> dN_dy;

	Eigen::Matrix<double, matrix_size, matrix_size> dNdx_dNdx_T_1pc;
	Eigen::Matrix<double, matrix_size, matrix_size> dNdx_dNdx_T_2pc;
	Eigen::Matrix<double, matrix_size, matrix_size> dNdx_dNdx_T_3pc;
	Eigen::Matrix<double, matrix_size, matrix_size> dNdx_dNdx_T_4pc;

	Eigen::Matrix<double, matrix_size, matrix_size> dNdy_dNdy_T_1pc;
	Eigen::Matrix<double, matrix_size, matrix_size> dNdy_dNdy_T_2pc;
	Eigen::Matrix<double, matrix_size, matrix_size> dNdy_dNdy_T_3pc;
	Eigen::Matrix<double, matrix_size, matrix_size> dNdy_dNdy_T_4pc;

	Eigen::Matrix<double, matrix_size, matrix_size> K_transX_transY_det_1pc;
	Eigen::Matrix<double, matrix_size, matrix_size> K_transX_transY_det_2pc;
	Eigen::Matrix<double, matrix_size, matrix_size> K_transX_transY_det_3pc;
	Eigen::Matrix<double, matrix_size, matrix_size> K_transX_transY_det_4pc;

	Eigen::Matrix<double, matrix_size, matrix_size> Matrix_H;

	Eigen::Matrix<double, matrix_size, matrix_size> c_ro_NN_det_1pc;
	Eigen::Matrix<double, matrix_size, matrix_size> c_ro_NN_det_2pc;
	Eigen::Matrix<double, matrix_size, matrix_size> c_ro_NN_det_3pc;
	Eigen::Matrix<double, matrix_size, matrix_size> c_ro_NN_det_4pc;

	Eigen::Matrix<double, matrix_size, matrix_size> Matrix_C;

	Eigen::Vector4d pow1_pc1_N;
	Eigen::Vector4d pow1_pc2_N;

	Eigen::Vector4d pow2_pc1_N;
	Eigen::Vector4d pow2_pc2_N;

	Eigen::Vector4d pow3_pc1_N;
	Eigen::Vector4d pow3_pc2_N;

	Eigen::Vector4d pow4_pc1_N;
	Eigen::Vector4d pow4_pc2_N;

	Eigen::Matrix<double, matrix_size, matrix_size> pow1_pc1;
	Eigen::Matrix<double, matrix_size, matrix_size> pow1_pc2;
	Eigen::Matrix<double, matrix_size, matrix_size> pow2_pc1;
	Eigen::Matrix<double, matrix_size, matrix_size> pow2_pc2;
	Eigen::Matrix<double, matrix_size, matrix_size> pow3_pc1;
	Eigen::Matrix<double, matrix_size, matrix_size> pow3_pc2;
	Eigen::Matrix<double, matrix_size, matrix_size> pow4_pc1;
	Eigen::Matrix<double, matrix_size, matrix_size> pow4_pc2;

	Eigen::Matrix<double, matrix_size, matrix_size> pow1_sum;
	Eigen::Matrix<double, matrix_size, matrix_size> pow2_sum;
	Eigen::Matrix<double, matrix_size, matrix_size> pow3_sum;
	Eigen::Matrix<double, matrix_size, matrix_size> pow4_sum;

	Eigen::Vector4d pow1_P;
	Eigen::Vector4d pow2_P;
	Eigen::Vector4d pow3_P;
	Eigen::Vector4d pow4_P;

	Eigen::Vector4d Vector_P;

	Eigen::Matrix<double, matrix_size, matrix_size> Matrix_H_BC;
	Eigen::Matrix<double, matrix_size, matrix_size> Matrix_H_Final;


	Element(Node *a, Node *b, Node *c, Node *d,
		std::vector<ShapeFunction*>&,
		Eigen::Matrix<double, matrix_size, matrix_size>&,
		Eigen::Matrix<double, matrix_size, matrix_size>&);
	~Element();

	//void calculateInterpolationForEachNode();

	void calculateJacobianMatrix();
	void calculateJacobianDet();
	void calculateInverseJacobianMatrix();

	void calculate_dn_dx();
	void calculate_dn_dy();
	void calculate_transposed();
	void multiply_transponed_by_det();
	void calculate_K_transX_transY_det();
	void calculate_Matrix_H();

	void calculate_c_ro_NN_det();
	void calculate_Matrix_C();

	void calculate_pow_pc_N();
	void calculate_pow_pc();
	void calculate_pow_sum_and_pow_P();
	void calculate_Vector_P();

	void calculate_Matrix_H_BC();
	void calculate_Matrix_H_Final();

	void print();
};

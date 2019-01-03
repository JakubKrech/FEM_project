#include "Element.h"
#include <iostream>
#include <cmath>
#include <iomanip>

#define  shape_and_interpolation_debug_logging 0
#define      jacobian_matrix_one_debug_logging 0
#define             jacobian_det_debug_logging 0
#define          jacobian_matrix_debug_logging 0

#define                    dn_dx_debug_logging 0
#define                    dn_dy_debug_logging 0
#define               transponed_debug_logging 0
#define      K_transX_transY_det_debug_logging 0
#define                 matrix_H_debug_logging 0

#define              c_ro_NN_det_debug_logging 0
#define                 matrix_C_debug_logging 0

#define       calculate_pow_pc_N_debug_logging 0
#define         calculate_pow_pc_debug_logging 0
#define        calculate_pow_sum_debug_logging 0
#define          calculate_pow_P_debug_logging 0
#define    calculate_Matrix_H_BC_debug_logging 0
#define       calculate_Matrix_P_debug_logging 0
#define calculate_Matrix_H_Final_debug_logging 0

#define ksi_global 1/sqrt(3)
#define eta_global 1/sqrt(3)
#define conductivity 25 //przewodnictwo 30
#define specific_heat 700 //cieplo wlasciwe
#define density 7800 //gestosc/density
#define convection 300 //konwekcja (alfa) 25
#define ambient_temp 1200

Element::Element(
	Node *a,
	Node *b,
	Node *c,
	Node *d,
	std::vector<ShapeFunction*> &sf,
	Eigen::Matrix<double, matrix_size, matrix_size> &toKsi,
	Eigen::Matrix<double, matrix_size, matrix_size> &toEta)
{
	Nodes.push_back(a);
	Nodes.push_back(b);
	Nodes.push_back(c);
	Nodes.push_back(d);

	ShapeFunctions = sf;
	ShapeFunctionsDerivativesToKsi = toKsi;
	ShapeFunctionsDerivativesToEta = toEta;

	id = global_id++;

	//calculateInterpolationForEachNode();
	calculateJacobianMatrix();
	calculateJacobianDet();
	calculateInverseJacobianMatrix();

	calculate_dn_dx();
	calculate_dn_dy();
	calculate_transposed();
	multiply_transponed_by_det();
	calculate_K_transX_transY_det();
	calculate_Matrix_H();

	calculate_c_ro_NN_det();
	calculate_Matrix_C();

	calculate_pow_pc_N();
	calculate_pow_pc();
	calculate_pow_sum_and_pow_P();
	calculate_Matrix_P();
	calculate_Matrix_H_BC();
	
	calculate_Matrix_H_Final();
}

Element::~Element()
{
	//std::cout << "Element[" << id << "] deleted" << std::endl;
}

int Element::global_id = 1;

void Element::print()
{
	std::cout << "Element[" << id << "] : (" <<
		Nodes[0]->id << "," <<
		Nodes[1]->id << "," <<
		Nodes[2]->id << "," <<
		Nodes[3]->id << ")  BC: ";
	if (Nodes[0]->BC && Nodes[1]->BC) std::cout << "pow1 ";
	if (Nodes[1]->BC && Nodes[2]->BC) std::cout << "pow2 ";
	if (Nodes[2]->BC && Nodes[3]->BC) std::cout << "pow3 ";
	if (Nodes[3]->BC && Nodes[0]->BC) std::cout << "pow4 ";
	std::cout << "\n\n";

	/*if (shape_and_interpolation_debug_logging) {	
		for (int i = 0; i < ShapeFunctions.size(); i++) {
			std::cout << "  --Node #" << Nodes[i]->id << "--\n" <<
				"   N1: " << ShapeFunctions[i]->N(0) <<
				" N2: " << ShapeFunctions[i]->N(1) <<
				" N3: " << ShapeFunctions[i]->N(2) <<
				" N4: " << ShapeFunctions[i]->N(3) <<
				"\n   InterpolationCoordinates: (" << InterpolCoord[i]->x << ", " << InterpolCoord[i]->y << ")\n\n";
		}
	}*/

	if (jacobian_matrix_one_debug_logging) {
		std::cout << "  --Jacobian Matrix-- \n" <<
			std::fixed << JacobianMatrix << "\n\n";
	}
	
	if (jacobian_det_debug_logging)
	{
		std::cout << "  --Jacobian Det-- \n";

		for (int i = 0; i < matrix_size; i++)
		{
			std::cout << "   " << std::setw(5) << JacobianDet[i] << " ";
		}
		std::cout << "\n\n";
	}
	
	if (jacobian_matrix_debug_logging)
	{
		std::cout << "  --Inverse Jacobian Matrix-- \n" << 
			std::fixed << InverseJacobianMatrix << "\n\n";
	}
	
	if (dn_dx_debug_logging)
	{
		std::cout << "  --dN/dx-- \n" <<
			std::fixed << dN_dx << "\n\n";
	}

	if (dn_dy_debug_logging)
	{
		std::cout << "  --dN/dy-- \n" <<
			std::fixed << dN_dy << "\n\n";
	}
	
	if (transponed_debug_logging)
	{
		std::cout << " -- TRANSPONED (Multiplied by det or not) --\n\n"; // multiplied by det or not
		
		std::cout << "  --dNdx_dNdx_T_1pc-- \n" <<
			std::fixed << dNdx_dNdx_T_1pc << "\n\n";

		std::cout << "  --dNdx_dNdx_T_2pc-- \n" <<
			std::fixed << dNdx_dNdx_T_2pc << "\n\n";

		std::cout << "  --dNdx_dNdx_T_3pc-- \n" <<
			std::fixed << dNdx_dNdx_T_3pc << "\n\n";

		std::cout << "  --dNdx_dNdx_T_4pc-- \n" <<
			std::fixed << dNdx_dNdx_T_4pc << "\n\n";


		std::cout << "  --dNdy_dNdy_T_1pc-- \n" <<
			std::fixed << dNdy_dNdy_T_1pc << "\n\n";

		std::cout << "  --dNdy_dNdy_T_2pc-- \n" <<
			std::fixed << dNdy_dNdy_T_2pc << "\n\n";

		std::cout << "  --dNdy_dNdy_T_3pc-- \n" <<
			std::fixed << dNdy_dNdy_T_3pc << "\n\n";

		std::cout << "  --dNdy_dNdy_T_4pc-- \n" <<
			std::fixed << dNdy_dNdy_T_4pc << "\n\n";
	}
	
	if (K_transX_transY_det_debug_logging)
	{
		std::cout << " -- K_transX_transY_det --\n\n";

		std::cout << "  --K_transX_transY_det_1pc-- \n" <<
			std::fixed << K_transX_transY_det_1pc << "\n\n";

		std::cout << "  --K_transX_transY_det_2pc-- \n" <<
			std::fixed << K_transX_transY_det_2pc << "\n\n";

		std::cout << "  --K_transX_transY_det_3pc-- \n" <<
			std::fixed << K_transX_transY_det_3pc << "\n\n";

		std::cout << "  --K_transX_transY_det_4pc-- \n" <<
			std::fixed << K_transX_transY_det_4pc << "\n\n";
	}
	
	if (matrix_H_debug_logging)
	{
		std::cout << "-------------------Matrix H------------------- \n" <<
			std::fixed << Matrix_H <<
			"\n----------------------------------------------\n\n";
	}

	if(c_ro_NN_det_debug_logging)
	{
		std::cout << " -- c_ro_NN_det --\n\n";

		std::cout << "  --c_ro_NN_det_1pc-- \n" <<
			std::fixed << c_ro_NN_det_1pc << "\n\n";

		std::cout << "  --c_ro_NN_det_2pc-- \n" <<
			std::fixed << c_ro_NN_det_2pc << "\n\n";

		std::cout << "  --c_ro_NN_det_3pc-- \n" <<
			std::fixed << c_ro_NN_det_3pc << "\n\n";

		std::cout << "  --c_ro_NN_det_4pc-- \n" <<
			std::fixed << c_ro_NN_det_4pc << "\n\n";
	}

	if (matrix_C_debug_logging)
	{
		std::cout << "-------------------Matrix C-------------------\n" <<
			std::fixed << Matrix_C <<
			"\n----------------------------------------------\n\n";
	}

	if (calculate_pow_pc_N_debug_logging)
	{
		std::cout << " -- pow_pc_N --\n\n";

		std::cout << "  --pow1_pc1_N-- \n" <<
			std::fixed << pow1_pc1_N << "\n\n";

		std::cout << "  --pow1_pc2_N-- \n" <<
			std::fixed << pow1_pc2_N << "\n\n";

		std::cout << "  --pow2_pc1_N-- \n" <<
			std::fixed << pow2_pc1_N << "\n\n";

		std::cout << "  --pow2_pc2_N-- \n" <<
			std::fixed << pow2_pc2_N << "\n\n";

		std::cout << "  --pow3_pc1_N-- \n" <<
			std::fixed << pow3_pc1_N << "\n\n";

		std::cout << "  --pow3_pc2_N-- \n" <<
			std::fixed << pow3_pc2_N << "\n\n";

		std::cout << "  --pow4_pc1_N-- \n" <<
			std::fixed << pow4_pc1_N << "\n\n";

		std::cout << "  --pow4_pc2_N-- \n" <<
			std::fixed << pow4_pc2_N << "\n\n";
	}

	if (calculate_pow_pc_debug_logging)
	{
		std::cout << " -- pow_pc --\n\n";

		std::cout << "  --pow1_pc1-- \n" <<
			std::fixed << pow1_pc1 << "\n\n";

		std::cout << "  --pow1_pc2-- \n" <<
			std::fixed << pow1_pc2 << "\n\n";

		std::cout << "  --pow2_pc1-- \n" <<
			std::fixed << pow2_pc1 << "\n\n";

		std::cout << "  --pow2_pc2-- \n" <<
			std::fixed << pow2_pc2 << "\n\n";

		std::cout << "  --pow3_pc1-- \n" <<
			std::fixed << pow3_pc1 << "\n\n";

		std::cout << "  --pow3_pc2-- \n" <<
			std::fixed << pow3_pc2 << "\n\n";

		std::cout << "  --pow4_pc1-- \n" <<
			std::fixed << pow4_pc1 << "\n\n";

		std::cout << "  --pow4_pc2-- \n" <<
			std::fixed << pow4_pc2 << "\n\n";
	}

	

	if (calculate_pow_sum_debug_logging)
	{
		std::cout << " -- pow_sum --\n\n";

		std::cout << "  --pow1_sum-- \n" <<
			std::fixed << pow1_sum << "\n\n";

		std::cout << "  --pow2_sum-- \n" <<
			std::fixed << pow2_sum << "\n\n";

		std::cout << "  --pow3_sum-- \n" <<
			std::fixed << pow3_sum << "\n\n";

		std::cout << "  --pow4_sum-- \n" <<
			std::fixed << pow4_sum << "\n\n";
	}

	if (calculate_Matrix_H_BC_debug_logging)
	{
		std::cout << "---------Matrix H Boundary Conditions--------- \n" <<
			std::fixed << Matrix_H_BC <<
			"\n----------------------------------------------\n\n";
	}

	if (calculate_pow_P_debug_logging)
	{
		std::cout << " -- pow_P --\n\n";

		std::cout << "  --pow1_P-- \n" <<
			std::fixed << pow1_P << "\n\n";

		std::cout << "  --pow2_P-- \n" <<
			std::fixed << pow2_P << "\n\n";

		std::cout << "  --pow3_P-- \n" <<
			std::fixed << pow3_P << "\n\n";

		std::cout << "  --pow4_P-- \n" <<
			std::fixed << pow4_P << "\n\n";
	}


	if (calculate_Matrix_P_debug_logging)
	{
		std::cout << "-------------------Matrix P-------------------\n" <<
			std::fixed << Matrix_P <<
			"\n----------------------------------------------\n\n";
	}

	if (calculate_Matrix_H_Final_debug_logging)
	{
		std::cout << "----------------Matrix H Final---------------- \n" <<
			std::fixed << Matrix_H_Final <<
			"\n----------------------------------------------\n\n";
	}
}

/*void Element::calculateInterpolationForEachNode()
{
	double Xp = ShapeFunctions[0]->N(0) * Nodes[0]->x +
		ShapeFunctions[0]->N(1) * Nodes[1]->x +
		ShapeFunctions[0]->N(2) * Nodes[2]->x +
		ShapeFunctions[0]->N(3) * Nodes[3]->x;
	double Yp = ShapeFunctions[0]->N(0) * Nodes[0]->y +
		ShapeFunctions[0]->N(1) * Nodes[1]->y +
		ShapeFunctions[0]->N(2) * Nodes[2]->y +
		ShapeFunctions[0]->N(3) * Nodes[3]->y;
	InterpolationCoordinates *intcoor0 = new InterpolationCoordinates(Xp, Yp);
	InterpolCoord.push_back(intcoor0);

	Xp = ShapeFunctions[1]->N(0) * Nodes[0]->x +
		ShapeFunctions[1]->N(1) * Nodes[1]->x +
		ShapeFunctions[1]->N(2) * Nodes[2]->x +
		ShapeFunctions[1]->N(3) * Nodes[3]->x;
	Yp = ShapeFunctions[1]->N(0) * Nodes[0]->y +
		ShapeFunctions[1]->N(1) * Nodes[1]->y +
		ShapeFunctions[1]->N(2) * Nodes[2]->y +
		ShapeFunctions[1]->N(3) * Nodes[3]->y;
	InterpolationCoordinates *intcoor1 = new InterpolationCoordinates(Xp, Yp);
	InterpolCoord.push_back(intcoor1);

	Xp = ShapeFunctions[2]->N(0) * Nodes[0]->x +
		ShapeFunctions[2]->N(1) * Nodes[1]->x +
		ShapeFunctions[2]->N(2) * Nodes[2]->x +
		ShapeFunctions[2]->N(3) * Nodes[3]->x;
	Yp = ShapeFunctions[2]->N(0) * Nodes[0]->y +
		ShapeFunctions[2]->N(1) * Nodes[1]->y +
		ShapeFunctions[2]->N(2) * Nodes[2]->y +
		ShapeFunctions[2]->N(3) * Nodes[3]->y;
	InterpolationCoordinates *intcoor2 = new InterpolationCoordinates(Xp, Yp);
	InterpolCoord.push_back(intcoor2);

	Xp = ShapeFunctions[3]->N(0) * Nodes[0]->x +
		ShapeFunctions[3]->N(1) * Nodes[1]->x +
		ShapeFunctions[3]->N(2) * Nodes[2]->x +
		ShapeFunctions[3]->N(3) * Nodes[3]->x;
	Yp = ShapeFunctions[3]->N(0) * Nodes[0]->y +
		ShapeFunctions[3]->N(1) * Nodes[1]->y +
		ShapeFunctions[3]->N(2) * Nodes[2]->y +
		ShapeFunctions[3]->N(3) * Nodes[3]->y;
	InterpolationCoordinates *intcoor3 = new InterpolationCoordinates(Xp, Yp);
	InterpolCoord.push_back(intcoor3);
}*/

void Element::calculateJacobianMatrix()
{
	for (int i = 0; i < matrix_size; i++)
	{
		JacobianMatrix(0,i) = ShapeFunctionsDerivativesToKsi(0, i) * Nodes[0]->x +
			ShapeFunctionsDerivativesToKsi(1, i) * Nodes[1]->x +
			ShapeFunctionsDerivativesToKsi(2, i) * Nodes[2]->x +
			ShapeFunctionsDerivativesToKsi(3, i) * Nodes[3]->x;

		JacobianMatrix(1, i) = ShapeFunctionsDerivativesToKsi(0, i) * Nodes[0]->y +
			ShapeFunctionsDerivativesToKsi(1, i) * Nodes[1]->y +
			ShapeFunctionsDerivativesToKsi(2, i) * Nodes[2]->y +
			ShapeFunctionsDerivativesToKsi(3, i) * Nodes[3]->y;

		JacobianMatrix(2, i) = ShapeFunctionsDerivativesToEta(0, i) * Nodes[0]->x +
			ShapeFunctionsDerivativesToEta(1, i) * Nodes[1]->x +
			ShapeFunctionsDerivativesToEta(2, i) * Nodes[2]->x +
			ShapeFunctionsDerivativesToEta(3, i) * Nodes[3]->x;

		JacobianMatrix(3, i) = ShapeFunctionsDerivativesToEta(0, i) * Nodes[0]->y +
			ShapeFunctionsDerivativesToEta(1, i) * Nodes[1]->y +
			ShapeFunctionsDerivativesToEta(2, i) * Nodes[2]->y +
			ShapeFunctionsDerivativesToEta(3, i) * Nodes[3]->y;
	}
}

void Element::calculateJacobianDet()
{
	for (int i = 0; i < matrix_size; i++)
	{
		JacobianDet(i) = JacobianMatrix(0, i) * JacobianMatrix(3, i) -
			JacobianMatrix(1, i) * JacobianMatrix(2, i);
	}
	
}

void Element::calculateInverseJacobianMatrix()
{
	for (int i = 0; i < matrix_size; i++)
	{
		InverseJacobianMatrix(0,i) = JacobianMatrix(3,i) / JacobianDet(i);
		InverseJacobianMatrix(1,i) = -(JacobianMatrix(2,i) / JacobianDet(i));
		InverseJacobianMatrix(2,i) = JacobianMatrix(1,i) / JacobianDet(i);
		InverseJacobianMatrix(3,i) = JacobianMatrix(0,i) / JacobianDet(i);   // 
	}
}

void Element::calculate_dn_dx()
{
	for (int i = 0; i < matrix_size; i++)
	{
		dN_dx(0,i) = InverseJacobianMatrix(0,0) * ShapeFunctionsDerivativesToKsi(i,0) +
			InverseJacobianMatrix(1,0) * ShapeFunctionsDerivativesToEta(i,0);

		dN_dx(1,i) = InverseJacobianMatrix(0,1) * ShapeFunctionsDerivativesToKsi(i,1) +
			InverseJacobianMatrix(1,1) * ShapeFunctionsDerivativesToEta(i,1);

		dN_dx(2,i) = InverseJacobianMatrix(0,2) * ShapeFunctionsDerivativesToKsi(i,2) +
			InverseJacobianMatrix(1,2) * ShapeFunctionsDerivativesToEta(i,2);

		dN_dx(3,i) = InverseJacobianMatrix(0,3) * ShapeFunctionsDerivativesToKsi(i,3) +
			InverseJacobianMatrix(1,3) * ShapeFunctionsDerivativesToEta(i,3);
	}
}

void Element::calculate_dn_dy()
{
	for (int i = 0; i < matrix_size; i++)
	{
		dN_dy(0,i) = InverseJacobianMatrix(2,0) * ShapeFunctionsDerivativesToKsi(i,0) +
			InverseJacobianMatrix(3,0) * ShapeFunctionsDerivativesToEta(i,0);

		dN_dy(1,i) = InverseJacobianMatrix(2,1) * ShapeFunctionsDerivativesToKsi(i,1) +
			InverseJacobianMatrix(3,1) * ShapeFunctionsDerivativesToEta(i,1);

		dN_dy(2,i) = InverseJacobianMatrix(2,2) * ShapeFunctionsDerivativesToKsi(i,2) +
			InverseJacobianMatrix(3,2) * ShapeFunctionsDerivativesToEta(i,2);

		dN_dy(3,i) = InverseJacobianMatrix(2,3) * ShapeFunctionsDerivativesToKsi(i,3) +
			InverseJacobianMatrix(3,3) * ShapeFunctionsDerivativesToEta(i,3);
	}
}

void Element::calculate_transposed()
{
	for (int i = 0; i < matrix_size; i++)
	{
		for (int j = 0; j < matrix_size; j++)
		{
			dNdx_dNdx_T_1pc(i,j) = dN_dx(0,i) * dN_dx(0,j);
			dNdx_dNdx_T_2pc(i,j) = dN_dx(1,i) * dN_dx(1,j);
			dNdx_dNdx_T_3pc(i,j) = dN_dx(2,i) * dN_dx(2,j);
			dNdx_dNdx_T_4pc(i,j) = dN_dx(3,i) * dN_dx(3,j);

			dNdy_dNdy_T_1pc(i,j) = dN_dy(0,i) * dN_dy(0,j);
			dNdy_dNdy_T_2pc(i,j) = dN_dy(1,i) * dN_dy(1,j);
			dNdy_dNdy_T_3pc(i,j) = dN_dy(2,i) * dN_dy(2,j);
			dNdy_dNdy_T_4pc(i,j) = dN_dy(3,i) * dN_dy(3,j);
		}
	}
}

void Element::multiply_transponed_by_det()
{
	for (int i = 0; i < matrix_size; i++)
	{
		for (int j = 0; j < matrix_size; j++)
		{
			dNdx_dNdx_T_1pc(i,j) *= JacobianDet(i);
			dNdx_dNdx_T_2pc(i,j) *= JacobianDet(i);
			dNdx_dNdx_T_3pc(i,j) *= JacobianDet(i);
			dNdx_dNdx_T_4pc(i,j) *= JacobianDet(i);

			dNdy_dNdy_T_1pc(i,j) *= JacobianDet(i);
			dNdy_dNdy_T_2pc(i,j) *= JacobianDet(i);
			dNdy_dNdy_T_3pc(i,j) *= JacobianDet(i);
			dNdy_dNdy_T_4pc(i,j) *= JacobianDet(i);
		}
	}
}

void Element::calculate_K_transX_transY_det()
{
	for (int i = 0; i < matrix_size; i++)
	{
		for (int j = 0; j < matrix_size; j++)
		{
			K_transX_transY_det_1pc(i,j) = conductivity *(dNdx_dNdx_T_1pc(i,j) + dNdy_dNdy_T_1pc(i,j));
			K_transX_transY_det_2pc(i,j) = conductivity *(dNdx_dNdx_T_2pc(i,j) + dNdy_dNdy_T_2pc(i,j));
			K_transX_transY_det_3pc(i,j) = conductivity *(dNdx_dNdx_T_3pc(i,j) + dNdy_dNdy_T_3pc(i,j));
			K_transX_transY_det_4pc(i,j) = conductivity *(dNdx_dNdx_T_4pc(i,j) + dNdy_dNdy_T_4pc(i,j));
		}
	}
}

void Element::calculate_Matrix_H()
{
	for (int i = 0; i < matrix_size; i++)
	{
		for (int j = 0; j < matrix_size; j++)
		{
			Matrix_H(i,j) = (K_transX_transY_det_1pc(i,j) +
				K_transX_transY_det_2pc(i,j) +
				K_transX_transY_det_3pc(i,j) +
				K_transX_transY_det_4pc(i,j));
		}
	}
}

void Element::calculate_c_ro_NN_det()
{
	for (int i = 0; i < matrix_size; i++)
	{
		for (int j = 0; j < matrix_size; j++)
		{
			c_ro_NN_det_1pc(i, j) = ShapeFunctions[0]->N(i) * ShapeFunctions[0]->N(j) * JacobianDet[0] * specific_heat * density;
			c_ro_NN_det_2pc(i, j) = ShapeFunctions[1]->N(i) * ShapeFunctions[1]->N(j) * JacobianDet[1] * specific_heat * density;
			c_ro_NN_det_3pc(i, j) = ShapeFunctions[2]->N(i) * ShapeFunctions[2]->N(j) * JacobianDet[2] * specific_heat * density;
			c_ro_NN_det_4pc(i, j) = ShapeFunctions[3]->N(i) * ShapeFunctions[3]->N(j) * JacobianDet[3] * specific_heat * density;
		}
	}
}

void Element::calculate_Matrix_C()
{
	for (int i = 0; i < matrix_size; i++)
	{
		for (int j = 0; j < matrix_size; j++)
		{
			Matrix_C(i, j) = c_ro_NN_det_1pc(i, j) + c_ro_NN_det_2pc(i, j) +
				c_ro_NN_det_3pc(i, j) + c_ro_NN_det_4pc(i, j);
		}
	}
}

void Element::calculate_pow_pc_N()
{
	double ksi, eta;
	if (Nodes[0]->BC == true && Nodes[1]->BC == true) {
		{
			ksi = -ksi_global;
			eta = -1;
			pow1_pc1_N(0) = 0.25 * (1 - ksi)*(1 - eta);
			pow1_pc1_N(1) = 0.25 * (1 + ksi)*(1 - eta);
			pow1_pc1_N(2) = 0.25 * (1 + ksi)*(1 + eta);
			pow1_pc1_N(3) = 0.25 * (1 - ksi)*(1 + eta);
		}
		{
			ksi = ksi_global;
			eta = -1;
			pow1_pc2_N(0) = 0.25 * (1 - ksi)*(1 - eta);
			pow1_pc2_N(1) = 0.25 * (1 + ksi)*(1 - eta);
			pow1_pc2_N(2) = 0.25 * (1 + ksi)*(1 + eta);
			pow1_pc2_N(3) = 0.25 * (1 - ksi)*(1 + eta);
		}
	}
	if (Nodes[1]->BC == true && Nodes[2]->BC == true) {
		{
			ksi = 1;
			eta = -eta_global;
			pow2_pc1_N(0) = 0.25 * (1 - ksi)*(1 - eta);
			pow2_pc1_N(1) = 0.25 * (1 + ksi)*(1 - eta);
			pow2_pc1_N(2) = 0.25 * (1 + ksi)*(1 + eta);
			pow2_pc1_N(3) = 0.25 * (1 - ksi)*(1 + eta);
		}
		{
			ksi = 1;
			eta = eta_global;
			pow2_pc2_N(0) = 0.25 * (1 - ksi)*(1 - eta);
			pow2_pc2_N(1) = 0.25 * (1 + ksi)*(1 - eta);
			pow2_pc2_N(2) = 0.25 * (1 + ksi)*(1 + eta);
			pow2_pc2_N(3) = 0.25 * (1 - ksi)*(1 + eta);
		}
	}
	if (Nodes[2]->BC == true && Nodes[3]->BC == true) {
		{
			ksi = ksi_global;
			eta = 1;
			pow3_pc1_N(0) = 0.25 * (1 - ksi)*(1 - eta);
			pow3_pc1_N(1) = 0.25 * (1 + ksi)*(1 - eta);
			pow3_pc1_N(2) = 0.25 * (1 + ksi)*(1 + eta);
			pow3_pc1_N(3) = 0.25 * (1 - ksi)*(1 + eta);
		}
		{
			ksi = -ksi_global;
			eta = 1;
			pow3_pc2_N(0) = 0.25 * (1 - ksi)*(1 - eta);
			pow3_pc2_N(1) = 0.25 * (1 + ksi)*(1 - eta);
			pow3_pc2_N(2) = 0.25 * (1 + ksi)*(1 + eta);
			pow3_pc2_N(3) = 0.25 * (1 - ksi)*(1 + eta);
		}
	}
	if (Nodes[3]->BC == true && Nodes[0]->BC == true) {
		{
			ksi = -1;
			eta = eta_global;
			pow4_pc1_N(0) = 0.25 * (1 - ksi)*(1 - eta);
			pow4_pc1_N(1) = 0.25 * (1 + ksi)*(1 - eta);
			pow4_pc1_N(2) = 0.25 * (1 + ksi)*(1 + eta);
			pow4_pc1_N(3) = 0.25 * (1 - ksi)*(1 + eta);
		}
		{
			ksi = -1;
			eta = -eta_global;
			pow4_pc2_N(0) = 0.25 * (1 - ksi)*(1 - eta);
			pow4_pc2_N(1) = 0.25 * (1 + ksi)*(1 - eta);
			pow4_pc2_N(2) = 0.25 * (1 + ksi)*(1 + eta);
			pow4_pc2_N(3) = 0.25 * (1 - ksi)*(1 + eta);
		}
	}
}

void Element::calculate_pow_pc()
{
	for (int i = 0; i < matrix_size; i++)
	{
		for (int j = 0; j < matrix_size; j++)
		{
			if (Nodes[0]->BC == true && Nodes[1]->BC == true) {
				pow1_pc1(i, j) = pow1_pc1_N(i) * pow1_pc1_N(j) * convection;
				pow1_pc2(i, j) = pow1_pc2_N(i) * pow1_pc2_N(j) * convection;
			}
			if (Nodes[1]->BC == true && Nodes[2]->BC == true) {
				pow2_pc1(i, j) = pow2_pc1_N(i) * pow2_pc1_N(j) * convection;
				pow2_pc2(i, j) = pow2_pc2_N(i) * pow2_pc2_N(j) * convection;
			}
			if (Nodes[2]->BC == true && Nodes[3]->BC == true) {
				pow3_pc1(i, j) = pow3_pc1_N(i) * pow3_pc1_N(j) * convection;
				pow3_pc2(i, j) = pow3_pc2_N(i) * pow3_pc2_N(j) * convection;
			}
			if (Nodes[3]->BC == true && Nodes[0]->BC == true) {
				pow4_pc1(i, j) = pow4_pc1_N(i) * pow4_pc1_N(j) * convection;
				pow4_pc2(i, j) = pow4_pc2_N(i) * pow4_pc2_N(j) * convection;
			}
		}
	}
}

void Element::calculate_pow_sum_and_pow_P()
{
	for (int i = 0; i < matrix_size; i++)
	{
		for (int j = 0; j < matrix_size; j++)
		{
			if (Nodes[0]->BC == true && Nodes[1]->BC == true) {
				double pow1_detJ = (Nodes[1]->x - Nodes[0]->x) / 2;
				pow1_sum(i, j) = (pow1_pc1(i, j) + pow1_pc2(i, j))*pow1_detJ;
				pow1_P(i) += ambient_temp * (pow1_pc1(i, j) + pow1_pc2(i, j)) * pow1_detJ;
			}
			if (Nodes[1]->BC == true && Nodes[2]->BC == true) {
				double pow2_detJ = (Nodes[2]->y - Nodes[1]->y) / 2;
				pow2_sum(i, j) = (pow2_pc1(i, j) + pow2_pc2(i, j))*pow2_detJ;
				pow2_P(i) += ambient_temp * (pow2_pc1(i, j) + pow2_pc2(i, j)) * pow2_detJ;
			}
			if (Nodes[2]->BC == true && Nodes[3]->BC == true) {
				double pow3_detJ = (Nodes[2]->x - Nodes[3]->x) / 2;
				pow3_sum(i, j) = (pow3_pc1(i, j) + pow3_pc2(i, j))*pow3_detJ;
				pow3_P(i) += ambient_temp * (pow3_pc1(i, j) + pow3_pc2(i, j)) * pow3_detJ;
			}
			if (Nodes[3]->BC == true && Nodes[0]->BC == true) {
				double pow4_detJ = (Nodes[3]->y - Nodes[0]->y) / 2;
				pow4_sum(i, j) = (pow4_pc1(i, j) + pow4_pc2(i, j))*pow4_detJ;
				pow4_P(i) += ambient_temp * (pow4_pc1(i, j) + pow4_pc2(i, j)) * pow4_detJ;
			}
		}
	}
}

void Element::calculate_Matrix_P()
{
	for (int i = 0; i < matrix_size; i++)
	{
		Matrix_P(i) = pow1_P(i) + pow2_P(i) +
			pow3_P(i) + pow4_P(i);
	}
}

void Element::calculate_Matrix_H_BC()
{
	for (int i = 0; i < matrix_size; i++)
	{
		for (int j = 0; j < matrix_size; j++)
		{
			Matrix_H_BC(i, j) = pow1_sum(i, j) + pow2_sum(i, j) +
				pow3_sum(i, j) + pow4_sum(i, j);
		}
	}
}

void Element::calculate_Matrix_H_Final()
{
	Matrix_H_Final = Matrix_H + Matrix_H_BC;
}
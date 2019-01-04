#include "Grid.h"
#include <fstream>
#include <iostream>

#define ksi_global 1/sqrt(3)
#define eta_global 1/sqrt(3)

#define shape_function_derivatives_logging 1

#define initial_temp 100
#define simulation_step_time 50
#define simulation_time 500

Grid::Grid(double H, double L, int nH, int nL) {

	double H_step = H / (nH - 1);
	double L_step = L / (nL - 1);

	for (double i = 0; i <= L; i += H_step) {
		for (double j = 0; j <= H; j += L_step) {
			if (i == 0 || i == H || j == 0 || j == L)
				nodes.push_back(new Node{ i, j, true, initial_temp });
			else
				nodes.push_back(new Node{ i, j, false, initial_temp });
		}
	}

	ShapeFunc = calculateShapeFunctions(nodes[0], nodes[nH], nodes[nH + 1], nodes[1]);
	ShFnDerivativesToKsi = calculateShapeFunctionDerivativesToKsi(nodes[0], nodes[nH], nodes[nH + 1], nodes[1]);
	ShFnDerivativesToEta = calculateShapeFunctionDerivativesToEta(nodes[0], nodes[nH], nodes[nH + 1], nodes[1]);

	for (int i = 0; i < nL - 1; i++)
	{
		for (int j = 0; j < nH - 1; j++)
		{

			Element *el = new Element(
				nodes[nH*i + j],
				nodes[nH*(i + 1) + j],
				nodes[nH*(i + 1) + j + 1],
				nodes[nH*i + j + 1],
				ShapeFunc,
				ShFnDerivativesToKsi,
				ShFnDerivativesToEta
			);

			elements.push_back(el);
		}
	}
	calculate_global_matrix_H();
	calculate_global_matrix_C();
	calculate_global_matrix_P();
}

Grid::~Grid()
{
}

void Grid::print_nodes()
{
	std::cout << "-- GRID: NODES --" << std::endl;

	for (auto x : nodes) {
		x->print();
	}

	if (shape_function_derivatives_logging) print_shape_function_derivatives();
}

void Grid::print_elements()
{
	std::cout << "\n-- GRID: ELEMENTS --" << std::endl;

	for (auto x : elements) {
		x->print();
	}
}

void Grid::print_element_by_index(int a)
{
	std::cout << "\n";
	elements[a]->print();
}

Grid* Grid::createGridFromDataFile(std::string fileName)
{
	/* File "data.txt" should contain 4 values - H, L, nH, nL - each in new line */
	double H, L;
	int nH, nL;

	std::fstream plik(fileName, std::ios_base::in);

	std::string line;

	getline(plik, line);
	H = atof(line.c_str());

	getline(plik, line);
	L = atof(line.c_str());

	getline(plik, line);
	nH = atoi(line.c_str());

	getline(plik, line);
	nL = atoi(line.c_str());

	plik.close();

	Grid *g = new Grid(H, L, nH, nL);
	return g;
}

std::vector<ShapeFunction*> Grid::calculateShapeFunctions(Node* a, Node* b, Node* c, Node* d)
{
	// wartosc ksi oraz eta zalezy od polozenia globalnego, ktore zostaje przeksztalcone na polozenie lokalne

	std::vector<ShapeFunction*> sh_func;

	double ksi = -ksi_global;
	double eta = -eta_global;
	double s1 = (0.25*(1 - ksi)*(1 - eta));
	double s2 = (0.25*(1 + ksi)*(1 - eta));
	double s3 = (0.25*(1 + ksi)*(1 + eta));
	double s4 = (0.25*(1 - ksi)*(1 + eta));

	ShapeFunction *sh0 = new ShapeFunction(s1, s2, s3, s4);
	sh_func.push_back(sh0);

	ksi = ksi_global;
	eta = -eta_global;
	s1 = (0.25*(1 - ksi)*(1 - eta));
	s2 = (0.25*(1 + ksi)*(1 - eta));
	s3 = (0.25*(1 + ksi)*(1 + eta));
	s4 = (0.25*(1 - ksi)*(1 + eta));

	ShapeFunction *sh1 = new ShapeFunction(s1, s2, s3, s4);
	sh_func.push_back(sh1);

	ksi = ksi_global;
	eta = eta_global;
	s1 = (0.25*(1 - ksi)*(1 - eta));
	s2 = (0.25*(1 + ksi)*(1 - eta));
	s3 = (0.25*(1 + ksi)*(1 + eta));
	s4 = (0.25*(1 - ksi)*(1 + eta));

	ShapeFunction *sh2 = new ShapeFunction(s1, s2, s3, s4);
	sh_func.push_back(sh2);

	ksi = -ksi_global;
	eta = eta_global;
	s1 = (0.25*(1 - ksi)*(1 - eta));
	s2 = (0.25*(1 + ksi)*(1 - eta));
	s3 = (0.25*(1 + ksi)*(1 + eta));
	s4 = (0.25*(1 - ksi)*(1 + eta));

	ShapeFunction *sh3 = new ShapeFunction(s1, s2, s3, s4);
	sh_func.push_back(sh3);

	return sh_func;
}

Eigen::Matrix<double, 4, 4> Grid::calculateShapeFunctionDerivativesToKsi(Node*, Node*, Node*, Node*)
{
	Eigen::Matrix<double, 4, 4> toKsi;


	std::vector<double> eta_local{ -eta_global, -eta_global, eta_global, eta_global };

	for (int i = 0; i < matrix_size; i++)
	{ //obliczane po kolei kolumnami
		toKsi(0,i) = -0.25*(1 - eta_local[i]);
		toKsi(1,i) = 0.25*(1 - eta_local[i]);
		toKsi(2,i) = 0.25*(1 + eta_local[i]);
		toKsi(3,i) = -0.25*(1 + eta_local[i]);
	}

	return toKsi;
}

Eigen::Matrix<double, 4, 4> Grid::calculateShapeFunctionDerivativesToEta(Node*, Node*, Node*, Node*)
{
	Eigen::Matrix<double, 4, 4> toEta;

	std::vector<double> ksi_local{ -ksi_global, ksi_global, ksi_global, -ksi_global };

	for (int i = 0; i < matrix_size; i++)
	{ //obliczane po kolei kolumnami
		toEta(0,i) = -0.25*(1 - ksi_local[i]);
		toEta(1,i) = -0.25*(1 + ksi_local[i]);
		toEta(2,i) = 0.25*(1 + ksi_local[i]);
		toEta(3,i) = 0.25*(1 - ksi_local[i]);
	}

	return toEta;
}

void Grid::print_shape_function_derivatives()
{
	std::cout << "\nShapeFunctionDerivativestoKsi\n" << ShFnDerivativesToKsi 
		<< "\nShapeFunctionDerivativestoEta\n" << ShFnDerivativesToEta << "\n";
}

void Grid::calculate_global_matrix_H()
{
	for (const auto &x : elements) {
		Eigen::Vector4i ID = { x->Nodes[0]->id, x->Nodes[1]->id, x->Nodes[2]->id, x->Nodes[3]->id };
		//std::cout << "TEST " << ID(0) << " " << ID(1) << " " << ID(2) << " " << ID(3) << " " << "\n";

		for (int i = 0; i < matrix_size; i++) {
			for (int j = 0; j < matrix_size; j++) {
				Global_Matrix_H(ID(i) - 1, ID(j) - 1) += x->Matrix_H_Final(i, j);
			}
		}
	}
}

void Grid::print_global_matrix_H()
{
	const Eigen::IOFormat fmt(3);
	std::cout << "----------------Global_Matrix_H--------------- \n" <<
		std::fixed << Global_Matrix_H.format(fmt) <<
		"\n----------------------------------------------\n\n";
}

void Grid::calculate_global_matrix_C()
{
	for (const auto &x : elements) {
		Eigen::Vector4i ID = { x->Nodes[0]->id, x->Nodes[1]->id, x->Nodes[2]->id, x->Nodes[3]->id };

		for (int i = 0; i < matrix_size; i++) {
			for (int j = 0; j < matrix_size; j++) {
				Global_Matrix_C(ID(i) - 1, ID(j) - 1) += x->Matrix_C(i, j);
			}
		}
	}
}

void Grid::print_global_matrix_C()
{
	const Eigen::IOFormat fmt(2);
	std::cout << "----------------Global_Matrix_C--------------- \n" <<
		std::fixed << Global_Matrix_C.format(fmt) <<
		"\n----------------------------------------------\n\n";
}

void Grid::calculate_global_matrix_P()
{
	for (const auto &x : elements) {
		Eigen::Vector4i ID = { x->Nodes[0]->id, x->Nodes[1]->id, x->Nodes[2]->id, x->Nodes[3]->id };
		for (int i = 0; i < matrix_size; i++) {
			Global_Vector_P(ID(i) - 1) += x->Vector_P(i);
		}
	}
}

void Grid::print_global_matrix_P()
{
	const Eigen::IOFormat fmt(3);
	std::cout << "----------------Global_Vector_P--------------- \n" <<
		std::fixed << Global_Vector_P.format(fmt) <<
		"\n----------------------------------------------\n\n";
}

void Grid::calculate_next_iterations(int iterations)
{
	for (int x = 1; x <= iterations; x++) {

		std::cout << "<<<<< ITERATION " << x - 1 << ", time " << x*simulation_step_time<< " >>>>>\n\n";

		Eigen::Matrix<double, matrix_size*matrix_size, matrix_size*matrix_size> temporaryH = Global_Matrix_H + Global_Matrix_C / (simulation_step_time);

		std::cout << "--------------H_Matrix ([H]+[C]/dT)------------- \n" <<
			std::fixed << temporaryH <<
			"\n----------------------------------------------\n\n";

		Eigen::Matrix<double, 16, 1> temporaryP = Global_Vector_P;
		for (int i = 0; i < matrix_size*matrix_size; i++) {
			for (int j = 0; j < matrix_size*matrix_size; j++) {
				temporaryP(i) += ((Global_Matrix_C(i, j) / (simulation_step_time)) * nodes[j]->temperature);
			}
		}

		std::cout << "-----------P_Vector ([{P}+{[C]/dT}*{T0})--------- \n" <<
			std::fixed << temporaryH <<
			"\n----------------------------------------------\n\n";

		Eigen::Matrix<double, 16, 1> new_temps = temporaryH.inverse() * temporaryP;

		std::cout << "TIME " << x * 50 << ", TEMPERATURES\n";
		for (int i = 0; i < matrix_size*matrix_size; i++) {
			nodes[i]->temperature = new_temps(i);
			if (i == 4 || i == 8 || i == 12) std::cout << "\n";
			std::cout << nodes[i]->temperature << " ";
		}
		std::cout << "\n\n\n";
	}
}
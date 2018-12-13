#include "Grid.h"
#include <fstream>
#include <iostream>

#define ksi_global 1/sqrt(3)
#define eta_global 1/sqrt(3)

#define shape_function_derivatives_logging 0

Grid::Grid(double H, double L, int nH, int nL) {

	for (int i = 0; i < nL; i++) {
		for (int j = 0; j < nH; j++) {
			if (i == 0 || i == nL - 1 || j == 0 || j == nH - 1)
				nodes.push_back(new Node{ L*i, H*j, true });
			else
				nodes.push_back(new Node{ L*i, H*j, false });
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
	std::cout << "\nShapeFunctionDerivativestoKsi\n" << ShFnDerivativesToKsi << "\nShapeFunctionDerivativestoEta\n" << ShFnDerivativesToEta << "\n";
}
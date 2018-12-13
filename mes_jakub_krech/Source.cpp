#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "Grid.h"

#ifndef EIGEN
#define EIGEN
#include <Eigen>
#endif

int main()
{
	Grid *g = Grid::createGridFromDataFile("data.txt");
	g->print_nodes();
	//g->print_elements();
	g->print_element_by_index(0);

	std::cout << "\n";
	system("pause");
	return 0;
}
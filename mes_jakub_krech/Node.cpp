#include "Node.h"
#include <iostream>

Node::Node(double xx, double yy, bool boundary_condition, double initial_temperature) :x{ xx }, y{ yy }, BC{ boundary_condition }, temperature{initial_temperature}
{
	id = global_id++;
}

Node::~Node()
{
	//std::cout << "Node[" << id << "] deleted" << std::endl;
}

int Node::global_id = 1;

void Node::print()
{
	std::cout << "Node[" << id << "] = (" << x << " , " << y << ")  BC: " << BC << "\n";
}
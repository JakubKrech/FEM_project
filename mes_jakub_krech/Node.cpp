#include "Node.h"
#include <iostream>

Node::Node(double xx, double yy, bool boundary_condition) :x{ xx }, y{ yy }, BC{ boundary_condition }
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
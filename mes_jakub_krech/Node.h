#pragma once

class Node
{
public:
	static int global_id;
	int id;
	double temperature;
	double x, y;
	bool BC;

	Node(double, double, bool, double);
	~Node();

	void print();
};
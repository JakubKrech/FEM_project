#pragma once

class Node
{
public:
	static int global_id;
	int id, temperature;
	double x, y;
	bool BC;

	Node(double, double, bool, int);
	~Node();

	void print();
};
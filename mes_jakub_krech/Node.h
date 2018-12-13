#pragma once

class Node
{
public:
	static int global_id;
	int id;
	double x, y, t0;
	bool BC;

	Node(double, double, bool);
	~Node();

	void print();
};
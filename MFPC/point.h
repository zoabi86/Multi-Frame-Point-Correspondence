
#ifndef POINT_H
#define POINT_H

#include <stdio.h>
#include <iostream>


using namespace std;


class Point
{
public:
	Point(){}
	Point(double _x,double _y):x(_x),y(_y){}
	double x;
	double y;

	//friend istream& operator>>(istream&, Point&);
	friend ostream& operator<<(ostream&, const Point&);

};
/*
istream& Point::operator>>(istream& in, Complex& c)
{
	double real, imag;
	in >> real >> imag;
	if (in.good())
	{
		c.myReal = real;
		c.myImag = imag;
	}
	return in;
}*/

ostream& operator<<(ostream& out, const Point& p)
{
	cout << "(" << p.x << "," << p.y << ")";
	return cout;
}

#endif
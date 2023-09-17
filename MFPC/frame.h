#ifndef FRAME_H
#define FRAME_H

#include <lemon/list_graph.h>
#include <lemon/core.h>
#include "point.h"
#include <stdio.h>
#include "def.h"

using namespace lemon;
using namespace std;
static int idzGen=1;

class Frame
{
public:
	int id,numOfvertices;
	double width,height;

	vector<DiNode> vertices;
	vector<Point> points;

	Frame(vector<Point> pts,int pointsNum,double dx,double dy): 
	numOfvertices(0),id(idzGen++),width(dx),height(dy)

	{
		for(int i=0; i<pointsNum; ++i)
			points.push_back(pts[i]);
	}

	void addNode(DiNode u)
	{
		vertices.push_back(u);
		numOfvertices++;
	}

	int size()
	{
		return numOfvertices;
	}

	DiNode& operator[] (int i) {return vertices[i];}

};


#endif
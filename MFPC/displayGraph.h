

#ifndef	DISPLAY_H
#define DISPLAY_H


#include<lemon/list_graph.h>
#include<lemon/graph_to_eps.h>
#include<lemon/math.h>
#include "def.h"

#include <iostream>
#include <sstream>
#include <string.h>


using namespace std;
using namespace lemon;
static int index=0;


//int DisplayGraph(DiGraph& g,int blockSize,Node2PointMap& pmp)
int DisplayGraph(MFPC& mfpc)

{
	Palette palette;
	Palette paletteW(true);
	
	typedef ListDigraph::NodeIt NodeIt;
	typedef dim2::Point<int> Point;


	ListDigraph::NodeMap<Point> coords(g);
	ListDigraph::NodeMap<double> sizes(g);
	ListDigraph::NodeMap<int> colors(g);
	ListDigraph::NodeMap<int> shapes(g);
	ListDigraph::ArcMap<int> acolors(g);
	ListDigraph::ArcMap<int> widths(g);

	Node2PointMap pmp = mfpc->matchingPoint;
	DiGraph g = mfpc->g;
	vector<Frame> Frames = mfpc->Frames;
	
	DiNode curr;
	for(int k=0; k < mfpc->numOfFrames ;++k)
	{
		for(int i=0;i<Frames[k].numOfvertices;++i)
		{
			curr  = Frames[k][i];
			coords[curr]=Point(10*pmp[curr].x ,10*pmp[curr].y); sizes[curr]=2; colors[curr]=1; shapes[curr]=0;
		}

	}

/*
	int i,t=0,counter;
	DiGraph::NodeIt curr(g);

	for(;curr!=INVALID;++t)
	{
		counter=0;
		i=0;
		for(; curr!=INVALID && counter< blockSize; ++curr,i++,counter++)
		{
			//coords[curr]=Point(100-20*t,50+i*10); sizes[curr]=2; colors[curr]=1; shapes[curr]=0;
			coords[curr]=Point(10*pmp[curr].x ,10*pmp[curr].y); sizes[curr]=2; colors[curr]=1; shapes[curr]=0;
		}
	}
		//5/countNodes(g);

		/*i=0;
		for(; curr!=INVALID; ++curr,i++)
		{
			coords[curr]=Point(50,50+i*10); sizes[curr]=1; colors[curr]=1; shapes[curr]=0;
		}*/




	i=0;
	for(DiGraph::ArcIt curr(g); curr!=INVALID; ++curr,++i)
	{
		acolors[curr]=0; widths[curr]=2;
	}

	IdMap<ListDigraph,DiNode> id(g);

	std::stringstream ss,sx,sy;
	ss << index;
	DiGraph::NodeMap<string> names(g);

	for(DiGraph::NodeIt curr(g) ; curr!=INVALID ; ++curr)
	{
		sx.precision(2);
		sy.precision(2);
		sx << pmp[curr].x;
		sy << pmp[curr].y;
		
		//cout << pmp[curr].x << "," << pmp[curr].y << endl; 
		cout << sx.str() << "," << endl << sy.str() << endl; 
		names[curr] = "(" + sx.str() + "," + sy.str() + ")" ;
		sx.str("");
		sy.str("");
		//cout << "-" << names[curr] << "-" << endl;
	}

	cout << "Create 'MFPCgraph" << index << ".eps'" << endl;
	graphToEps(g,"MFPCgraph" +  ss.str() + ".eps").scale(1000).
		title("Sample .eps figure (with arrowheads)").
		copyright("(C) 2003-2009 LEMON Project").
		absoluteNodeSizes().absoluteArcWidths().
		nodeColors(composeMap(palette,colors)).
		coords(coords).
		nodeScale(5).nodeSizes(sizes).
		nodeShapes(shapes).
		arcColors(composeMap(palette,acolors)).
		arcWidthScale(.5).arcWidths(widths).
		nodeTexts(names).nodeTextSize(18).
		drawArrows().arrowWidth(0.4).arrowLength(1).
		enableParallel().parArcDist(0.2).
		run();
	index++;

	return 0;
}

#endif	
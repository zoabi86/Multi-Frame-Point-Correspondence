#ifndef DEF_H
#define DEF_H

#include <lemon/list_graph.h>
#include <lemon/core.h>
#include <lemon/adaptors.h>
#include "point.h"

using namespace lemon;
using namespace std;


typedef ListDigraph			DiGraph;
typedef	ListDigraph::Node	DiNode;
typedef	ListDigraph::Arc	DiEdge;


typedef	ListGraph			uGraph;
typedef	ListGraph::Node		uNode;
typedef	ListGraph::Edge		uEdge;



typedef enum {OLD=0,EXT,CORR} EdgeType;

typedef DiGraph::ArcMap<double>		Edge2doubleMap;
typedef DiGraph::NodeMap<Point>		Node2PointMap;
typedef DiGraph::ArcMap<EdgeType>	Edge2TypeMap;
typedef DiGraph::ArcMap<bool>		Edge2BoolMap;
typedef DiGraph::ArcMap<uEdge>		Arc2uEdgeMap;


typedef DiGraph::NodeMap<uNode>		DiNode2uNodeMap;
typedef DiGraph::NodeMap<int>		DiNode2intMap;
typedef DiGraph::NodeMap<double>	DiNode2doubleMap;
typedef DiGraph::NodeMap<DiNode>	DiNode2DiNodeMap;
typedef uGraph::NodeMap<DiNode>		uNode2DiNodeMap;
typedef DiGraph::NodeMap<bool>		DiNode2boolMap;


#endif
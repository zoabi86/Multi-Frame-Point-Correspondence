

//#include <lemon/list_graph.h>
//#include <lemon/core.h>

//#include "frame.h

#include "parser.h"
#include "def.h"
#include "frame.h"
#include <lemon/math.h>
#include <lemon/random.h>
#include <lemon/matching.h>
#include <lemon/graph_to_eps.h>
#include "myGainFunc.h"
#include <cmath> 
#include "parser_converter.h"

#include "MFPC.h"


using namespace lemon;
using namespace std;

int main(int argc, char *argv[])
{
	vector<Frame> vec;
	//cout << argv[0] << endl;
	//cout << argv[1] << endl;
	//system("pause");

	//return 0;
	
	//int res = Parse("argv[1]",&vec);

	//int res = Parse(argv[1],&vec);
	//	MFPC* mfpc = new MFPC(vec,res,atoi(argv[2]));
	Convertor("newout3.txt");
	int res = Parse("newTracks.txt",&vec);
	MFPC* mfpc = new MFPC(vec,res,5);

	mfpc->MultiframeCorrespondence();	
	delete mfpc;	

	//system("pause");

	return 0;
}



//
//
//int main()
//{
//	uGraph g3;
//	DiGraph g,g2,graph;
//	DiNode u = g.addNode();
//	DiNode v = g.addNode();
//
////	extendGraph(g2);
//	cout<< "number of nodes:" << countNodes(g2) << endl;
//
//	
//	
//	// an example of mapping double value to each node, i.e : weight
//	// NOTE: accessing a mapped value is done in O(1) as said in LEMON tutorial.
//
//	typedef DiGraph::ArcMap<double> DoubleArcMap;
//	DoubleArcMap weight(g);
//
//	// An example of mapping names to nodes.
//	ListDigraph::NodeMap<std::string> name(g);
//	name[u] = "DiNode A";
//
//	// An example of adding 10 edges to the graph 'g' and mapping for them a weight.
//	for(int i=0; i<10; ++i){
//		DiEdge arci = g.addArc(u,v);
//		cout << " A new arc: " << g.id(arci) << " has been added to the graph" << endl;
//	}
//
//	cout << "-------------------------------" << endl;
//	int j=1;
//	for(DiGraph::ArcIt curr(g); curr!=INVALID; ++curr,++j)
//		weight[curr] = rnd(100.0);//(j*PI/3);
//
//
//	for(DiGraph::ArcIt curr(g); curr!=INVALID; ++curr)
//		cout << "the weight of arc " << g.id(curr) <<" is: " << weight[curr] << endl;
//
//
//	cout << "******************************" << endl;
//	
//
//	uNode u3 = g3.addNode();
//	uNode v3 = g3.addNode();
//
//	for(int i=0; i<10; ++i){
//		uEdge edgei = g3.addEdge(u3,v3);
//	}
//
//	uGraph::EdgeMap<double> wm(g3);
//	
//	for(uGraph::EdgeIt curr(g3); curr!=INVALID; ++curr)
//		wm[curr] = rnd(100.0);
//
//	//typedef  uGraph::NodeMap<DiEdge> MatchingMap;
//
//	MaxWeightedMatching<uGraph,uGraph::EdgeMap<double> > m = MaxWeightedMatching<uGraph,uGraph::EdgeMap<double> >(g3,wm);
//	m.run();
////	MaxWeightedMatching<uGraph,uGraph::EdgeMap<double> >::MatchingMap x = m.matchingMap();
//
//	for(uGraph::NodeIt curr(g3); curr!=INVALID; ++curr)
//		cout << 	"id of edge is: ---- " << g3.id(m.matchingMap()[curr])  << endl;
//	
//
//	for(DiGraph::NodeIt curr(g); curr!=INVALID; ++curr)
//		cout << g.id(v) << endl;
//	
//	//DisplayGraph(g);
//
//	cout << "hello world, graphs world" << endl;
//
//
//	/*################################################################333
//	########################################################################
//	#############################################################################
//	##############################################################################*/
//	
//	vector<Frame> vec;
//
//	int res = Parse("tracks.txt",&vec);
//	cout << "res= " << res << endl;
//		//system("pause");
//	MFPC mfpc(vec,res);
//	mfpc.MultiframeCorrespondence();	
//
//	system("pause");
//
//	return 0;
//	/*vector<Point> vp1,vp2,vp3,vp4,vp5;
//	for(int i=0;i<10;++i){
//		
//		Point* p1 = new Point(0,i*10 + rnd(10.0) );
//		vp1.push_back(*p1);
//		
//		Point* p2 = new Point(10,i*10 + rnd(6.0) );
//		vp2.push_back(*p2);
//		
//		Point* p3 = new Point(20,i*10 + rnd(5.0)*(-1) );
//		vp3.push_back(*p3);
//		
//		Point* p4 = new Point(30,i*10 + rnd(14.0) );
//		vp4.push_back(*p4);
//		
//		Point* p5 = new Point(40,i*10 + rnd(2.0) );
//		vp5.push_back(*p5);
//	}
//	
//	
//
//	Frame* f1  = new Frame(vp1,10);
//	Frame* f2  = new Frame(vp2,10);
//	Frame* f3  = new Frame(vp3,10);
//	Frame* f4  = new Frame(vp4,10);
//	Frame* f5  = new Frame(vp5,10);
//	vec.push_back(*f1);
//	vec.push_back(*f2);
//	vec.push_back(*f3);
//	vec.push_back(*f4);
//	vec.push_back(*f5);
//*/
//
//}
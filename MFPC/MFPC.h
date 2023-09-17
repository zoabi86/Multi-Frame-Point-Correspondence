#ifndef MFPC_H
#define MFPC_H


#include <lemon/list_graph.h>
#include <lemon/core.h>
#include <lemon/graph_to_eps.h>
#include <lemon/math.h>
#include <lemon/random.h>

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string.h>

#include "frame.h"
#include "point.h"
#include "def.h"
#include "myGainFunc.h"

//#include "displayGraph.h"


static int index=0;

using namespace lemon;
using namespace std;


class MFPC
{
	DiGraph g;
	Edge2doubleMap weight;
	Edge2BoolMap  isValid;
	Edge2BoolMap  isFH; // this means: if this edge is False Hypothesis.
	Node2PointMap matchingPoint;
	Edge2TypeMap edgeType;
	vector<Frame> Frames;
	int windowSize;
	int totalPoints,numOfFrames,currCreatedFrames;
	DiNode2intMap numPointsInTrack;
	DiNode2intMap indexOfIncludingFrame;
	DiNode2doubleMap sumOfSlopesInTrack;

	//*********************************//
	DiNode2intMap bt_numPointsInTrack;
	DiNode2doubleMap bt_sumOfSlopesInTrack;
	//*******************************//

public:
	
	MFPC(vector<Frame> vecOfFrames,int framesNum,int wSize):weight(g),
		matchingPoint(g),edgeType(g),
		isValid(g),isFH(g),numOfFrames(framesNum),windowSize(wSize),numPointsInTrack(g),
		sumOfSlopesInTrack(g),bt_numPointsInTrack(g),bt_sumOfSlopesInTrack(g),indexOfIncludingFrame(g)
	{
		Frames = vecOfFrames;
		/*for(int i=0;i<numOfFrames;++i)
		{
			Frame* f = new Frame(); // create all the empty frames, they will be filled on demand.
			*f  = vec[i];
			Frames.push_back(*f);		
		}*/

	}

	void MultiframeCorrespondence()
	{
		// solve the 2-frames problem
		create2FramesGraph();
		calcWeights(true,false);
		maxPathCover();
		gainFuncAux(1,false);

		int r=2;
		for(int i=2; i<numOfFrames; ++i)
		{
			extendGraph(Frames[i].points,i-r);
			//DisplayGraph(g);
			calcWeights(false,false);
			printEdgesWeights();
			DisplayGraph(g);
			maxPathCover();
			DisplayGraph(g);
			falseHypothesisReplacement(1,0,currCreatedFrames-1,false);
			gainFuncAux(i,false);
			if(i == windowSize-1)
			{	
				cout << "backtracking is started.." << endl;
				backTrack();
				reverseAllArcs();

			}
			++r;
			if(r >= windowSize)
				r = windowSize-1;

		}
		//DisplayGraph(g);
	}
	void backTrack()
	{
		reverseAllArcs();
		int farestFrame = findMostFarFrame(windowSize-1);
		for(int i=0; i<Frames[windowSize-1].size(); ++i)
		{
			bt_numPointsInTrack[Frames[windowSize-1][i]] = 0;
			bt_sumOfSlopesInTrack[Frames[windowSize-1][i]] = 0;
		}
		
		for(int i=windowSize-2;i>=farestFrame; --i)
		{
			gainFuncAux(i,true);
		}
		//DisplayGraph(g);
		bt_clearEdges(farestFrame);
		//DisplayGraph(g);
		//DisplayGraph(g);
		for(int i=farestFrame; i>0; --i)
		{
			gainFuncAux(i,true);
			bt_Extend(i-1,farestFrame);
			//DisplayGraph(g);
			if(i==1)
			{
				for(int j=0; j<Frames[2].size() ; ++j)
				{
					cout << "|||||||||||||||||" << endl;
					cout << "bef- SUM directly after update = " << bt_sumOfSlopesInTrack[Frames[2][j]] << endl;
					cout << "bef- NUM directly after update = " << bt_numPointsInTrack[Frames[2][j]] << endl;
					cout << "|||||||||||||||||" << endl;
				}
			}
			
			calcWeights(false,true);
			if(i==1)
			{
				for(int j=0; j<Frames[2].size() ; ++j)
				{
					cout << "|||||||||||||||||" << endl;
					cout << "SUM directly after update = " << bt_sumOfSlopesInTrack[Frames[2][j]] << endl;
					cout << "NUM directly after update = " << bt_numPointsInTrack[Frames[2][j]] << endl;
					cout << "|||||||||||||||||" << endl;
				}
			}
			maxPathCover();

			MonitorTracks(windowSize-1);
			//DisplayGraph(g);
						
			falseHypothesisReplacement(-1,farestFrame,i,true);
			
			

			//DisplayGraph(g);
		}
		MonitorTracks(farestFrame);
		//DisplayGraph(g);

	}

	void create2FramesGraph()
	{
		DiNode u;
		for(int i=0;i<Frames[0].points.size();++i)
		{
			u = g.addNode();
			matchingPoint[u] = Frames[0].points[i];
			numPointsInTrack[u] = 1;
			sumOfSlopesInTrack[u] = 0;
			Frames[0].addNode(u);
			indexOfIncludingFrame[u] = 0;
		}

		for(int i=0;i<Frames[0].size();++i)
			cout << g.id(Frames[0][i]) << endl;

		for(int i=0;i<Frames[1].points.size();++i)
		{
			u = g.addNode();
			matchingPoint[u] = Frames[1].points[i];
			numPointsInTrack[u] = 2;
			sumOfSlopesInTrack[u] = 0;
			Frames[1].addNode(u);
			indexOfIncludingFrame[u] = 1;
		}
	
		// create all the edges from Frame1 to Frame2
		for(int i=0;i<Frames[0].size();++i)
		{
			for(int j=0;j<Frames[1].size();++j)
			{
				DiEdge e = g.addArc(Frames[0][i],Frames[1][j]);
				edgeType[e]=OLD;
				isValid[e]=true;
			}
		}

		currCreatedFrames=2;
	}

	void extendGraph(vector<Point> nVec,int idx){
		
		DiNode u;
		for(int i=0;i<nVec.size();++i)
		{
			u = g.addNode();
			matchingPoint[u] = nVec[i];
			numPointsInTrack[u] = 0;
			sumOfSlopesInTrack[u] = 0;
			Frames[currCreatedFrames].addNode(u);
			indexOfIncludingFrame[u] = currCreatedFrames;
		}
		// create all the edges from Frame-k to the new created Frame.

		for(int k=idx;k<currCreatedFrames;++k)
		{
			for(int i=0;i<Frames[k].size();++i)
			{
				for(int j=0;j<Frames[currCreatedFrames].size();++j)
				{
					DiEdge e;
					if ( isTerminal(Frames[k][i]) ) // It's important to check isTerminal before adding the new edge.
					{
						e = g.addArc(Frames[k][i],Frames[currCreatedFrames][j]);
						edgeType[e]=EXT;
					}
					else
					{
						e = g.addArc(Frames[k][i],Frames[currCreatedFrames][j]);
						edgeType[e]=CORR;
					}
					isValid[e]=true;
				}
			}
		}
	
		currCreatedFrames++;

	}

	void bt_Extend(int idx,int farest)
	{

		for(int k=farest;k>idx;--k)
		{
			for(int i=0;i<Frames[k].size();++i)
			{
				for(int j=0;j<Frames[idx].size();++j)
				{
					DiEdge e;
					if ( isTerminal(Frames[k][i]) ) // It's important to check isTerminal before adding the new edge.
					{
						e = g.addArc(Frames[k][i],Frames[idx][j]);
						edgeType[e]=EXT;
					}
					else
					{
						e = g.addArc(Frames[k][i],Frames[idx][j]);
						edgeType[e]=CORR;
					}
					isValid[e]=true;
				}
			}
		}
	}
	bool isTerminal(DiNode u)
	{
		DiGraph::OutArcIt a(g, u);
		return (a==INVALID);
	}

	bool isSource(DiNode u)
	{
		DiGraph::InArcIt a(g, u);
		return (a==INVALID);
	}
	void markFalseHypothesisEdges()
	{
		for(DiGraph::ArcIt curr(g); curr!=INVALID; ++curr)
		{
			if(!isValid[curr] && edgeType[curr]==OLD)
			{
				//Attention: we assume here that for each node there is only one outer arc, because this method must be
				// called after calling : maximum path cover.
				if( !isTerminal(g.target(curr)) )
					cout << "yeah" << endl;
				for(DiGraph::OutArcIt outE(g,g.target(curr)) ; outE!=INVALID ; outE = DiGraph::OutArcIt(g,g.target(outE)) )
					isFH[outE] = true;
			}
		}

	}

	bool checkCond(bool isBackTrack, int a , int b)
	{
		if(isBackTrack)
			return a > b;
		return a < b;
	}
	void falseHypothesisReplacement(int step,int a,int b,bool isBackTrack)
		// the non-recursive version.
	{


		for(int k=a;checkCond(isBackTrack,k,b);k+=step)
		{

			for(int i=0;i<Frames[k].size();++i)
			{
				for(DiGraph::OutArcIt outE(g,Frames[k][i]); outE!=INVALID ; ++outE)
				{	
					if(isFH[outE]==true)
						g.erase(outE);
				}
				
			}
		

			DiGraph tmpG;
			DiNode2boolMap isInF1(tmpG);
			DiNode2DiNodeMap mp(tmpG),revmp(g);
			Node2PointMap n2pmp(tmpG);
			// create a temp graph that has only the vertices from Fk  and Fk+step
			for(int i=0;i<Frames[k].size();++i)
			{
					
				if(isTerminal(Frames[k][i]))
				{
					DiNode newNode = tmpG.addNode();
					isInF1[newNode]=true;
					mp[newNode] = Frames[k][i];
					//revmp[Frames[k][i]] = newNode;
					n2pmp[newNode] = matchingPoint[Frames[k][i]];
				}
				
			}
			for(int i=0;i<Frames[k+step].size();++i)
			{
				if(isSource(Frames[k+step][i]))
				{
					DiNode newNode = tmpG.addNode();
					isInF1[newNode] = false;
					mp[newNode] = Frames[k+step][i];
					//revmp[Frames[k+step][i]] = newNode;
					n2pmp[newNode] = matchingPoint[Frames[k+step][i]];
				}
			}
			// create all the edges from Frame k to Frame k+step

			if(countNodes(tmpG)==0) continue;

			for(DiGraph::NodeIt curr(tmpG); curr!=INVALID; ++curr)
			{
				if(!isInF1[curr])
					continue;				
				for(DiGraph::NodeIt curr2(tmpG); curr2!=INVALID; ++curr2)
				{
					if(isInF1[curr2])
						continue;
					DiEdge e = tmpG.addArc(curr,curr2);
					edgeType[e]=OLD;
					isValid[e]=true;
				}
			}

			if(countArcs(tmpG)<1)
				continue;
			giveWeights(tmpG,n2pmp,false,isBackTrack);
			//DisplayGraph(tmpG);
			maxPathCoverAux(tmpG);

			for(DiGraph::ArcIt curr(tmpG); curr!=INVALID; ++curr)
			{
				DiEdge e = g.addArc(mp[tmpG.source(curr)],mp[tmpG.target(curr)]);
				edgeType[e] = OLD;
				isValid[e] = true;
				weight[e] = weight[curr];
			}
			
		}
	}



	void maxPathCover()
	{
		maxPathCoverAux(this->g);
	}

	void maxPathCoverAux(DiGraph& gr)
	{
		uGraph ug;
		DiNode2uNodeMap mpIn(gr);
		DiNode2uNodeMap mpOut(gr);
		uNode2DiNodeMap revMpIn(ug);
		Arc2uEdgeMap arcMp(gr);

		uNode newNodeOut,newNodeIn;

		for(DiGraph::NodeIt curr(gr); curr!=INVALID; ++curr)
		{
			newNodeIn = ug.addNode();
			newNodeOut = ug.addNode();
			mpIn[curr] = newNodeIn;
			revMpIn[newNodeIn] = curr;
			mpOut[curr] = newNodeOut;
		}
				
		uGraph::EdgeMap<double> wm(ug);

		for(DiGraph::ArcIt currArc(gr); currArc!=INVALID; ++currArc)
		{
			DiNode u = gr.source(currArc);
			DiNode v = gr.target(currArc);
			uEdge theNewEdge  = ug.addEdge(mpOut[u],mpIn[v]);
			arcMp[currArc] = theNewEdge;
			wm[theNewEdge] = weight[currArc];

		}
		//DisplayGraph(gr);
		MaxWeightedMatching<uGraph,uGraph::EdgeMap<double> > x = MaxWeightedMatching<uGraph,uGraph::EdgeMap<double> >(ug,wm);
		x.run();

		for(DiGraph::ArcIt curr(gr); curr!=INVALID ; ++curr)
		{
			cout << " current arc is: ";
			printEdgeType(curr);
			if (isValid[curr]) 
				cout << " isValid = true" << endl;
			else 
				cout << " isValid = false" << endl;

		}


		for(DiGraph::ArcIt curr(gr); curr!=INVALID ; ++curr)
		{
				if(x.matching(arcMp[curr])==false)
				{
					isValid[curr]=false;
				}
				else
				{
					edgeType[curr]=OLD;
				}	
		}

		for(DiGraph::ArcIt curr(gr); curr!=INVALID ; ++curr)
		{
			cout << "********* current arc is: ";
			printEdgeType(curr);
			if (isValid[curr]) 
				cout << " isValid = true" << endl;
			else 
				cout << " isValid = false" << endl;

		}
		//DisplayGraph(g);
		markFalseHypothesisEdges();

		for(DiGraph::ArcIt curr(gr); curr!=INVALID ; ++curr)
		{
			if(!isValid[curr])
				gr.erase(curr);

		}
		
	}


	void DBGgiveWeights()
	{
		int i=0;
		for(DiGraph::ArcIt curr(g); curr!=INVALID ; ++curr,++i)
		{
			if(i%2==0)
				weight[curr] = 1000;
			else
				weight[curr] = rnd(100.0);
		}


	}

	void calcWeights(bool isInit,bool isBackTrack)
	{
		giveWeights(this->g,matchingPoint,isInit,isBackTrack);
	}

	void giveWeights(DiGraph& gr,Node2PointMap& mp,bool isInit,bool isBackTrack)
	{
		int i=0;
		for(DiGraph::ArcIt curr(gr); curr!=INVALID ; ++curr,++i)
		{
			DiNode u = gr.source(curr);
			DiNode v = gr.target(curr);

			if (indexOfIncludingFrame[g.source(curr)] == currCreatedFrames-1 && isBackTrack)
				continue;
			Matrix observed(2,1),predicted(2,1);
				
			observed(0,0) = mp[v].x - mp[u].x;
			observed(1,0) = mp[v].y - mp[u].y;
			if(isInit)
				initial_predict(u,predicted);
			predict(u,predicted,isBackTrack);
			

			double advantage=0;
			double sx,sy;
			//calcDims(sx,sy);
			int deltaFrames = indexOfIncludingFrame[v]-indexOfIncludingFrame[u];
			if(deltaFrames > 1)
				advantage=0.9;

			weight[curr] = gainFunc(predicted,observed) + advantage;

			//if ((matchingPoint[u].x == 11.0 && matchingPoint[u].y == 17.19) &&
			//	((matchingPoint[v].x == 8.26 && matchingPoint[v].y == 18.19) ||
			//	(matchingPoint[v].x == 11.6 && matchingPoint[v].y == 19.17))
			//	)
			//{
				cout << "----------------------------------BEGIN-----------------------------------" << endl;
				cout << "curr edge is: " << "(" << matchingPoint[u].x << "," << matchingPoint[u].y << ")" <<
					" -> (" << matchingPoint[v].x << "," << matchingPoint[v].y << ")" << endl;
				cout << "observed is: " << "(" << observed(0,0) << "," << observed(1,0) << ")" << endl;
				cout << "predicted is: " << "(" << predicted(0,0) << "," << predicted(1,0) << ")" << endl;
				cout << "gain returned: " << weight[curr] << endl;
				cout << "-----------------------------------END------------------------------------" << endl;

			
			//weight[curr] = (rand()%10) / 10 + 0.1;
			//if(matchingPoint[u].x == 8.14 && matchingPoint[v].x == 9.25){
			//	weight[curr] = 0.9;
			//}
			//else weight[curr]=0.00001;
			//cout << "The weight of: " << matchingPoint[u] << "--->" << matchingPoint[v] << " is:" << weight[curr] << endl;
		
		}

	}
	void initial_predict(DiNode n,Matrix& m)
	{
		m(0,0) = 0; // this should be changed
		m(1,0) = 1; // this means the vector: (1,0)
	}
	void predict(DiNode n,Matrix& m,bool isBackTrack)
	{

		double slopesAvg;
		if(isBackTrack){
			if(!bt_numPointsInTrack[n])
				cout<<endl;
			slopesAvg = bt_sumOfSlopesInTrack[n] / bt_numPointsInTrack[n];
			cout << "^^^^^^^^^^^^^^^^^^^^ PREDICT BEGIN ^^^^^^^^^^^^^^^^^^^" << endl;
			cout << "U:= " << "(" << matchingPoint[n].x << "," << matchingPoint[n].y << ")" << "in frame: " << indexOfIncludingFrame[n] << endl;
			cout << "Slopes avg :" << slopesAvg<<"=" <<bt_sumOfSlopesInTrack[n] <<"/" <<bt_numPointsInTrack[n] << endl;
			cout << "^^^^^^^^^^^^^^^^^^^^ PREDICT END ^^^^^^^^^^^^^^^^^^^^^" << endl;
		}
		else
			slopesAvg = sumOfSlopesInTrack[n] / numPointsInTrack[n];
		
		double vectorLen = 2;

		m(1,0) = -sqrt( vectorLen /(1+pow(slopesAvg,2)) ); // check here is kevon hetkadmot, and put sign=+-1 correspondingly
		m(0,0) = -slopesAvg * m(1,0);
		
		//m(1,0) = slopesAvg;
		//m(0,0) = 1;


		if(isBackTrack)
		{
			m(0,0) = -1 * m(0,0); 
			m(1,0) = -1 * m(1,0);
		}
		//else
		//{
		//	m(1,0) = slopesAvg; 
		//	m(0,0) = 1.0;
		//}



		//m(0,0) = 8.98; 
		//m(1,0) = 6.73;

		//cout << "Slopes avg:" << slopesAvg << endl;
		// should be changed to vector: Dy/Dx : the average slope of a point in the track that includes this point.
	}

	void printEdgesWeights()
	{
		for(DiGraph::ArcIt curr(g); curr!=INVALID ; ++curr)
		{
			DiNode u = g.source(curr);
			DiNode v = g.target(curr);
			//cout << "The weight of: " << matchingPoint[u] << "--->" << matchingPoint[v] << " is:" << weight[curr] << endl;
		}

	}
	

	int DisplayGraph(DiGraph& g)

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
		
		double nodeSize = 1.2; // should be calculated properly -  as a function of :
		//num of nodes in g, the density of the nodes in the graph ( by coordinations ).					
		double arcWidth = 1.3;

		int numOfLowRecs=0,numOfHighRecs=0;

		for(int k=0 ; k < numOfFrames ; ++k)
		{
			if(Frames[k].width > Frames[k].height)
				numOfLowRecs++;
			else
				numOfHighRecs++;
		}

		int printUpwards = (numOfLowRecs > numOfHighRecs) ? 1 : 0;

		DiNode curr;
		int offset=0;
		for(int k=0 ; k < numOfFrames ; ++k)
		{
			int acc=0;
			for(int i=0;i<Frames[k].numOfvertices;++i)
			{
				curr  = Frames[k][i];
				
				sizes[curr]=nodeSize;
				colors[curr]=1;
				shapes[curr]=0;

				if(printUpwards)
					coords[curr]=Point(11*(matchingPoint[curr].x) ,offset + 8*matchingPoint[curr].y);
				else
					coords[curr]=Point(offset + 10*matchingPoint[curr].x ,10*matchingPoint[curr].y);
			}
			
			if(printUpwards)
				offset += 10*Frames[k].height;
			else
				offset += 10*Frames[k].width;

		}



		
		int i=0;

		for(DiGraph::ArcIt curr(g); curr!=INVALID; ++curr,++i)
		{
			acolors[curr]=38; 
			widths[curr]=arcWidth;
		}

		int currentColor=1;
		for(int k=0; k<Frames[0].size();++k)
		{
			for(DiGraph::OutArcIt outE(g,Frames[0][k]) ; outE!=INVALID ; outE = DiGraph::OutArcIt(g,g.target(outE)) )
			{
				acolors[outE]=currentColor;
				widths[outE]=arcWidth;
			}
			currentColor+=1;
		}



		IdMap<ListDigraph,DiNode> id(g);

		std::stringstream ss,sx,sy;
		ss << index;
		DiGraph::NodeMap<string> names(g);

		for(DiGraph::NodeIt curr(g) ; curr!=INVALID ; ++curr)
		{
			sx.precision(4);
			sy.precision(4);
			sx << matchingPoint[curr].x;
			sy << matchingPoint[curr].y;
			
			//cout << pmp[curr].x << "," << pmp[curr].y << endl; 
			cout << sx.str() << "," << sy.str() << endl; 
			names[curr] = "(" + sx.str() + "," + sy.str() + ")" ;
			sx.str("");
			sy.str("");
			//cout << "-" << names[curr] << "-" << endl;
		}

		cout << "Create 'MFPCgraph" << index << ".eps'" << endl;
		graphToEps(g,"Plots/MFPCgraph" +  ss.str() + ".eps").scale(1000).
			title("Sample .eps figure (with arrowheads)").
			copyright("(C) 2003-2009 LEMON Project").
			absoluteNodeSizes().absoluteArcWidths().
			nodeColors(composeMap(palette,colors)).
			coords(coords).
			nodeScale(nodeSize).nodeSizes(sizes).
			nodeShapes(shapes).
			arcColors(composeMap(palette,acolors)).
			arcWidthScale(.2).arcWidths(widths).
			nodeTexts(names).nodeTextSize(0.9*nodeSize).
			drawArrows().arrowWidth(0.4).arrowLength(1).
			enableParallel().parArcDist(0.2).
			run();
		index++;

		return 0;
		
	}

	void gainFuncAux(int i,bool isBackTrack)
	{	
	
		for (int j=0; j< Frames[i].size(); ++j)
		{
			int count=0;
			for(DiGraph::InArcIt InE(g,Frames[i][j]) ; InE!=INVALID ; InE++,count++ )
			{
				if(count==1)
				{
					cout << "exit(0)" << endl;
					system ("pause");
					exit(0);
				}
				DiNode u = g.source(InE);
				DiNode v = Frames[i][j];
				

				double slope = (matchingPoint[v].x - matchingPoint[u].x ) / (matchingPoint[v].y - matchingPoint[u].y ) ;

				

				cout << "calcSlop = " << bt_sumOfSlopesInTrack[u]<< " + " << slope << endl;
				
				if(isBackTrack)
				{
					cout << "++++++++++++++++++++++++START+++++++++++++++++++++++" << endl;
					cout << "+++++++++++++++++++(" << i << ")+++++++++++++++++++++" << endl;

					cout << "numOfPOintsINtrack: " << bt_numPointsInTrack[u] << endl;
					cout << "someOfSlopes: " << bt_sumOfSlopesInTrack[u] << endl;

					bt_numPointsInTrack[v] = bt_numPointsInTrack[u] + 1;
					bt_sumOfSlopesInTrack[v] = bt_sumOfSlopesInTrack[u] + slope;

					cout << "NEW SUM = " << bt_sumOfSlopesInTrack[v] << endl;
					cout << "NEW NUM = " << bt_numPointsInTrack[v] << endl;
					cout << "+++++++++++++++++++++++++END+++++++++++++++++++++++" << endl;
				}
				else
				{
					numPointsInTrack[v] = numPointsInTrack[u] + 1;
					sumOfSlopesInTrack[v] = sumOfSlopesInTrack[u] + slope;
				}

				


			}
		}
	} 

	//++++++++++++++++++++++++++++++++++++++++++++++++++

	void reverseAllArcs()
	{
		//int cnt=countArcs(g);
		//Edge2BoolMap isReversedEdge(g);
		//for(DiGraph::ArcIt curr(g); curr!=INVALID ; ++curr)
		//	isReversedEdge[curr]=false;
		
		typedef struct arc_atrr_t
		{
			bool _isFH;
			bool _isValid;
			double _weight;
			EdgeType _type;
			DiNode target;
			DiNode source;
		}ArcAtrr;
		
		vector<ArcAtrr> tmpvec;

		
		for(DiGraph::ArcIt curr(g); curr!=INVALID ; ++curr)
		{	

			ArcAtrr* at = new ArcAtrr();

			at->_isFH=isFH[curr];
			at->_isValid=isValid[curr];
			at->_type=edgeType[curr];
			at->_weight=weight[curr];
			at->source=g.source(curr);
			at->target=g.target(curr);

			tmpvec.push_back(*at);
			
			g.erase(curr);
			//if(isReversedEdge[curr]==false)
			//{
			//	g.erase(curr);
			//	cnt--;
			//}
			//DiEdge newE = g.addArc(target,source);
			//isReversedEdge[newE]=true;
			//
			//isFH[newE]=_isFH;
			//isValid[newE]=_isValid;
			//weight[newE]=_weight;
			//edgeType[newE]=_type;
			
		}		
		for(int i=0; i < tmpvec.size(); i++)
		{
			DiEdge newE = g.addArc(tmpvec[i].target,tmpvec[i].source);
			isFH[newE]=tmpvec[i]._isFH;
			isValid[newE]=tmpvec[i]._isValid;
			weight[newE]=tmpvec[i]._weight;
			edgeType[newE]=tmpvec[i]._type;
		}


		

	}

	int findMostFarFrame(int k)
	{
		int farest=k;
		for (int j=0; j< Frames[k].size(); ++j)
		{
			int count=0;
			for(DiGraph::OutArcIt outE(g,Frames[k][j]) ; outE!=INVALID ; outE++,count++ )
			{
				if(indexOfIncludingFrame[g.target(outE)] < farest) // indexOfIncludingFrame is a map from
					farest = indexOfIncludingFrame[g.target(outE)];
			}
		}
		return farest;
		
	}

	void bt_clearEdges(int farest)
	{
		for(int k=farest;k>0;--k)
		{
			for (int j=0; j< Frames[k].size(); ++j)
			{
				for(DiGraph::OutArcIt outE(g,Frames[k][j]) ; outE!=INVALID ; outE++)
				{
					g.erase(outE);
				}
			}
		}
	}


	void MonitorTracks(int i)
	{
		cout << "######################################################################" << endl;
		for(int k=0; k<Frames[i].size();++k)
		{
			cout << "Track number: " << k+1 ;
			for(DiGraph::OutArcIt outE(g,Frames[i][k]) ; outE!=INVALID ; outE = DiGraph::OutArcIt(g,g.target(outE)) )
			{
				cout << "  %% num = " << bt_numPointsInTrack[g.source(outE)] << " , sum = " << bt_sumOfSlopesInTrack[g.source(outE)] << " | ";
			}
			cout << endl;
		}
		cout << "######################################################################" << endl;
	}
	private:
		void printEdgeType(DiEdge e)
		{
			switch (edgeType[e])
			{
			case OLD:
				cout << "OLD";
				break;

			case EXT:
				cout << "EXT";
				break;

			case CORR:
				cout << "CORR";
				break;
			}
		}

		
		void calcDims(double& sx,double& sy)
		{
			double maxWidth=0,maxHeight=0;

			for(int i=0; i<numOfFrames; ++i)
			{
				if(Frames[i].width > maxWidth)
					maxWidth = Frames[i].width;

				if(Frames[i].height > maxHeight)
					maxHeight = Frames[i].height;
			}
			sx = maxWidth;
			sy = maxHeight;
		}
};


#endif

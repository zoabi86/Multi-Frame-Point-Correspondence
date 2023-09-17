
#ifndef CON_H
#define CON_H

#define T_SPC "   "
#define D_TAB "\t\t"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream> 
#include <fstream>
#include "point.h"

using namespace std;


int Convertor(char* ijFileName)
{

	//Open file for reading
	ifstream inFile(ijFileName);
	ofstream tracksFile("newTracks.txt", ios::out);

	//Temporary buffer for line read from file
	string line;

	//Skipping first line
//	getline(inFile,line);

	int counter=0,frameIndex=1;
	vector<Point> vec;
	//Reading rest of the lines
	while(!inFile.eof() && inFile.good())
	{
		string  trash;
		double xPoint,yPoint;
		inFile >> trash;

		if(trash[0]=='*')
		{
				tracksFile << frameIndex << D_TAB << counter << D_TAB;
				for(int i=0;i<counter;++i)
					tracksFile << vec[i].x << T_SPC << vec[i].y << D_TAB;
				tracksFile << endl;
				
				vec.clear();
				counter=0;
				frameIndex++;
				getline(inFile,line);
				//getline(inFile,line);
				continue;
		}
		else if( (trash[0]-'0')<0 || (trash[0]-'0')>9)
		{
			getline(inFile,line);
			continue;
		}

		inFile >> trash >> xPoint >> yPoint;
		vec.push_back(Point(xPoint,yPoint));
		counter++;
		
	}

	return 0;

}

#endif



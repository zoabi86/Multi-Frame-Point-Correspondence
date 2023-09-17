
#ifndef PARS_H
#define PARS_H


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream> 
#include "frame.h"

//#include "aux_functions.h"
//#include "input_functions.h"

//#define delim '\t '

using namespace std;


int Parse(char* filename,vector<Frame>* allFrames)
{
	FILE *ifd;
	//char current_line[MAX_LEN],temp_line[MAX_LEN],the_command[MAX_LEN],*strtok_result;

	int frameNum,numOfPoints;
	double currPointX,currPointY;

	ifd=fopen(filename,"r");
	
	if(ifd==NULL)
		return -1;


	//vector<Frame> allFrames;

	int numOfFrames=0;

	while(feof(ifd)==0)
	{
		if(fscanf(ifd, "%d", &frameNum)==EOF) break;
		fscanf(ifd, "%d", &numOfPoints);
		cout << "frame: " << frameNum << " , numOfPoints: " << numOfPoints << endl;
		

		double maxBoundX,minBoundX,maxBoundY,minBoundY;
		vector<Point> vp;

		for(int i=0 ; i < numOfPoints ; ++i)
		{
			fscanf(ifd, "%lf", &currPointX);
			fscanf(ifd, "%lf", &currPointY);		
			if(i==0)
			{
				maxBoundX = currPointX;
				minBoundX = currPointX;

				maxBoundY = currPointY;
				minBoundY = currPointY;
			}
			else
			{
				maxBoundX = (currPointX > maxBoundX) ? currPointX : maxBoundX;
				minBoundX = (currPointX < minBoundX) ? currPointX : minBoundX;
				
				maxBoundY = (currPointY > maxBoundY) ? currPointY : maxBoundY;
				minBoundY = (currPointY < minBoundY) ? currPointY : minBoundY;
			}
			//cout << "(" << currPointX << "," << currPointY << ")	";

			Point* p = new Point(currPointX,currPointY);
			vp.push_back(*p);

		}
		//RectangleSizeX = maxBoundX - minBoundX;
		//RectangleSizeY = maxBoundY - minBoundY;

		Frame* f  = new Frame(vp,numOfPoints,maxBoundX - minBoundX,maxBoundY - minBoundY);
		(*allFrames).push_back(*f);
		numOfFrames++;

		//cout << endl;

	}

	if(fclose(ifd)==EOF)
		return -1;

	//system("pause");
	return numOfFrames;
}

#endif

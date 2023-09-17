
#ifndef GAIN_H
#define GAIN_H

#include <time.h>  
#include <math.h>
#include "matrix.h"  


#define SX	10// BEWARE this is wrong -  should be the actual frame size i.e width and height.
#define SY	10
#define A	0.9
#define N	2.0

using namespace std;
using namespace math;

typedef matrix<double> Matrix;


double gainFunc( Matrix& mPre, Matrix& mObs);
double pNorm(Matrix v);
double scalarMult(const Matrix& m1,const Matrix& m2);
//double Norm(Matrix v);

/* calculate gain function as proposed by the article*/
double gainFunc(  Matrix& mPre,  Matrix& mObs)
{
	Matrix res(mPre.RowNo(),1);
	cout << "Observed = " << "(" << mObs(0,0) << "," << mObs(1,0) << ")"  << endl;
	cout << "Predicted = " << "(" << mPre(0,0) << "," << mPre(1,0) << ")"  << endl;

	double result =  A * ( 0.5 +  (scalarMult(mPre,mObs) /  (2*pNorm(mPre) *pNorm(mObs))) ) + 
		(1-A) * ( 1 -  (pNorm(mPre - mObs) /  (sqrt(pow(SX,N) + pow(SY,N)))) );
	cout  << "gain returned: " << result << endl;

	
	return result;

}

/* calculate p-norm for vectors*/
double pNorm(Matrix v)
{
	double sum=0;
	for(int i=0;i<v.RowNo();i++){
		sum += pow( abs( v(i,0) ) , N );	
	}
	sum = pow(sum,1/N);
	return sum;
}

//double Norm(Matrix v)
//{
//	double sum=0;
//	for(int i=0;i<v.RowNo();i++){
//		sum += abs(v(i,0)) * abs(v(i,0))  ;	
//	}
//	sum = sqrt(sum);
//	return sum;
//}

/*scalar multiply for vectors*/
double scalarMult(const Matrix& m1,const Matrix& m2)
{
	double sum=0;
	for(int i=0;i<2;i++){
		sum+=m1(i,0) * m2(i,0);

	}
	return sum;

}

#endif
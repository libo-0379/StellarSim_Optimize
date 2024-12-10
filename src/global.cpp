/*
 * =====================================================================================
 *
 *       Filename:  global.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  02/16/22 15:47:15
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  monkey
 *          email:  monkey@icloud.com
 *   Organization:  
 *
 * =====================================================================================
 */
#include <stdlib.h>

int N = 10;
double n = 1.0;
double t;
double tEnd;
double dt;
double M = 2.0;
double R = 0.75;
double m = 0.0;
	

double** pos = NULL;
double** vel = NULL;
double** acc = NULL;
double*  rho = NULL;    //密度
double*    P = NULL;    
double lmbda = 0.0;

double** dx = NULL;
double** dy = NULL;
double** dz = NULL;
double**  r = NULL;
double**  W = NULL;
double** wx = NULL;
double** wy = NULL;
double** wz = NULL;
double** gradPara = NULL;

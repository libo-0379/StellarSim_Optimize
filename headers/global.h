/*
 * =====================================================================================
 *
 *       Filename:  global.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  02/16/22 15:29:15
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  monkey
 *   Organization:  
 *          email:  monkey@icloud.com
 *
 * =====================================================================================
 */

#ifndef _GLOBAL_H
#define _GLOBAL_H
extern int N;
extern double n;
extern double t;
extern double tEnd;
extern double dt;
extern double M;
extern double R;
extern double** pos;
extern double** vel;
extern double** acc;
extern double*  rho;
extern double*    P;
extern double lmbda;
extern double m;
extern double** dx;
extern double** dy;
extern double** dz;
extern double**  r;
extern double**  W;
extern double** wx;
extern double** wy;
extern double** wz;
extern double** gradPara;

const double h  = 0.1;
const double k  = 0.1;
const double nu = 1.0;
const bool	 plotRealTime = true;
const double pi = 3.1415926;

// #define OPT_BASE //基础优化
#define OPT_OMP //omp优化
// #define OPT_SIMD //simd优化

#ifdef OPT_SIMD
#include <arm_neon.h>
#endif
#ifdef OPT_OMP
#include <omp.h>
#endif
#endif

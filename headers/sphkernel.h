/*
 * =====================================================================================
 *
 *       Filename:  sphkernel.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  02/27/22 16:24:41
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  monkey
 *   Organization:  
        *   email:  monkey@icloud.com
 *
 * =====================================================================================
 */

#ifndef _SPHKERNEL_H
#define _SPHKERNEL_H
#include "global.h"

extern void getPairwiseSeparations(double** &pos);
extern void getW(double** &dx, double** &dy, double** &dz, const double h);
extern void getGradW(double** &dx, double** &dy, double** &dz, const double h);
extern void getDensity(double** &pos, double &m, const double h);
extern void getPressure(double* &rho, const double k, double &n);
extern void getAcc(double** &pos, double** &vel, double &m, const double h, const double k, double &n, double lmbda, const double nu);

#ifdef OPT_SIMD
#define EPSILON 1e-6
extern float64_t exp_(float64_t x);
extern float64_t pow_(float64_t a,float64_t b);
#endif

#endif

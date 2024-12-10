/*
 * =====================================================================================
 *
 *       Filename:  initial.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  02/27/22 02:47:42
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  monkey
 *          email:  monkey@icloud.com
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <assert.h>
#include "../headers/global.h"
using namespace std;

double* assignParticleMem(int Nparticle)
{
  double* ptr = new double[Nparticle];
  return ptr;

}

void allocMEM()
{
  pos = new double*[3]; assert(pos!=NULL);
  vel = new double*[3]; assert(vel!=NULL);
  acc = new double*[3]; assert(acc!=NULL);
  rho = new double [N]; assert(rho!=NULL);
    P = new double [N]; assert(  P!=NULL);
   dx = new double*[N]; assert( dx!=NULL);
   dy = new double*[N]; assert( dy!=NULL);
   dz = new double*[N]; assert( dz!=NULL);
    r = new double*[N]; assert(  r!=NULL);
    W = new double*[N]; assert(  W!=NULL);
   wx = new double*[N]; assert( wx!=NULL);
   wy = new double*[N]; assert( wy!=NULL);
   wz = new double*[N]; assert( wz!=NULL);
  gradPara = new double*[N]; assert( gradPara!=NULL);
  for (int i = 0; i < 3; i++) {
    pos[i] = assignParticleMem(N); assert(pos[i]!=NULL);
    vel[i] = assignParticleMem(N); assert(vel[i]!=NULL);
    acc[i] = assignParticleMem(N); assert(acc[i]!=NULL);
  }
  // #pragma omp paralle for
  for (int i = 0; i < N; i++) {
    dx[i] = assignParticleMem(N); assert( dx[i]!=NULL);
    dy[i] = assignParticleMem(N); assert( dy[i]!=NULL);
    dz[i] = assignParticleMem(N); assert( dz[i]!=NULL);
     r[i] = assignParticleMem(N); assert(  r[i]!=NULL);
     W[i] = assignParticleMem(N); assert(  W[i]!=NULL);
    wx[i] = assignParticleMem(N); assert( wx[i]!=NULL);
    wy[i] = assignParticleMem(N); assert( wy[i]!=NULL);
    wz[i] = assignParticleMem(N); assert( wz[i]!=NULL);
    gradPara[i] = assignParticleMem(N); assert( gradPara[i]!=NULL);
  }
}

void initialize()
{
  srand(42);
  int randS, randE;
  randS = 0; randE = 2;
  // #pragma omp paralle for
  for (int i = 0; i < N; i++) {
    pos[0][i] = pow(-1, (rand()%(randS - randE) + randS +1))*rand() / double(RAND_MAX);
    pos[1][i] = pow(-1, (rand()%(randS - randE) + randS +1))*rand() / double(RAND_MAX);
    pos[2][i] = pow(-1, (rand()%(randS - randE) + randS +1))*rand() / double(RAND_MAX);
    vel[0][i] = 0.0;
    vel[1][i] = 0.0;
    vel[2][i] = 0.0;
    rho[i]    = 0.0;
      P[i]    = 0.0;
  // #pragma omp paralle for
    for (int j = 0; j < N; j++) {
      dx[i][j] = 0.0;
      dy[i][j] = 0.0;
      dz[i][j] = 0.0;
       r[i][j] = 0.0;
       W[i][j] = 0.0;
      wx[i][j] = 0.0;
      wy[i][j] = 0.0;
      wz[i][j] = 0.0;
      gradPara[i][j] = 0.0;
    }
  }
}

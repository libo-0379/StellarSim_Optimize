/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:
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
#include <iomanip>
#include <math.h>
#include <cmath>
#include "../headers/global.h"
#include "../headers/initial.h"
#include "../headers/sphkernel.h"
#include <sys/time.h>
#include <fstream>
#include <iomanip>

using namespace std;

int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    cout << "SPHexe [number of particles:int]" << endl;
    exit(0);
  }

  #if 0
    printf("%f\n",exp_(5.0));

    return 0;
  #endif


  //粒子个数
  N = atoi(argv[1]);
  cout << "StellarSim performs on " << N << "particles" << endl;
  allocMEM();
  lmbda = 2.0 * k * (1.0 + n) * pow(pi, (-3.0 / (2.0 * n))) * (M * tgamma(5.0 / 2.0 + n) / pow(R, 3.0) / pow(tgamma(1.0 + n), (1.0 / n))) / pow(R, 2.0);
  cout << "lmbda is " << lmbda << endl;
  m = M / N;
  #ifdef OPT_OMP
  omp_set_num_threads(32);
  #endif

  initialize();
  getAcc(pos, vel, m, h, k, n, lmbda, nu);

  t = 0.0;
  tEnd = 12;
  dt = 0.04;
  //迭代次数
  int Nt = int(ceil(tEnd / dt));
  
  struct timeval tpstart, tpend;
  gettimeofday(&tpstart, NULL);
  for (int i = 0; i <= Nt; i++)
  {
    for (int j = 0; j < N; j++)
    {
      vel[0][j] += acc[0][j] * dt / 2;
      vel[1][j] += acc[1][j] * dt / 2;
      vel[2][j] += acc[2][j] * dt / 2;
      pos[0][j] += vel[0][j] * dt;
      pos[1][j] += vel[1][j] * dt;
      pos[2][j] += vel[2][j] * dt;
    }

    getAcc(pos, vel, m, h, k, n, lmbda, nu);

    for (int j = 0; j < N; j++)
    {
      vel[0][j] += acc[0][j] * dt / 2;
      vel[1][j] += acc[1][j] * dt / 2;
      vel[2][j] += acc[2][j] * dt / 2;
    }

    cout << "Time Step is: " << setw(4) << t << endl;
    
    t += dt;
    getDensity(pos, m, h);
  }

  gettimeofday(&tpend, NULL);
  double timeuse = 1000000 * (tpend.tv_sec - tpstart.tv_sec) + tpend.tv_usec - tpstart.tv_usec;
  timeuse /= 1000000;
  cout << "elapsed time: " << setprecision(3) << setw(6) << timeuse << " s" << endl;

  // print acc[][]
  ofstream file("acc_value.txt");
  file.setf(ios::scientific);
  for (int j = 0; j < N; j++)
  {
    file << setprecision(8) << setw(18) << acc[0][j] << ","
         << setw(18) << acc[1][j] << "," << setw(18) << acc[2][j] << endl;
  }
  file.close();

  return 0;
}

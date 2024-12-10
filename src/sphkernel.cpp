/*
 * =====================================================================================
 *
 *       Filename:  sphkernel.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  02/27/22 14:55:57
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  monkey
 *   Organization:  
 *          email:  monkey@icloud.com
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "../headers/global.h"
#include "../headers/sphkernel.h"

using namespace std;
 
//计算距离
void getPairwiseSeparations(double** &pos)
{
  // 1. 提出循环不变量
  // 2. todo simd
  // 3. 对 j 进行循环分块 影响dx dy dz的命中，不适用
  // 4. omp
  #if defined(OPT_BASE) && (defined(OPT_SIMD)||defined(OPT_OMP))
    #ifdef OPT_OMP
    #pragma omp parallel for schedule(guided) proc_bind(close) 
    #endif
    for (int i = 0; i < N; i++) 
    {
      #ifndef OPT_SIMD
        double temp1 = pos[0][i];
        double temp2 = pos[1][i];
        double temp3 = pos[2][i];
        for (int j = 0; j < N; j++) 
        {
          // dx[i][j] = -dx[j][i] 粒子彼此计算相对距离
          dx[i][j] = pos[0][j] - temp1;
          dy[i][j] = pos[1][j] - temp2;
          dz[i][j] = pos[2][j] - temp3;
        }
      #else
        float64x2_t v0 = vld1q_dup_f64(&pos[0][i]);
        float64x2_t v1 = vld1q_dup_f64(&pos[1][i]);
        float64x2_t v2 = vld1q_dup_f64(&pos[2][i]);
        for(int j=0;j<N/2*2;j+=2)
        {
            float64x2_t v0_0 = vld1q_f64(&pos[0][j]);
            float64x2_t v1_0 = vld1q_f64(&pos[1][j]);
            float64x2_t v2_0 = vld1q_f64(&pos[2][j]);
            vst1q_f64(&dx[i][j],vsubq_f64(v0_0,v0));
            vst1q_f64(&dy[i][j],vsubq_f64(v1_0,v1));
            vst1q_f64(&dz[i][j],vsubq_f64(v2_0,v2));
        }
        for (int j = N/2*2; j < N; j++) 
        {
          dx[i][j] = pos[0][j] - pos[0][i];
          dy[i][j] = pos[1][j] - pos[1][i];
          dz[i][j] = pos[2][j] - pos[2][i];
        }
      #endif
    }
  #else
    #ifdef OPT_OMP
    #pragma omp parallel for schedule(guided) proc_bind(close) 
    #endif
    for (int i = 0; i < N; i++) 
    {
      for (int j = 0; j < N; j++) 
      {
        // dx[i][j] = -dx[j][i] 粒子彼此计算相对距离
        dx[i][j] = pos[0][j] - pos[0][i];
        dy[i][j] = pos[1][j] - pos[1][i];
        dz[i][j] = pos[2][j] - pos[2][i];
        //fprintf(stdout, "%12.6f", dz[i][j]);
        //fflush(stdout);
      }
    //fprintf(stdout,"\n");
    }
  #endif
}

void getW(double** &dx, double** &dy, double** &dz, const double h)
{
  // 1. 循环不变量提出
  // 2. omp
  #if defined(OPT_OMP) || defined(OPT_BASE)
    double value1 = pow((1.0 / (h*sqrt(pi))), 3.0);
    double value2 = pow(h,2);
    #ifdef OPT_OMP
    #pragma omp parallel for schedule(guided) proc_bind(close)
    #endif
    for (int i = 0; i < N; i++) 
    {
      for (int j = 0; j < N; j++) 
      {
        r[i][j] = sqrt(pow(dx[i][j],2.0) + pow(dy[i][j],2.0) + pow(dz[i][j],2.0));
        W[i][j] = value1 * exp((-pow(r[i][j],2) / value2)); 
      }
    }
  #else   
    #ifdef OPT_OMP
    #pragma omp parallel for schedule(guided) proc_bind(close) 
    #endif
    for (int i = 0; i < N; i++) 
    {
      for (int j = 0; j < N; j++)
      {
        r[i][j] = sqrt(pow(dx[i][j],2.0) + pow(dy[i][j],2.0) + pow(dz[i][j],2.0));
        W[i][j] = pow((1.0 / (h*sqrt(pi))), 3.0) * exp((-pow(r[i][j],2) / pow(h,2))); 
        //fprintf(stdout, "%12.6f", r[i][j]);
        //fprintf(stdout, "%12.6f", W[i][j]);
        //fflush(stdout);
      }
    }
  #endif
    //fprintf(stdout,"\n");
}

void getGradW(double** &dx, double** &dy, double** &dz, const double h)
{
  // 1. 循环不变量提出
  // 2. omp
  #if defined(OPT_OMP) || defined(OPT_BASE)
    double value1 = pow(h,2);
    double value2 = -2/pow(h,5)/pow(pi,(3/2));
    #ifdef OPT_OMP
    #pragma omp parallel for schedule(guided) proc_bind(close)
    #endif
    for (int i = 0; i < N; i++) 
    {
      for (int j = 0; j < N; j++) 
      {
        r[i][j]  = sqrt(pow(dx[i][j],2.0) + pow(dy[i][j],2.0) + pow(dz[i][j],2.0));
        gradPara[i][j] = exp(-pow(r[i][j],2) / value1) * value2;
        wx[i][j] = gradPara[i][j]*dx[i][j];
        wy[i][j] = gradPara[i][j]*dy[i][j];
        wz[i][j] = gradPara[i][j]*dz[i][j];
      }
    }
  #else
    #ifdef OPT_OMP
    #pragma omp parallel for schedule(guided) proc_bind(close) 
    #endif
    for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      // r[i][j] = r[j][i] 
      r[i][j]  = sqrt(pow(dx[i][j],2.0) + pow(dy[i][j],2.0) + pow(dz[i][j],2.0));
      // gradPara[i][j] = gradPara[j][i]
      gradPara[i][j] = -2 * exp(-pow(r[i][j],2) / pow(h,2)) / pow(h,5) / pow(pi,(3/2));
      // wx[i][j] = -wx[j][i]
      wx[i][j] = gradPara[i][j]*dx[i][j];
      wy[i][j] = gradPara[i][j]*dy[i][j];
      wz[i][j] = gradPara[i][j]*dz[i][j];
      //fprintf(stdout, "%12.6f", wy[i][j]);
      //fflush(stdout);
    }
    //fprintf(stdout,"\n");
  }
  #endif
}

void getDensity(double** &pos, double &m, const double h)
{
  getPairwiseSeparations(pos);
  getW(dx, dy, dz, h);
  // 1. todo 访问顺序 
  // 2. 内层 simd
  // 3. 外层omp有数据竞争,内层omp可能伪共享, omp 不适用
  // 4. 简化计算公式 rho[j] += W[i][j] rho[i] *= m；W 每一列计算出一个 rho[j] 此处不适用
  #ifdef OPT_BASE
    for(int i = 0; i < N; i++) 
    {
      for(int j = 0; j < N; j++) 
        rho[j] += m * W[i][j];
    }
  #else
    for (int j = 0; j < N; j++) 
    {
      for (int i = 0; i < N; i++) 
        rho[j] += m * W[i][j];
      //fprintf(stdout, "%12.6f", rho[j]);
      //fflush(stdout);
      //fprintf(stdout,"\n");
    }
  #endif
}

void getPressure(double* &rho, const double k, double &n)
{
  // 1. 提出循环不变量,循环展开
  // 2. omp 线程调度开销大，不合适
  // 3. simd pow较复杂
  #ifdef OPT_BASE
    double value = 1+1/n;
    for (int j = 0; j < N; j++) 
      P[j] = k * pow(rho[j], value);
  #else
    for (int j = 0; j < N; j++) 
    {
      P[j] = k * pow(rho[j], (1+1/n));
      //fprintf(stdout, "%12.6f\n", P[j]);
    }
  #endif
  
}

void getAcc(double** &pos, double** &vel, double &m, const double h, const double k, double &n, double lmbda, const double nu)
{
  getDensity(pos, m, h);
  getPressure(rho, k, n);
  getPairwiseSeparations(pos);
  getGradW(dx, dy, dz, h);
  #if defined(OPT_BASE)
  #ifdef OPT_OMP
  #pragma omp parallel for schedule(guided) proc_bind(close)
  #endif
  for (int j = 0; j < N; j++) 
  {
    // 1. wx[i][j] = -wx[j][i] 访问 w[j][i] 可以增加缓存命中，并且可以向量化
    // 如果for循环交换 j i访问顺序，先访问i 后 j，内层循环无法做向量化优化(内循环每次计算不同的目标元素)
    // 并且，如果for循环先访问j,计算 acc[0][j], P[j]/pow(rho[j],2)是常量可以提出最后计算
    // 2. 简化计算 m* 放在求和之后
    // 3. 循环不变量提出
    double temp1 = P[j]/pow(rho[j],2);
    double temp3 =0.0,temp4=0.0,temp5=0.0;
    for (int i = 0; i < N; i++)
    {
      double temp2 = pow(P[i]/rho[i],2);
      temp3 += (temp1 + temp2) * wx[j][i];
      temp4 += (temp1 + temp2) * wy[j][i];
      temp5 += (temp1 + temp2) * wz[j][i];
    }
    acc[0][j] += (temp3 *=m);
    acc[1][j] += (temp4 *=m);
    acc[2][j] += (temp5 *=m);
  }
  #else
  #ifdef OPT_OMP
  #pragma omp parallel for schedule(guided) proc_bind(close)
  #endif
  for (int j = 0; j < N; j++) 
  {
    for (int i = 0; i < N; i++) 
    {
      acc[0][j] -= m * ( P[j]/pow(rho[j],2) + pow(P[i]/rho[i],2)  ) * wx[i][j];
      acc[1][j] -= m * ( P[j]/pow(rho[j],2) + pow(P[i]/rho[i],2)  ) * wy[i][j];
      acc[2][j] -= m * ( P[j]/pow(rho[j],2) + pow(P[i]/rho[i],2)  ) * wz[i][j];
    }
  }
  #endif
  // 1. simd
  // 2. 循环合并
  #ifdef OPT_BASE
  for (int j = 0; j < N; j++) 
  {
    acc[0][j] -= (lmbda * pos[0][j] + nu * vel[0][j]); 
    acc[1][j] -= (lmbda * pos[1][j] + nu * vel[1][j]); 
    acc[2][j] -= (lmbda * pos[2][j] + nu * vel[2][j]); 
  }
  #else
  for (int j = 0; j < N; j++) 
  {
    acc[0][j] -= lmbda * pos[0][j]; 
    acc[1][j] -= lmbda * pos[1][j];
    acc[2][j] -= lmbda * pos[2][j];
  }
  for (int j = 0; j < N; j++) 
  {
    acc[0][j] -= nu * vel[0][j];
    acc[1][j] -= nu * vel[1][j];
    acc[2][j] -= nu * vel[2][j];
  }
  #endif
}

#ifdef OPT_SIMD
// 需要用到 泰勒展开
float64_t exp_(float64_t x)
{
  //初始化第一个值
  int n = 0;
  double prior = 1.0;
  double sum = prior; //求和保存结果
  while(1)
  {
    double cur = prior * x /++n;
    sum += cur;
    prior = cur;
    if(cur<=EPSILON)
      break;
  }

  return sum;
}

// a^b = e^(b*ln(a)); neon 未提供ln，设想采用 cmath ln函数,向量化对每个元素的计算用 omp task
// float64_t pow_(float64_t a,float64_t b)
// {
//   logf()
// }

#endif
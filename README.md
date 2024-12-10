# StellarSim_Optimize
使用循环、openmp、simd对粒子移动状态进行模拟的计算优化。

1. 文件后缀说明
original: 源码
BASE：循环优化
OMP：openmp 多线程优化
SIMD：使用arm neon 指令集向量化优化

|项目           |耗时 s|相比源码加速比|
|-------------|----|-------|
|original(源码) |141 |       |
|BASE(循环优化)   |58.7|2.4    |
|BASE+SIMD    |51  |2.8    |
|BASE+OMP     |10  |14.1   |
|OMP          |23.1|6.1    |
|BASE+OMP+SIMD|6.47|21.8   |

2. 结果分析
2.1 程序内计时
   (1) 单核串行耗时从 141s 优化为 51s，加速比为 2.8；
   (2) 多核32线程耗时 23.1s(实测16线程耗时一致)，加速比为 6.1；
   (3) 综合优化后耗时 6.47，加速比为 21.8。
2.2 其他方面
   (1) 多线程优化方面 schedule(guided) 策略拥有最高的效率，比 dynamic 略优。
   (2) 基础优化和多线程有较大的加速比提升。
   (3) 16线程与32线程的效率一致。
   (4) 注意多线程数据竞争问题，可能由于此问题造成多线程效率下降。
   (5) 注意cacheline 竞争问题，避免多线程效率下降和计算错误。
2.3 gprof分析
   (1) 使用openmp 会明显增加线程框架开销，每增加一个线程使用会增加开销；
   (2) getAcc 函数从源码的时间从 36.9s 最终将为几乎为 0;
   (3) 函数计算的耗时比例从 100% 降为 15.9%，其余 83.78% 为多线程带来的管理开销。


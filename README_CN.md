# driftScalarDyFoam: 基于OpenFOAM的风致雪漂求解器

## 介绍

[English](./README.md) | [简体中文](./README_CN.md)

本求解器使用标量输运方程模拟雪相在湍流大气层中的输运[<sup>[1]](#refer-1)[<sup>[2]](#refer-2)[<sup>[3]](#refer-3)

![\frac{\partial \phi}{\partial t}+\frac{\partial \phi u_{j}}{\partial x_{j}}+\frac{\partial \phi w_{f}}{\partial x_{3}}=\frac{\partial}{\partial x_{j}} D_{t} \frac{\partial \phi}{\partial x_{j}}](http://www.sciweavers.org/tex2img.php?eq=1%2Bsin%28mc%5E2%29&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=)

并基于侵蚀/沉积方程调整积雪表面网格

$$q_{\text{dep}} = -\phi_{p} w_f \left(1-\frac{u_{*}^{2}}{u_{*t}^{2}}\right),\quad u_*<u_{*t}$$

$$q_{\text{ero}}=-A_\text{ero} u_{*}^{2}\left(1-\frac{u_{*t}^{2}}{u_*^{2}}\right), \quad u_*>u_{*t}$$


## 相关论文

如您在您的研究中使用或借鉴了该项目，请引用：(相关论文正在投稿)

## 鸣谢

本项目受到以下项目的资助：

1. 大跨度钢结构屋盖的风雪流作用及响应（国家自然科学基金（面上项目），51378428）

2. 建筑风雪流相界自适应模型与分阶段准动态耦合分析理论（国家自然科学基金（面上项目），52178506）
## 参考文献
<div id="refer-1"></div>
[1] Tominaga Y, Okaze T, Mochida A. CFD modeling of snowdrift around a building: An overview of models and evaluation of a new approach[J]. Building and Environment, 2011, 46(4): 899-910.

<div id="refer-2"></div>
[2] Zhou X, Kang L, Gu M, et al. Numerical simulation and wind tunnel test for redistribution of snow on a flat roof[J]. Journal of Wind Engineering and Industrial Aerodynamics, 2016, 153: 92-105.

<div id="refer-3"></div>
[3] Zhu F, Yu Z, Zhao L, et al. Adaptive-mesh method using RBF interpolation: A time-marching analysis of steady snow drifting on stepped flat roofs[J]. Journal of Wind Engineering and Industrial Aerodynamics, 2017, 171: 1-11.
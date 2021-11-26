# driftScalarDyFoam: 基于OpenFOAM的风致雪漂求解器

## [简体中文](./README_CN.md) | [English](./README.md)
## 介绍

本求解器使用标量输运方程模拟雪相在湍流大气层中的输运[<sup>[1]](#refer-1)[<sup>[2]](#refer-2)[<sup>[3]](#refer-3)

<a href="https://www.codecogs.com/eqnedit.php?latex=\bg_white&space;\color{Red}\frac{\partial&space;\phi}{\partial&space;t}&plus;\frac{\partial&space;\phi&space;u_{j}}{\partial&space;x_{j}}&plus;\frac{\partial&space;\phi&space;w_{f}}{\partial&space;x_{3}}=\frac{\partial}{\partial&space;x_{j}}&space;D_{t}&space;\frac{\partial&space;\phi}{\partial&space;x_{j}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\bg_white&space;\color{Red}\frac{\partial&space;\phi}{\partial&space;t}&plus;\frac{\partial&space;\phi&space;u_{j}}{\partial&space;x_{j}}&plus;\frac{\partial&space;\phi&space;w_{f}}{\partial&space;x_{3}}=\frac{\partial}{\partial&space;x_{j}}&space;D_{t}&space;\frac{\partial&space;\phi}{\partial&space;x_{j}}" title="\color{Red}\frac{\partial \phi}{\partial t}+\frac{\partial \phi u_{j}}{\partial x_{j}}+\frac{\partial \phi w_{f}}{\partial x_{3}}=\frac{\partial}{\partial x_{j}} D_{t} \frac{\partial \phi}{\partial x_{j}}" /></a>

并基于侵蚀/沉积方程调整积雪表面网格

<a href="https://www.codecogs.com/eqnedit.php?latex=\bg_white&space;\color{Red}q_{\text{dep}}&space;=&space;-\phi_{p}&space;w_f&space;\left(1-\frac{u_{*}^{2}}{u_{*t}^{2}}\right),\quad&space;u_*<u_{*t}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\bg_white&space;\color{Red}q_{\text{dep}}&space;=&space;-\phi_{p}&space;w_f&space;\left(1-\frac{u_{*}^{2}}{u_{*t}^{2}}\right),\quad&space;u_*<u_{*t}" title="\color{Red}q_{\text{dep}} = -\phi_{p} w_f \left(1-\frac{u_{*}^{2}}{u_{*t}^{2}}\right),\quad u_*<u_{*t}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\bg_white&space;\color{Red}q_{\text{ero}}=-A_\text{ero}&space;u_{*}^{2}\left(1-\frac{u_{*t}^{2}}{u_*^{2}}\right),&space;\quad&space;u_*>u_{*t}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\bg_white&space;\color{Red}q_{\text{ero}}=-A_\text{ero}&space;u_{*}^{2}\left(1-\frac{u_{*t}^{2}}{u_*^{2}}\right),&space;\quad&space;u_*>u_{*t}" title="\color{Red}q_{\text{ero}}=-A_\text{ero} u_{*}^{2}\left(1-\frac{u_{*t}^{2}}{u_*^{2}}\right), \quad u_*>u_{*t}" /></a>

## 版本适配
本项目的主分支面向OpenFOAM的ESI分支开发([https://www.openfoam.com](https://www.openfoam.com))，基金会分支尚在开发中。

目前已验证的版本适配包括：

driftScalarDyFoam-master - [OpenFOAM-v2012](https://www.openfoam.com/download/release-history), [OpenFOAM-v2106](https://www.openfoam.com/download/release-history)

## 下载与安装

```shell
cd $WM_PROJECT_DIR  #进入OpenFOAM根目录
mkdir extend        #建立一个任意名称的空文件夹
cd extend

git clone https://github.com/fightingxiaoxiao/driftScalarDyFoam.git

cd driftScalarDyFoam
wmake #编译项目
```

## 参考算例

[平屋面积雪](./tutorials/flatRoof3D)

[建筑周边积雪](./tutorials/snowAroundBuilding)

## 快速开始

driftScalarDyFoam参考scalarTransportFoam、simpleFoam及moveDynamicMesh开发。其算例框架与标准的OpenFOAM算例是一致的。需要额外补充的内容包括：

### erosionDepositionProperties

该字典应当放置在constant文件夹下，需要声明的参数包括：

```cpp
ca              7e-4;   // 即A_ero，控制侵蚀方程的常数
rhoSnow         150;    // 雪的堆积密度         
rhoAir          1;      // 空气密度，在不可压求解器中恒为1
Uthreshold      0.2;    // 阈值剪切风速

UResidual       (2e-4 5e-4 5e-4);   // 子循环中的速度残差阈值
pResidual       5e-4;               // 子循环中的压力残差阈值
TResidual       5e-6;               // 子循环中的浓度残差阈值
nSubCycles      1000;               // 每个计算阶段内的最大子循环数
```

### 雪漂浓度边界

### 动网格

## 相关论文

如您在研究中使用或借鉴了该项目，请引用：(相关论文正在投稿)

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
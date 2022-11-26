# driftScalarDyFoam: 基于OpenFOAM的风致雪漂求解器

## [简体中文](./README.md) | [English](./README_EN.md)
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

积雪表面的浓度根据当前的侵蚀量给出一个Neumann边界条件[<sup>[1]](#refer-1)：

<a href="https://www.codecogs.com/eqnedit.php?latex=\bg_white&space;\color{Red}&space;\left.\left(\frac{\partial&space;\phi}{\partial&space;\mathbf{n}}\right)\right|_{\text&space;{surface}}=&space;-\frac{1}{D_{t}}&space;q_\text{ero}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\bg_white&space;\color{Red}&space;\left.\left(\frac{\partial&space;\phi}{\partial&space;\mathbf{n}}\right)\right|_{\text&space;{surface}}=&space;-\frac{1}{D_{t}}&space;q_\text{ero}" title="\bg_white \color{Red} \left.\left(\frac{\partial \phi}{\partial \mathbf{n}}\right)\right|_{\text {surface}}= -\frac{1}{D_{t}} q_\text{ero}" /></a>

在数值计算中，该Neumann边界条件被离散为如下形式：

<a href="https://www.codecogs.com/eqnedit.php?latex=\color{Red}&space;\phi_{face}&space;=&space;\phi_{cell}&space;-&space;\Delta&space;\frac{q_\text{ero}}{D_{t}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\color{Red}&space;\phi_{face}&space;=&space;\phi_{cell}&space;-&space;\Delta&space;\frac{q_\text{ero}}{D_{t}}" title="\color{Red} \phi_{face} = \phi_{cell} - \Delta \frac{q_\text{ero}}{D_{t}}" /></a>

由于原生的OpenFOAM没有类似的边界条件，因此现有的方案是采用`codedFixedValue`边界条件在计算时进行动态编译,相关的内容可见算例中的`0.orig/T`及`system/codeDict`的`erosionFlux`代码段。

此外，为了确保边界上的质量交换率能被正确更新，雪面的边界名称应当附加".snow"后缀，如"roof.snow"。
### 动网格

根据每个面网格的雪质量交换率，driftScalarDyFoam创建了一个名为`deltaH`的面标量场来监控雪面的高度变化。在`0.orig/pointMotionU`及`system/codeDict`中，我们利用OpenFOAM原有的场映射组件，将面标量场映射为节点向量场，并同样采用`codedFixedValue`边界指定每个时间步（或称阶段）中的边界节点位移速度。因此，使用者应当结合自己的模型特征对`0.orig/pointMotionU`及`system/codeDict`中的`erosionDeposition`代码段进行改动。

## 相关论文

如您对我们的项目感兴趣，可以在GitHub中点击Star以进行支持和持续关注。

如您在研究中使用或借鉴了该项目，请引用：

[1] Chen X and Yu Z (2022) DriftScalarDyFoam: An OpenFOAM-Based Multistage Solver for Drifting Snow and Its Distribution Around Buildings. Front. Earth Sci. 10:822140. doi: 10.3389/feart.2022.822140


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
本篇主要介绍去膜的基本思路和流程。详细流程请见：
* [1-得到膜中心位于图像中心的2D average](Tutorials/Obtain-2D-Averages_zh-CN.md)
* [2-使用Radonfit方法去除膜信号](Tutorials/Radonfit_tutorials/index.md)
* [3-使用Bezierfit方法去除膜信号](Tutorial/Bezierfit_tutorials/index.md)

# 基本思路

1. 使用cryoSPARC挑出生物膜位于图像中心的particles以及对应的2D average；
2. 对于每一个2D average，使用Radon变换、贝塞尔曲线拟合等方法，得到该类生物膜的函数表达式，包含了膜中心、角度、曲率等信息；
3. 根据该类生物膜的函数表达式以及cryoSPARC中的alignment信息，对该类生物膜中所有的raw image进行去膜操作，得到去膜后的particle stack；
4. 将去膜后的particle stack放回原来的micrograph中，得到去膜后的micrograph，进行后续分析。


# 基本流程

## 1 得到膜中心位于图像中心的2D average
目的主要是使得第二步的膜分析更加准确和容易，因为生物膜的膜中心如果不在图像中心，特别是与图像中心偏差较大时，Radon变换有可能分析不出来。

### 1.1 生成一系列膜中心位于图像中心的细胞膜2D average模版

从cryoSPRAC中选择一个信噪比比较好的、生物膜信号占主导的2D average，根据它生成一系列2D average模版，这些模版中的生物膜的膜中心都是在图像中心，方向也是竖直方向，只有曲率不同。

<center><img src="../img/kappa-templates-image.gif" alt="Kappa templates"></center>


用它生成的一系列模版，可以再作为templates，再次pick particles，这样得到的particles的中心就都能准确定位在细胞膜中心了。



### 1.2 Extract particles & 2D classification
将上述particles提取出来，进行2D classification，得到的2D averages，应该都是膜中心位于图像中心的。在使用这些2D averages进行后续分析。

## 2 分析2D average，得到相应的函数表达式

### 2.1 Radonfit

这种方法主要是利用Radon变换、cross-correlation等方法，以**简单的直线和圆弧为模型**，拟合2D average，得到函数表达式，包含膜中心、角度、曲率等信息。适用于**较为简单的生物膜模型**，比如病毒的包膜。

#### 2.1.1 Radon Analysis Blinking

逐个分析确定每个2D average能够进行Radon变换的参数，包括crop rate、阈值、Radon变换角度范围，最后将参数记录到JSON文件中。这些参数将用在下一步对膜的分析。

#### 2.1.2 Membrane Analysis
进行2D averages的分析，对于2D average，每一个得到膜中心、角度、脂双层距离、单层厚度（高斯拟合后的sigma大小）、以及曲率。这些信息也就对应了函数表达式，利用表达式可以得到平均后的生物膜以及对应的蒙版。

### 2.2 Bezierfit

这种方法主要是利用结合了蒙特卡洛方法、遗传算法等方法的贝塞尔曲线来对2D average中的生物膜进行拟合，以**更加复杂的不规则曲线为模型**，拟合2D average，得到若干个控制点以及对应的函数解析式。适用于**形状较为复杂的生物膜模型**，例如S型、W型等，比如线粒体膜。

#### 2.2.1 使用蒙特卡洛方法在膜区域随机生成若干个点
对于每一个2D average，首先使用最大值滤波器，将2D average中大致的生物膜区域提取出来，然后在这个区域随机生成若干个点，这些点将作为首次拟合贝塞尔曲线的参考点。

#### 2.2.2 使用遗传算法进行初步拟合
对上一步生成的若干点，使用遗传算法进行初步拟合，得到若干个控制点，以及对应的函数解析式。通过这一步，根据上一步若干点的位置信息，可以得到膜的大致走向及形状。

#### 2.2.3 使用遗传算法调整控制点，优化拟合结果
对上一步初步拟合得到的若干个控制点，再次使用遗传算法进行优化。通过最大化cross-correlation的值，来调整控制点的位置，得到拟合2D average中膜的形状的若干个最优控制点，这些最优的控制点就能对应一个解析式，描绘膜的形状。

## 3 Membrane Subtraction
一类2D average对应了若干个raw particle stacks。通过之前对于2D average的分析，我们能够得到该类生物膜的函数表达式，也就得到了每个raw particle的生物膜信息。

对于每一个之前extract出来的raw particle，都进行去膜操作，得到对应的去膜后的particle，组成另一个对应的particle stack。这个particle stack就是去膜后的particle stack。

## 4 将去膜后的particle放回micrograph中
将相应的去膜后的particles取代micrographs中原来的particles，得到新的micrographs，其中膜信号就能得到减弱甚至去除。

## 5 使用膜蛋白模版在新的micrograph中pick particle
在新的去膜后的micrograph中，使用膜蛋白模版pick particles，进行后续计算。


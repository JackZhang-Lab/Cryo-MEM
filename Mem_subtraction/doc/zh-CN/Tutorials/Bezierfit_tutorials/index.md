# Membrane Subtraction Workflow Using Bezierfit

## 方法介绍

这种方法主要结合了蒙特卡洛方法、遗传算法等方法的贝塞尔曲线来对2D average中的生物膜进行拟合，以**更加复杂的不规则曲线为模型**，拟合2D average，得到若干个控制点以及对应的函数解析式。适用于**形状较为复杂的生物膜模型**，例如S型、W型等，比如线粒体膜。

## 具体流程

* [使用Bezierfit进行膜分析](./Bezierfit-Membrane-Analysis_zh-CN.md)；
* [使用Bezierfit对particles进行去膜](./Bezierfit-Membrane-Subtraction_zh-CN.md)；
* [将去膜后的particles放回原来的micrograph中，得到去膜后的micrograph](./Bezierfit-Micrograph-Membrane-Subtraction_zh-CN.md)。


# ax-b-interpolation-fitting
解决AX=B，插值以及拟合的一些程序

解决AX=B，插值以及拟合的一些程序，代码可见

1 拟合
1.1 法方程
由函数是对的最佳平方逼近函数的充分必要条件是与所有的正交，
1.2 超定方程
当待求变量个数少于方程等式数量时便会形成超定方程组
2 插值
插值的任务就是由已知的观测点(xi,yi)为物理量(未知量),建立一个简单的、连续的解析模型g(x) ，以便能根据该模型推测该物理量在非观测点处的特性。即由实验或测量的方法得到所求函数 y=f(x) 在互异点x0 , x1 , ... , xn 处的值 y0 , y1 , … , yn , 构造一个简单函数 F(x) 作为函数 y=f(x) 的近似表达式 y= f(x)  F(x) 使 F(x0)=y0 , F(x1)=y1 , , F(xn)=yn 插值条件成立，这类问题称为插值问题。 f(x) 称为被插值函数，F(x) 称为插值函数， x0 , x1 , ... , xn 称为插值节点。
2.1 高次Lagrange插值
基于Lagrange插值方法编写程序对测井曲线插值，。
2.2 高次Newton插值
值得注意的是，每增加一个结点，Newton插值多项式只增加一项，克服了 Lagrange插值的缺点。对差商的求取，更一般地是使用差商表：
使用牛顿插值方法较Lagrang方法快，但是由于是高次插值，其插值结果在区间边界仍然存在较大的浮动。
插值的目的就是数值逼近的一种手段，而数值逼近是为得到一个数学问题的精确解或足够精确的解。插值节点的增多,尽管使插值多项式在更多的插值节点上与函数 f(x) 的值相等,但在两个节点之间Pn(x)不一定能很好地逼近 f(x) , 有时误差会大得惊人。著名的龙格(Runge)现象证实了这个观点。从计算的数值运算误差看,对于等距节点的差分形式,由于高阶差分的误差传播,函数值的微小变化都将使插值产生很大的误差。龙格(Runge)现象表明插值多项式序列不收敛，实际上，严格的理论分析可知插值多项式序列确是不收敛的，而且高阶插值还是不稳定的。以上两种高次插值方法计算时间复杂度高，而且插值效果在边界表现地并不理想。
2.3 分段Lagrang插值
2.4 分段Newton插值

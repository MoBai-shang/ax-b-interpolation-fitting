X=[2 4 5 7 8 12 13 15]';%;3 5 6 9 11 16 17 19
Y=[5 17 8 20 14 12 10 8]';
%[data,text]  = xlsread('单井数据统计.xlsx');
%Y=str2double(text);
point_num=200;
x=zeros(point_num,size(X,2));
for k=1:size(X,2)
    x(:,k)=linspace(min(X(:,k)),max(X(:,k)),point_num)';
end
N=3;
[y0,formula]=fit_mutil_var_N(X,Y,N,X);
[y,formula]=fit_mutil_var_N(X,Y,N,x);
loss=sqrt(sum(power(Y-y0,2)));
if rank(X)==1
    plot(X,Y,'*');
    hold on
    plot(x,y);
    gtext(formula)
else
    %hold on
    plot3(X(:,1),X(:,2),Y,'r*')
    hold on
    plot3(x(:,1),x(:,2),y,'b*')
    %gtext(formula)
    zlabel('砂厚');
end
title( {[num2str(N),'次平方逼近'],['均方根误差 = ', num2str(loss)]})
legend('true','predict');
xlabel('h^{1/2}');
ylabel('砂厚');



%一元多次函数拟合（基于法方程GC=F）
%返回基于数据X,Y对x的拟合结果及拟合函数公式
function [res,formula]=fit_one_var_N(X,Y,N,x)
%{
parameters:
    X:abscissa of fit point
    Y:ordinate of fit point
    x:abscissa of point to be solved
%}
assert (sum(size(X)-size(Y))==0,'The dimensions of matrices X and Y must be the same');
assert (rank(X)==1,'X and Y must be one-dimensional matrices');
n=max(size(X));
N=N+1;
assert (N>0 & N <=n,'the interpolation times must bigger than 0 and smaller than length of matrix X')
f=zeros(N,n);
for k=1:N
    f(k,:)=power(X,k-1);
end
G=zeros(N,N);
for i=1:N
    for j=1:i
        G(i,j)=f(i,:)*f(j,:)';
        G(j,i)=G(i,j);
    end
end
F=zeros(N,1);
for k=1:N
    F(k)=f(k,:)*Y;
end
C=G\F;
res=ones(size(x,1),1)*C(1);
formula=num2str(C(1));
for k=2:N
    res=res+power(x,k-1)*C(k);
    formula=[formula,'+',num2str(C(k)),'x^{',num2str(k-1),'}'];
end
end

%多元N次函数拟合（基于法方程GC=F）
%返回基于数据X,Y对x的拟合结果及拟合函数公式
function [res,formula]=fit_mutil_var_N(X,Y,N,x)
%{
parameters:
    X:Matrix, each column represents a variable, and each row represents a sample
    Y:Fitting target value,Matrix with multiple rows and one column
    x:abscissa of point to be solved
%}
n=size(X,1);
assert (n==size(Y,1),'The first dimensions of matrices X and Y must be the same');
mutils=size(X,2);
assert (N>0 & N*mutils <n,'The number of rows of matrix X must be greater than the product of its number of columns and N')
vars=N*mutils+1;
f=zeros(vars,n);
f(1,:)=ones(1,n);
count=1;
for k=1:N
    for m=1:mutils
        count=count+1;
        f(count,:)=power(X(:,m),k);
    end
end
G=zeros(vars,vars);
for i=1:vars
    for j=1:i
        G(i,j)=f(i,:)*f(j,:)';
        G(j,i)=G(i,j);
    end
end
F=zeros(vars,1);
for k=1:vars
    F(k)=f(k,:)*Y;
end
C=G\F;
res=ones(size(x,1),1)*C(1);
formula=num2str(C(1));
count=1;
for k=1:N
    for m=1:mutils
        count=count+1;
        res=res+power(x(:,m),k)*C(count);
        formula=[formula,'+',num2str(C(count)),'x_{',num2str(m),'}^{',num2str(k),'}'];
    end
end
end

%多元一次函数拟合（基于超定方程A'AC=A'Y）
%返回基于数据X,Y对x的拟合结果及拟合函数公式
function [res,formula]=fit_mutil_var_1(X,Y,x)
%{
methord is to resolve XC=Y,and find the answer of C
parameters:
    X:Matrix, each column represents a variable, and each row represents a sample
    Y:Fitting target value,Matrix with multiple rows and one column
    x:abscissa of point to be solved
%}
X=[X ones(size(X,1),1)];
C=(X'*X)\X'*Y;
x=[x ones(size(x,1),1)];
res=x*C;
formula=num2str(C(1));
for k=2:size(C,1)
    formula=[formula,'+',num2str(C(k)),'x^{',num2str(k-1),'}'];
end
end
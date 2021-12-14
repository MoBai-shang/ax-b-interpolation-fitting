X=[2 5 4 7];
Y=[5 7 10 4];
Piecewise_Linear(X,Y,2.6)
Piecewise_Newton(X,Y,2.6,2)
Lagrange(X,Y,2.6)
Newton(X,Y,2.6)


function [s]=Lagrange(X,Y,x)
%{
parameters:
    X:abscissa of interpolation point
    Y:ordinate of interpolateion point
    x:abscissa of point to be solved
%}
assert (sum(size(X)-size(Y))==0,'The dimensions of matrices X and Y must be the same');
assert (rank(X)==1,'X and Y must be one-dimensional matrices');
n=max(size(X));
s=sum(arrayfun(@(j)Y(j)*(prod(arrayfun(@(i)(x-X(i))/(X(j)-X(i)),1:j-1))*prod(arrayfun(@(i)(x-X(i))/(X(j)-X(i)),j+1:n))),1:n));
end

function [res]=Newton(X,Y,x)
%{
Calculate the polynomial coefficients based on the difference quotient table
parameters:
    X:abscissa of interpolation point
    Y:ordinate of interpolateion point
    x:abscissa of point to be solved
%}
assert (sum(size(X)-size(Y))==0,'The dimensions of matrices X and Y must be the same');
assert (rank(X)==1,'X and Y must be one-dimensional matrices');
n=max(size(X));
res=Y(1);
w=x-X(1);
f_dif=zeros(n,n);%initialize the difference quotient table
f_dif(:,1)=Y;
for i=2:n
    for j=2:i
        f_dif(i,j)=(f_dif(i,j-1)-f_dif(i-1,j-1))/(X(i)-X(i-j+1));
    end
    res=res+w*f_dif(i,i);
    w=w*(x-X(i));
end
end

function [res]=Piecewise_Linear(X,Y,x)
%{
parameters:
    X:abscissa of interpolation point
    Y:ordinate of interpolateion point
    x:abscissa of point to be solved
%}
assert (sum(size(X)-size(Y))==0,'The dimensions of matrices X and Y must be the same');
assert (rank(X)==1,'X and Y must be one-dimensional matrices');
assert (x>=min(X) & x<=max(X),'The abscissa of the point to be evaluated must be within the abscissa of the interpolation point')
n=max(size(X));
sort_X=sort(X);
x1=sort_X(n-1);x2=sort_X(n);
for i=1:n
    if sort_X(i)>x
        x1=sort_X(i-1);
        x2=sort_X(i);
        break;
    end
end
res=(x-x2)/(x1-x2)*Y(X==x1)+(x-x1)/(x2-x1)*Y(X==x2);
end

function [res]=Piecewise_Cubic(X,Y,x)
%{
parameters:
    X:abscissa of interpolation point
    Y:ordinate of interpolateion point
    x:abscissa of point to be solved
    N:the number of interpolation points,defined as 4
%}
N=2;
assert (sum(size(X)-size(Y))==0,'The dimensions of matrices X and Y must be the same');
assert (rank(X)==1,'X and Y must be one-dimensional matrices');
assert (x>=min(X) & x<=max(X),'The abscissa of the point to be evaluated must be within the abscissa of the interpolation point')
n=max(size(X));
sort_X=sort(X);
points=sort_X(n-N+1:n);
for i=1:n
    if sort_X(i)>x
        lim_1=i-N+round(N/2);lim_2=i+round(N/2)-1;
        if lim_1<1
            points=sort_X(1:N);
        elseif lim_2<=n
            points=sort_X(lim_1:lim_2);
        end
        break;
    end
end
res=sum(arrayfun(@(j)Y(X==points(j))*(prod(arrayfun(@(i)(x-points(i))/(points(j)-points(i)),1:j-1))*prod(arrayfun(@(i)(x-points(i))/(points(j)-points(i)),j+1:N))),1:N));
end

function [res]=Piecewise_Lagrange(X,Y,x,N)
%{
parameters:
    X:abscissa of interpolation point
    Y:ordinate of interpolateion point
    x:abscissa of point to be solved
    N:Interpolation times
%}
N=N+1;%the number of interpolation points
assert (sum(size(X)-size(Y))==0,'The dimensions of matrices X and Y must be the same');
assert (rank(X)==1,'X and Y must be one-dimensional matrices');
assert (x>=min(X) & x<=max(X),'The abscissa of the point to be evaluated must be within the abscissa of the interpolation point')
n=max(size(X));
assert (N>0 & N <=n,'the interpolation times must between 0 and length of matrix X')
sort_X=sort(X);
points=sort_X(n-N+1:n);
for i=1:n
    if sort_X(i)>x
        lim_1=i-N+round(N/2);lim_2=i+round(N/2)-1;
        if lim_1<1
            points=sort_X(1:N);
        elseif lim_2<=n
            points=sort_X(lim_1:lim_2);
        end
        break;
    end
end
res=sum(arrayfun(@(j)Y(X==points(j))*(prod(arrayfun(@(i)(x-points(i))/(points(j)-points(i)),1:j-1))*prod(arrayfun(@(i)(x-points(i))/(points(j)-points(i)),j+1:N))),1:N));
end

function [res]=Piecewise_Newton(X,Y,x,N)
%{
parameters:
    X:abscissa of interpolation point
    Y:ordinate of interpolateion point
    x:abscissa of point to be solved
    N:Interpolation times
%}
N=N+1;%the number of interpolation points
assert (sum(size(X)-size(Y))==0,'The dimensions of matrices X and Y must be the same');
assert (rank(X)==1,'X and Y must be one-dimensional matrices');
assert (x>=min(X) & x<=max(X),'The abscissa of the point to be evaluated must be within the abscissa of the interpolation point')
n=max(size(X));
assert (N>0 & N <=n,'the interpolation times must between 0 and length of matrix X')
sort_X=sort(X);
points=sort_X(n-N+1:n);
for i=1:n
    if sort_X(i)>x
        lim_1=i-N+round(N/2);lim_2=i+round(N/2)-1;
        if lim_1<1
            points=sort_X(1:N);
        elseif lim_2<=n
            points=sort_X(lim_1:lim_2);
        end
        break;
    end
end
point_y=zeros(1,N);
for k=1:N
    point_y(k)=Y(X==points(k));
end
%methord first:Direct call Newton method
%res=Newton(points,point_y,x);

%methord second:
res=point_y(1);
w=x-points(1);
f_dif=zeros(N,N);%initialize the difference quotient table
f_dif(:,1)=point_y;
for i=2:N
    for j=2:i
        f_dif(i,j)=(f_dif(i,j-1)-f_dif(i-1,j-1))/(points(i)-points(i-j+1));
    end
    res=res+w*f_dif(i,i);
    w=w*(x-points(i));
end
end

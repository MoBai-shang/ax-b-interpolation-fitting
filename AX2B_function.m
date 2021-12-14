%Define variable coefficient matrix A with shape that is equal to (n,n) 
fprintf('Matrix A is:')
A=[12 1 2 4 4;
    1 9 -1 2 -3;
    2 -1 7 3 -5;
    3 2 3 12 -1;
    4 -3 -5 -1 15]

%Define target value matrix B whose shape is equal to (n,1)
fprintf('Matrix B is:')
B=[12 -27 14 -17 12]
rand('seed',1);

epoch=16
%BSOR methord
fprintf('BSOR methord:')
X=BSOR(A,B,epoch)'
%Jacobi method
fprintf('Jacobi methord:')
X=Jacobi(A,B,epoch)'
%GaussSeidel method
fprintf('GaussSeidel methord:')
X=GaussSeidel(A,B,epoch)'
%SteepestDescent method
fprintf('SteepestDescent methord:')
X=SteepestDescent(A,B,epoch)'
%ConjugateGradient method
fprintf('ConjugateGradient methord:')
X=ConjugateGradient(A,B,epoch)'

legend('BSOR','Jacobi','Gauss-Seidel','SteepestDescent','Conjugate Gradient');

%parameter:
%A :coefficient matrix
%B:Target value matrix
%epoch:Number of iterations
function [X]=BSOR(A,B,epoch)
    if size(B,2)>1
        B=B.';
    end
    verify(A,B);
    D=diag(diag(A));
    L=-tril(A,-1);
    U=-triu(A,1);
    W=-triu(A,0);
    Bj=diag(1./diag(A))*(L+U);
    w=2/(1+sqrt(1-power(vrho(Bj),2)));
    X=randn(size(A,1),1);%Initialize x
    errors=zeros(1,epoch);%Storage loss
    for k=1:epoch
        %Matrix iteration formula of BSOR method
        X=inv(D-w*L)* (w*B+(w*W+D)*X);% DX'=DX+w(B+LX'+WX)
        errors(k)=sqrt(mean(power(A*X-B,2)));%Euclidean distance loss
    end
    hold on
    %legend('BSOR');
    plot(errors);
    title('Loss variation in iterative process')
    xlabel('epoch');
    ylabel('loss');
end

function [X]=Jacobi(A,B,epoch)
    if size(B,2)>1
        B=B.';
    end
    verify(A,B);
    D=diag(diag(A));
    L=-tril(A,-1);
    U=-triu(A,1);
    Bj=diag(1./diag(A))*(L+U);
    g=inv(D)*B;
    X=randn(size(A,1),1);%Initialize x
    errors=zeros(1,epoch);%Storage loss
    for k=1:epoch
        %Matrix iteration formula of Jacobi method
        X=Bj*X+g;% X'=BjX+g
        errors(k)=sqrt(mean(power(A*X-B,2)));%Euclidean distance loss
    end
    hold on
    %legend('Jacobi');
    plot(errors);
end

function [X]=GaussSeidel(A,B,epoch)
    if size(B,2)>1
        B=B.';
    end
    verify(A,B);
    D=diag(diag(A));
    L=-tril(A,-1);
    U=-triu(A,1);
    Bj=diag(1./diag(A))*(L+U);
    Bg=inv(D-L)*U;
    g=inv(D-L)*B;
    X=randn(size(A,1),1);%Initialize x
    errors=zeros(1,epoch);%Storage loss
    for k=1:epoch
        %Matrix iteration formula of Jacobi method
        X=Bg*X+g;% X'=BjX+g
        errors(k)=sqrt(mean(power(A*X-B,2)));%Euclidean distance loss
    end
    hold on
    %legend('Gauss-Seidel');
    plot(errors);
end

function [X]=SteepestDescent(A,B,epoch)
    if size(B,2)>1
        B=B.';
    end
    verify(A,B);
    X=rand(size(A,1),1);%Initialize x
    errors=zeros(1,epoch);%Storage loss
    for k=1:epoch
        r=B-A*X;
        alpa=dot(r,r)/dot(A*r,r);
        X=X+alpa*r;
        errors(k)=sqrt(mean(power(A*X-B,2)));%Euclidean distance loss
    end
    hold on
    plot(errors)
end

function [X]=ConjugateGradient(A,B,epoch)
    if size(B,2)>1
        B=B.';
    end
    verify(A,B);
    X=rand(size(A,1),1);%Initialize x
    r0=B-A*X;
    p=r0;
    errors=zeros(1,epoch);%Storage loss
    for k=1:epoch
        lr=dot(r0,r0)/dot(A*p,p);
        %lr=r0'*r0/(p'*A*p);
        X=X+lr*p;
        r1=r0-lr*A*p;
        beta=dot(r1,r1)/dot(r0,r0);
        %beta=r1'*r1/(r0'*r0);
        p=r1+beta*p;
        r0=r1;
        errors(k)=sqrt(mean(power(A*X-B,2)));%Euclidean distance loss
    end
    hold on
    %legend('Conjugate Gradient')
    plot(errors)
end

function []=verify(A,B)
    assert(size(A,1)==size(A,2),['Matrix A must be square with shape (n,n)'])
    assert(1==size(B,2),['The shape of matrix B must be (n,1) or (1,n)'])
    assert(size(A,1)==size(B,1),['Matrix A and B must have the same first or second dimension'])
end

%{
/**********************************************************/
/*************** VIVEK PAL ********************************/
/*************** 17MA20048 ********************************/
/**********************************************************/

The problem statement is to solve the boundary value problem using finite
difference method (fictious point at the boundary conditions) and to print
nodes, value of unknowns as well as absolute error.

The Equation is... 
u"(x) = -(pi^2)cos(pi*x); u(0)=1, u'(0.5)=-pi;
Take n=30. (exact solution: cos(pi*x))

And you have to find Y(x+k*h), k = 0,1,...,n; h=(b-a)/n;
%}

a=0;        %% Initial node/point
b=0.5;      %% Final node/point
N=30;
dx=(b-a)/N;
Xk=zeros(1,N);      %% Node Values
for i=1:N
    Xk(1,i)=a+i*dx;
end
Ck=zeros(1,N);      %% cos(x) at different nodes
for i=1:N
    Ck(1,i)=cos(pi*Xk(1,i));
end

mat=zeros(N,N);     %%% It is being used to form Thomos algo matrix.
                    %%% In System of Equations AX=B, here mat is A.
for i=1:N
    if i==1
        mat(i,i)=-2;
        mat(i,i+1)=1;
    elseif i==N
        mat(i,i-1)=2;
        mat(i,i)=-2;
    else
        mat(i,i-1)=1;
        mat(i,i)=-2;
        mat(i,i+1)=1;
    end
end

%%% There is somewhere multiplication or division of 3600
%%% Due to some h^2, because h=1/60, it is according to this only
%%% particular question.
D=zeros(1,N);     %%% In System of Equations AX=B, here B is D.
for i=1:N
    if i==1
        D(1,i)=(-pi*pi*Ck(1,i))/3600-1;
    elseif i==N
        D(1,i)=(-pi*pi*cos(pi*0.5))/3600+12*pi/360;
    else
        D(1,i)=(-pi*pi*Ck(1,i))/3600;
    end
end
mat=3600*mat;
D=3600*D;
for i=1:N
    for j=1:N
        fprintf('%d ', mat(i,j));
    end
    fprintf('%d\n', D(1,i));
end

%%% Operating Thomas algorithm
%%% One may do inverse of mat
%%% as in AX=B, for solution of X; X=A^(-1)B
for i=1:N
   if i==1
       mat(i,i+1)=mat(i,i+1)/mat(i,i);
       D(1,i)=D(1,i)/mat(i,i);
       mat(i,i)=1;
   elseif i==N
       mat(i,i)=mat(i,i)-mat(i-1,i)*mat(i,i-1);
       D(1,i)=D(1,i)-D(1,i-1)*mat(i,i-1);
       D(1,i)=D(1,i)/mat(i,i);
       mat(i,i)=1;
       mat(i,i-1)=0;
   else
       mat(i,i)=mat(i,i)-mat(i-1,i)*mat(i,i-1);
       D(1,i)=D(1,i)-D(1,i-1)*mat(i,i-1);
       mat(i,i+1)=mat(i,i+1)/mat(i,i);
       D(1,i)=D(1,i)/mat(i,i);
       mat(i,i)=1;
       mat(i,i-1)=0;
   end
end

for i=1:N
    for j=1:N
        fprintf('%d ', mat(i,j));
    end
    fprintf('%d\n', D(1,i));
end
anst=zeros(1,N);
%%% Back-substitution
for i=N:-1:1
    if i==N
        anst(1,i)=D(1,i);
    else
        anst(1,i)=D(1,i)-mat(i,i+1)*anst(1,i+1);
    end
end

fprintf('Point Xk   ---  Expected Value -- Exact Value -- Abs Error\n');
for i=0:N
    if i==0
        fprintf('%d            **    %d         **       %d      **      %d\n', a+i*dx, 1, cos(pi*a), abs(cos(pi*a)-1));
    else
        fprintf('%d ** %d ** %d ** %d\n', a+i*dx, anst(1,i), cos(pi*(a+i*dx)), abs(cos(pi*(a+i*dx))-anst(1,i)));
    end
end
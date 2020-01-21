%{
/**********************************************************/
/*************** VIVEK PAL ********************************/
/*************** 17MA20048 ********************************/
/**********************************************************/

The problem statement is to solve the boundary value problem using finite
difference method (fictious point at the boundary conditions) and to print
nodes, value of unknowns as well as absolute error.

The Equation is... 
u"(x) + u(x) = sin(3x); u(0)+u'(0)=-1, u'(pi/2)=1;
Take n=30. (exact solution: -cos(x)+(3/8)sin(x)-(1/8)sin(3x) )

And you have to find Y(x+k*h), k = 0,1,...,n; h=(b-a)/n;
%}

a=0;        %% Initial node/point
b=pi/2;     %% Final node/point
N=30;
dx=(b-a)/N;
Xk=zeros(1,N+1);      %% Node Values
for i=1:N+1
    Xk(1,i)=a+(i-1)*dx;
end

mat=zeros(N+1,N+1);     %%% It is being used to form Thomos algo matrix.
                        %%% In System of Equations AX=B, here mat is A.
for i=1:N+1
    if i==1
        mat(i,i)=dx*dx+2*dx-2;
        mat(i,i+1)=2;
    elseif i==N+1
        mat(i,i-1)=2;
        mat(i,i)=dx*dx-2;
    else
        mat(i,i-1)=1;
        mat(i,i)=dx*dx-2;
        mat(i,i+1)=1;
    end
end

D=zeros(1,N+1);     %%% In System of Equations AX=B, here B is D.
for i=1:N+1
    if i==1
        D(1,i)=dx*dx*sin(3*Xk(1,i))-2*dx;
    elseif i==N+1
        D(1,i)=dx*dx*sin(3*Xk(1,i))-2*dx;
    else
        D(1,i)=dx*dx*sin(3*Xk(1,i));
    end
end
for i=1:N+1
    for j=1:N+1
        fprintf('%d ', mat(i,j));
    end
    fprintf('%d\n', D(1,i));
end

%%% Operating Thomas algorithm
%%% One may do inverse of mat
%%% as in AX=B, for solution of X; X=A^(-1)B
for i=1:N+1
   if i==1
       mat(i,i+1)=mat(i,i+1)/mat(i,i);
       D(1,i)=D(1,i)/mat(i,i);
       mat(i,i)=1;
   elseif i==N+1
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

for i=1:N+1
    for j=1:N+1
        fprintf('%d ', mat(i,j));
    end
    fprintf('%d\n', D(1,i));
end
anst=zeros(1,N+1);
%%% Back-substitution
for i=N+1:-1:1
    if i==N+1
        anst(1,i)=D(1,i);
    else
        anst(1,i)=D(1,i)-mat(i,i+1)*anst(1,i+1);
    end
end
%%% Exk is to store exact solution at nodes
Exk=zeros(1,N+1);
for i=1:N+1
    Exk(1,i)=-cos(Xk(1,i))+(3/8)*sin(Xk(1,i))-(1/8)*sin(3*Xk(1,i));
end
E=-cos(0)+(3/8)*sin(0)-(1/8)*sin(3*0); %%% Exact solution at initial node
fprintf('Point Xk   ---  Expected Value --- Exact Value --- Abs Error\n');
for i=0:N
    if i==0
        fprintf('%d     **        %d **     %d    **      %d\n', a+i*dx, anst(1,i+1), E, abs(E-anst(1,i+1)));
    else
        fprintf('%d ** %d ** %d ** %d\n', a+i*dx, anst(1,i+1), Exk(1,i+1), abs(Exk(1,i+1)-anst(1,i+1)));
    end
end
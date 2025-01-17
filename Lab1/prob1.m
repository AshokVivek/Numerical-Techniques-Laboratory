%{
/**********************************************************/
/*************** VIVEK PAL ********************************/
/*************** 17MA20048 ********************************/
/**********************************************************/

The Problem statement is to solve Linear Second Order BVP
by Finite Difference Method using Thomos Algorithm.

The Equation would be like this... 
Y" + A(x)Y' + B(x)Y = C(x), a<x<b ; Y(a) = ya, Y(b) = yb.

And you have to find Y(x+k*h), k = 0,1,...,n; h=(b-a)/n;
%}

%{ Taking the parameters of the equation as input...}
a=input('Please enter initial boundary values of x, a :- ');
b=input('Please enter final boundary values of x, b :- ');
ya=input("Please enter value y(a) :- ");
yb=input("Please enter value y(b) :- ");
dx=input("Please input value of delta x :- ");
Na=input("Please enter the degree of polynomial A(x) :- ");
fprintf('Please enter the coefficient of A(x) in increasing order of degree :- \n');
A=zeros(1,Na+1);
for i=1:Na+1
    fprintf("Coeff. of A(x^%d) :- ", i-1);
    A(1,i)=input('');
end
Nb=input("Please enter the degree of B(x) :- ");
fprintf('Please enter the coefficient of B(x) in increasing order of degree :- \n');
B=zeros(1,Nb+1);
for i=1:Nb+1
    fprintf("Coeff. of B(x^%d) :- ", i-1);
    B(1,i)=input('');
end
Nc=input("Please enter the degree of C(x) :- ");
fprintf('Please enter the coefficient of C(x) in increasing order of degree :- \n');
C=zeros(1,Nc+1);
for i=1:Nc+1
    fprintf("Coeff. of C(x^%d) :- ", i-1);
    C(1,i)=input('');
end

N=(b-a)/dx;         %%% Calculating total no. of steps.
Xk=zeros(1,N+1);    %%% Calculating value of 'x' at each step
for i=1:N+1
    Xk(1,i)=a+(i-1)*dx;
end

Ak=zeros(1,N+1);    %%% Calcualting value of polynomial A(x) at each step
for i=1:N+1
    x=1;
    for j=1:Na+1
        Ak(1,i)=Ak(1,i)+x*A(1,j);
        x=x*Xk(1,i);
    end
end
Bk=zeros(1,N+1);    %%% Calcualting value of polynomial B(x) at each step
for i=1:N+1
    x=1;
    for j=1:Nb+1
        Bk(1,i)=Bk(1,i)+x*B(1,j);
        x=x*Xk(1,i);
    end
end
Ck=zeros(1,N+1);    %%% Calcualting value of polynomial C(x) at each step
for i=1:N+1
    x=1;
    for j=1:Nc+1
        Ck(1,i)=Ck(1,i)+x*C(1,j);
        x=x*Xk(1,i);
    end
end

disp(Xk);
disp(Ak);
disp(Bk);
disp(Ck);

%{
In Finite Difference Method:
The general form of equation comes after discretization of differential eqn as
Y(k-1)[1/(h*h) - Ak/(2*h)] + Y(k)[Bk-2/(h*h)] + Y(k+1)[1/(h*h) + Ak/(2*h)] = 0.
Here Y(k) represent value of Y at Xk, and Ak represent value of polynomial
A at value Xk.
%}
aak=zeros(1,N-1);
bbk=zeros(1,N-1);
cck=zeros(1,N-1);

for i=1:N-1
    aak(1,i)=(1/(dx*dx)-Ak(1,i+1)/dx);
    bbk(1,i)=Bk(1,i+1)-2/(dx*dx);
    cck(1,i)=(1/(dx*dx)+Ak(1,i+1)/dx);
end

fprintf("Now aak, bbk, cck, time :- \n");
disp(aak);
disp(bbk);
disp(cck);

mat=zeros(N-1,N-1);     %%% It is being used to form Thomos algo matrix.
                        %%% In System of Equations AX=B, here mat is A.
for i=1:N-1
    if i==1
        mat(i,i)=bbk(1,i);
        mat(i,i+1)=cck(1,i);
    elseif i==N-1
        mat(i,i-1)=aak(1,i);
        mat(i,i)=bbk(1,i);
    else
        mat(i,i-1)=aak(1,i);
        mat(i,i)=bbk(1,i);
        mat(i,i+1)=cck(1,i);
    end
end
D=zeros(1,N-1);     %%% In System of Equations AX=B, here B is D.
for i=1:N-1
    if i==1
        D(1,i)=Ck(1,i+1)-aak(1,i)*ya;
    elseif i==N-1
        D(1,i)=Ck(1,i+1)-cck(1,i)*yb;
    else
        D(1,i)=Ck(1,i+1);
    end
end

fprintf("Now intermidiate mat and D time :-\n");
mat1=mat;
disp(mat);
D1=D;
disp(D);
for i=1:N-1
   if i==1
       mat(i,i+1)=mat(i,i+1)/mat(i,i);
       D(1,i)=D(1,i)/mat(i,i);
       mat(i,i)=1;
   elseif i==N-1
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

fprintf("Now final mat and D time :- \n");
disp(mat);
disp(D);

anst=zeros(1,N-1);
for i=N-1:-1:1
    if i==N-1
        anst(1,i)=D(1,i);
    else
        anst(1,i)=D(1,i)-mat(i,i+1)*anst(1,i+1);
    end
end
for i=0:N
    if i==0
        fprintf("%d -- %d\n", a+i*dx, ya);
    elseif i==N
        fprintf("%d -- %d\n", a+i*dx, yb);
    else
        fprintf("%d -- %d\n", a+i*dx, anst(1,i));
    end
end
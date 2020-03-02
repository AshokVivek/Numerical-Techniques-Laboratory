%*********************************************
%******* VIVEK PAL ***************************
%******* 17MA20048 ***************************
%*********************************************

%{
Solve the Boundary Value Problem using finite difference Method.
Print nodes, numerical value of unknown y, exact values, as well as
absolute error.
    Problem: Y(IV)+4Y=1; Y(-1)=Y(1)=0=Y"(-1)=Y"(1), h=0.05.
    Exact Solution Y(x)=(1/4)*[1-2*{sin1*sinh1*sinx*sinhx + cos1*cosh1*cosx*coshx}/(cos2+cosh2)].
%}

h=0.05;
N=2/h;
%N=2;

%A=zeros(N-1);
A=zeros(2,2,N-1,N-1);
for i=1:N-1
   if i==1
       A(:,:,i,i)=[-2/(h*h),-1;4,-2/(h*h)];
       A(:,:,i,i+1)=[1/(h*h),0;0,1/(h*h)];
   elseif i==N-1
       A(:,:,i,i-1)=[1/(h*h),0;0,1/(h*h)];
       A(:,:,i,i)=[-2/(h*h),-1;4,-2/(h*h)];
   else
       A(:,:,i,i-1)=[1/(h*h),0;0,1/(h*h)];
       A(:,:,i,i)=[-2/(h*h),-1;4,-2/(h*h)];
       A(:,:,i,i+1)=[1/(h*h),0;0,1/(h*h)];
   end
end
%disp(A);
D=zeros(2,1,N-1);
for i=1:N-1
   D(:,:,i)=[0;1];
end
%disp(D);

for i=1:N-1
   sun=A(:,:,i,i);
   %disp(sun);
   if i==1
        D(:,:,i)=A(:,:,i,i)\D(:,:,i);
        A(:,:,i,i+1)=A(:,:,i,i)\A(:,:,i,i+1);
        A(:,:,i,i)=eye(2);
   elseif i==N-1
        A(:,:,i,i)=A(:,:,i,i)-A(:,:,i,i-1)*A(:,:,i-1,i); % Bi'=Bi-AiC'(i-1)
        D(:,:,i)=A(:,:,i,i)\(D(:,:,i)-A(:,:,i,i-1)*D(:,:,i-1));
        A(:,:,i,i)=eye(2);
        A(:,:,i,i-1)=zeros(2);
   else
        A(:,:,i,i)=A(:,:,i,i)-A(:,:,i,i-1)*A(:,:,i-1,i); % Bi'=Bi-AiC'(i-1)
        A(:,:,i,i+1)=A(:,:,i,i)\A(:,:,i,i+1);
        D(:,:,i)=A(:,:,i,i)\(D(:,:,i)-A(:,:,i,i-1)*D(:,:,i-1));
        A(:,:,i,i)=eye(2);
        A(:,:,i,i-1)=zeros(2);
   end
end
%disp(A);
%disp(D);
ans1=zeros(2,1,N-1);
for i=1:N-1
   ans1(:,:,i)=zeros(2,1);
end

for i=N-1:-1:1
   if i==N-1
       ans1(:,:,i)=D(:,:,i);
   else
       ans1(:,:,i)=D(:,:,i)-A(:,:,i,i+1)*ans1(:,:,i+1);
   end
end
%disp(ans1);

a=-1;
b=1;
Xk=zeros(1,N+1);
Yk=zeros(1,N+1);    %%% Calculating value of 'x' at each step
for i=1:N+1
    Xk(1,i)=a+(i-1)*h;
end
for i=1:N+1
    Yk(1,i)=(1-2*(sin(1)*sinh(1)*sin(Xk(1,i))*sinh(Xk(1,i))+cos(1)*cosh(1)*cos(Xk(1,i))*cosh(Xk(1,i)))/(cos(2)+cosh(2)))/4;
end

z=0;
fprintf("Node Point x       Y(x)       Exact Soln   Absolute Error\n");
for i=2:N
    fprintf("%d  %d  %d  %d\n", Xk(1,i), ans1(1,1,i-1), Yk(1,i), abs(ans1(1,1,i-1)-Yk(1,i)));
end
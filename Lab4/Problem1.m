%*********************************************
%******* VIVEK PAL ***************************
%******* 17MA20048 ***************************
%*********************************************

%{
Solve the Boundary Value Problem using finite difference Method.
Print nodes, numerical value of unknown u, exact values, as well as
absolute error.
    Problem: u"(x)=-(1/8)[32+2x^3-uu']; u(1)=17, u(3)=43/3, h=0.05.
    Exact Solution u(x)=x^2+(16/x).
%}

% After discretizing the given eqation by finite difference method
% We get a function G(xi,u(i-1), u(i), u(i+1))=0  % here u(i) is ui.
% We have to get the equation in terms for delta(yi), we denote it as dyi

function Problem1()

    %-----SOME INITIAL VALUES-----
    h=0.05;
    a=1;
    b=3;
    ua=17;
    ub=43/3;

    %-----Total Nodes Calculation-----
    N=(b-a)/h;

    %-----Nodes Value Storage-----
    Xk=zeros(1,N);
    for i=1:N
       Xk(1,i)=a+i*h; 
    end
    %print(Xk,N,0);
    
    %-----Initial guess for ui-----
    %   Initial guess is done by: u(x)=(x-a)(ub/(b-a))+(b-x)(ua/(b-a)).
    Uk=zeros(1,N);
    for i=1:N
       Uk(1,i)=(Xk(1,i)-a)*(ub/(b-a))+(b-Xk(1,i))*(ua/(b-a));
    end
    %print(Uk, N, 0);
    
    %-----Now start iterative Optimization-----
    A=zeros(N-1,N-1);
    D=zeros(N-1,1);
    maxm=1;
    while(maxm>(10.^(-10)))
        %{
        fprintf('--%d----------------------\n', cnt);
        print(Uk,N-1,0);
        fprintf('--------------------------\n');
        cnt=cnt-1;
        %}
        for i=1:N-1
           if i==1
               A(i,i)=-ua/(16*h)+Uk(1,i+1)/(16*h)-2/(h*h);
               A(i,i+1)=1/(h*h)+Uk(1,i)/(16*h);
           elseif i==N-1
               A(i,i-1)=1/(h*h)-Uk(1,i)/(16*h);
               A(i,i)=-Uk(1,i-1)/(16*h)+ub/(16*h)-2/(h*h);
           else
               A(i,i-1)=1/(h*h)-Uk(1,i)/(16*h);
               A(i,i)=-Uk(1,i-1)/(16*h)+Uk(1,i+1)/(16*h)-2/(h*h);
               A(i,i+1)=1/(h*h)+Uk(1,i)/(16*h);
           end
        end
        %print2D(A,N-1,N-1);
        
        
        for i=1:N-1
            if i==1
                D(i,1)=ua*(-1/(h*h)+Uk(1,i)/(16*h))+2*Uk(1,i)/(h*h)-Uk(1,i+1)*(1/(h*h)+Uk(1,i)/(16*h))+4+(Xk(1,i)*Xk(1,i)*Xk(1,i))/4;
            elseif i==N-1
                D(i,1)=Uk(1,i-1)*(-1/(h*h)+Uk(1,i)/(16*h))+2*Uk(1,i)/(h*h)-ub*(1/(h*h)+Uk(1,i)/(16*h))+4+(Xk(1,i)*Xk(1,i)*Xk(1,i))/4;
            else
                D(i,1)=Uk(1,i-1)*(-1/(h*h)+Uk(1,i)/(16*h))+2*Uk(1,i)/(h*h)-Uk(1,i+1)*(1/(h*h)+Uk(1,i)/(16*h))+4+(Xk(1,i)*Xk(1,i)*Xk(1,i))/4;
            end
        end
        %print(D,N-1,1);
        
        dUk=A\D;
        %print(dUk,N-1,1)
        maxm=maxt(dUk,N-1,1);
        %fprintf('---**-- %d --**---\n', maxm);
        
        for i=1:N-1
           Uk(1,i)=Uk(1,i)+dUk(i,1); 
        end
    end
    %print(Uk,N,0);
    
    eUk=zeros(1,N);
    for i=1:N
       eUk(1,i)=Xk(1,i)*Xk(1,i)+16/Xk(1,i);
    end
    
    fprintf('Node Xk        Appx. Value    Exact Value    Absolute error\n');
    for i=1:N-1
       fprintf('%d   %d   %d   %d\n', Xk(1,i), Uk(1,i), eUk(1,i), abs(Uk(1,i)-eUk(1,i))); 
    end
end

function print(A,N,k)
    if k==0
        for i=1:N
           fprintf('%d ',A(1,i)); 
        end
    else
        for i=1:N
           fprintf('%d ',A(i,1)); 
        end
    end
    fprintf('\n');
end

function print2D(A, N, M)
    for i=1:N
       for j=1:M
          fprintf('%d ', A(i,j)); 
       end
       fprintf('\n');
    end
end

function maxm=maxt(A,N,k)
    maxm=0;
    if k==0
        for i=1:N
           if A(1,i)<0
               A(1,i)=-A(1,i);
           end
        end
        for i=1:N
           if maxm<A(1,i)
               maxm=A(1,i);
           end
        end
    else
        for i=1:N
           if A(i,1)<0
               A(i,1)=-A(i,1);
           end
        end
        for i=1:N
           if maxm<A(i,1)
               maxm=A(i,1);
           end
        end
    end
end
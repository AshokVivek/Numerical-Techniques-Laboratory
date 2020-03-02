%*********************************************
%******* VIVEK PAL ***************************
%******* 17MA20048 ***************************
%*********************************************

%{
Solve the Boundary Value Problem using finite difference Method.
Print nodes, numerical value of unknown u, exact values, as well as
absolute error.
    Problem: u"(x)=-u'^2-u+ln(x) ; u(1)=0, u(2)=ln(2); h=0.025.
    Exact Solution u(x)=ln(x).
%}

% After discretizing the given eqation by finite difference method
% We get a function G(xi,u(i-1), u(i), u(i+1))=0  % here u(i) is ui.
% We have to get the equation in terms for delta(yi), we denote it as dyi

function Problem2()

    %-----SOME INITIAL VALUES-----
    h=0.025;
    a=1;
    b=2;
    ua=0;
    ub=log(2);

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
    %for k=1:4
        %print(Uk,N-1,0);
        %fprintf('--------------------------\n');
        
        for i=1:N-1
           if i==1
               A(i,i)=1-2/(h*h);
               A(i,i+1)=Uk(1,i+1)/(2*h*h)+1/(h*h)-ua/(2*h*h);
           elseif i==N-1
               A(i,i-1)=Uk(1,i-1)/(2*h*h)+1/(h*h)-ub/(2*h*h);
               A(i,i)=1-2/(h*h);
           else
               A(i,i-1)=Uk(1,i-1)/(2*h*h)+1/(h*h)-Uk(1,i+1)/(2*h*h);
               A(i,i)=1-2/(h*h);
               A(i,i+1)=Uk(1,i+1)/(2*h*h)+1/(h*h)-Uk(1,i-1)/(2*h*h);
           end
        end
        %print2D(A,N-1,N-1);
        
        
        for i=1:N-1
            if i==1
                D(i,1)=-(Uk(1,i+1)*Uk(1,i+1)/(4*h*h)+Uk(1,i+1)/(h*h)+ua*ua/(4*h*h)+ua/(h*h)-Uk(1,i+1)*(ua/(2*h*h))+Uk(1,i)*(1-2/(h*h))-log(Xk(1,i)));
            elseif i==N-1
                D(i,1)=-(ub*ub/(4*h*h)+ub/(h*h)+Uk(1,i-1)*Uk(1,i-1)/(4*h*h)+Uk(1,i-1)/(h*h)-ub*(Uk(1,i-1)/(2*h*h))+Uk(1,i)*(1-2/(h*h))-log(Xk(1,i)));
            else
                D(i,1)=-(Uk(1,i+1)*Uk(1,i+1)/(4*h*h)+Uk(1,i+1)/(h*h)+Uk(1,i-1)*Uk(1,i-1)/(4*h*h)+Uk(1,i-1)/(h*h)-Uk(1,i+1)*(Uk(1,i-1)/(2*h*h))+Uk(1,i)*(1-2/(h*h))-log(Xk(1,i)));
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
       eUk(1,i)=log(Xk(1,i));
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
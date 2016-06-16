%THIS IS THE AUXILIARY FUNCTION OF THE MAIN CODE "periodic_benchmark"



function[Z,extra,n1,B,L,A,h] = SN_per_bench_solver(T,m1,m2,n,N,Es1,Es2,Et1,Et2,yo,y_,Q1,Q2,u,wt,a)

%Periodic Medium........
interval=T/n;
m12=zeros(2,1);
m12(1)=m1;m12(2)=m2;
mm=0;
i=0;
s=0;
if a>1
    if a<m2/interval+2
		i=i+1;
		x1=(a-1)*interval;
		s=s+x1;
		x(i,1)=s;
		x(i,2)=2;
	else
    	i=i+1;
		x1=(a-m2/interval-1)*interval;
		s=s+x1;
		x(i,1)=s;
		x(i,2)=mod(mm,2)+1;
       mm=mm+1;
    end
end

while s<T
    i=i+1;
    x1=m12(mod(mm,2)+1);
    s=s+x1;
    if s<=T
        x(i,1)=s;
    else
        x(i,1)=T;
    end
    x(i,2)=mod(mm,2)+1;
    mm=mm+1;
end

H=T/n;n1=1;i=1;j=1;
extra=zeros(n,1);
t1=i*H; %%%%%%WHATS THE POINT OF THIS???
dlmwrite('xsn_m.txt',x,'delimiter',',','precision',18);
display('saved xsn');
L(1,1)=0;L(1,2)=x(j,2);
if t1==x(j,1)
    display('hi');
    extra(i)=1;
    h(n1)=x(j,1)/2;
    L(n1+1,1)=h(n1);
    L(n1+1,2)=x(j,2);
    h(n1+1)=h(n1);
    L(n1+2,1)=x(j,1);
    L(n1+2,2)=3;
    n1=n1+2;
     i=i+1;
     t1=i*H;
else
    while t1<=x(j,1)
        h(n1)=H;
        display(H);
        L(n1+1,1)=i*H;
        if L(n1+1,1)==x(j,1)                 % checking is point is interface point?
            L(n1+1,2)=3;
        else
            L(n1+1,2)=x(j,2);
        end
        n1=n1+1;
        i=i+1;
        t1=i*H;
    end
end
j=2;
if x(1,1)~=T
    while i<=n
        while t1<=x(j,1)
                L(n1+1,1)=i*H;
                if L(n1+1,1)==x(j,1)          %checking is point is interface point
                    L(n1+1,2)=3;
                else
                    L(n1+1,2)=x(j,2);
                end
                h(n1)=L(n1+1,1)-L(n1,1);
                n1=n1+1;
                i=i+1;
                t1=i*H;
        end
        j=j+1;
        if L(n1,1)==T && L(n1,2)==3 && L(n1-1,2)==3
            display('boop');
            display(n1);
            i=i-1; 
            extra(i)=extra(i)+1;
            L(n1,1)=(x(j-1,1)+x(j-2,1))/2;
            L(n1,2)=x(j-1,2);
            h(n1-1)=L(n1,1)-L(n1-1,1);
            n1=n1+1;
            display(n1);
            L(n1,1)=x(j-1,1);
            L(n1,2)=3;
            
            h(n1-1)=L(n1,1)-L(n1-1,1);
            i=i+1;
        end       
    end
end
n1=n1-1;
A=zeros(N*n1,N*n1); B=zeros(N*n1,1);
Et12=zeros(2,1); Es12=zeros(2,1); Q12=zeros(2,1);
Et12(1)=Et1;Et12(2)=Et2;
Es12(1)=Es1;Es12(2)=Es2;
Q12(1)=Q1;Q12(2)=Q2;


% Diagonal Block of matrix up to N/2.................................
for t=1:N/2
    s=(t-1)*n1;
    A(s+1,s+1)=-u(t)*(1/h(1))+Et12(L(1,2))-Es12(L(1,2))*wt(t)/2;
    A(s+1,s+2)=u(t)*(1/h(1));
    B(s+1)=Q12(L(1,2))/2;
    for i=2:n1-1
        if L(i,2)==3
            if i==n1-1
                A(s+i,s+i)=-u(t)*(1/h(i)+1/(h(i)+h(i+1)))+Et12(L(i+1,2))-Es12(L(i+1,2))*wt(t)/2;
                A(s+i,s+i+1)=u(t)*(1/h(i)+1/h(i+1));
                B(s+i)=u(t)*y_*(h(i)/(h(i+1)*(h(i)+h(i+1))))+Q12(L(i+1,2))/2;
            else
                A(s+i,s+i)=-u(t)*(1/h(i)+1/(h(i)+h(i+1)))+Et12(L(i+1,2))-Es12(L(i+1,2))*wt(t)/2;
                A(s+i,s+i+1)=u(t)*(1/h(i)+1/h(i+1));
                A(s+i,s+i+2)=-u(t)*(h(i)/(h(i+1)*(h(i)+h(i+1))));
                B(s+i)=Q12(L(i+1,2))/2;
        
            end
        else
            A(s+i,s+i-1)=-u(t)*h(i)/(h(i-1)*(h(i-1)+h(i)));
            A(s+i,s+i)=u(t)*(h(i)-h(i-1))/(h(i)*h(i-1))+Et12(L(i,2))-Es12(L(i,2))*wt(t)/2;
            A(s+i,s+i+1)=u(t)*h(i-1)/(h(i)*(h(i-1)+h(i)));
            	B(s+i)=Q12(L(i,2))/2;
        end
    end
    A(s+n1,s+n1-1)=-u(t)*h(n1)/(h(n1-1)*(h(n1-1)+h(n1)));
    A(s+n1,s+n1)=u(t)*(h(n1)-h(n1-1))/(h(n1)*h(n1-1))+Et12(L(n1,2))-Es12(L(n1,2))*wt(t)/2;
    B(s+n1)=-u(t)*y_*h(n1-1)/(h(n1)*(h(n1-1)+h(n1)))+Q12(L(n1,2))/2;
    
    % Blocks from N/2 to N........................................
    a=0;
    for p=1:N/2
        S=(N/2+p-1)*n1;
        for i=2:n1
            if L(i,2)==3
                A(s+i,S+i-1)=-Es12(L(i+1,2))*wt(N/2+p)/2;
            else
                A(s+i,S+i-1)=-Es12(L(i,2))*wt(N/2+p)/2;
            end
        end
        a=a+(Es12(L(1,2))*wt(N/2+p)*yo/2);
    end
    B(s+1)=B(s+1)+a;
end
    
% Diagonal Block of matrix from N/2+1 to N.........................
for t=N/2+1:N
    s=(t-1)*n1;
    A(s+1,s+1)=u(t)*(h(2)-h(1))/(h(2)*h(1))+Et12(L(1,2))-Es12(L(1,2))*wt(t)/2;
    A(s+1,s+1+1) = u(t)*h(1)/(h(2)*(h(1)+h(2)));
    B(s+1)=u(t)*yo*h(2)/(h(1)*(h(1)+h(2)))+Q12(L(1,2))/2;
    for i=2:n1-1
        if L(i+1,2)==3
            if i==2
                A(s+i,s+i)=u(t)*(1/h(i)+1/(h(i)+h(i-1)))+Et12(L(i,2))-Es12(L(i,2))*wt(t)/2;
                A(s+i,s+i-1)=-u(t)*(1/h(i)+1/h(i-1));
                B(s+i)=-u(t)*yo*(h(i)/(h(i-1)*(h(i)+h(i-1))))+Q12(L(i,2))/2;
            else
                A(s+i,s+i)=u(t)*(1/h(i)+1/(h(i)+h(i-1)))+Et12(L(i,2))-Es12(L(i,2))*wt(t)/2;
                A(s+i,s+i-1)=-u(t)*(1/h(i)+1/h(i-1));
                A(s+i,s+i-2)=u(t)*(h(i)/(h(i-1)*(h(i)+h(i-1))));
                B(s+i)=Q12(L(i,2))/2;
                if s+i ==105
                    fprintf('105');
                end
	        end
        else
            A(s+i,s+i-1)=-u(t)*h(i+1)/(h(i)*(h(i)+h(i+1)));
            A(s+i,s+i)=u(t)*(h(i+1)-h(i))/(h(i+1)*h(i))+Et12(L(i+1,2))-Es12(L(i+1,2))*wt(t)/2;
            A(s+i,s+i+1)=u(t)*h(i)/(h(i+1)*(h(i)+h(i+1)));
            B(s+i)=Q12(L(i+1,2))/2;
	    end
    end
    A(s+n1,s+n1)=u(t)*(1/h(n1))+Et12(L(n1,2))-Es12(L(n1,2))*wt(t)/2;
    A(s+n1,s+n1-1)=-u(t)*(1/h(n1));
    B(s+n1)=Q12(L(n1,2))/2;
    
    % Blocks from 1 to N/2........................................
    a=0;
    for p=1:N/2
        S=(p-1)*n1;
        for i=1:n1-1
            if L(i+1,2)==3
                A(s+i,S+i+1)=-Es12(L(i,2))*wt(p)/2;
            else
                A(s+i,S+i+1)=-Es12(L(i+1,2))*wt(p)/2;
            end
        end
        a=a+(Es12(L(n1,2))*wt(p)*y_/2);
    end
    B(s+n1)=B(s+n1)+a;
end

Z=A\B; % Solving for angular flux.

return
end

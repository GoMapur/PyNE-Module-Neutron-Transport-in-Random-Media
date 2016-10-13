%THIS IS THE AUXILIARY FUNCTION OF THE MAIN CODE "New_SN_bench"

%RICHARD VASQUES & NITIN KUMAR YADAV


function[Z,extra,n1,A,B] = New_SN_bench_solver_MC(T,m1,m2,n,N,Es1,Es2,Et1,Et2,yo,y_,Q1,Q2,u,wt,med,a)

%Building realization----------
%EXP. Random Medium....
if med==0
    m12=zeros(2,1);
    m12(1)=m1;m12(2)=m2;
    xx=rand(1);
    if xx<=m1/(m1+m2)
        mm=0;
    else
        mm=1;
    end
    s=0;i=0;
    while s<T
        i=i+1;
        x1=exprnd(m12(mod(mm,2)+1),1,1);
        s=s+x1;
        if s<=T
            x(i,1)=s;
        else
            x(i,1)=T;
        end
        x(i,2)=mod(mm,2)+1;
        mm=mm+1;
    end
else
    interval=T/n;
    %Periodic Medium........
    m12=zeros(2,1);
    m12(1)=m1;m12(2)=m2;
    mm=0;
    i=0;
    s=0;
    if a>1
	if a<m1/interval+2
		i=i+1;
		x1=(a-1)*interval;
		s=s+x1;
		x(i,1)=s;
		x(i,2)=2;
	else
    		i=i+1;
		x1=(a-m1/interval-1)*interval;
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
end
%Adding interfaces and extra points inside layers--------------

%L(*,1) is the spatial point
%L(*,2) indicates in which material is the spatial point (1 or 2).
%If L(*,2)=3, the point is an interface.
H=T/n;n1=1;i=1;j=1;
extra=zeros(n,1);
t1=i*H;
L(1,1)=0;L(1,2)=x(j,2);
if t1>x(j,1)
    extra(i)=2;
    h(n1)=x(j,1)/2;
    L(n1+1,1)=h(n1);
    L(n1+1,2)=x(j,2);
    h(n1+1)=h(n1);
    L(n1+2,1)=x(j,1);
    L(n1+2,2)=3;
    n1=n1+2;
else
    while t1<=x(j,1)
        h(n1)=H;
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
    if x(j,1)<T && L(n1)~=x(j,1)
        extra(i)=1;
        L(n1+1,1)=x(j,1);
        L(n1+1,2)=3;
        h(n1)=L(n1+1,1)-L(n1,1);
        n1=n1+1;
    end
    if t1==T && (x(j+1,1)-x(j,1))<H
        extra(i)=extra(i)+1;
        j=j+1;
        L(n1+1,1)=(x(j,1)+x(j-1,1))/2;
        L(n1+1,2)=x(j-1,2);
        h(n1)=L(n1+1,1)-L(n1,1);
        n1=n1+1;
        L(n1+1,1)=x(j,1);
        L(n1+1,2)=3;
        h(n1)=L(n1+1,1)-L(n1,1);
        n1=n1+1;
        i=i+1;
    end
end
j=2;
if x(1,1)~=T
    while i<=n
        if t1>x(j,1)
            extra(i)=extra(i)+2;
            L(n1+1,1)=(x(j,1)+x(j-1,1))/2;
            L(n1+1,2)=x(j,2);
            h(n1)=L(n1+1,1)-L(n1,1);
            n1=n1+1;
            L(n1+1,1)=x(j,1);
            L(n1+1,2)=3;
            h(n1)=L(n1+1,1)-L(n1,1);
            n1=n1+1;
        else
            while t1<=x(j,1)
                L(n1+1,1)=i*H;
                if L(n1+1,1)==x(j,1)          %checking is point is interface point?
                    L(n1+1,2)=3;
                else
                    L(n1+1,2)=x(j,2);
                end
                h(n1)=L(n1+1,1)-L(n1,1);
                n1=n1+1;
                i=i+1;
                t1=i*H;
            end
            if x(j,1)<T && L(n1)~=x(j,1)
                extra(i)=extra(i)+1;
                L(n1+1,1)=x(j,1);
                L(n1+1,2)=3;
                h(n1)=L(n1+1,1)-L(n1,1);
                n1=n1+1;
            end
        end
        j=j+1;

        if t1==T && (x(j,1)-x(j-1,1))<H
            extra(i)=extra(i)+1;
            L(n1+1,1)=(x(j,1)+x(j-1,1))/2;
            L(n1+1,2)=x(j-1,2);
            h(n1)=L(n1+1,1)-L(n1,1);
            n1=n1+1;
            L(n1+1,1)=x(j,1);
            L(n1+1,2)=3;
            h(n1)=L(n1+1,1)-L(n1,1);
            n1=n1+1;
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
minleft=T/2-m1/2;
maxright=T/2+m1/2;

% Diagonal Block of matrix up to N/2.................................
for t=1:N/2
    s=(t-1)*n1;
    A(s+1,s+1)=-u(t)*(1/h(1)+1/(h(1)+h(2)))+Et12(L(1,2))-Es12(L(1,2))*wt(t)/2;
    A(s+1,s+2)=u(t)*(1/h(1)+1/h(2));
    A(s+1,s+3)=-u(t)*(h(1)/(h(2)*(h(1)+h(2))));
    %B(s+1)=Q12(L(1,2))/2;
    B(s+1)=0.;
    for i=2:n1-1
        if L(i,2)==3
            if i==n1-1
                A(s+i,s+i)=-u(t)*(1/h(i)+1/(h(i)+h(i+1)))+Et12(L(i+1,2))-Es12(L(i+1,2))*wt(t)/2;
                A(s+i,s+i+1)=u(t)*(1/h(i)+1/h(i+1));
		                if L(i+1,1)>=minleft && L(i+1,1)<=maxright
                	      B(s+i)=u(t)*y_*(h(i)/(h(i+1)*(h(i)+h(i+1))))+Q12(L(i+1,2))/2;
		                else
			                  B(s+i)=u(t)*y_*(h(i)/(h(i+1)*(h(i)+h(i+1))));
		                end
            else
                A(s+i,s+i)=-u(t)*(1/h(i)+1/(h(i)+h(i+1)))+Et12(L(i+1,2))-Es12(L(i+1,2))*wt(t)/2;
                A(s+i,s+i+1)=u(t)*(1/h(i)+1/h(i+1));
                A(s+i,s+i+2)=-u(t)*(h(i)/(h(i+1)*(h(i)+h(i+1))));
		                if L(i+1,1)>=minleft && L(i+1,1)<=maxright
                	      B(s+i)=Q12(L(i+1,2))/2;
                		else
                			  B(s+i)=0.;
                		end
            end
        else
            A(s+i,s+i-1)=-u(t)*h(i)/(h(i-1)*(h(i-1)+h(i)));
            A(s+i,s+i)=u(t)*(h(i)-h(i-1))/(h(i)*h(i-1))+Et12(L(i,2))-Es12(L(i,2))*wt(t)/2;
            A(s+i,s+i+1)=u(t)*h(i-1)/(h(i)*(h(i-1)+h(i)));
  	    if L(i,1)>=minleft && L(i,1)<=maxright
              	B(s+i)=Q12(L(i,2))/2;
  	    else
  		B(s+i)=0.;
  	    end
        end
    end
    A(s+n1,s+n1-1)=-u(t)*h(n1)/(h(n1-1)*(h(n1-1)+h(n1)));
    A(s+n1,s+n1)=u(t)*(h(n1)-h(n1-1))/(h(n1)*h(n1-1))+Et12(L(n1,2))-Es12(L(n1,2))*wt(t)/2;
    %B(s+n1)=-u(t)*y_*h(n1-1)/(h(n1)*(h(n1-1)+h(n1)))+Q12(L(n1,2))/2;
    B(s+n1)=-u(t)*y_*h(n1-1)/(h(n1)*(h(n1-1)+h(n1)));

    % Remaining Blocks in same direction up to N/2..............
    l=t;
    if l==1 && N>2
        for p=l+1:N/2
            S=(p-1)*n1;
            for i=1:n1
                if L(i,2)==3
                    A(s+i,S+i)=-Es12(L(i+1,2))*wt(p)/2;
                else
                    A(s+i,S+i)=-Es12(L(i,2))*wt(p)/2;
                end
            end
        end
    elseif l>1 && N>2
        for p=1:l-1
            S=(p-1)*n1;
            for i=1:n1
                if L(i,2)==3
                    A(s+i,S+i)=-Es12(L(i+1,2))*wt(p)/2;
                else
                    A(s+i,S+i)=-Es12(L(i,2))*wt(p)/2;
                end
            end
        end
        for p=l+1:N/2
            S=(p-1)*n1;
            for i=1:n1
                if L(i,2)==3
                    A(s+i,S+i)=-Es12(L(i+1,2))*wt(p)/2;
                else
                    A(s+i,S+i)=-Es12(L(i,2))*wt(p)/2;
                end;
            end
        end

    end

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
    %B(s+1)=u(t)*yo*h(2)/(h(1)*(h(1)+h(2)))+Q12(L(1,2))/2;
    B(s+1)=u(t)*yo*h(2)/(h(1)*(h(1)+h(2)));
    for i=2:n1-1
        if L(i+1,2)==3
            if i==2
                A(s+i,s+i)=u(t)*(1/h(i)+1/(h(i)+h(i-1)))+Et12(L(i,2))-Es12(L(i,2))*wt(t)/2;
                A(s+i,s+i-1)=-u(t)*(1/h(i)+1/h(i-1));
		if L(i,1)>=minleft && L(i,1)<=maxright
                	B(s+i)=-u(t)*yo*(h(i)/(h(i-1)*(h(i)+h(i-1))))+Q12(L(i,2))/2;
		else
			B(s+i)=-u(t)*yo*(h(i)/(h(i-1)*(h(i)+h(i-1))));
		end
            else
                A(s+i,s+i)=u(t)*(1/h(i)+1/(h(i)+h(i-1)))+Et12(L(i,2))-Es12(L(i,2))*wt(t)/2;
                A(s+i,s+i-1)=-u(t)*(1/h(i)+1/h(i-1));
                A(s+i,s+i-2)=u(t)*(h(i)/(h(i-1)*(h(i)+h(i-1))));
		if L(i,1)>=minleft && L(i,1)<=maxright
                	B(s+i)=Q12(L(i,2))/2;
		else
			B(s+i)=0.;
		end
            end
        else
            A(s+i,s+i-1)=-u(t)*h(i+1)/(h(i)*(h(i)+h(i+1)));
            A(s+i,s+i)=u(t)*(h(i+1)-h(i))/(h(i+1)*h(i))+Et12(L(i+1,2))-Es12(L(i+1,2))*wt(t)/2;
            A(s+i,s+i+1)=u(t)*h(i)/(h(i+1)*(h(i)+h(i+1)));
	    if L(i+1,1)>=minleft && L(i+1,1)<=maxright
            	B(s+i)=Q12(L(i+1,2))/2;
	    else
		B(s+i)=0.;
	    end
        end
    end
    A(s+n1,s+n1)=u(t)*(1/h(n1)+1/(h(n1)+h(n1-1)))+Et12(L(n1,2))-Es12(L(n1,2))*wt(t)/2;
    A(s+n1,s+n1-1)=-u(t)*(1/h(n1)+1/h(n1-1));
    A(s+n1,s+n1-2)=u(t)*(h(n1)/(h(n1-1)*(h(n1)+h(n1-1))));
    %B(s+n1)=Q12(L(n1,2))/2;
    B(s+n1)=0.;

    % Remaining Blocks in same direction up to N..............
    l=t;
    if l==N/2+1 && N>2
        for p=l+1:N
            S=(p-1)*n1;
            for i=1:n1
                if L(i+1,2)==3
                    A(s+i,S+i)=-Es12(L(i,2))*wt(p)/2;
                else
                    A(s+i,S+i)=-Es12(L(i+1,2))*wt(p)/2;
                end
            end
        end
    elseif l>N/2+1 && N>2
        for p=N/2+1:l-1
            S=(p-1)*n1;
            for i=1:n1
                if L(i+1,2)==3
                    A(s+i,S+i)=-Es12(L(i,2))*wt(p)/2;
                else
                    A(s+i,S+i)=-Es12(L(i+1,2))*wt(p)/2;
                end
            end
        end
        for p=l+1:N
            S=(p-1)*n1;
            for i=1:n1
                if L(i+1,2)==3
                    A(s+i,S+i)=-Es12(L(i,2))*wt(p)/2;
                else
                    A(s+i,S+i)=-Es12(L(i+1,2))*wt(p)/2;
                end
            end
        end

    end

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

%If you want to see Flux outgoing in positive directin see Y(n) or Z(M+1)
%If you want see Flux scattering back see Y(n+1) or Z(M+2)

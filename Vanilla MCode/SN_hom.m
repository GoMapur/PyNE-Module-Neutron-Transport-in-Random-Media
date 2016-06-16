%THIS CODE SOLVES THE STEADY-STATE, MONOENERGETIC TRANSPORT EQUATION 
%IN A HOMOGENEOUS MEDIUM WITH ISOTROPIC SCATTERING AND ISOTROPIC SOURCE
%IN ROD AND SLAB (SN) GEOMETRIES. IT USES 2-POINT CENTRAL DIFFERENCING
%(ORDER H^2) WITH 3-POINT FORWARD/BACKWARD AT THE BOUNDARIES.

%MAIN OUTPUTS: transmission "TR"; reflection "RL"
%              scalar flux "SCAL"


clear all
clc
% INPUTS--------------------------------
display('Default is ROD geometry!!');
display('  ');
rod_slab=input('Enter 1 if you want to change to SLAB Geometry:   ');
if rod_slab==1
    N=1;
    while floor(N/2)~=ceil(N/2) 
        N=input(   'Enter the number of Discrete ordinate directions: ');
    end
    clc
    display('SOLVING PROBLEM IN SLAB GEOMETRY!!');
    display('  ')
else
    rod_slab=0;
    clc
    display('SOLVING PROBLEM IN ROD GEOMETRY!!');
    display('  ');
    N=2;
end
T=input('Enter the total length of the system:               ');
n=1;
while n<4
    n=input('Enter at how many points you want to calculate:     ');
end
Et=input('Enter the total cross section Sigma_t :             ');
cc=10;
while cc>1 || cc<0
    cc=input('Enter the Scattering Ratio c (between 0 and 1):     ');
end
Es = cc*Et;
yo=input('Enter the boundary value in the positive direction: ');
y_=input('Enter the boundary value in the negative direction: ');
Q=input('Enter the homogeneous isotropic Source:             ');

M=n*N;  h=T/n;
A=zeros(M,M);      B(1:M)=Q/2;

% GAUSS-LEGENDRE QUADRATURE............................
beta = (1:N-1)./sqrt(4*((1:N-1)).^2-1);
[w,x] = eig(diag(beta,-1)+diag(beta,1));
u = diag(x);
wt = 2*w(1,:)'.^2;
if rod_slab~=1
    u(1)=-1;u(2)=1;
end

% Diagonal Block of matrix up to N/2.................................
for t=1:N/2
    s=(t-1)*n;
    A(s+1,s+1)=(-11*u(t)/(6*h)+Et-Es*wt(t)/2);
    A(s+1,s+2) = 3*u(t)/h;
    A(s+1,s+3) = -3*u(t)/(2*h);
    A(s+1,s+4) = u(t)/(3*h);
    for i=2:n-1
        A(s+i,s+i-1) = -u(t)/(2*h); 
        A(s+i,s+i)=(Et-Es*wt(t)/2);
        A(s+i,s+i+1)= u(t)/(2*h);
    end
    A(s+n,s+n-1) = -u(t)/(2*h);
    A(s+n,s+n)=(Et-Es*wt(t)/2);
    B(s+n)=-u(t)*y_/(2*h)+Q/2;
    % Remaining Blocks in same direction up to N/2..............
    l=t;
    if l==1 && N>2
        for p=l+1:N/2
            S=(p-1)*n;
            for i=1:n
                A(s+i,S+i)=-Es*wt(p)/2;
            end
        end
    elseif l>1 && N>2
        for p=1:l-1
            S=(p-1)*n;
            for i=1:n
                A(s+i,S+i)=-Es*wt(p)/2;
            end
        end
        for p=l+1:N/2
            S=(p-1)*n;
            for i=1:n
                A(s+i,S+i)=-Es*wt(p)/2;
            end
        end
      
    end
    % Blocks from N/2 to N........................................
    a=0;
    for p=1:N/2
        S=(N/2+p-1)*n;
        for i=2:n
            A(s+i,S+i-1)=-Es*wt(N/2+p)/2;
        end
        a=a+(Es*wt(N/2+p)*yo/2);
    end
    B(s+1)=a+Q/2;
end


% Diagonal Block of matrix from N/2+1 to N.........................
for t=N/2+1:N
    s=(t-1)*n;
    A(s+1,s+1)=(Et-Es*wt(t)/2);
    A(s+1,s+1+1) = u(t)/(2*h);
    B(s+1)=u(t)*yo/(2*h)+Q/2;
    for i=2:n-1
        A(s+i,s+i-1) = -u(t)/(2*h); 
        A(s+i,s+i)=(Et-Es*wt(t)/2);
        A(s+i,s+i+1)= u(t)/(2*h);
    end
    A(s+n,s+n-3)=-u(t)/(3*h);
    A(s+n,s+n-2)=3*u(t)/(2*h);
    A(s+n,s+n-1)=-3*u(t)/h;
    A(s+n,s+n) =(11*u(t)/(6*h)+Et-Es*wt(t)/2);
    
    % Remaining Blocks in same direction up to N..............
    l=t;
    if l==N/2+1 && N>2
        for p=l+1:N
            S=(p-1)*n;
            for i=1:n
                A(s+i,S+i)=-Es*wt(p)/2;
            end
        end
    elseif l>N/2+1 && N>2
        for p=N/2+1:l-1
            S=(p-1)*n;
            for i=1:n
                A(s+i,S+i)=-Es*wt(p)/2;
            end
        end
        for p=l+1:N
            S=(p-1)*n;
            for i=1:n
                A(s+i,S+i)=-Es*wt(p)/2;
            end
        end
      
    end
    
    % Blocks from 1 to N/2........................................
    a=0;
    for p=1:N/2
        S=(p-1)*n;
        for i=1:n-1
            A(s+i,S+i+1)=-Es*wt(p)/2;
        end
        a=a+(Es*wt(p)*y_/2);
    end
    B(s+n)=a+Q/2;
end

X=A\B';

Y=zeros(M+N,1);

% Adding boundary conditions to the array...............
i=1;
for t=1:N/2
    s=(t-1);
    for j=i:t*n+s
        Y(j)=X(j-s);
    end
    Y(j+1)=y_;
    i=j+2;
end
i=i-1;
for t=N/2+1:N
    Y(i+1)=yo;
    i=i+1;
    for j=i+1:t*n+t
        Y(j)=X(j-t);
    end 
    i=j;
end

% Calculating Reflection, Transmission and Scalar Flux...........
RL=0;
for t=1:N/2
    s=(t-1)*n;
    RL=RL+abs(wt(t)*u(t)*Y(s+t));
end
TR=0;
for t=N/2+1:N
 s=(t-1)*n;
 S=(t-N/2-1)*n;
 k=N/2+1;
    TR=TR+wt(t)*u(t)*Y(s+n+t);
end
m=n+1;
SCAL=zeros(m,1);
for t=1:N
    s=(t-1)*m;
    for i=1:m
        SCAL(i)=SCAL(i)+wt(t)*Y(s+i);
    end 
end

x=0:h:T;
%plot(x,SCAL,'--');
plot(x,SCAL,'b'); hold on

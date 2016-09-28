%THIS CODE GENERATES ENSEMBLE-AVERAGED BENCHMARK RESULTS FOR THE
%STEADY-STATE, MONOENERGETIC TRANSPORT EQUATION IN A PERIODIC MEDIUM
%COMPOSED OF TWO MATERIALS WITH ISOTROPIC SCATTERING AND ISOTROPIC SOURCE
%IN ROD GEOMETRY. THE TOTAL NUMBER OF REALIZATIONS IS DEFINED BY THE
%SPATIAL DISCRETIZATION.

%THIS CODE CALLS THE SOLVER ROUTINE "SN_per_bench_solver", USING 2-POINT CENTRAL
%DIFFERENCING (ORDER H^2).

clear all
clc

%------------------------------------------------------------
%INPUTS
T= 20;       % total length of the system
n= 2560;     % # of points (discretization)

yo=0.;      % left boundary condition
y_=0.;      % right boundary condition

m1= 1;    % thickness of material 1 layers
m2= 1;    % thickness of material 2 layers

Et1= 1.0;    % total cross section Sigma_t1 of material 1
cc= 0.5;    % Scattering Ratio c1 of material 1
Q1= 1.0;     % homogeneous isotropic Source of material 1

Et2=0.;      % total cross section Sigma_t2 of material 2
cc2= 0.;    % Scattering Ratio c2 of material 2
Q2= 0.;     % homogeneous isotropic Source of material 2
%------------------------------------------------------------
%------------------------------------------------------------

Es1=cc*Et1;
Ea1= Et1-Es1;
Es2=cc2*Et2;
Ea2= Et2-Es2;

%weights and directions
N=2;
wt(1)= 1;
wt(2)=1;
u(1)=-1;
u(2)=1;

a=1;
reflec=0;reflec2=0;
transm=0;transm2=0;
SF=zeros(n+1,1);
SF2=zeros(n+1,1);
cond=0;

%MAIN LOOP
total = (m1+m2)/(T/n);
while(a<total+1)
    fprintf('Problem %d of %d\n', a, total);
    [Z,extra,n1,B] = SN_per_bench_solver(T,m1,m2,n,N,Es1,Es2,Et1,Et2,yo,y_,Q1,Q2,u,wt,a);
    
            % Adjusting points..........................
            X=zeros(n*N,1);
            for i=1:N/2
                X((i-1)*n+1)=Z((i-1)*n1+1);
                k=2;
                for j=2:n
                    k=k+extra(j-1);
                    X((i-1)*n+j)=Z((i-1)*n1+k);
                    k=k+1;
                end
            end            
            for i=N/2+1:N
                k=1;
                for j=1:n
                    k=k+extra(j);
                    X((i-1)*n+j)=Z((i-1)*n1+k);
                    k=k+1;
                end
            end
            
            Y=zeros(N*(n+1),1);
            % Adding boundary conditions...............
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
 
            % Calculating Reflection, Transmission and Scalar Flux......     
            RL=0;
            for t=1:N/2
                s=(t-1)*n;
                RL=RL+abs(wt(t)*u(t)*Y(s+t));
            end
            TR=0;
            for t=N/2+1:N
                s=(t-1)*n;
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
            
                reflec=reflec+RL;
                reflec2=reflec2+(RL)^2;
                transm=transm+TR;
                transm2=transm2+(TR)^2;
                for i=1:n+1
                    SF(i)=SF(i)+SCAL(i);
                    SF2(i)=SF2(i)+(SCAL(i))^2;
                end
            
            
            a=a+1;
          
end
a=a-1;
SF=SF./a;
SF2=SF2./a;

kk=T/n;
p=0:kk:T;                   % p is line along x-axis.
plot(p,SF,'r'); hold on
plot(p,SF2,'b');

save per_bench.mat


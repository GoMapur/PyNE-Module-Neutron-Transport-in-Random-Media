%THIS CODE GENERATES ENSEMBLE-AVERAGED BENCHMARK RESULTS FOR THE
%STEADY-STATE, MONOENERGETIC TRANSPORT EQUATION IN A RANDOM & PERIODIC MEDIUM
%COMPOSED OF TWO MATERIALS WITH ISOTROPIC SCATTERING AND ISOTROPIC SOURCE
%IN ROD AND SLAB (SN) GEOMETRIES. IT (1)GENERATES A RANDOM REALIZATION;
%(2)SOLVES THE DETERMINISTIC PROBLEM; (3)STORES THE RESULTS; (4)REPEAT. THE
%CONVERGENCE CRITERIA IS GIVEN BY A DYNAMICAL CHECK OF THE STATISTICAL
%ERROR OF THE OUTBOUND FLUXES AT THE BOUNDARIES, USING THE CENTRAL LIMIT
%THEOREM.

%THIS CODE CALLS THE SOLVER ROUTINE "BENCH_SOLVER", USING 2-POINT CENTRAL
%DIFFERENCING (ORDER H^2) WITH 2-POINT FORWARD/BACKWARD AT THE BOUNDARIES.

%MAIN OUTPUTS: average transmission "MTR"; average reflection "MRL"
%              average scalar flux "SF"

%RICHARD VASQUES & NITIN KUMAR YADAV


clear all
clc
% INPUTS---------------------------
%display('Default is ROD geometry & Random Medium!!');
%display('  ');
rod_slab=input('Enter 1 if you want to change to SLAB Geometry:   ');
med=input('Enter 1 if you want to change to Periodic Medium:      ');
if rod_slab==1
    N=1;
    while floor(N/2)~=ceil(N/2)
        N=input(   'Enter the number of Discrete ordinate directions: ');
    end
    if med==1
        clc
        %display('SOLVING PROBLEM IN SLAB GEOMETRY & PERIODIC Medium!!');
        %display('  ')
    else
        med=0;
        clc
        %display('SOLVING PROBLEM IN SLAB GEOMETRY & RANDOM Medium!!');
        %display('  ')
    end

else
    rod_slab=0;
    if med==1
        clc
        %display('SOLVING PROBLEM IN ROD GEOMETRY & PERIODIC Medium!!');
        %display('  ')
    else
        med=0;
        clc
        %display('SOLVING PROBLEM IN ROD GEOMETRY & RANDOM Medium!!');
        %display('  ')
    end
    N=2;
end

T=input('Enter the total length of the system:                 ');
n=1;
while n<3
    n=input('Enter at how many points you want to calculate:       ');
end
yo=input('Enter the boundary value in the positive direction:   ');
y_=input('Enter the boundary value in the negative direction:   ');

%display('  ');
m1=input('Enter average thickness of material 1:                ');
Et1=input('Enter the total cross section Sigma_t1 of material 1: ');
% cc=10;
% while cc>1 || cc<0
%     cc=input('Enter the Scattering Ratio c1 of material 1:          ');
% end
%Es1=cc*Et1;
Ea1=input('Enter the Absorbtion cross section of material 1:   ');
Es1 = Et1-Ea1;
Q1=input('Enter the homogeneous isotropic Source of material 1: ');

%display('  ');
m2=input('Enter average thickness of material 2:                ');
Et2=input('Enter the total cross section Sigma_t2 of material 2: ');
% cc=10;
% while cc>1 || cc<0
%     cc=input('Enter the Scattering Ratio c2 of material 2:          ');
% end
%Es2=cc*Et2;
Ea2=input('Enter the Absorbtion cross section of material 2:   ');
Es2 = Et2-Ea2;
Q2=input('Enter the homogeneous isotropic Source of material 2: ');


% GAUSS-LEGENDRE QUADRATURE............................
beta = (1:N-1)./sqrt(4*((1:N-1)).^2-1);
[w,v] = eig(diag(beta,-1)+diag(beta,1));
u = diag(v);
wt = 2*w(1,:)'.^2;
if rod_slab~=1
    u(1)=-1;u(2)=1;
end

%SETTING THE RANDOM # GENERATOR
%rng('shuffle');
%rng('default');
%randseed=rng;
randseed=1

a=1;
reflec=0;reflec2=0;
transm=0;transm2=0;
SF=zeros(n+1,1);
SF2=zeros(n+1,1);
cond=0;

%MAIN LOOP
rr=0;tt=0; b=0;
while(cond==0)
    [Z,extra,n1,randseed] = New_SN_bench_solver(T,m1,m2,n,N,Es1,Es2,Et1,Et2,yo,y_,Q1,Q2,u,wt,randseed,med);

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

            %ADDING VALUES
            maximum=T*max(Q1,Q2)+yo+y_;
            if abs(RL+TR)<=maximum
                reflec=reflec+RL;
                reflec2=reflec2+(RL)^2;
                transm=transm+TR;
                transm2=transm2+(TR)^2;
                for i=1:n+1
                    SF(i)=SF(i)+SCAL(i);
                    SF2(i)=SF2(i)+(SCAL(i))^2;
                end
            else
                a=a-1;
                b=b+1;
            end

            %CHECKING STATISTICAL ERROR
            if a>5000
                test1=transm/a-MTR;
                test2=reflec/a-MRL;
                if test1<0.00001
                    tt=tt+1;
                else
                    tt=0;
                end
                if test2<0.00001
                    rr=rr+1;
                else
                    rr=0;
                end
            end
            if tt>200 && rr>200
                cond=1;
            end
            MTR=transm/a;
            MTRSQ=transm2/a;
            MRL=reflec/a;
            MRLSQ=reflec2/a;
            SgTR=sqrt(abs(MTRSQ-MTR^2));
            SgRL=sqrt(abs(MRLSQ-MRL^2));
            TRERR=2*SgTR/(MTR*sqrt(a));
            RLERR=2*SgRL/(MRL*sqrt(a));
            if a>1000 && TRERR<.01 && RLERR<.01
              cond=1;
            end
            if a>100000
                cond=1;
            end
            if med==1
                cond=1;
            end
            a
            a=a+1;

end
a=a-1;
SF=SF./a;
SF2=SF2./a;
kk=T/n;
p=0:kk:T;                   % p is line along x-axis.
plot(p,SF,'r'); hold on
plot(p,SF2,'b');

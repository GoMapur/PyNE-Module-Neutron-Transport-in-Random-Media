import numpy as np

def non_homogenous(T=20, n=1000, yo=0.0, y_=0.0, m1=1, m2=1, Et1=1.0, cc=0.5, Q1=1.0, Et2=0.0, cc2=0.0, Q2=0.0, med=0):
    Es1=cc * Et1
    Ea1=Et1 - Es1
    Es2=cc2 * Et2
    Ea2=Et2 - Es2
    N=2
    beta=(arange(1,N - 1)) / sqrt(4 * ((arange(1,N - 1))) ** 2 - 1)
    w,v=eig(diag(beta,- 1) + diag(beta,1))
    u=diag(v)
    wt=2 * w(1,arange()).T ** 2
    randseed=1
    a=1
    reflec=0
    reflec2=0
    transm=0
    transm2=0
    SF=zeros(n + 1,1)
    SF2=zeros(n + 1,1)
    cond=0
    rr=0
    tt=0
    b=0
    while (cond == 0):

        Z,extra,n1,randseed = non_homogenous_solver(T,m1,m2,n,N,Es1,Es2,Et1,Et2,yo,y_,Q1,Q2,u,wt,randseed,med,a)
        X=zeros(n * N,1)
        for i in arange(1,N / 2).reshape(-1):
            X[(i - 1) * n + 1]=Z((i - 1) * n1 + 1)
            k=2
            for j in arange(2,n).reshape(-1):
                k=k + extra(j - 1)
                X[(i - 1) * n + j]=Z((i - 1) * n1 + k)
                k=k + 1
        for i in arange(N / 2 + 1,N).reshape(-1):
            k=1
            for j in arange(1,n).reshape(-1):
                k=k + extra(j)
                X[(i - 1) * n + j]=Z((i - 1) * n1 + k)
                k=k + 1
        Y=zeros(N * (n + 1),1)
        i=1
        for t in arange(1,N / 2).reshape(-1):
            s=(t - 1)
            for j in arange(i,t * n + s).reshape(-1):
                Y[j]=X(j - s)
            Y[j + 1]=y_
            i=j + 2
        i=i - 1
        for t in arange(N / 2 + 1,N).reshape(-1):
            Y[i + 1]=yo
            i=i + 1
            for j in arange(i + 1,t * n + t).reshape(-1):
                Y[j]=X(j - t)
            i=copy(j)
        RL=0
        for t in arange(1,N / 2).reshape(-1):
            s=(t - 1) * n
            RL=RL + abs_(wt(t) * u(t) * Y(s + t))
        TR=0
        for t in arange(N / 2 + 1,N).reshape(-1):
            s=(t - 1) * n
            TR=TR + wt(t) * u(t) * Y(s + n + t)
        m=n + 1
        SCAL=zeros(m,1)
        for t in arange(1,N).reshape(-1):
            s=(t - 1) * m
            for i in arange(1,m).reshape(-1):
                SCAL[i]=SCAL(i) + wt(t) * Y(s + i)
        maximum=T * max(Q1,Q2) + yo + y_
        if abs_(RL + TR) <= maximum:
            reflec=reflec + RL
            reflec2=reflec2 + (RL) ** 2
            transm=transm + TR
            transm2=transm2 + (TR) ** 2
            for i in arange(1,n + 1).reshape(-1):
                SF[i]=SF(i) + SCAL(i)
                SF2[i]=SF2(i) + (SCAL(i)) ** 2
        else:
            a=a - 1
            b=b + 1
        if a > 5000:
            test1=transm / a - MTR
            test2=reflec / a - MRL
            if test1 < 1e-05:
                tt=tt + 1
            else:
                tt=0
            if test2 < 1e-05:
                rr=rr + 1
            else:
                rr=0
        if tt > 200 and rr > 200:
            cond=1
        MTR=transm / a
        MTRSQ=transm2 / a
        MRL=reflec / a
        MRLSQ=reflec2 / a
        SgTR=sqrt(abs_(MTRSQ - MTR ** 2))
        SgRL=sqrt(abs_(MRLSQ - MRL ** 2))
        TRERR=2 * SgTR / (MTR * sqrt(a))
        RLERR=2 * SgRL / (MRL * sqrt(a))
        if a > 1000 and TRERR < 0.01 and RLERR < 0.01:
            cond=1
        if a > 100000:
            cond=1
        if med == 1:
            cond=1
        a
        a=a + 1

    a=a - 1
    SF=SF / a
    SF2=SF2 / a
    kk=T / n
    p=arange(0,T,kk)
    plot(p,SF,char('r'))
    hold(char('on'))
    plot(p,SF2,char('b'))


def non_homogenous_solver(rod_slab, N, T, m1, m2, n, N, Es1, Es2, Et1, Et2, yo, y_, Q1, Q2, u, wt, randseed, med, a):

    # %Building realization----------
    # %EXP. Random Medium....

    # random python docs: https://docs.python.org/2/library/random.html
    if not med:
        m12 = [0, m1, m2]
        if random.random() <= m1/(m1+m2):
            mm = 0
        else:
            mm = 1
        s = 0
        x = []  # is this a guaranteed loop?
        while s < T:
            s += random.expovariate(1 / m12[mm % 2 + 1])
            x += [[min(s, T), mm % 2 + 1]]
            mm += 1
    else:
        interval = T / n
        m12 = [0, m1, m2]
        mm = 0
        s = 0
        if a > 1:
            if a < m1/interval + 2:  # check order of op here
                x1 = (a-1) * interval
                s += x1
                x = [[s, 2]]
            else:  # not sure if this is inner or outer if CHECK
                x1 = (a - m1/interval - 1) * interval
                s += x1
                x = [s, mm % 2 + 1]
                mm += 1
        while s < T:
            x1 = m12[mm % 2 + 1]
            s += x1
            x += [[min(s, T), mm % 2 + 1]]
            mm += 1
    # %Adding interfaces and extra points inside layers--------------
    H = T / n
    i = 1
    j = n1 = 0
    extra = []  # using a list will be easier for extra follow periodic code
    t1 = i * H
    L = [[0, x[j, 1]]]
    if t1 > x[j, 0]:
        extra = [2]
        h = [x[j,0]/2, x[j,0]/2]
        L += [[h[n1], x[j, 1]]]
        L += [[x[j, 0], 3]]
        n1 += 2
    else:
        # Why does this begin with i = 2?
        h = []
        while t1 <= x[j, 0]:
            h += [H]
            # Begin translation
            if t1 == x(j, 0):
                # checking is point is interface point?
                L += [[t1, 3]]
            else:
                L += [[t1, x[j, 1]]]
            n1 += 1
            i += 1
            t1 = i * H
        if x[j, 0] < T and L[n1] != x[j, 0]:
            extra += [0]*(i-1) + [1]
            L += [[x[j, 0], 3]]
            h += [L[n1+1, 0] - L[n1, 0]]
            n1 += 1
        # Question: Is this guarenteed to run sequentially?
        if t1 == T and (x[j+1, 0] - x[j, 0]) < H:
            extra[-1] = extra[-1] + 1
            j = j + 1
            L += [[(x[j,0]+x[j-1,0])/2, x[j-1,1]]]
            h += [L[n1+1, 0] - L[n1, 0]]
            n1 += 1
            L += [[x[j,0], 3]]
            h += [L[n1+1, 0] - L[n1, 0]]
            n1 += 1
            i += 1
    j = 2
    if x[0, 0] != T:
        while i <= n:
            if t1 > x[j, 0]:
                extra[i-1] = extra[i-1] + 2
                L += [[(x[j,0]+x[j-1,0])/2, x[j, 1]]]
                h += [L[n1+1, 0] - L[n1, 0]]
                n1 += 1
                L += [[x[j, 0], 3]]
                h += [L[n1+1, 0] - L[n1, 0]]
                n1 += 1
            else:
                while t1 <= x[j, 0]:
                    # Checking if the point is a interface
                    if L[n1+1, 0] == x[j, 0]:
                        L += [[i*H, 3]]
                    else:
                        L += [[i*H, x[j, 1]]]
                    h += [L[n1, 0] - L[n1-1, 0]]
                    n1 += 1
                    i += 1
                    t1 = i*H
                if x[j, 0] < T and L[n1] != x[j, 0]:
                    extra[i-1] += 1
                    L += [[x[j-1, 0], 3]]
                    h += [L[n1, 0] - L[n1-1, 0]]
                    n1 += 1
            j += 1

            if t1 == T and (x[j, 0] - x[j-1, 0]) < H:
                extra[i-1] += 1
                L += [[(x[j,0]+x[j-1,0])/2, x[j-1, 1]]]
                h += [L[n1+1, 0] - L[n1, 0]]
                n1 += 1
                L += [[x[j, 0], 3]]
                h += [L[n1+1, 0] - L[n1, 0]]
                n1 += 1
                i += 1

    A = np.zeros((N * (n1+1), N * (n1+1)))
    B = np.zeros((N * (n1+1), 1))
    Et12 = [Et1, Et2]
    Es12 = [Es1, Es2]
    Q12 = [Q1, Q2]
    minleft = T / 2 - m1 / 2
    maxright = T / 2 + m1 / 2

    # Diagonal Block of matrix up to N/2.................................
    for t in range(N/2+1):
        s = t * (n1+1)
        A[s, s] = -u[t] * (1/h[0] + 1/(h[0]+h[1])) + Et12[L[0, 1]-1] - Es12[L[0,1]-1] * wt[t] / 2
        A[s, s+1] = u[t] * (1/h[0] + 1/h[1])
        A[s, s+2] = -u[t] * h[0] / (h[1] * (h[0]+h[1]))
        B[s] = 0.0
        for i in range(1, n1):
            if L[i, 1] == 3:
                if i == n1-1:
                    A[s+i-1, s+i-1] = -u[t] * (1/h[i] + 1/(h[i] + h[i+1])) + Et12[L[i+1, 1]-1] - Es12[L[i+1, 1]-1] * wt[t] / 2
                    A[s+i-1, s+i-1] = u[t] * (1/h[i] + 1/h[i+1])
                    if L[i+1, 0] >= minleft and L[i+1, 0] <= maxright:
                        B[s+i-1] = u[t] * y_ * (h[i] / (h[i+1] * (h[i] + h[i+1]))) + Q12[L[i+1, 1]-1] / 2
                    else:
                        B[s+i-1] = u[t] * y_ * (h[i] / (h[i+1] * (h[i] + h[i+1])))
                else:
                    A[s+i-1, s+i-1] = -u[t] * (1/h[i] + 1/(h[i] + h[i+1])) + Et12[L[i+1, 1]-1] - Es12[L[i+1, 1]-1] * wt[t] / 2
                    A[s+i-1, s+i] = u[t] * (1/h[i] + 1/h[i+1])
                    A[s+i-1, s+i+1] = -u[t] * h[i] / (h[i+1] * (h[i] + h[i+1]))
                    if L[i+1, 0] >= minleft and L[i+1, 0] <= maxright:
                        B[s+i-1] = Q12[L[i+1, 1]-1] / 2
                    else:
                        B[s+i-1] = 0.0
            else:
                A[s+i-1, s+i-2] = -u[t] * h[i] / (h[i-1] * (h[i-1] + h[i]))
                A[s+i-1, s+i-1] = u[t] * (h[i] - h[i-1]) / (h[i] * h[i - 1]) + Et12[L[i, 1]-1] - Es12[L[i, 1]-1] * wt[t] / 2
                A[s+i-1, s+i] = u[t] * h[i-1] / (h[i] * (h[i] + h[i-1]))
                if L[i, 0] >= minleft and L[i, 0] <= maxright:
                    B[s+i-1] = Q12[L[i, 1]-1] / 2
                else:
                    B[s+i-1] = 0.0
        A[s+n1-1, s+n1-2] = -u[t] * h[n1] / (h[n1-1] * (h[n1-1] + h[n1]))
        A[s+n1-1, s+n1-1] = u[t] * (h[n1] - h[n1-1]) / (h[n1] * h[n1-1]) + Et12[L[n1, 1]-1] - Es12[L[n1, 1]-1] * wt[t] / 2
        B[s+n1-1] = -u[t] * y_ * h[n1-1] / (h[n1] * (h[n1-1] + h[n1]))
        # Will this cause type unmatching
        l = t
        if l == 0 and N > 2:
            for p in range(l, N/2):
                S = p * (n1+1)
                for i in range(n1+1):
                    if L[i, 1] == 3:
                        A[s+i-1, S+i-1] = -Es12[L[i+1, 1]-1] * wt[p] / 2
                    else:
                        A[s+i-1, S+i-1] = -Es12[L[i, 1]-1] * wt[p] / 2
        else:
            if l > 0 and N > 2:
                for p in range(l):
                    S = p * (n1+1)
                    for i in range(n1+1):
                        if L[i, 1] == 3:
                            A[s+i-1, S+i-1] = -Es12[L[i+1, 1]-1] * wt[p] / 2
                        else:
                            A[s+i-1, S+i-1] = -Es12[L[i, 1]-1] * wt[p] / 2
                for p in range(l, N/2):
                    S = p * (n1+1)
                    for i in range(n1+1):
                        if L[i, 1] == 3:
                            A[s+i-1, S+i-1] = -Es12[L[i+1, 1]-1] * wt[p] / 2
                        else:
                            A[s+i-1, S+i-1] = -Es12[L[i, 1]-1] * wt[p] / 2
        # Blocks from N/2 to N........................................
        a = 0
        for p in range(N/2):
            S = (N/2 + p) * (n1 + 1)
            for i in range(1, n1+1):
                if L[i, 1] == 3:
                    A[s+i-1, S+i-1] = -Es12[L[i+1, 1]-1] * wt[N/2 + p] / 2
                else:
                    A[s+i-1, S+i] = -Es12[L[i-1, 1]-1] * wt[N/2 + p] / 2
            a += (Es12[L[0, 1]-1] * wt[N/2 + p] * yo / 2)
        B[s] += a

    # Diagonal Block of matrix from N/2+1 to N.........................
    for t in range(N/2, N):
        s =  * n1
        A[s, s] = u[t-1] * (h[1] - h[0]) / (h[1] * h[0]) + Et12[L[0, 1]-1] - Es12[L[0, 1]-1] * wt[t-1] / 2
        A[s, s+1] = u[t-1] * h[0] / (h[1] * (h[0] + h[1]))
        B[s] = u[t-1] * yo * h[1] / (h[0] * (h[0] + h[1]))
        for i in range(1, n1):
            if L[i, 1] == 3:
                if i == 2:
                    A[s+i-1, s+i-1] = u[t-1] * (1/h[i-1] + 1 / (h[i-1] + h[i])) + Et12[L[i-1, 1]-1] - Es12[L[i-1, 1]-1] * wt[t-1] / 2
                    A[s+i-1, s+i] = -u[t-1] * (1/h[i-1] + 1/h[i])
                    if L[i-1, 0] >= minleft and L[i-1, 0] <= maxright:
                        B[s+i-1] = -u[t-1] * yo * (h[i-1] / (h[i] * (h[i-1] + h[i]))) + Q12[L[i-1, 1]-1] / 2
                    else:
                        B[s+i-1] = -u[t-1] * yo * (h[i-1] / (h[i] * (h[i-1] + h[i])))
                else:
                    A[s+i-1, s+i-1] = u[t-1] * (1/h[i-1] + 1 / (h[i-1] + h[i])) + Et12[L[i-1, 1]-1] - Es12[L[i-1, 1]-1] * wt[t-1] / 2
                    A[s+i-1, s+i] = -u[t-1] * (1/h[i-1] + 1/h[i])
                    A[s+i-1, s+i-3] = u[t-1] * (h[i-1] / (h[i] * (h[i-1] + h[i])))
                    if L[i-1, 0] >= minleft and L[i-1, 0] <= maxright:
                        B[s+i-1] = Q12[L[i-1, 1]-1] / 2
                    else:
                        B[s+i-1] = 0.0
            else:
                A[s+i-1, s+i] = -u[t-1] * h[i] / (h[i-1] * (h[i-1] + h[i]))
                A[s+i-1, s+i-1] = u[t-1] * (h[i] - h[i-1]) / (h[i] * h[i-1]) + Et12[L[i, 1]-1] - Es12[L[i, 1]-1] * wt[t-1] / 2
                A[s+i-1, s+i] = u[t-1] * h[i-1] / (h[i] * (h[i-1] + h[i]))
                if L[i, 0] >= minleft and L[i, 0] <= maxright:
                    B[s+i-1] = Q12[L[i, 1]-1] / 2
                else:
                    B[s+i-1] = 0.0
        A[s+n1-1, s+n1-1] = u[t-1] * (1 / h[n1-1] + 1 / (h[n1-1] + h[n1])) + Et12[L[n1-1, 1]-1] - Es12[L[n1-1, 1]-1] * wt[t-1] / 2
        A[s+n1-1, s+n1] = -u[t-1] * (1 / h[n1-1] + 1 / h[n1])
        A[s+n1-1, s+n1-3] = u[t-1] * (h[n1-1] / (h[n1] * (h[n1-1] + h[n1])))
        B[s+n1-1] = 0.0
        l = copy(t)
        if l == N/2 + 1 and N > 2:
            for p in range(l, N+1):
                S = (p-1) * n1
                for i in range(n1+1):
                    if L[i, 1] == 3:
                        A[s+i-1, S+i-1] = -Es12[L[i-1, 1]-1] * wt[p-1] / 2
                    else:
                        A[s+i-1, S+i-1] = -Es12[L[i, 1]-1] * wt[p-1] / 2
        else:
            if l > N/2 + 1 and N > 2:
                for p in range(N/2, l):
                    S = (p-1) * n1
                    for i in range(n1+1):
                        if L[i, 1] == 3:
                            A[s+i-1, S+i-1] = -Es12[L[i-1, 1]-1] * wt[p-1] / 2
                        else:
                            A[s+i-1, S+i-1] = -Es12[L[i, 1]-1] * wt[p-1] / 2
                for p in range(l, N+1):
                    S = (p-1) * n1
                    for i in range(n1+1):
                        if L[i, 1] == 3:
                            A[s+i-1, S+i-1] = -Es12[L[i-1, 1]-1] * wt[p-1] / 2
                        else:
                            A[s+i-1, S+i-1] = -Es12[L[i, 1]-1] * wt[p-1] / 2
        a = 0
        for p in range(N/2+1):
            S = (p-1) * n1
            for i in range(n1):
                if L[i, 1] == 3:
                    A[s+i-1, S+i] = -Es12[L[i-1, 1]-1] * wt[p-1] / 2
                else:
                    A[s+i-1, S+i] = - Es12[L[i, 1]-1] * wt[p-1] / 2
            a = a + (Es12[L[n1-1, 1]-1] * wt[p-1] * y_ / 2)
        B[s+n1-1] = B[s+n1-1] + a
    Z = numpy.linalg.solve(A, B)
    return

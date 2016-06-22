import numpy as np

def non_homogenous(var):
    rod_slab = var[0]
    N = var[1]
    T = var[2]
    m1 = var[3]
    m2 = var[4]
    n = var[5]
    N = var[6]
    Es1 = var[7]
    Es2 = var[8]
    Et1 = var[9]
    Et2 = var[10]
    yo = var[11]
    y_ = var[12]
    Q1 = var[13]
    Q2 = var[14]
    u = var[15]
    wt = var[16]
    randseed = var[17]
    med = var[18]
    a = var[19]

    # %Building realization----------
    # %EXP. Random Medium....

    # random python docs: https://docs.python.org/2/library/random.html
    if not med:
        random.seed(randseed)
        # refers to which material, pad with 0 so index 1 == material 1
        # (slightly less confusing this way)
        m12 = [0, m1, m2]
        xx = random.random()  # python random gives a float in [0.0, 1.0)
        # turns out both Matlab and python use deterministic randomgenerators
        if xx <= m1/(m1+m2):
            mm = 0
        else:
            mm = 1
        s = 0
        i = 0
        x = None  # is this a guaranteed loop?
        while s < T:
            i += 1
            # exprnd returns a 1x1
            x1 = random.expovariate(1 / m12[mm % 2 + 1])
            # why does x1 matter? its not being used in the loop unless it
            # changes the previous value for random need to check matlab docs
            s += x1
            tmp = []
            if s <= T:
                tmp.append(s)
            else:
                tmp.append(T)
            tmp.append(mm % 2 + 1)
            if not x:
                x = np.asarray([tmp])
            else:
                x = np.vstack((x, np.asarray([tmp])))
            mm += 1
        randseed = random.getstate()  # use setstate() later
    else:
        interval = T / n
        m12 = [0, m1, m2]
        mm = 0
        i = 0
        s = 0
        if a > 1:
            if a < m1/interval + 2:  # check order of op here
                i += 1
                x1 = (a-1) * interval
                s += x1
                x = np.asarray([[s, 2]])
            else:  # not sure if this is inner or outer if CHECK
                i += 1
                x1 = (a - m1/interval - 1) * interval
                s += x1
                x = np.asarray([s, mm % 2 + 1])
                mm += 1
        while s < T:
            i += 1
            x1 = m12[mm % 2 + 1]
            s = s + x1
            tmp = []
            if s <= T:
                tmp.append(s)
            else:
                tmp.append(T)
            tmp.append(mm % 2 + 1)
            x = np.vstack((x, np.asarray([tmp])))
            mm += 1
    # %Adding interfaces and extra points inside layers--------------
    H = T / n
    n1 = 1
    i = 1
    j = 1
    extra = []  # using a list will be easier for extra follow periodic code
    t1 = i * H
    L = np.asarray([[0, x[j-1, 1]]])
    if t1 > x[j-1, 0]:
        extra[i-1] = 2
        h[n1-1] = x[j-1, 0] / 2
        L[n1, 0] = h[n1-1]
        L[n1, 1] = x[j-1, 1]
        h[n1] = h[n1 - 1]
        L[n1+1, 0] = x[j-1, 0]
        L[n1+1, 1] = 3
        n1 += 2
    else:
        while t1 <= x[j-1, 0]:
            h[n1-1] = H
            # Begin translation
            L[n1, 0] = i * H
            if L[n1, 0] == x(j-1, 0):
                # checking is point is interface point?
                L[n1, 1] = 3
            else:
                L[n1, 1] = x[j-1, 1]
            n1 = n1 + 1
            i = i + 1
            t1 = i * H
        if x[j-1, 0] < T and L[n1-1] != x[j-1, 0]:
            extra[i-1] = 1
            L[n1, 0] = x[j-1, 0]
            L[n1, 1] = 3
            h[n1-1] = L[n1, 0] - L[n1-1, 0]
            n1 += 1
        if t1 == T and (x[j, 0] - x[j-1, 0]) < H:
            extra[i-1] = extra[i-1] + 1
            j = j + 1
            L[n1, 0] = (x[j-1, 0] + x[j-2, 0]) / 2
            L[n1, 1] = x[j-2, 1]
            h[n1-1] = L[n1, 0] - L[n1-1, 0]
            n1 = n1 + 1
            L[n1, 0] = x[j-1, 0]
            L[n1, 1] = 3
            h[n1-1] = L[n1, 0] - L[n1-1, 0]
            n1 = n1 + 1
            i = i + 1
    j = 2
    if x[0, 0] != T:
        while i <= n:
            if t1 > x[j-1, 0]:
                extra[i-1] = extra[i-1] + 2
                L[n1, 0] = (x[j-1, 0] + x[j-2, 0]) / 2
                L[n1, 1] = x[j - 1, 1]
                h[n1-1] = L[n1, 0] - L[n1-1, 0]
                n1 = n1 + 1
                L[n1, 0] = x[j-1, 0]
                L[n1, 1] = 3
                h[n1-1] = L[n1, 0] - L[n1-1, 0]
                n1 = n1 + 1
            else:
                while t1 <= x[j-1, 0]:

                    L[n1, 0] = i * H
                    if L[n1, 0] == x[j-1, 0]:
                        L[n1, 1] = 3
                    else:
                        L[n1, 1] = x[j-1, 1]
                    h[n1-1] = L[n1, 0] - L[n1-1, 0]
                    n1 = n1 + 1
                    i = i + 1
                    t1 = i * H

                if x[j-1, 0] < T and L[n1-1] != x[j-1, 0]:
                    extra[i-1] = extra[i-1] + 1
                    L[n1, 0] = x[j-1, 0]
                    L[n1, 1] = 3
                    h[n1] = L[n1, 0] - L[n1-1, 0]
                    n1 = n1 + 1
            j = j + 1
            if t1 == T and (x[j-1, 0] - x[j-2, 0]) < H:
                extra[i-1] = extra[i-1] + 1
                L[n1, 0] = (x[j-1, 0] + x[j-2, 0]) / 2
                L[n1, 1] = x[j-2, 1]
                h[n1-1] = L[n1, 0] - L[n1-1, 0]
                n1 = n1 + 1
                L[n1, 0] = x[j-1, 0]
                L[n1, 1] = 3
                h[n1] = L[n1, 0] - L[n1-1, 0]
                n1 = n1 + 1
                i = i + 1

    n1 = n1 - 1
    A = zeros(N * n1, N * n1)
    B = zeros(N * n1, 1)
    Et12 = [Et1, Et2]
    Es12 = [Es1, Es2]
    Q12 = [Q1, Q2]
    minleft = T / 2 - m1 / 2
    maxright = T / 2 + m1 / 2
    for t in range(1, N/2+1):
        s = (t-1) * n1
        A[s, s] = - u[t-1] * (1/h[0] + 1/(h[0]+h[1])) + Et12[L[0, 1]-1] - Es12[L[0,1]-1] * wt[t-1] / 2
        A[s, s+1] = u[t-1] * (1/h[0] + 1/h[1])
        A[s, s+2] = -u[t-1] * h[0] / (h[1] * (h[0]+h[1]))
        B[s] = 0.0
        for i in range(2, n1):
            if L[i-1, 1] == 3:
                if i == n1-1:
                    A[s+i-1, s+i-1] = -u[t-1] * (1/h[i-1] + 1/(h[i-1] + h[i])) + Et12[L[i, 1]-1] - Es12[L[i, 1]-1] * wt[t-1] / 2
                    A[s+i-1, s+i] = u[t-1] * (1/h[i-1] + 1/h[i])
                    if L[i, 0] >= minleft and L[i, 0] <= maxright:
                        B[s+i-1] = u[t-1] * y_ * (h[i-1] / (h[i] * (h[i-1] + h[i]))) + Q12[L[i, 1]-1] / 2
                    else:
                        B[s+i-1] = u[t-1] * y_ * (h[i-1] / (h[i] * (h[i-1] + h[i])))
                else:
                    A[s+i-1, s+i-1] = -u[t-1] * (1/h[i-1] + 1/(h[i-1] + h[i])) + Et12[L[i, 1]-1] - Es12[L[i, 1]-1] * wt[t-1] / 2
                    A[s+i-1, s+i] = u[t-1] * (1/h[i-1] + 1/h[i])
                    A[s+i-1, s+i+1] = -u[t-1] * h[i-1] / (h[i] * (h[i-1] + h[i]))
                    if L[i, 0] >= minleft and L[i, 0] <= maxright:
                        B[s+i-1] = Q12[L[i, 1]-1] / 2
                    else:
                        B[s+i-1] = 0.0
            else:
                A[s+i-1, s+i] = -u[t-1] * h[i-1] / (h[i-0] * (h[i-0] + h[i-1]))
                A[s+i-1, s+i-1] = u[t-1] * (h[i-1] - h[i-0]) / (h[i - 1] * h[i - 0]) + Et12[L[i-1, 1]-1] - Es12[L[i-1, 1]-1] * wt[t-1] / 2
                A[s+i-1, s+i] = u[t-1] * h[i] / (h[i-1] * (h[i] + h[i-1]))
                if L[i-1, 0] >= minleft and L[i-1, 0] <= maxright:
                    B[s+i-1] = Q12[L[i-1, 1]-1] / 2
                else:
                    B[s+i-1] = 0.0
        A[s+n1-1, s+n1] = -u[t-1] * h[n1-1] / (h[n1] * (h[n1] + h[n1-1]))
        A[s+n1-1, s+n1-1] = u[t-1] * (h[n1-1] - h[n1]) / (h[n1-1] * h[n1]) + Et12[L[n1-1, 1]-1] - Es12[L[n1-1, 1]-1] * wt[t-1] / 2
        B[s+n1-1] = -u[t-1] * y_ * h[n1-0] / (h[n1-1] * (h[n1] + h[n1-1]))
        l = copy(t)
        if l == 1 and N > 2:
            for p in range(l+1, N/2+1):
                S = (p-1) * n1
                for i in range(1, n1+1):
                    if L[i-1, 1] == 3:
                        A[s+i-1, S+i-1] = -Es12[L[i, 1]-1] * wt[p-1] / 2
                    else:
                        A[s+i-1, S+i-1] = -Es12[L[i - 1, 1]-1] * wt[p-1] / 2
        else:
            if l > 1 and N > 2:
                for p in range(1, l):
                    S = (p-1) * n1
                    for i in range(1, n1+1):
                        if L[i-1, 1] == 3:
                            A[s+i-1, S+i-1] = -Es12[L[i, 1]-1] * wt[p-1] / 2
                        else:
                            A[s+i-1, S+i-1] = -Es12[L[i-1, 1]-1] * wt[p-1] / 2
                for p in range(l+1, N/2+1):
                    S = (p-1) * n1
                    for i in range(1, n1+1):
                        if L[i-1, 1] == 3:
                            A[s+i-1, S+i-1] = -Es12[L[i, 1]-1] * wt[p-1] / 2
                        else:
                            A[s+i-1, S+i-1] = -Es12[L[i-1, 1]-1] * wt[p-1] / 2
        a = 0
        for p in range(1, N/2+1):
            S = (N/2 + p - 1) * n1
            for i in range(2, n1+1):
                if L[i-1, 1] == 3:
                    A[s+i-1, S+i] = -Es12[L[i, 1]-1] * wt[N/2 + p - 1] / 2
                else:
                    A[s+i-1, S+i] = -Es12[L[i-1, 1]-1] * wt[N/2 + p - 1] / 2
            a = a + (Es12[L[0, 1]-1] * wt[N/2 + p - 1] * yo / 2)
        B[s] = B[s] + a
    for t in range(N/2 + 1, N+1):
        s = (t-1) * n1
        A[s, s] = u[t-1] * (h[1] - h[0]) / (h[1] * h[0]) + Et12[L[0, 1]-1] - Es12[L[0, 1]-1] * wt[t-1] / 2
        A[s, s+1] = u[t-1] * h[0] / (h[1] * (h[0] + h[1]))
        B[s] = u[t-1] * yo * h[1] / (h[0] * (h[0] + h[1]))
        for i in range(2, n1):
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
            for p in range(l + 1, N+1):
                S = (p-1) * n1
                for i in range(1, n1+1):
                    if L[i, 1] == 3:
                        A[s+i-1, S+i-1] = -Es12[L[i-1, 1]-1] * wt[p-1] / 2
                    else:
                        A[s+i-1, S+i-1] = -Es12[L[i, 1]-1] * wt[p-1] / 2
        else:
            if l > N/2 + 1 and N > 2:
                for p in range(N/2 + 1, l):
                    S = (p-1) * n1
                    for i in range(1, n1+1):
                        if L[i, 1] == 3:
                            A[s+i-1, S+i-1] = -Es12[L[i-1, 1]-1] * wt[p-1] / 2
                        else:
                            A[s+i-1, S+i-1] = -Es12[L[i, 1]-1] * wt[p-1] / 2
                for p in range(l+1, N+1):
                    S = (p-1) * n1
                    for i in range(1, n1+1):
                        if L[i, 1] == 3:
                            A[s+i-1, S+i-1] = -Es12[L[i-1, 1]-1] * wt[p-1] / 2
                        else:
                            A[s+i-1, S+i-1] = -Es12[L[i, 1]-1] * wt[p-1] / 2
        a = 0
        for p in range(1, N/2+1):
            S = (p-1) * n1
            for i in range(1, n1):
                if L[i, 1] == 3:
                    A[s+i-1, S+i] = -Es12[L[i-1, 1]-1] * wt[p-1] / 2
                else:
                    A[s+i-1, S+i] = - Es12[L[i, 1]-1] * wt[p-1] / 2
            a = a + (Es12[L[n1-1, 1]-1] * wt[p-1] * y_ / 2)
        B[s+n1-1] = B[s+n1-1] + a
    Z = numpy.linalg.solve(A, B)
    return

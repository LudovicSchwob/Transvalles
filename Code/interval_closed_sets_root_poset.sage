from collections import defaultdict


def ICS_covers(x, down, up, m):
    l = [0] if down else []
    l.extend(x)
    if up:
        l.append(2)
    if m:
        C1, C2 = [], [[]]
    else:
        C1, C2 = [[]], []
    for i in range(len(l)-1):
        c1, c2 = [], []
        for j in range(l[i], l[i+1]+1):
            if l[i] == 1 or l[i+1] == 1:
                for t in C1:
                    c1.append(t + [j])
                if j == 1:
                    for t in C2:
                        c1.append(t + [j])
                else:
                    for t in C2:
                        c2.append(t + [j])
            else:
                if j == 1:
                    for t in C1:
                        if i == 0 or t[i-1] != 1:
                            c2.append(t + [j])
                        else:
                            c1.append(t + [j])
                    for t in C2:
                        c1.append(t + [j])
                else:
                    for t in C1:
                        c1.append(t + [j])
                    for t in C2:
                        c2.append(t + [j])
        C1, C2 = c1, c2
    return [tuple(t) for t in C1], [tuple(t) for t in C2]


"""
Dm keeps intervals such that in the current slice, there only one entry equal to 1,
and which is incomparable to other 1s
"""


def Dyck(n):
    R = PolynomialRing(QQ, 'X')
    X = R.gen()
    D, Dm = {(): 1}, {}
    for k in range(2*n):
        down = not bool(k % 2)
        up = bool(k < n)
        D2, Dm2 = defaultdict(int), defaultdict(int)
        for x in D:
            C1, C2 = ICS_covers(x, down, up, False)
            for y in C1:
                D2[y] += D[x]
            for y in C2:
                Dm2[y] += D[x]
        for x in Dm:
            C1, C2 = ICS_covers(x, down, up, True)
            for y in C1:
                D2[y] += Dm[x]
            for y in C2:
                D2[y] += Dm[x]*X
        D, Dm = D2, Dm2
    return D[()]


def TSP(n):
    R = PolynomialRing(QQ, 'X')
    X = R.gen()
    D, Dm = {(): 1}, {}
    for k in range(n):
        D2, Dm2 = defaultdict(int), defaultdict(int)
        for x in D:
            C1, C2 = ICS_covers(x, True, True, False)
            for y in C1:
                D2[y] += D[x]
            for y in C2:
                Dm2[y] += D[x]
        for x in Dm:
            C1, C2 = ICS_covers(x, True, True, True)
            for y in C1:
                D2[y] += Dm[x]
            for y in C2:
                D2[y] += Dm[x]*X
        D, Dm = D2, Dm2
    return sum(D.values()) + X * sum(Dm.values())

def Rect(m, n):
    R = PolynomialRing(QQ, 'X')
    X = R.gen()
    D, Dm = {(): 1}, {}
    for k in range(n + m):
        down = bool(k < m)
        up = bool(k < n)
        D2, Dm2 = defaultdict(int), defaultdict(int)
        for x in D:
            C1, C2 = ICS_covers(x, down, up, False)
            for y in C1:
                D2[y] += D[x]
            for y in C2:
                Dm2[y] += D[x]
        for x in Dm:
            C1, C2 = ICS_covers(x, down, up, True)
            for y in C1:
                D2[y] += Dm[x]
            for y in C2:
                D2[y] += Dm[x]*X
        D, Dm = D2, Dm2
    return D[()]
"""
number of preposets:
1, 4, 29, 355, 6942... (A000798)

number of connected Tamari interval preposets (to add to the OEIS):
1, 3, 14, 82, 546, 3945


"""

def Preorders(n):
    LP = [DiGraph({k:[] for k in range(1, n+1)})]
    L = LP
    while LP!=[]:
        L2 = []
        for P in LP:
            for i in range(1, n+1):
                for j in range(1, n+1):
                    if i!=j and not P.has_edge(i,j):
                        P2 = P.copy()
                        P2.add_edge(i,j)
                        P2 = P2.transitive_closure()
                        if P2 not in L and P2 not in L2:
                            L2.append(P2)
        LP = L2
        L.extend(LP)
    return [G.copy(immutable=True) for G in L]

def preposet_order(p, q):
    for x, y, _ in p.edges():
        if x < y and not q.has_edge(x, y):
            return False
    for x, y, _ in q.edges():
        if x > y and not p.has_edge(x, y):
            return False
    return True

def is_WO_preposet(G):
    for a, c, _ in G.edges():
        if a < c:
            for b in range(a+1, c):
                if not (G.has_edge(a, b) or G.has_edge(b, c)):
                    return False
        else:
            for b in range(c+1, a):
                if not (G.has_edge(a, b) or G.has_edge(b, c)):
                    return False
    return True


def is_Tamari_preposet(G):
    for a, c, _ in G.edges():
        if a < c:
            for b in range(a+1, c):
                if not G.has_edge(a, b):
                    return False
        else:
            for b in range(c+1, a):
                if not G.has_edge(a, b):
                    return False
    return True

class WeakOrderEdge():
    def __init__(self,main,size, type = 'default'):
        self.main = main
        self.size = size
        self.type = type
    def __str__(self):
        if self.type == 'default':
            return str(self.main)
        elif self.type == 'quotient':
            S1,i,j,S2 = self.main
            s = ''
            for k in S1:
                s += S1[k]*str(k)
            s += '('+str(i)+str(j)+')'
            for k in S2:
                s += S2[k]*str(k)
            return s
        else:
            raise Exception('type must be default of quotient')
    def __repr__(self):
        return str(self.main)
    def __hash__(self):
        return str(self.main).__hash__()
    def __eq__(self, other):
        if type(other)!=WeakOrderEdge:
            return False
        return self.main == other.main and self.size == other.size
    def __le__(self,other):
        E,F = self.main,other.main
        if self.type == 'default':
            return E[1]<=F[1] and E[2]>=F[2] and F[0].issubset(E[0]) and F[3].issubset(E[3])
        elif self.type == 'quotient':
            if E[1]>F[1] or E[2]<F[2]:
                return False
            for k in E[0]:
                if k in F[0] and F[0][k] > E[0][k]:
                    return False
            for k in E[3]:
                if k in F[3] and F[3][k] > E[3][k]:
                    return False
            return True
        else:
            raise Exception('type must be default of quotient')

def WeakOrderEdges(n):
    E = []
    for i in range(1,n):
        for j in range(i+1,n+1):
            for k in range(2**(j-i-1)):
                S1,S2 = [],[]
                for a in range(i+1,j):
                    if (k>>(a-i-1))&1:
                        S1.append(a)
                    else:
                        S2.append(a)
                E.append(WeakOrderEdge((Set(S1),i,j,Set(S2)),size = n))
    return E

def is_WO_quotient_min(I, p):
    D = defaultdict(list)
    for x in I:
        S1, i, j, S2 = x.main
        D[(i,j)].append((S1, S2))
    l = []
    for i in range(1, len(p)):
        if p[i-1] > p[i]:
            for S1, S2 in D[(p[i], p[i-1])]:
                if all(k not in l for k in S1) and all(k in l for k in S2):
                    return False
        l.append(p[i-1])
    return True

def is_WO_quotient_max(I, p):
    D = defaultdict(list)
    for x in I:
        S1, i, j, S2 = x.main
        D[(i,j)].append((S1, S2))
    l = []
    for i in range(1, len(p)):
        if p[i-1] < p[i]:
            for S1, S2 in D[(p[i-1], p[i])]:
                if all(k not in l for k in S1) and all(k in l for k in S2):
                    return False
        l.append(p[i-1])
    return True

def WO_quotient_max(I, p):
    p = list(p)
    D = defaultdict(list)
    for x in I:
        S1, i, j, S2 = x.main
        D[(i,j)].append((S1, S2))
    while True:
        l, m = [], True
        for i in range(1, len(p)):
            if p[i-1] < p[i]:
                for S1, S2 in D[(p[i-1], p[i])]:
                    if all(k not in l for k in S1) and all(k in l for k in S2):
                        p[i-1], p[i] = p[i], p[i-1]
                        m = False
                        break
            l.append(p[i-1])
        if m:
            return Permutation(p)

def WO_quotient_min_max(I, lattice = True):
    n = I[0].size
    mM = []
    for p in Permutations(n):
        if is_WO_quotient_min(I, p):
            mM.append((p, WO_quotient_max(I, p)))
    if lattice:
        return LatticePoset((mM, lambda p, q: p[0].permutohedron_lequal(q[0])))
    return mM

# P = Poset((WeakOrderEdges(4), lambda p,q:p<=q))

def is_WOQ_preposet(I, G):
    """
    I = Ideal of the forcing order on S_n
    """
    D = defaultdict(list)
    for x in I:
        S1, i, j, S2 = x.main
        D[(i,j)].append((S1, S2))
    for a, c, _ in G.edges():
        if a > c:
            a2, c2 = c, a
        else:
            a2, c2 = a, c
        for b in range(a2+1, c2):
            if not (G.has_edge(a, b) or G.has_edge(b, c)):
                return False
        for S1, S2 in D[(a2, c2)]:
            t = False
            for b in S1:
                if G.has_edge(a, b):
                    t = True
                    break
            for b in S2:
                if G.has_edge(b, c):
                    t = True
                    break
            if not t:
                return False
    return True

# les préposets correspondent bien aux transvalles !
# testé jusqu'à n = 4
def Test_WO_quotients(n):
    P = Poset((WeakOrderEdges(n), lambda p,q:p<=q))
    preorders = Preorders(n)
    for I in P.order_ideals_lattice():
        if len(I) > 0:
            LG = LatticePoset(([G for G in preorders if is_WOQ_preposet(I, G)], preposet_order))
            L = WO_quotient_min_max(I)
            L2 = semidistributive_transvals(L)
            if not L2.is_isomorphic(LG):
                print('ERREUR', I)
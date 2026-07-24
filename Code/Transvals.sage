from collections import defaultdict

def join_kappa(L, j):
    if j not in L.join_irreducibles():
        raise ValueError("element is not join-irreducible")
    js = L.lower_covers(j)[0]
    while True:
        if L.is_lequal(j, L.upper_covers(js)[0]):
            if len(L.upper_covers(js)) == 1:
                return js
            js = L.upper_covers(js)[1]
        else:
            js = L.upper_covers(js)[0]

def semidistributive_transvals(L, lattice = True):
    ST = []
    for x,y in L.intervals_poset():
        I = LatticePoset(L.subposet(L.interval(x,y)))
        K = {j: join_kappa(I, j) for j in I.join_irreducibles()}
        C = set(I.coatoms())
        cc = [(c,K[c]) for c in I.atoms() if K[c] in C]
        for S in Subsets(cc):
            ST.append((I.join([s[0] for s in S]),I.meet([s[1] for s in S])))
    if lattice:
        return LatticePoset((ST,lambda p,q:L.is_lequal(p[0],q[0]) and L.is_lequal(p[1],q[1])))
    return ST

def semidistributive_face_transvals(L, lattice = True):
    SF = []
    for x in L:
        for J in Subsets(L.upper_covers(x)):
            y = L.join([x]+list(J))
            I = LatticePoset(L.subposet(L.interval(x,y)))
            K = {j: join_kappa(I, j) for j in I.join_irreducibles()}
            C = set(I.coatoms())
            cc = [(c, K[c]) for c in I.atoms() if K[c] in C]
            for S in Subsets(cc):
                SF.append((I.join([s[0] for s in S]), I.meet([s[1] for s in S])))
    if lattice:
        return LatticePoset((SF,lambda p,q:L.is_lequal(p[0],q[0]) and L.is_lequal(p[1],q[1])))
    return SF


def transvals_canonical_joinands(L):
    L2 = semidistributive_transvals(L)
    m = L.minimal_elements()[0]
    D = defaultdict(list)
    for x in L2:
        y = defaultdict(list)
        for j in L2.canonical_joinands(x):
            if j[0] == m:
                y[j[1]].append(0)
            else:
                y[j[0]].append(1)
        t = tuple(sorted(y, key = lambda t: t.__hash__()))
        D[t].append(tuple(y[j] for j in t))
    return D

def face_transvals_canonical_joinands(L):
    L2 = semidistributive_face_transvals(L)
    m = L.minimal_elements()[0]
    D = defaultdict(list)
    for x in L2:
        y = defaultdict(list)
        for j in L2.canonical_joinands(x):
            if L.is_lequal(j[0],j[1]):
                y[j[1]].append(0)
            else:
                y[j[0]].append(1)
        t = tuple(sorted(y, key = lambda t: t.__hash__()))
        D[t].append(tuple(y[j] for j in t))
    return D


"""
Canonical joinands of the transvals of Tamari give
the following statistic when restricted to Tamari:

1
1	1
3	3	1
12	12	6	1
51	62	32	10	1
238	330	200	70	15	1
1183	1833	1275	515	135	21	1
6148	10556	8218	3815	1155	238	28	1
33092	62352	53760	28000	9744	2352	392	36	1
183185	375741	355428	205044	80262	22302	4452	612	45	1
1037501	2301305	2369115	1499370	650580	203868	47040	7950	915	55	1
5988996	14284116	15897090	10947915	5212350	1812888	473088	93060	13530	1320	66	1

p = sum(int(k)*x^i for i,k in enumerate('...'.split('\t')))

row sums : 
1, 2, 7, 31, 156, 854, 4963, 30159, 189729, 1227072, 8117700, 54724420
p(2) : intervals of Tamari (A000260)
1, 3, 13, 68, 399, 2530, 16965, 118668, 857956, 6369883, 48336171, 373537388
p(3) : transvals of Tamari (A234268)
1, 4, 21, 129, 876, 6376, 48829, 388771, 3191849, 26864936, 230807084, 2017470636


Same thing for the weak order on S_n:

1
1	1
5	4	1
45	27	11	1
515	360	110	26	1

row sums : 1, 2, 10, 84, 1012
p(2) : intervals of weak order (A007767)
p(3) : transvals of weak order (A213438)

Same thing for the Stanley lattice:

1
1	1
4	3	1
24	14	6	1
164	97	35	10	1

1, 2, 8, 45, 307 (A367316)
1, 3, 14, 84, 594
1, 4, 22, 147, 1121
"""

############### INTERVAL-CLOSED SUBSETS ###############

# for posets (distributive lattices) :

def disjoint_antichains(P):
    L = []
    for a in P.antichains():
        P2 = P.subposet([x for x in P.order_ideal(a) if x not in a])
        for b in P2.antichains():
            L.append((tuple(a), tuple(b)))
    return Poset((L, lambda p, q: Set(P.order_ideal(p[0])).issubset(Set(P.order_ideal(q[0]))) and Set(P.order_ideal(p[1])).issubset(Set(P.order_ideal(q[1])))))

def disjoint_antichains2(P):
    L = []
    for a in P.antichains():
        P2 = P.subposet([x for x in P.order_ideal(a) if x not in a])
        for b in P2.antichains():
            a2, c = [], []
            for x in a:
                if any(P.is_lequal(y, x) for y in b):
                    a2.append(x)
                else:
                    c.append(x)
            L.append((tuple(a2), tuple(b), tuple(c)))
    return Poset((L, lambda p, q: Set(P.order_ideal(list(p[0])+list(p[2]))).issubset(Set(P.order_ideal(list(q[0])+list(q[2])))) and Set(P.order_filter(list(p[1])+list(p[2]))).issubset(Set(P.order_filter(list(q[1])+list(q[2]))))))


# for two-acyclic factorization systems (semidistributive lattices)

def min_max_subsets(G):
    """
    returns subsets S of G such that for all x in S, 
    there cannot exist y,z in S such that y -> x -> z
    """
    lG = list(G)
    D = {x: i for i, x in enumerate(lG)}
    n = len(G)
    # 1: mins, 2: max, 3: nor min or max
    l, L = [([],[],[])], []
    def min_max_aux(t, k):
        a, b, c = t
        a2, b2, c2 = list(a), list(b), []
        m, M = True, True
        for x in a:
            if G.has_edge(k, x):
                return
            if G.has_edge(x, k):
                m = False
        for x in b:
            if G.has_edge(x, k):
                return
            if G.has_edge(k, x):
                if not m:
                    return
                M = False
        for x in c:
            if G.has_edge(k, x):
                if not m:
                    return
                M = False
                if G.has_edge(x, k):
                    return
                b2.append(x)
            else:
                if G.has_edge(x, k):
                    if not M:
                        return
                    m = False
                    a2.append(x)
                else:
                    c2.append(x)
        if m and M:
            c2.append(k)
        elif M:
            b2.append(k)
        else:
            a2.append(k)
        return a2, b2, c2
    while len(l) > 0:
        l2 = []
        for t in l:
            m = max(-1 if len(a)==0 else max([D[x] for x in a]) for a in t) + 1
            for k in range(m, n):
                t2 = min_max_aux(t, lG[k])
                if t2 != None:
                    l2.append(t2)
        L.extend(l)
        l = l2
    return L

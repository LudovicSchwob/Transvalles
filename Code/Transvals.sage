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

def intervals_canonical_joinands(L):
    L2 = semidistributive_transvals(L)
    D = defaultdict(list)
    m = L.minimal_elements()[0]
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

def transvals_canonical_joinands(L):
    L2 = semidistributive_transvals(L)
    D = defaultdict(list)
    m = L.minimal_elements()[0]
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

1
0	1
1	1	1
5	3	3	1
12	24	8	6	1
52	85	70	20	10	1
225	387	345	165	45	15	1
939	1974	1722	1050	350	91	21	1
4249	9508	9856	5908	2674	700	168	28	1
19672	47700	54000	36960	17136	6048	1344	288	36	1
92604	244780	298140	226440	115920	44016	12600	2490	465	45	1
445820	1264604	1674860	1368180	781440	321090	103026	24750	4455	715	55	1

####################

Same thing for the weak order on S_n:

1
1	1
5	4	1
45	27	11	1
515	360	110	26	1

row sums : 1, 2, 10, 84, 1012
p(2) : intervals of weak order (A007767)
p(3) : transvals of weak order (A213438)

1
0	1
2	2	1
28	8	8	1
240	214	38	22	1

##############################

Same thing for the Stanley lattice:

1
1	1
4	3	1
24	14	6	1
164	97	35	10	1
1277	747	271	74	15	1
10926	6324	2313	637	140	21	1
100392	57718	21128	5956	1345	244	28	1

1, 2, 8, 45, 307 (A367316)
1, 3, 14, 84, 594 (A005700)
1, 4, 22, 147, 1121


1
0	1
2	1	1
15	5	3	1
93	53	11	6	1
741	372	129	24	10	1
6398	3148	1047	267	50	15	1
58974	29009	9289	2491	510	97	21	1
575362	280623	90391	23220	5329	935	175	28	1
"""

######################## MIN-MAX SUBSETS ########################


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


def min_max_subsets2(G):
    G = G.canonical_label()
    n = len(G)
    # a = indicatrice des éléments maximaux
    # b = indicatrice des éléments minimaux
    l = [(0, 0)]
    for k in range(n):
        l2 = []
        x_in = sum(1<<x for x in G.neighbors_in(k) if x < k)
        x_out = sum(1<<x for x in G.neighbors_out(k) if x < k)
        for a, b in l:
            l2.append((a, b))
            t_in = x_in & (a | b)
            t_out = x_out & (a | b)
            if t_in == 0:
                if t_out == 0:
                    l2.append((a|(1<<k), b|(1<<k)))
                elif t_out & a == t_out:
                    l2.append((a, (b|(1<<k)) - (t_out & b)))
            elif t_in & b == t_in:
                if t_out == 0:
                    l2.append(((a|(1<<k)) - (t_in & a), b))
        l = l2
    def minmax_stat(a, b):
        x, s = a & b, 0
        while x != 0:
            s += x & 1
            x >>= 1
        return s
    R = PolynomialRing(QQ, 'x')
    x = R.gen()
    return sum(x^minmax_stat(a, b) for a, b in l)


def interval_labels(L):
    E = {}
    for x, y in L.cover_relations():
        E[(x, y)] = Edge_to_JoinIrr(L, x, y)
    S = set()
    for x, y in L.intervals_poset():
        I = L.subposet(L.interval(x, y))
        S.add(Set([E[(a, b)] for a, b in I.cover_relations()]))
    return list(S)

"""
Number of sets of labels of edges of Tamari intervals :
7, 33, 182, 1104

"""

def interval_labels2(L):
    E = {}
    for x, y in L.cover_relations():
        E[(x, y)] = Edge_to_JoinIrr(L, x, y)
    S = set()
    for x, y in L.intervals_poset():
        I = L.subposet(L.interval(x, y))
        S1 = Set([E[(x, a)] for a in I.upper_covers(x)])
        S2 = Set([E[(b, y)] for b in I.lower_covers(y)])
        S.add((S1,S2))
    return list(S)

def labelled_ideals(P, L):
    """
    P = Poset, L = dictionary {element: label}
    """
    if len(P) == 0:
        return 1
    k = L[P[0]]
    l, l2 = [], []
    for x in P:
        if L[x] == k:
            l.append(x)
        else:
            l2.append(x)
    P2 = P.subposet(l2)
    P3 = P2.subposet([x for x in P2 if not any(P.is_lequal(y, x) or P.is_lequal(x, y) for y in l)])
    r = labelled_antichains(P2, L) + labelled_antichains(P3, L)
    print(P2.hasse_diagram().to_dictionary())
    print(P3.hasse_diagram().to_dictionary())
    print(P.hasse_diagram().to_dictionary())
    print(r)
    return r
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


"""
Canonical joinands of the transvals of Tamari give
the following statistic when restricted to Tamari:

1
1	1
3	3	1
12	12	6	1
51	62	32	10	1
238	330	200	70	15	1

row sums : 1, 2, 7, 31, 156, 854 (quid ?)
p(2) : intervals of Tamari (A000260)
p(3) : transvals of Tamari (A390811)


"""





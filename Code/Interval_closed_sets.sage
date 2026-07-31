def join_dis_low(P, C):
    """
    C = coloring of the poset P
    """
    return LatticePoset(([Set([C[i] for i in I]) for I in P.order_ideals_lattice()], lambda p, q: p.issubset(q)))

def join_dis_coloring(L):
    maxs = {}
    for j in L.join_irreducibles_poset():
        M = []
        l = [L.lower_covers(j)[0]]
        while l != []:
            l2 = []
            for x in l:
                m = True
                for y in L.upper_covers(x):
                    if not L.is_lequal(j, y):
                        l2.append(y)
                        m = False
                if m:
                    M.append(x)
            l = l2
        maxs[j] = Set(M)
    return L.join_irreducibles_poset(), {x: maxs[x][0] for x in maxs}

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

def interval_closed_sets(P):
    L = []
    for a in P.antichains():
        P2 = P.subposet([x for x in P.order_ideal(a) if x not in a])
        for b in P2.antichains():
            L.append(Set(list(a)+P2.order_filter(b)))
    return LatticePoset((L, lambda p, q: p.issubset(q)))
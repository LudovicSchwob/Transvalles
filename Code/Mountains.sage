from collections import defaultdict

"""
Number of nested mountains: 
2, 7, 34, 199, 1308, 9300, 69978, 549559 (A393920)
 = number of subsets of the AR quiver closed by extension

what about other AR quivers ?
It seems that the number of extension-closed subsets only depend on the type
(True: extensions can be obtained from the root poset)

D_4 : 496
D_5 : 9884

Nested mountains with umbrella:

same sequence as extension-closed subset of the AR quiver,
but with different offset

"""


def nested_mountains(W):
    S = dict(RootSystem(W).root_space().simple_roots())
    Phi = RootSystem(W).root_space().root_poset()
    def join(a, b):
        a = {s: c for s, c in a.element}
        b = {s: c for s, c in b.element}
        j = sum(max(a[s],b[s])*S[s] if s in b else a[s]*S[s] for s in a)
        for s in b:
            if s not in a:
                j += b[s]*S[s]
        return j
    def meet(a, b):
        a = {s: c for s, c in a.element}
        b = {s: c for s, c in b.element}
        return sum(min(a[s],b[s])*S[s] if s in b else 0 for s in a)
    l = []
    for s in Subsets(Phi):
        t = True
        for a in s:
            for b in s:
                j, m = join(a, b), meet(a, b)
                if (j in Phi and Phi(j) not in s) or (m in Phi and Phi(m) not in s):
                    t = False
                    break
        if t:
            l.append(s)
    return l

def nested_mountains(W):
    S = dict(RootSystem(W).root_space().simple_roots())
    Phi = RootSystem(W).root_space().root_poset()
    def join(a, b):
        return sum(max(a.element[s], b.element[s])*S[s] for s in S)
    def meet(a, b):
        return sum(min(a.element[s], b.element[s])*S[s] for s in S)
    implications = {}
    for a in Phi:
        for b in Phi:
            if not (Phi.is_lequal(a, b) or Phi.is_lequal(b, a)):
                j, m = join(a, b), meet(a, b)
                l = []
                if j in Phi:
                    l.append(Phi(j))
                if m in Phi:
                    l.append(Phi(m))
                implications[(a, b)] = l
    def closure(S, a):
        l = [a]
        S = list(S)
        S.append(a)
        while len(l) > 0:
            l2 = []
            for x in S:
                for y in l:
                    if (x, y) in implications:
                        for z in implications[(x, y)]:
                            if z not in S and z not in l2:
                                l2.append(z)
            l = l2
            S.extend(l2)
        return Set(S)
    L, l = set([Set()]), [Set()]
    while len(l) > 0:
        l2 = []
        for S in l:
            for a in Phi:
                if a not in S:
                    S2 = closure(S, a)
                    if S2 not in L and S2 not in l2:
                        l2.append(S2)
        L = L.union(l)
        l = l2
    return list(L)

def coclosed_nested_mountains(W):
    Phi = RootSystem(W).root_space().root_poset()
    sums = defaultdict(list)
    for a in Phi:
        for b in Phi:
            if not (Phi.is_lequal(a, b) or Phi.is_lequal(b, a)):
                c = a.element + b.element
                if c in Phi:
                    sums[Phi(c)].append((a, b))
    L = []
    for S in nested_mountains(W):
        if all(a in S or b in S for x in S for a, b in sums[x]):
            L.append(S)
    return L
                 

"""
1	2	1
3	10	7	1
12	48	52	16	1
52	248	360	185	30	1
241	1334	2439	1792	519	50	1

cf. A000256

2	4	1
6	17	10	1
22	80	76	20	1
91	403	533	245	35	1
408	2128	3619	2444	644	56	1

cf. A000139
"""


########### EXTENSION-CLOSED SUBSETS OF AR QUIVERS #########


def extension_closed_subsets(G):
    n = len(G)
    def dimension_vector(x):
        l = n*[0]
        for i, c in x.dimension_vector():
            l[i-1] = c
        return tuple(l)
    AR = G.auslander_reiten_quiver()
    AR = AR.digraph().copy(immutable = False)
    AR.relabel(dimension_vector)
    extensions = defaultdict(list)
    for a in AR:
        for b in AR:
            if not (all(i <= j for i, j in zip(a, b)) or all(i >= j for i, j in zip(a, b))):
                j, m = tuple(max(i, j) for i, j in zip(a, b)), tuple(min(i, j) for i, j in zip(a, b))
                l = []
                if j in AR:
                    l.append(j)
                if m in AR:
                    l.append(m)
                extensions[(a, b)] = l
    def closure(S, a):
        l = [a]
        S = list(S)
        S.append(a)
        while len(l) > 0:
            l2 = []
            for x in S:
                for y in l:
                    for z in extensions[(x, y)]:
                        if z not in S and z not in l2:
                            l2.append(z)
            l = l2
            S.extend(l2)
        return Set(S)
    L, l = set([Set()]), [Set()]
    while len(l) > 0:
        l2 = []
        for S in l:
            for a in AR:
                if a not in S:
                    S2 = closure(S, a)
                    if S2 not in L and S2 not in l2:
                        l2.append(S2)
        L = L.union(l)
        l = l2
    return list(L)
            

def extension_closed_subsets_with_max(G):
    l = extension_closed_subsets(G)
    n = len(G)
    m = n*[0]
    for s in l:
        for x in s:
            m = [max(i, j) for i, j in zip(m, x)]
    m = tuple(m)
    return [s for s in l if m in s]
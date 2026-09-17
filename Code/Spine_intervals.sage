def Spine_intervals_graph(G):
    """
    G must be a acyclic & two-acylic factorization system
    (corresponds to an extremal and semidistributive lattice)
  
    For Tamari : Pairs(G) gives a lattice of size A001764
    """
    G.remove_loops()
    P = Poset(G)
    G2 = DiGraph()
    for x in G:
        G2.add_vertex((x,0))
        G2.add_vertex((x,1))
    for x, y, _ in G.edges():
        G2.add_edge((x,0),(y,0))
        G2.add_edge((x,1),(y,1))
    for x in G:
        for y in P.order_filter([x]):
            G2.add_edge((x,0),(y,1))
    return G2

def Spine_intervals(L, lattice = False):
    S = Spine(L)
    l = []
    for x, y in L.intervals_poset():
        if any(z in S for z in L.interval(x, y)):
            l.append((x,y))
    if lattice:
        return LatticePoset((l, lambda x, y: L.is_lequal(x[0], y[0]) and L.is_lequal(x[1], y[1])))
    return l

"""
L = LatticePoset({0: [1,2], 1:[3], 2:[4,5], 3:[4], 4:[6,7], 5:[7], 6:[8], 7:[9], 8:[9]})
G = SDLGraph(L)
G2 = Spine_intervals_graph(G)
L2 = Spine_intervals(L, True)
L3 = MOPLattice(G2)
L2.is_isomorphic(L3)
"""
    

def Spine_transvals_graph(G):
    """
    G must be a acyclic & two-acylic factorization system
    (corresponds to an extremal and semidistributive lattice)
  
    For Tamari : 1, 4, 17, 81, 412, 2192, 12049 (A121545 ????)
    """
    G.remove_loops()
    P = Poset(G)
    G2 = DiGraph()
    for x in G:
        G2.add_vertex((x,0))
        G2.add_vertex((x,1))
    for x, y, _ in G.edges():
        G2.add_edge((x,0),(y,0))
        G2.add_edge((x,1),(y,1))
    for x in G:
        for y in P.order_filter([x]):
            if y != x:
                G2.add_edge((x,0),(y,1))
    return G2


# C'est pas encore ça !!!!!
def Spine_transvals(L, lattice = False):
    S = Spine(L)
    l = []
    for x, y in L.intervals_poset():
        I = LatticePoset(L.subposet(L.interval(x, y)))
        K = {j: join_kappa(I, j) for j in I.join_irreducibles()}
        C = set(I.coatoms())
        cc = [(c,K[c]) for c in I.atoms() if K[c] in C]
        for Ss in Subsets(cc):
            x2, y2 = I.join([s[0] for s in Ss]), I.meet([s[1] for s in Ss])
            if any(z in S for z in L.interval(x2, y)) and any(z in S for z in L.interval(x, y2)):
                l.append((x2, y2))
    if lattice:
        return LatticePoset((l, lambda x, y: L.is_lequal(x[0], y[0]) and L.is_lequal(x[1], y[1])))
    return l


# ne donne pas un TAFS, à corriger pour obtenir le treillis des nested topologies
# il semblerait qu'il faille rajouter des arêtes dans les deux copies de G, mais lesquelles ?
def Spine_faces_graph(G):
    """
    G must be a acyclic & two-acylic factorization system
    (corresponds to an extremal and semidistributive lattice)
  
    For Tamari :
    """
    G.remove_loops()
    P = Poset(G)
    G2 = DiGraph()
    for x in G:
        G2.add_vertex((x,0))
        G2.add_vertex((x,1))
    for x, y, _ in G.edges():
        G2.add_edge((x,0),(y,0))
        G2.add_edge((x,1),(y,1))
        G2.add_edge((x,1),(y,0))
    for x in G:
        for y in P.order_filter([x]):
            G2.add_edge((x,0),(y,1))
    return G2


def test_spine_faces(L):
    """
    échoue pour L = LatticePoset({0: [1,2], 1:[3,4], 2:[4,7], 3:[5], 4:[6], 5:[6], 6:[8], 7:[8]})
    """
    G = SDLGraph(L)
    L2 = Spine_faces(L, True)
    G2 = SDLGraph(L2)
    G0, G1 = [], []
    for j in G2:
        jd = L2.lower_covers(j)[0]
        if jd[0] == j[0]:
            G0.append(j)
        elif jd[1] == j[1]:
            G1.append(j)
        else:
            print('ERREUR 404', j)
    return G, G2.subgraph(G0), G2.subgraph(G1)

def test_spine_faces_Tamari(n):
    L = NestedLattice(n)
    vertex_colors = {'red': [], 'green': []}
    G = SDLGraph(L)
    for j in G:
        jstar = L.lower_covers(j)[0]
        for i in range(n):
            if j[i] != jstar[i]:
                if j[i][0] == 0:
                    vertex_colors['green'].append(j)
                else:
                    vertex_colors['red'].append(j)
    return G, vertex_colors

def spine_face_relabel1(x):
    """
    x = join-irréductible de NestedLattice(n) de type 'red'
    renvoie un élément de TamariGraph(n)
    """
    n = len(x)
    for i in range(n - 1):
        if x[i+1][0] == x[i+1][1]:
            t = n * [0]
            t[i] = x[i][1] - i
            return tuple(t)

def spine_face_relabel2(x):
    """
    x = join-irréductible de NestedLattice(n) de type 'green'
    renvoie un élément de TamariGraph(n)
    """
    n = len(x)
    for i in range(n-1, 0, -1):
        for j in range(i-1, -1, -1):
            if x[i] == x[j]:
                t = n * [0]
                t[j] = i - j
                return tuple(t)

def Spine_faces_bv1(n):
    G, vertex_colors = test_spine_faces_Tamari(n)
    G2 = G.subgraph(vertex_colors['red'])
    G2.relabel(spine_face_relabel1)
    L = MOPLattice(G2)
    def Max(X):
        t = n*[0]
        for x in X:
            t = [max(i, j) for i, j in zip(t, x)]
        return tuple(t)
    return L.relabel(Max)

def Spine_faces_bv2(n):
    G, vertex_colors = test_spine_faces_Tamari(n)
    G2 = G.subgraph(vertex_colors['green'])
    G2.relabel(spine_face_relabel2)
    L = MOPLattice(G2)
    def Max(X):
        t = n*[0]
        for x in X:
            t = [max(i, j) for i, j in zip(t, x)]
        return tuple(t)
    return L.relabel(Max)
            
def bvs(n):
    L = []
    for t in cartesian_product([range(n-k) for k in range(n)]):
        if all(t[i]-t[j]>=j-i for i in range(n-1) for j in range(i+1,i+t[i]+1)):
            L.append(t)
    return L

"""
G, vertex_colors = test_spine_faces_Tamari(5)
G2 = G.subgraph(vertex_colors['green'])
G2.relabel(spine_face_relabel2)
E = set(TamariGraph(4).edges())
for e in G2.edges():
    if e not in E:
        print(e)

L = Spine_faces_bv1(5)
for t in bvs(5):
    if t not in L:
        print(t)

G, vertex_colors = test_spine_faces_Tamari(5)
G.relabel( lambda x: (spine_face_relabel1(x),0) if x in vertex_colors['red'] else (spine_face_relabel2(x),1))
for e in G.edges():
    if e[0][0] == e[1][0]:
        G.delete_edge(e)

15, 56, 209, 780, 2911
"""

def Tamari_Ribs(n, lattice = True):
    G, vertex_colors = test_spine_faces_Tamari(n)
    G2 = G.subgraph(vertex_colors['red'])
    G3 = G.subgraph(vertex_colors['green'])
    G2.relabel(spine_face_relabel1)
    G3.relabel(spine_face_relabel2)
    for e in G3.edges():
        G2.add_edge(e)
    if not is_TAFS(G2):
        print('NOT A 2-ACYCLIC FACTORIZATION SYSTEM')
    return MOPLattice(G2, lattice)

"""
size of Tamari_Ribs(n) :
1, 2, 5, 12, 28, 64, 144, 320 (A045623)
"""

####################################

def lattice_popdown(L, x):
    return L.meet([x, L.meet(L.lower_covers(x))])

def lattice_popup(L, x):
    return L.join([x, L.join(L.upper_covers(x))])

def test_spine_faces2(L):
    """
    La construction des treillis ne semble pas tout à fait bonne non plus
    plutôt enlever les éléments x tels qu'il n'y aie pas de face (x,y) minimale intersectant S ?
    """
    G,G0,G1 = test_spine_faces(L)
    S = Spine(L)
    L0, L1 = MOPLattice(G0), MOPLattice(G1)
    S0 = L.subposet([x for x in L if any(y in S for y in L.interval(lattice_popdown(L, x), x))])
    S1 = L.subposet([x for x in L if any(y in S for y in L.interval(x, lattice_popup(L, x)))])
    print(len(L), len(L0), len(L1))
    return S0.is_isomorphic(L0), S1.is_isomorphic(L1)

# !! ne donne pas un treillis semidistributif pour L = JurassianLattice([4,1,2,3])
def Spine_faces(L, lattice = False):
    S = Spine(L)
    l = []
    for x in L:
        for s in Subsets(L.upper_covers(x)):
            y = L.join([x]+list(s))
            if any(z in S for z in L.interval(x, y)):
                l.append((x,y))
    if lattice:
        return LatticePoset((l, lambda x, y: L.is_lequal(x[0], y[0]) and L.is_lequal(x[1], y[1])))
    return l

"""
L = TamariLattice(5)
L2 = Spine_faces(L, True)
L3 = NestedLattice(5)
L2.is_isomorphic(L3)
"""



def Spine_transvals_restriction(L):
    G = SDLGraph(L)
    L = MOPLattice(G)
    G2 = Spine_transvals_graph(G)
    rest = []
    for X in MOPLattice(G2, False):
        s1, s2 = [x[0] for x in X if x[1]==0], [x[0] for x in right_orthogonal(G2, X) if x[1]==1]
        rest.append((Set(s1), Set(left_orthogonal(G, s2))))
    L2 = LatticePoset((rest, lambda p,q: p[0].issubset(q[0]) and p[1].issubset(q[1])))
    return L2

# pas nécessaire ! La restriction suffit
def Spine_transvals_cloture_restriction(L):
    G = SDLGraph(L)
    L = MOPLattice(G)
    G2 = Spine_transvals_graph(G)
    rest = []
    for X in MOPLattice(G2, False):
        s1, s2 = [x[0] for x in X if x[1]==0], [x[0] for x in right_orthogonal(G2, X) if x[1]==1]
        rest.append((Set(left_orthogonal(G, right_orthogonal(G, s1))),Set(left_orthogonal(G, s2))))
    L2 = LatticePoset((rest, lambda p,q: p[0].issubset(q[0]) and p[1].issubset(q[1])))
    return L2


def Spine_transvals(L):
    G = SDLGraph(L)
    D = MOPLattice(G).is_isomorphic(L, certificate= True)[1]
    L2 = Spine_transvals_restriction(L).relabel(lambda x: (D[x[0]],D[x[1]]))
    L3 = semidistributive_transvals(L)
    if not L2.is_sublattice(L3):
        print('spine transvals do not form a sublattice of transvals')
    if not all(x in L3 for x in L2):
        print('spine transvals are not transvals')
    return list(L2), [x for x in L3 if x not in L2]
    
def Test_spine_transvals(L):
    l1, l2 = Spine_transvals(L)
    S = Spine(L)
    def test(x, y):
        a, b = L.meet([x, y]), L.join([x, y])
        return (any(z in S for z in L.interval(a, y)) or any(z in S for z in L.interval(x, b))) and (all(z in S for z in L.interval(a, x)) or all(z in S for z in L.interval(y, b)))
    for x, y in l1:
        if not test(x, y):
            print(f'le transvalle {(x, y)} ne vérifie pas le test')
    for x, y in l2:
        if test(x, y):
            print(f"{(x, y)} vérifie pas le test mais n'est pas un transvalle")
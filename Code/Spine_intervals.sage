def Nested_intervals_graph(G):
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
G2 = Nested_intervals_graph(G)
L2 = Spine_intervals(L, True)
L3 = MOPLattice(G2)
L2.is_isomorphic(L3)
"""
    

def Nested_transvals_graph(G):
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
# il semblerait qu'il faille rajouter des arêtes dans la copie (G, 1), mais lesquelles ?
def Nested_faces_graph(G):
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

# marche pas - pas trouvé de description des irréductibles 
def test_nested_faces(L):
    G = SDLGraph(L)
    G2 = Nested_faces_graph(G)
    L2 = Spine_faces(L, True)
    G3 = SDLGraph(L2)
    J = L.join_irreducibles()
    def relabel(j):
        if j[0] == j[1] or j[0] not in J:
            return (j[1], 0)
        else:
            return (j[0], 1)
    for x in G3:
        print(x, relabel(x))
    G3.relabel(relabel)
    return G3

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




from collections import defaultdict

def ConvexTopologies(n):
    LT = [((0,0),)]
    for k in range(1,n):
        LT2 = []
        for D in LT:
            C,a = [k],k
            while a!=0:
                a = D[a-1][0]
                C.append(a)
            DS,a = [],k
            for x,I in enumerate(D):
                if k-1==I[1]:
                    if I[0]==a:
                        DS[-1].append(x)
                    else:
                        DS.append([x])
                        a = I[0]
            for i in range(len(DS)+1):
                D2 = list(D)
                for j in range(i):
                    for x in DS[j]:
                        D2[x] = (D2[x][0],k)
                for c in C:
                    if i==0 or D[DS[i-1][0]][0]<=c:
                        D3 = D2.copy()
                        D3.append((c,k))
                        LT2.append(tuple(D3))
        LT = LT2
    return LT

def CT_order(T1,T2):
    for A,B in zip(T1,T2):
        if A[0]>B[0] or A[1]>B[1]:
            return False
    return True

def is_nested(T):
    for i in range(len(T)):
        a,b = T[i]
        for j in range(1,len(T)):
            a2,b2 = T[j]
            if not (a<=a2<=b2<=b or a2<=a<=b<=b2):
                return False
    return True

def ConvexLattice(n):
    return LatticePoset((ConvexTopologies(n),lambda p,q:CT_order(p,q)))

def NestedLattice(n):
    return LatticePoset(([T for T in ConvexTopologies(n) if is_nested(T)],lambda p,q:CT_order(p,q)))

def is_interval(T):
    MN = set()
    for s in T:
        if s in MN:
            return False
        MN.add(s)
    return True

#marche pas encore
def is_face(T):
    if not is_interval(T):
        return False
    DI = {}
    for i in range(len(T)):
        a,b = min(T[i]),max(T[i])
        DI[i] = set()
        for j in range(len(T)):
            a2,b2 = min(T[j]),max(T[j])
            if i!=j and a2<=a and b<=b2:
                DI[i].add(j)
    for i in range(len(T)-1):
        a,b = min(T[i]),max(T[i])
        for j in range(1,len(T)):
            a2,b2 = min(T[j]),max(T[j])
            if b+1==a2 and DI[i]!=DI[j]:
                return False
    return True

def is_face2(T):
    """
    T : topologie convexe sur une chaîne
    """
    for i in range(len(T)-1):
        a,b = T[i]
        for j in range(1,len(T)):
            a2,b2 = T[j]
            if b+1==a2 or (b>=a2 and a<a2 and b<b2):
                return False
    return True

def is_singleton(T):
    """
    T : topologie convexe sur une chaîne
    """
    if not is_interval(T):
        return False
    for i in range(len(T)-1):
        a,b = T[i]
        for j in range(1,len(T)):
            a2,b2 = T[j]
            if b+1==a2 or (b>=a2 and a<a2 and b<b2):
                return False
    return True

#vérifie si pour intervalle I, l'ensemble des k tels que MN(k)=I est un intervalle.
def is_connected(T):
    MN = {}
    for k in range(len(T)):
        if T[k] in MN:
            if k-1 not in MN[T[k]]:
                return False
            MN[T[k]].append(k)
        else:
            MN[T[k]] = [k]
    return True

#donne le treillis booléen de dim 2(n-1)
def is_connected2(T):
    """
    T : topologie convexe sur une chaîne
    """
    for i in range(len(T)-1):
        a,b = T[i]
        for j in range(1,len(T)):
            a2,b2 = T[j]
            if (a<=a2<=b2<=b or a2<=a<=b<=b2) and not (a==a2 or b==b2):
                return False
    return True

#convex topology -> binary tree
def CT_to_BT_down(T):
    if len(T)==0:
        return ()
    k = 0
    while T[k][1]+1<len(T):
        k = T[k][1]+1
    T1 = [(T[i][0],min(T[i][1],k-1)) for i in range(k)]
    T2 = [(max(T[i][0],k)-k-1,T[i][1]-k-1) for i in range(k+1,len(T))]
    return (CT_to_BT_down(T1), CT_to_BT_down(T2))

def CT_to_BT_up(T):
    if len(T)==0:
        return ()
    k = len(T)-1
    while T[k][0]>0:
        k = T[k][0]-1
    T1 = [(T[i][0],min(T[i][1],k)) for i in range(k)]
    T2 = [(max(T[i][0],k+1)-k-1,T[i][1]-k-1) for i in range(k+1,len(T))]
    return (CT_to_BT_up(T1), CT_to_BT_up(T2))

def CT_to_SkInt(T):
    return CT_to_BT_up(T), CT_to_BT_down(T)


def TamariLattice(n, lattice = True):
    T = ()
    for k in range(n):
        T = (T,())
    LT,G = [Binary_Tree(T)],{}
    while len(LT)!=0:
        LT2 = set()
        for T in LT:
            R = T.upper_covers()
            G[T] = R
            for T2 in R:
                if T2 not in G:
                    LT2.add(T2)
        LT = LT2
    if lattice:
        return LatticePoset(G)
    return DiGraph(G)

def TamariDual(T):
    if T==():
        return ()
    return (TamariDual(T[1]),TamariDual(T[0]))

def ArcBinaryTrees(n):
    L = [[[]]]
    for k in range(1,n+1):
        l = []
        for i in range(k):
            for t in L[i]:
                for u in L[k-i-1]:
                    l.append([(1,k)]+t+[(a+i+1,b+i+1) for a,b in u])
        L.append(l)
    return L[-1]

def BT_to_ArcBT(T,size = False):
    if T==():
        l,s = [],0
    else:
        T1,s1 = BT_to_ArcBT(T[0],True)
        T2,s2 = BT_to_ArcBT(T[1],True)
        s = s1+s2+1
        l = [(1,s)]+T1+[(a+s1+1,b+s1+1) for a,b in T2]
    if size:
        return l,s
    return l
    

def ArcBT_is_interval(T,U):
    for a,b in T:
        for c,d in U:
            if a<c<=b+1<d+1:
                return False
    return True

def ArcBT_is_skew_interval(T,U):
    for a,b in T:
        for c,d in U:
            if a<c<=b<d:
                return False
    return True
    
#################### symmetric intervals, skew intervals ... ####################

def SymCT(n):
    return [T for T in ConvexTopologies(n) if all(T[k][0]==n-T[-k-1][1]-1 and T[k][1]==n-T[-k-1][0]-1 for k in range((n+1)//2))]

"""
1, 	0, 	1, 	0,  	2, 	0,  	5, 	0, 	14	self-dual binary trees
1,  	1, 	3, 	3,	11,	11,	45,	45,	197	self-dual faces
1,  	1, 	3, 	4,  	15, 	22, 	91  	140,	612	self-dual intervals (A369082)
1,  	2, 	5, 	11, 	32, 	74,	233,	555, 	1833	self-dual skew intervals
"""



def bracket_vector_aux(T):
    if T == ():
        return [], 0
    b1, n1 = bracket_vector_aux(T[0])
    b2, n2 = bracket_vector_aux(T[1])
    return b1 + [n2] + b2, n1 + n2 + 1

def to_bracket_vector(T):
    return tuple(bracket_vector_aux(T)[0])

def dual_bracket_vector_aux(T):
    if T == ():
        return [], 0
    b1, n1 = dual_bracket_vector_aux(T[0])
    b2, n2 = dual_bracket_vector_aux(T[1])
    return b1 + [n1] + b2, n1 + n2 + 1

def to_dual_bracket_vector(T):
    return tuple(dual_bracket_vector_aux(T)[0])

class Binary_Tree():
    def __init__(self, T):
        self.T = T
        self.bv = to_bracket_vector(T)
        self.dual_bv = to_dual_bracket_vector(T)
    def main(self):
        return self.T
    def __str__(self):
        return str(self.T)
    def __repr__(self):
        return str(self)
    def __hash__(self):
        return str(self).__hash__()
    def __eq__(self, other):
        if type(other)!=Binary_Tree:
            try:
                return self == Binary_Tree(other)
            except:
                return False
        return self.T == other.T
    def __le__(self,other):
        return all(i<=j for i,j in zip(self.bv, other.bv))
    def left_child(self):
        if self.T == ():
            raise Exception('T is an empty tree')
        return Binary_Tree(self.T[0])
    def right_child(self):
        if self.T == ():
            raise Exception('T is an empty tree')
        return Binary_Tree(self.T[1])
    def tree_graph(self):
        T = self.T
        if T == ():
            return [], 0
        G1, n1 = self.left_child().tree_graph()
        G2, n2 = self.right_child().tree_graph()
        n = n1 + n2 + 1
        G = [([n1/2, n/2], [n1/2, n/2]), ([n- n2/2, n/2], [n2/2, n/2])]
        G.extend(G1)
        G.extend([([x + n1 + 1 for x in X], Y) for X, Y in G2])
        return G, n
    def draw(self, a = None, size = 70,style = None):
        if a==None:
            fig, ax = plt.subplots(1,1, figsize=(5,5),subplot_kw={'aspect': 'equal'})
            ax.axis('off')
        else:
            ax = a
        G, n = self.tree_graph()
#        ax.plot([n,0],[0,0],color='black', linewidth=.5)
        for X, Y in G:
            ax.plot(X, Y, color='blue', linewidth=1)
        if a==None:
            plt.show()
    def lower_covers(self):
        T = self.T
        if T==():
            return []
        LR = []
        if T[1]!=():
            LR.append(Binary_Tree(((T[0], T[1][0]),T[1][1])))
        for R in Binary_Tree(T[0]).lower_covers():
            LR.append(Binary_Tree((R.T,T[1])))
        for R in Binary_Tree((T[1])).lower_covers():
            LR.append(Binary_Tree((T[0],R.T)))
        return LR
    def upper_covers(self):
        T = self.T
        if T==():
            return []
        LR = []
        if T[0]!=():
            LR.append(Binary_Tree((T[0][0],(T[0][1],T[1]))))
        for R in Binary_Tree(T[0]).upper_covers():
            LR.append(Binary_Tree((R.T,T[1])))
        for R in Binary_Tree((T[1])).upper_covers():
            LR.append(Binary_Tree((T[0],R.T)))
        return LR
    def popdown(self):
        bv = self.bv
        for t in self.lower_covers():
            bv = [min(i, j) for i, j in zip(bv, t.bv)]
        return to_binary_tree(bv)
    def popup(self):
        dual_bv = self.dual_bv
        for t in self.upper_covers():
            dual_bv = [min(i, j) for i, j in zip(dual_bv, t.dual_bv)]
        return to_binary_tree_dual(dual_bv)
    def canopy(self):
        T = self.T
        if T == ():
            return (None,)
        T1, T2 = T
        C1 = (1,) if T1 == () else Binary_Tree(T1).canopy()
        C2 = (0,) if T2 == () else Binary_Tree(T2).canopy()
        C = list(C1)
        C.extend(C2)
        return tuple(C)

def to_binary_tree(bv):
    if len(bv) == 0:
        return Binary_Tree(())
    for i in range(len(bv)):
        if bv[-i-1] == i:
            k = i
    if k == 0:
        return Binary_Tree((to_binary_tree(bv[:-k-1]).main(), ()))
    return Binary_Tree((to_binary_tree(bv[:-k-1]).main(), to_binary_tree(bv[-k:]).main()))

def to_binary_tree_dual(bv):
    if len(bv) == 0:
        return Binary_Tree(())
    for i in range(len(bv)):
        if bv[i] == i:
            k = i
    return Binary_Tree((to_binary_tree_dual(bv[:k]).main(), to_binary_tree_dual(bv[k+1:]).main()))


def random_binary_tree(n):
    if n==0:
        return Binary_Tree(())
    r = random.random()*catalan_number(n)
    k, s =0, 0
    while s<r:
        s += catalan_number(k)*catalan_number(n-k-1)
        k += 1
    return Binary_Tree((random_binary_tree(k-1).main(), random_binary_tree(n-k).main()))

def Draw_BT_pair(T1, T2):
    if all(i <= j for i, j in zip(T1.bv, T2.bv)):
        print('T1 <= T2')
    else:
        print('T1 <!= T2')
    print('T1', T1.bv)
    print('T2', T2.bv)
    def graphs(T):
        if T == ():
            return [], [], 0
        Gl1, Gr1, n1 = graphs(T[0])
        Gl2, Gr2, n2 = graphs(T[1])
        n = n1 + n2 + 1
        Gl, Gr = [([n1/2, n/2], [n1/2, n/2])], [([n- n2/2, n/2], [n2/2, n/2])]
        Gl.extend(Gl1)
        Gl.extend([([x + n1 + 1 for x in X], Y) for X, Y in Gl2])
        Gr.extend(Gr1)
        Gr.extend([([x + n1 + 1 for x in X], Y) for X, Y in Gr2])
        return Gl, Gr, n

    T1 = T1.main()
    T2 = T2.main()
    Gl1, Gr1, n1 = graphs(T1)
    Gl2, Gr2, n2 = graphs(T2)

    if n1 != n2:
        raise Exception('T1 and T2 must have the same size')

    fig, ax = plt.subplots(1,1, figsize=(5,5),subplot_kw={'aspect': 'equal'})
    ax.axis('off')
    for X, Y in Gl2:
        ax.plot(X, Y, color='red', linewidth=5, alpha = 0.5)
    for X, Y in Gr1:
        ax.plot(X, Y, color='blue', linewidth=5, alpha = 0.5)
    plt.show()


def canopy_preimages(n):
    L = TamariLattice(n)
    D = defaultdict(list)
    for T in L:
        D[T.canopy()].append(T)
    n_int, n_trv = 0, 0
    for c in D:
        L2 = LatticePoset(L.subposet(D[c]))
        n_int += len(L2.intervals_poset())
        n_trv += len(semidistributive_transvals(L2, False))
    return n_int, n_trv

# nombre de canopy transvals 1, 2, 7, 28, 125, 598, 3011, 15760
# à comparer avec 2,7,34,203,1394
def canopy_transvals(n):
    L = TamariLattice(n)
    return [(x, y) for x,y in semidistributive_transvals(L) if x.canopy() == y.canopy()]


############# all topologies ############

# to be more efficient (eliminate doublons)
# code a function to generate labelled posets
def Topologies(n):
    """
    return topologies on {1..n}, as pairs of minimal neighborhoods of k for 1<=k<=n
      and minimal neighborhoods of the dual topology
    """
    L = set()
    for k in range(1, n+1):
        for P in Posets(k):
            ideals = [P.order_ideal([i]) for i in range(k)]
            filters = [P.order_filter([i]) for i in range(k)]
            for S in SetPartitions(n, k):
                for p in Permutations(S):
                    T, T2 = n*[None], n*[None]
                    for i, I in enumerate(ideals):
                        t = Set([])
                        for j in I:
                            t = t.union(Set(p[j]))
                        for j in p[i]:
                            T[j-1] = t
                    for i, I in enumerate(filters):
                        t = Set([])
                        for j in I:
                            t = t.union(Set(p[j]))
                        for j in p[i]:
                            T2[j-1] = t
                    L.add((tuple(T), tuple(T2)))
    return L

def preposet_to_topology(G):
    T, T2 = [], []
    for x in range(1, len(G)+1):
        ideal, filter = [], []
        for y in G:
            if x == y or G.has_edge(x, y):
                ideal.append(y)
            if x == y or G.has_edge(y, x):
                filter.append(y)
        T.append(Set(ideal))
        T2.append(Set(filter))
    return tuple(T), tuple(T2)

def test_topologies(n, I):
    S = set()
    lP = Preorders(n)
    for G in lP:
        if is_WOQ_preposet(I, G):
            T, T2 = preposet_to_topology(G)
            for s in T:
                for s2 in T2:
                    S.add((s, s2))
    for s in Subsets(n):
        for s2 in Subsets(n):
            if len(s)>0 and len(s2)>0 and (s, s2) not in S:
                print((s,s2))
    for G in lP:
        T, T2 = preposet_to_topology(G)
        if all((s, s2) in S for s in T for s2 in T2) != is_WOQ_preposet(I, G):
            print('erreur', T, T2, S1, S2, G.edges())

# pas au point
def is_WO_topology(T):
    n = len(T)
    for k in range(1, n+1):
        t = T[k-1]
        for a in t:
            if a < k:
                for i in range(a+1, k):
                    if i not in t and a not in T[i-1]:
                        return False
            elif k < a:
                for i in range(k+1, a):
                    if i not in t and a not in T[i-1]:
                        return False
    return True
        
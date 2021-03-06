# base packages
import time
import pickle
from itertools import combinations, product, izip, count

# 3rd party
import numpy as np
from scipy import sparse

# gabe's tool
from timer import show_progress, get_time_str

# OLD -- works, but makes matrix dense. this wont work for
# genus>4. Useful for testing.
def dense_mod2rank(M, in_place=False):
    M = np.array(M, copy=not in_place)
    # fuck = sparse.lil_matrix(M)

    for index, column in enumerate(M.T):
        nonzero_rows = np.flatnonzero(column)
        # shit = fuck.T[index].tocsr()
        # shit.eliminate_zeros()
        # shit = shit.indices
        # assert (fuck[shit,index].A==1).all()
        # assert ((shit==nonzero_rows).all())

    if any(nonzero_rows):
        first_one, other_ones = ( nonzero_rows[0], 
                                  nonzero_rows[1:] )
        M[other_ones] = (M[other_ones]+M[first_one])%2
        M[first_one, index+1:] = 0

        # for row in other_ones:
        #   fuck[row] = fuck[row]+fuck[first_one]
        #   fuck.data[row] = [i%2 for i in fuck[row].data[0]]
        #   assert (np.array(fuck[row].data[0])<2).all()
        # fuck[first_one,index+1:]=0
        # if not ((M == fuck.toarray()).all()):
        #   import pdb; pdb.set_trace()
        return M.shape[1] - M.sum()
        
def mod2rank(M, in_place=False):
    # fuck = M.toarray()
    print "computing rank"
    for column in show_progress(xrange(M.shape[1])):
        # this conversion is linear in num_cols each time. not great but
        # not terrible... might be able to forgo completely actually by
        # using csr throughout...
        shit = M.T[column].tocsr()

        shit.eliminate_zeros() 
        # necessary to make sure nonzero_rows are what they are

        nonzero_rows = shit.indices
        if any(nonzero_rows):
            first_one, other_ones = ( nonzero_rows[0], 
                                      nonzero_rows[1:] )
            for row in other_ones:
                M[row] = M[row]+M[first_one]
                M.data[row] = [i%2 for i in M[row].data[0]]
                #assert (np.array(M[row].data[0])<2).all()
            # fuck[other_ones] = (fuck[other_ones]+fuck[first_one])%2
            # assert not (fuck == M.toarray()).all()
            M[first_one, column+1:] = 0
            # fuck[first_one,column+1:]=0

    # assert (fuck==M.toarray()).all()
    return M.shape[1] - M.sum()

class Space:
    def __init__(self,M):
        self.length = len(M[0])           # fails for empty input
        self.rank = 0
        self.basis = np.array([np.zeros(self.length)],dtype=int) #include an extra row of zeros to detect zero space
        self.pivots = np.array([self.length]) #add an extra pivot for the row of zeros
        for vector in M:
            self.includevector(vector)

    def includevector(self,vector):
        # if self.rank == 3 and any(vector[self.pivots[:self.rank]]): 
        #     import pdb; pdb.set_trace()
        newvector = (vector - sum(self.basis[i] * vector[self.pivots[i]]
                                  for i in xrange(self.rank) ))%2
        newpivot = pivot(newvector)
        if newpivot == -1:
            return
        for i in range(self.rank):
            self.basis[i] = (self.basis[i] + self.basis[i][newpivot]*newvector)%2
        insertposition = np.where(self.pivots >= newpivot)[0][0]
        self.basis = np.insert(self.basis,insertposition,newvector,axis=0)
        self.pivots = np.insert(self.pivots,insertposition,newpivot)
        self.rank += 1

    def copy(self):
        T = zerospace(self.length)
        T.length = self.length
        T.rank = self.rank
        T.basis = self.basis.copy()
        T.pivots = self.pivots
        return T
       
    def __contains__(self,vector):
        # guess = sum(self.basis[i]*vector[self.pivots[i]] 
        #                 for i in xrange(self.rank) )%2
        guess = (vector[self.pivots[:self.rank]] * 
                 self.basis[:self.rank].T).T.sum(0)%2
        # assert np.array_equal(guess_old,guess)
        return np.array_equal(vector,guess)
 
    def orthogonal_to(self,vector):
        return all(orthogonal(vector,x) for x in self.basis)

    def __gt__(self, space):
        return all(v in self for v in space.basis)      
 
    def __lt__(self, space):
        return space > self
 
    def __eq__(self, N):
        return np.array_equal(self.basis,N.basis)

    def __hash__(self):
        return hash(tuple(map(tuple,self.basis)))

def pivot(vector):
    for i,v in enumerate(vector):
        if v!=0:
            return i
    return -1
    
def orthogonal(x,y):
    return (x[::-1].dot(y) %2) == 0

def is_orthogonal(M,N):     #M,N must by np.arrays
    return all(orthogonal(*x) for x in product(M,N))

def is_isotropic(M):
    return all(orthogonal(*x) for x in combinations(M,2))

def span(S,vector):
    T = S.copy()
    T.includevector(vector)
    return T

def itervectors(length):
    return (np.array(m) for m in product(*( [range(2)]*length ))
            if pivot(m)!=-1)

def zerospace(k):
    return Space(np.array([[0]*k]))

def perp(S):
    if S == zerospace(S.length):
        return itervectors(S.length)
    elif S.length < 2*GENUS :
        # cut down onthe search space for lower dimensions by
        # restricting with pivot
        p = S.pivots[-2]+1
        tmp = (np.array([0]*p+list(v)) for v in itervectors(S.length-p))
        return (i for i in tmp if S.orthogonal_to(i))
    else:
        return (i for i in itervectors(2*GENUS) 
                if S.orthogonal_to(i) and not i in S)

def get_data():
    start = time.time()

    lagrangians = set([zerospace(2*GENUS)])
    containment= set()
  
    for r in range(0,GENUS):
        triangles = lagrangians
        lagrangians = set()
        print "computing dim",r
        for S in show_progress(triangles):
            #       print S.basis, map(list,perp(S))
            for v in perp(S):
                T = span(S,v)
                lagrangians.add(T)
                if r==GENUS-1:
                    containment.add((S,T))

    print "done with data", get_time_str(time.time()-start)
    return containment,lagrangians,triangles

def build_matrix(containment,lagrangians,triangles):
    start = time.time()
    t,l = count(),count()
    triangle2row   = dict(izip(triangles,t))
    lagrangian2col = dict(izip(lagrangians,l))
    shape = t.next(),l.next()
    print "triangles,lagrangians = ",shape

    num_ones = len(containment)
    rows,cols = np.array([(triangle2row[i],
                           lagrangian2col[j])
                          for i,j in containment]).reshape(num_ones,2).T
    matrix = sparse.coo_matrix((np.ones(num_ones),(rows,cols)), shape=shape,dtype=int)
    print "building matrix took", get_time_str(time.time()-start)

    with open("isotropic_subspaces_%s.pkl"%GENUS,'w') as f:
        pickle.dump(matrix,f)

    return matrix

def tests(matrix,containment):
    assert (matrix.sum(1)==3).all()
    assert (matrix.sum(0)==2**GENUS-1).all()

    x,y = next(iter(containment))
    assert x<y

GENUS=3

if __name__=="__main__":
    data = get_data()
    matrix = build_matrix(*data)
    result = mod2rank(matrix.tolil())
    print "result", result

    tests(matrix,data[0])
    if GENUS==3: assert result==15
    if GENUS==4: assert result==51

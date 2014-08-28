from itertools import combinations, product
import numpy as np
import csv
import time
from scipy import sparse
import pickle

def mod2rank(M, in_place=False):
  for index, column in enumerate(M.T):
    nonzero_rows = column.tocoo().col
    assert all(column.T[nonzero_rows] == 1)
    if any(nonzero_rows):
      first_one, other_ones = ( nonzero_rows[0], 
                                nonzero_rows[1:] )
      for row in other_ones:
        M[row] = M[row]+M[first_one]
        M[row].data = [[i%2 for i in M[row].data[0]]]
      M[first_one, index+1:] = 0
  return M.sum()

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
    for i in range(len(vector)):
        if vector[i]!=0:
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

def makevectors(length):
    return (np.array(m) for m in product(*( [range(2)]*length ))
            if pivot(m)!=-1)

def zerospace(k):
    return Space(np.array([[0]*k]))

def perp(S):
    if S == zerospace(S.length):
        return makevectors(S.length)
    else:
        p = S.pivots[-2]+1
        return (np.array([0]*p+list(v)) for v in makevectors(S.length-p)
                if S.orthogonal_to(np.array([0]*p+list(v))))

def big_perp(S):
    return (i for i in makevectors(2*GENUS) 
            if S.orthogonal_to(i) and not i in S)

def test_one(matrix):
    assert all(matrix.sum(1)==3)
def test_two(containment):
    x,y = next(iter(containment))
    assert x<y

GENUS=3

def get_data():
  isotropics = [ set([zerospace(2*GENUS)]) ]
  containment= set()

  start = time.time()
  
  for r in range(0,GENUS):
    print "computing dim",r
    isotropics.append(set())
    for S in isotropics[r]:
      #       print S.basis, map(list,perp(S))
      if r < GENUS-1:
        for v in perp(S):
          T = span(S,np.array(v))
          isotropics[r+1].add(T)
      else:
        for v in big_perp(S):
          T = span(S,v)
          isotropics[r+1].add(T)
          containment.add((S,T))

          
  lagrangians = list(isotropics[-1])
  triangles = list(isotropics[-2])
  with open("isotropic_subspaces_%s.pkl"%GENUS,'w') as f:
    pickle.dump({"containment":containment,
                 "lagrangians":lagrangians,
                 "triangles":triangles},f)
  print "done with data", time.time()-start
  return containment,lagrangians,triangles

def compute_rank(containment,lagrangians,triangles):
  start = time.time()
  l = len(lagrangians)
  t = len(triangles)

  # matrix = sparse.coo_matrix((t,l),dtype=int)
  
  lagrangian2index = dict(zip(lagrangians,xrange(l)))
  triangle2index   = dict(zip(triangles,xrange(t)))
  # for triangle,lagrangian in containment:
  #     matrix[triangle2index  [triangle],
  #            lagrangian2index[lagrangian]] = 1
  num_ones = len(containment)
  rows,cols = np.array([(triangle2index[i],
                         lagrangian2index[j])
                        for i,j in containment]).reshape(num_ones,2).T
  matrix = sparse.coo_matrix(
    (np.ones(num_ones),(rows,cols)), 
    shape=(t,l),dtype=int)
  

  # print matrix
  print "result", len(lagrangians) - mod2rank(matrix.tolil())
  print "that took", time.time()-start

  test_two(containment)
  test_one(matrix)

if __name__=="__main__":
  data = get_data()
  compute_rank(*data)

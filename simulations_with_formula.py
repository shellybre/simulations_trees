from numpy import random
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.special import rel_entr, kl_div
from scipy.stats import wasserstein_distance


#-------------------for tree :
class Node:
  def __init__(self):
    global count
    self.changes = 0
    self.name = ""
    self.children = []
    self.lengthS = ""
    self.length = 0
    self.parent = None
    self.genom = []
    self.distances = []
    self.genomLength = []

class Tree:
    def __init__(self):
        self.leaves = []
        self.root = None
    def getLeaves (self):
        for x in range(len(self.leaves)):
            print(self.leaves[x].name)
    def getRandomLeaf (self):
         if (len(self.leaves) > 0):
           x = random.randint(len(self.leaves))
           return x
         return 0

#function to calculate bin
#get=lam=length to check,and limits of the bin
def calculateWIthK(lam,first,second):
    x=0
    for k in range (first,second+1):
        p_nj = math.exp(-lam) # Probability of no jump
        p_j = 1 - p_nj  # Probability of a jump
        z=p_j + 4*(p_nj**2)*(1-p_nj)/(1+p_nj)
        if (k==1):
            p=(p_j+4*(p_nj**2)*(p_j**2))/z
        else:
            p=(4*(p_nj**(2*k))*(p_j**2))/z
        x+=p
    return x

#function to create new genome from existent one
#@get source_gen = get existent genome (parent genome), changed_gen=the genome to be changed, change=number of jumps
#return the new genome
def change_vector (source_gen, changed_gen, change):
    x = 0
    while x < change :
        y = random.randint(len(source_gen))
        z = random.randint(len(source_gen))
        element = changed_gen.pop(y)
        changed_gen.insert (z,element)
        x+=1
    return changed_gen

#recursive function to create recursive tree (void function)
# @get node, n number of leaves, scaleexpL=length of left branch, scaleexpR=length of right branch, sizeofGenom=length of genome
def create_recr_tree (tree, n, scaleexpL,scaleexpR ,sizeofGenom):                                 #create tree
    if (len(tree.leaves) == n):
        return
    if (len(tree.leaves) == 0):
        root = Node()
        size = sizeofGenom
        for x in range(size):
            root.genom.append(x)
        tree.root = root
        tree.leaves = [root]
        create_recr_tree (tree, n, scaleexpL,scaleexpR, sizeofGenom)
    else:
        numLeaf = tree.getRandomLeaf()
        chosenLeaf = tree.leaves[numLeaf]
        left = Node()
        left.length = scaleexpL
        left.parent = chosenLeaf
        change = random.poisson(lam=sizeofGenom*left.length, size=1)
        left.genom= change_vector (left.parent.genom, left.parent.genom.copy(), change)
        right = Node()
        right.length = scaleexpR
        right.parent = chosenLeaf
        change = random.poisson(lam=sizeofGenom*right.length, size=1)
        right.genom = change_vector (right.parent.genom, right.parent.genom.copy(), change)
        chosenLeaf.children = [left,right]
        tree.leaves.pop(numLeaf)
        tree.leaves.append(left)
        tree.leaves.append(right)
        create_recr_tree (tree, n, scaleexpL,scaleexpR ,sizeofGenom)

#find relevant nodes for parents
#get the tree and the array to save the nodes
def arrParents (root, arr):
  if (root != None):
    if (len(root.children) > 0):
      arr.append (root)
      arrParents(root.children[0],arr)
      arrParents(root.children[1], arr)
      return arr

#find relevant nodes for grandparents
#get the tree and the array to save the nodes
def arrGrand (root, arr):
  if (root != None):
    if (len(root.children) > 0 and (len(root.children[0].children) > 0 or len(root.children[1].children) > 0)):
      arr.append (root)
      arrGrand(root.children[0],arr)
      arrGrand(root.children[1], arr)
      return arr

#function to find distance on tree
#get two genomes
def findLCA (first, second):                                    #find LCA
    root = []
    ancestor = first.parent
    while (ancestor != None):
        root.append (ancestor)
        ancestor = ancestor.parent
    length = 0
    while second != None and second not in root:
        length += second.length
        lca = second
        second = second.parent
        if (second == first):
            return length
    while first != None and first != second:
        length += first.length
        first = first.parent
    return length

#function to create bins between two genomes
#get the island distribution
def CreateArrReal (sample):
   size = np.sum(sample)
   xarr = []
   yarr = []
   x = 0
   while (sample[x] == 0 and x < len(sample)):
      x+=1
   y = x
   lim = x
   while (y < len(sample)):
     if (sample[y]> 0):
       lim = y
     y+=1
   while (x <= lim):
     count = sample[x]
     first = x
     second = x
     x+=1
     while (x <= lim and count < size/10):
       if (sample[x] + count > size/10):
         remainder = size/10 - count
         if remainder >= sample[x]/2 or x + 1 == len(sample) or count == 0:
           count+= sample[x]
           second = x
           x+=1
           break
         else:
           second = x -1
           break
       else:
         count+=sample[x]
         second = x
         x+=1
     yarr.append(count / size)
     if (first == second):
       xarr.append(str(first))
     else:
      xarr.append (str(first) + " to " + str(second))
   return xarr, yarr

#function to find islands
#get two genomes
def find_islands(first, second):                        #find islands for two leaves
    x = 0
    root = [0]*(len(first) +1)
    while x < len(first):
        count = 0;
        val = first[x]
        y = second.index(val)
        while (y < len(second) and x < len(first) and first[x] == second[y]):
            x+=1
            y+=1
            count+=1
        root[count]+=1
        count = 0
    return root

#create bins for simulation
#get = length to check, realarr=given bins
def CreateArrSimulation (lam,realArr):
  #size of bin
  arrSim = []
  arrSimX=[]
  x=-1
  #name of bin
  for x in range (len (realArr)-1):
    liml= int(realArr[x].split()[0])
    limh = int(realArr[x].split()[len(realArr[x].split())-1])
    if (liml == limh):
      arrSimX.append(str(liml))
    else:
      arrSimX.append(str(liml) + " to " + str(limh))
    #liml = first island in bin, limh = last island in bin
    if(lam==0):
    	p=0
    else:
        p=calculateWIthK(lam,liml,limh)
    arrSim.append(p)
  liml= int(realArr[x+1].split()[0])
  arrSim.append(1-sum(arrSim))
  arrSimX.append(str(liml) + " to inf")
  return arrSimX,arrSim

#function to get the closest simulation for two given genomes
#get two genomes and the array of sizes of simulations
def mainfun1 (first, second, temp):       #Q = model, P = real
  length = findLCA(first, second)
  sizes = temp.copy()
  #create bins from the two genomes
  realArr, yreal = CreateArrReal(find_islands(first.genom, second.genom))
  results = []
  for x in range (len(sizes)):
    #create bins for simulation
    arrSimX, arrSim = CreateArrSimulation(sizes[x], realArr)
    if (arrSimX == None):
        print("None")
    res = wasserstein_distance(yreal, arrSim)
    results.append(res)
  Z = [x for _,x in sorted(zip(results,sizes))]
  return Z[0], length

#function to find the closest simulation in tree
#get=scaleL=length of left branch, scaleR=length of right branch,sizes=simulations to compare to,relationType=parent-child/brothers/grandparent-grandson
def findClosestSimulation (scaleL,scaleR, sizes, relationType):
    #create tree with genomes of length 5000
    tree = Tree()
    create_recr_tree(tree, 20, scaleL, scaleR, 5000)
    arr = []
    # get relevant nodes to check
    if relationType=='grandparent':
        arr = arrGrand(tree.root, arr)
    else:
        arr = arrParents(tree.root, arr)
    y = 0
    results=[]
    #get the closest simulation for the wanted relationship
    for y in range (len(arr)):
        if relationType=='parent':
            res, real =mainfun1(arr[y], arr[y].children[0], sizes)
            results.append(res)
        elif  relationType=='brothers':
            res, real = mainfun1(arr[y].children[0], arr[y].children[1], sizes)
            results.append(res)
        else:
            if (len(arr[y].children[0].children) > 0):
                res, real  = mainfun1(arr[y], arr[y].children[0].children[1], sizes)
                results.append(res)
            if (len(arr[y].children[1].children) > 0):
                res, real = mainfun1(arr[y], arr[y].children[1].children[0], sizes)
                results.append(res)
    return results

#main function to find the closest simulation to the real distance
#get=amount=number of simulations, scaleL=length of left branch, scaleR=length of right branch,sizes=simulations to compare to,relationType=parent-child/brothers/grandparent-grandson
#returns percentage
def main (amount, scaleL,scaleR, sizes,relationType):
    count_sizes=[0]*(len(sizes))
    for x in range (amount):
        results=findClosestSimulation(scaleL,scaleR, sizes,relationType)
        for y in range(len(results)):
            index = sizes.index(results[y])
            count_sizes[index]+=1
    print(count_sizes)
    percentage=[x / sum(count_sizes) * 100 for x in count_sizes]
    print(percentage)
    return percentage



if __name__ == "__main__":
    sizes = [0.005, 0.01, 0.03, 0.05, 0.07, 0.1, 0.12, 0.14, 0.16, 0.2, 0.22, 0.24, 0.26, 0.3]
    relationType = ['parent', 'grandparent', 'brothers']
    l=0.05
    r=0.05
    amount=10
    percentage1=main(amount, l, r, sizes, relationType[0])
    percentage2=main(amount, l, r, sizes, relationType[1])
    percentage3=main(amount, l, r, sizes, relationType[2])
    plt.bar(sizes,percentage1 , width=0.01, color='blue', edgecolor='black')
    plt.show()
    plt.bar(sizes, percentage2, width=0.01, color='blue', edgecolor='black')
    plt.show()
    plt.bar(sizes, percentage3, width=0.01, color='blue', edgecolor='black')
    plt.show()


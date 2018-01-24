import numpy as np
import scipy as sp
import math as ma
import sys
import csv

def Dijkst(ist,isp,wei):
    # Dijkstra algorithm for shortest path in a graph
    #    ist: index of starting node
    #    isp: index of stopping node
    #    wei: weight matrix

    # exception handling (start = stop)
    if (ist == isp):
        shpath = [ist]
        return shpath

    # initialization
    N         =  len(wei)
    Inf       =  sys.maxint
    UnVisited =  np.ones(N,int)
    cost      =  np.ones(N)*1.e6
    par       = -np.ones(N,int)*Inf

    # set the source point and get its (unvisited) neighbors
    jj            = ist
    cost[jj]      = 0
    UnVisited[jj] = 0
    tmp           = UnVisited*wei[jj,:]
    ineigh        = np.array(tmp.nonzero()).flatten()
    L             = np.array(UnVisited.nonzero()).flatten().size

    # start Dijkstra algorithm
    while (L != 0):
        # step 1: update cost of unvisited neighbors,
        #         compare and (maybe) update
        for k in ineigh:
            newcost = cost[jj] + wei[jj,k]
            if ( newcost < cost[k] ):
                cost[k] = newcost
                par[k]  = jj

        # step 2: determine minimum-cost point among UnVisited
        #         vertices and make this point the new point
        icnsdr     = np.array(UnVisited.nonzero()).flatten()
        cmin,icmin = cost[icnsdr].min(0),cost[icnsdr].argmin(0)
        jj         = icnsdr[icmin]

        # step 3: update "visited"-status and determine neighbors of new point
        UnVisited[jj] = 0
        tmp           = UnVisited*wei[jj,:]
        ineigh        = np.array(tmp.nonzero()).flatten()
        L             = np.array(UnVisited.nonzero()).flatten().size

    # determine the shortest path
    shpath = [isp]
    while par[isp] != ist:
        shpath.append(par[isp])
        isp = par[isp]
    shpath.append(ist)

    return shpath[::-1]

def calcWei(RX,RY,RA,RB,RV):
    # calculate the weight matrix between the points

    n    = len(RX)
    wei = np.zeros((n,n),dtype=float)
    m    = len(RA)
    for i in range(m):
        xa = RX[RA[i]-1]
        ya = RY[RA[i]-1]
        xb = RX[RB[i]-1]
        yb = RY[RB[i]-1]
        dd = ma.sqrt((xb-xa)**2 + (yb-ya)**2)
        tt = dd/RV[i]
        wei[RA[i]-1,RB[i]-1] = tt
    return wei
    
def traffic(ist,isp,weio,it,k):
    
    n = len(weio)
    
    # to start with, there are no cars and
    # no edges have been used
    wei  = np.matrix.copy(weio)
    used = np.matrix.copy(weio)
    
    # we create a matrix C with 58 rows, 
    # the number of nodes, and it columns,
    # the number of iterations
    C        = np.zeros((58,it))
    C[ist,0] = 20
    
    # loop over time steps
    for t in range(1,it):
        
        # update the weights matrix wei 
        for i in range(n):
            for j in range(n):
                if (weio[i,j] != 0):
                    if (C[i,t-1]==C[j,t-1]==0):
                        break
                    else:
                        wei[i,j] = weio[i,j] + k*(C[i,t-1]+C[j,t-1])/2
                    
        # look at where cars were in the previous timestep
        L = np.array(C[:,t-1].nonzero()).flatten()
        
        for c in L:
            if (c==isp):
                # update the cars vector for the node 52
                # by making 40% of the cars leave the network
                C[c,t] += C[c,t-1] - ma.floor(0.3*C[c,t-1])
            else:
                # compute shortest path to Coliseum for cars at each node c
                shpath = Dijkst(c,isp,wei)
                # d is the first step to destination for cars in node c
                d = shpath[1]
            
                # update the cars vector by sending 70%
                # of cars in node c to node d 
                C[c,t] += C[c,t-1] - ma.floor(0.7*C[c,t-1])
                C[d,t] += ma.floor(0.7*C[c,t-1])
                
                # mark the visited nodes with a zero in the used matrix
                used[c,d] = 0
        
        # inject 20 cars at node 13 for first 180 iterations
        if (t<180):
            C[ist,t] += 20
        
    return C,used

if __name__ == '__main__':
    
    print 'M3SC Class Project 1'
    print 'Mathilde Duverger'
    
    # PROJECT 1 (path through Rome)
    RomeX = np.empty(0,dtype=float)
    RomeY = np.empty(0,dtype=float)
    with open('RomeVertices','r') as file:
        AAA = csv.reader(file)
        for row in AAA:
            RomeX = np.concatenate((RomeX,[float(row[1])]))
            RomeY = np.concatenate((RomeY,[float(row[2])]))
    file.close()

    RomeA = np.empty(0,dtype=int)
    RomeB = np.empty(0,dtype=int)
    RomeV = np.empty(0,dtype=float)
    with open('RomeEdgesProject','r') as file:
        AAA = csv.reader(file)
        for row in AAA:
            RomeA = np.concatenate((RomeA,[int(row[0])]))
            RomeB = np.concatenate((RomeB,[int(row[1])]))
            RomeV = np.concatenate((RomeV,[float(row[2])]))
    file.close()
    
    # we compute the original weights matrix weio
    weio = calcWei(RomeX,RomeY,RomeA,RomeB,RomeV)
    
    ist = 12 # Saint Paul's
    isp = 51 # Coliseum
    it = 200 # number of iterations
    k = 0.01 # parameter for the weights
    
    # run the traffic function
    C,used = traffic(ist,isp,weio,it,k)
    
    # compute the maximum load for each node  
    maxnode =  np.amax(C,axis=1)  
    print 'The maximum load for each node over the 200 iterations is'
    print maxnode
    
    
    # find the five most congested nodes
    maxfive  = (-maxnode).argsort()[:5]
    maxfive += 1
    print 'The five nodes with the highest maximums are in order'
    print maxfive
    
    avenode = np.mean(C,axis=1)
    maxavefive  = (-avenode).argsort()[:5]
    maxavefive += 1
    print 'The five nodes with maximum averages are in order'
    print maxavefive
    
    # find all unvisited nodes
    print 'The nodes',
    for i in range(58):
        if maxnode[i]==0:
            print i+1,',',
    print 'remain unvisited after 200 iterations.'  
    
    # find all unused edges
    unused1 = np.array(used.nonzero()).flatten()
    nbun    = len(unused1)
    unused  = np.zeros((nbun/2,2))
    
    for i in range(nbun/2):
        unused[i,0] = unused1[i] +1
        unused[i,1] = unused1[i+nbun/2] +1
    
    print 'The following edges remain unused after 200 iterations:' 
    print unused  

    
    # flow pattern for parameter k=0
    print 'Flow pattern for parameter k=0'
    Ck,usedk = traffic(ist,isp,weio,it,0)
    vn       = np.array(np.nonzero(np.amax(Ck,axis=1))).flatten()
    print 'The only visited nodes are'
    print vn+1

    print 'The normal Dijstra algorithm gives the fastest route between nodes 13 and 52' 
    shpath = Dijkst(ist,isp,weio)
    print 'Dijkstra:',ist+1,' -> ',isp+1,' is ',np.array(shpath)+1
    
    print 'An accident occurs at node 30 which blocks any route to or from node 30.'
    weio[29,:] = 0
    weio[:,29] = 0
    
    C30,used30 = traffic(ist,isp,weio,it,k)
    
    # compute the maximum load for each node  
    maxnode30 =  np.amax(C30,axis=1)  
    print 'The maximum load for each node over the 200 iterations is'
    print maxnode30
    
    # find the five most congested nodes
    maxfive30  = (-maxnode30).argsort()[:5]
    maxfive30 += 1
    print 'The five nodes with the highest maximums are in order'
    print maxfive30
    
    avenode30     = np.mean(C30,axis=1)
    maxavefive30  = (-avenode30).argsort()[:5]
    maxavefive30 += 1
    print 'The five nodes with maximum averages are in order'
    print maxavefive30
    
    # find the nodes with biggest decrease and increase in peak
    dif    = np.delete(maxnode-maxnode30,[29])
    difmax = np.argmax(dif)
    difmin = np.argmin(dif)
    print 'The node which has decreased the most in peak value is',
    if difmax>28 :
        print 'node', difmax +2, '.'
    else:
        print 'node', difmax +1, '.'
    print 'The node which has increased the most in peak value is',
    if difmin>28 :
        print 'node', difmin +2, '.'
    else:
        print 'node', difmin +1, '.'

    
   

    
      
        
        
    
  

    

  

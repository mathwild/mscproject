import numpy as np
import scipy as sp
import sys

dur = np.array([41,51,50,36,38,45,21,32,32,49,30,19,26])
dep = np.array([[0,1],
               [0,7],
               [0,10],
               [1,4],
               [1,12],
               [2,3],
               [5,7],
               [6,5],
               [6,9],
               [9,11],
               [10,12]])

def BellmanFord(ist,isp,wei):
    #----------------------------------
    #  ist:    index of starting node
    #  isp:    index of stopping node
    #  wei:    adjacency matrix (V x V)
    #
    #  shpath: shortest path
    #----------------------------------

    V = wei.shape[1]

    # step 1: initialization
    Inf    = sys.maxint
    d      = np.ones((V),float)*np.inf
    p      = np.zeros((V),int)*Inf
    d[ist] = 0

    # step 2: iterative relaxation
    for i in range(0,V-1):
        for u in range(0,V):
            for v in range(0,V):
                w = wei[u,v]
                if (w != 0):
                    if (d[u]+w < d[v]):
                        d[v] = d[u] + w
                        p[v] = u

    # step 3: check for negative-weight cycles
    for u in range(0,V):
        for v in range(0,V):
            w = wei[u,v]
            if (w != 0):
                if (d[u]+w < d[v]):
                    print('graph contains a negative-weight cycle')

    # step 4: determine the shortest path
    shpath = [isp]
    while p[isp] != ist:
        shpath.append(p[isp])
        isp = p[isp]
    shpath.append(ist)

    return shpath[::-1]

def initweights(dur,dep):
    
    # get non dependent nodes and finish nodes
    depn = np.ones(len(dur))
    depf = np.ones(len(dur))
    for i in range(len(dep)):
        depf[dep[i,0]] = 0
        depn[dep[i,1]] = 0
    nondep = depn.nonzero()
    fin = depf.nonzero()
    
    # create adjacency matrix
    adj = np.zeros([len(dur)+2,len(dur)+2])
    
    # connect virtual start to non dependent nodes
    adj[0, [x+1 for x in nondep]] = 1
    
    # connect dependent nodes
    adj[dep[:, 0] + 1, dep[:, 1] + 1] = dur[dep[:,0]]
    
    # connect fin nodes to virtual finish
    adj[[x+1 for x in fin],len(dur)+1] = dur[fin]
    
    wei = -adj 
    return wei
    
def Treestrings(dur,dep):
    
    wei = initweights(dur,dep)
    
    # compute iteratively the longest strings of jobs
    strings = []
    index = [0, len(dur)+1]
    for i in range(len(dur)):
        shpath = BellmanFord(0,len(dur)+1,wei)
        sh = np.asarray(shpath)
        lpath = sh[1:len(sh)-1] - 1
        strings.append(lpath)
        s = lpath[len(lpath)-1] 
        wei[s+1,:] = 0
        wei[:,s+1] = 0
        index.extend([s+1])
        up = np.where(wei.any(axis=1))[0]
        up = np.setdiff1d(up,index)
        for j in up:
            wei[j,len(dur)+1] = - dur[j-1]
    
    # get earliest start times and finish times
    Gantt = np.zeros([len(dur),2])
    for i in range(len(dur)):
        subs = strings[i]
        task = subs[len(subs)-1]
        for j in subs[0:len(subs)-1]:
            Gantt[task,0] += dur[j]
    Gantt[:,1] = Gantt[:,0] + dur
    
    
    # get substrings of jobs
    
    # copy of the strings list that will be modified
    copstrings = strings[:]
    
    # start strings are the jobs workers are going to start with
    ststrings = []
    # other strings are the jobs that will be distributed to workers
    ostrings = []
    
    # get the different strings
    while len(copstrings)!=0:
        w1 = copstrings[0]
        ststrings.append(w1)
        
        index2 = w1.flatten()
    
        dels = []
        
        for i in range(len(copstrings)):
            for j in range(len(w1)):
                if w1[j] in copstrings[i]:
                    dels.append(i)
                    tw = np.setdiff1d(copstrings[i],index2)
                    index2 = np.append(index2,tw)
                    if len(tw)!=0:
                        ostrings.append(tw)
        
        copstrings = np.delete(copstrings,dels,0)
    
    # all the strings and remaining substrings of jobs 
    substrings = np.append(ststrings,ostrings)
    
    # the first column is the earliest start time and the second
    # column is the duration of the substring
    subtimes = np.zeros([len(substrings),2])
    
    for i in range(len(substrings)):
        u = substrings[i]
        subtimes[i,0] = Gantt[u[0],0]
        subtimes[i,1] = Gantt[u[len(u)-1],1] - Gantt[u[0],0]

     
    return strings, Gantt, substrings, subtimes

if __name__ == '__main__':    
    
    print 'M3SC Class Project 2'
    print 'Mathilde Duverger'
        
    s,G,ss,st = Treestrings(dur,dep)
    print 'The longest strings are in order ', s
    print 'The earliest start and finish times for each job are ', G
    print 'I have decided to look for the longest substrings and the remaining substrings.'
    print 'The substrings of jobs are ', ss
    print 'The earliest start time and the duration of each substring are ', st
    print 'We can notice that the three last substrings of jobs would fit after the three above without lenghtening the time taken to finish all the jobs.'
    
    wjobs = ss[0:4]
    wtime = st[0:4]
    for i in [1,2,3]:
        u = wjobs[i]
        v = ss[i+3]
        wjobs[i] = [u,v]
        wtime[i,1] += st[i+3,1]

    
    print 'The workers should execute the following tasks in parallel ', wjobs
    print 'This is the duration for each parallel string of jobs ', wtime

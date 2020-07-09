#! /usr/bin/env python
import sys
import cvxopt
import re
from cvxopt import matrix, solvers

from os import listdir


def myblock(l):
    # numpy didn't work on MareNostrum!!!
    newm = []
    for blockrow in l:
        assert len(blockrow) >= 1
        nrows,ncol0 = blockrow[0].size
        for row in range(0,nrows):
            r = []
            for m in blockrow:
                r.extend( list(m[row,:]))
            newm.append(r)
    return matrix(newm).trans()

def zeros(rows,cols):
    return matrix(0, (rows,cols))


def read_current_alloc():
    fnames = listdir('.balance')
    max_instance = 0
    max_node = 0
    ranks = {}
    allocs = {}
    loads = {}
    for fname in fnames:
        m = re.match(r'load-([0-9*])-([0-9]*)', fname)
        if m:
            instance = int(m.group(1))
            node = int(m.group(2))
            max_instance = max(instance, max_instance)
            max_node = max(node, max_node)
            f = open('.balance/' + fname)
            l = f.readline().strip().split(' ')
            print l
            ranks[ (instance,node) ] = int(l[0])
            allocs[ (instance,node) ] = float(l[1])
            loads[ (instance,node) ] = float(l[2])
            f.close()
        
    return max_instance+1, max_node+1, ranks, allocs, loads
    
def write_new_alloc(ni, nn, B, opt_allocs):
    for instance in range(0,ni):
        for node in range(0,nn):
            if B[instance,node] != 0:
                f = open('.balance/alloc-%d-%d' % (instance,node), 'w')
                print >>f, opt_allocs[(instance,node)]
                f.close()

    
    

def printout(ni,nn,ranks,allocs,L):
    print '           ' + ' ' * 5 * nn + 'NUM CORES ' + ' ' * 5 * nn + 'Total_cores   Load   Load/Total_cores'
    print '           ' + ' ' * 5 * nn + '   node   ' + ' ' * 5 * nn + '    '
    print '           ' + ('%10d' * nn) % tuple(range(0,nn))

    for instance in range(0, ni):
        print 'Instance %2d' % instance,
        for node in range(0, nn):
            print '%10.2f' % allocs[(instance,node)],
        load = L[instance,0]
        total_c = 0
        for node in range(0, nn):
            total_c += allocs[(instance,node)]
        print '%10.2f' % total_c,
        print '%10.2f' % load,
        print '%10.2f' % (load/total_c)


    used = []
    for node in range(0,nn):
        u = 0
        for instance in range(0,ni):
            u += allocs[(instance,node)]
        used.append(u)
    print '           ' + ('%10s' % '---') * nn
    print '           ' + ('%10.2f' * nn) % tuple(used)

def optimize(ni, nn, L, B, C):
    # Minimise      -t                                        (i.e. maximise worst-case cores per load)
    #
    #                 N
    # subject to    Sum  a               <= C     all i       (available cores)
    #               j=1   ij                 i
    #
    #                   N
    #               - Sum  a     +  L  t   <= 0     all j       (Sum cores >= load * t)
    #                 i=1   ij       j

    # Variables are t, node0instance0, node0instance1, ... , node1instance0, ...

    # Objective function to minimise: t
    c = matrix( [[-1.0]] + [[0.0]] * ni * nn).trans()

    # LHS of constraints: one per node and one per instance
    Ali = []
    for node in range(0,nn):
        # Sum up the instances on this node
        row = [0] * ni * node + [1] * ni + [0] * ni * (nn-node-1)
        Ali.append(row)
    Aln = []
    for instance in range(0,ni):
        # Sum up the nodes on this instance
        row = [0.0] * nn * ni
        for node in range(0,nn):
            if B[instance,node] != 0:
                row[ node*ni + instance] = -1.0
        Aln.append(row)
    Alp = []
    for node in range(0,nn):
        for instance in range(0,ni):
            row = [0.0] * nn * ni
            row[node*ni + instance] = -1.0
            Alp.append(row)

    A = myblock( [ [zeros(nn,1),    matrix(Ali).trans()],
                   [L,              matrix(Aln).trans()],
                   [zeros(nn*ni,1), matrix(Alp).trans()]])
    # print 'A', A.size
    # print(A)

    # RHS of constraints
    b = 1.0 * myblock( [ [C],
                       [zeros(ni,1)],
                       [zeros(nn*ni,1)]] )
    # print 'b', b.size
    # print (b)
    # print 'c', c.size
    # print(c)

    sol = solvers.lp( c, A, b)

    opt_allocs = {}

    for instance in range(0, ni):
        for node in range(0, nn):
            opt_allocs[(instance,node)] = sol['x'][1 + node*ni + instance]
    return opt_allocs


def make_integer(ni, nn, allocs, L, B, C):
    int_allocs = {}
    for node in range(0,nn):
        indices = [instance for instance in range(0,ni) if B[(instance,node)]>0 and int(allocs[(instance,node)]) != allocs[(instance,node)] ]

        ncores      = [allocs[(instance,node)] for instance in indices]
        ncores_int  = [int(allocs[(instance,node)]) for instance in indices]
        total_cores = [sum([allocs[(instance,n)] for n in range(0,nn) if B[(instance,n)]>0]) for instance in indices]
        # load_per_cores = [L[indices[j]] / total_cores[j] for j in range(0,len(total_cores))]
        # print 'node', node, 'indices', indices, 'cores', ncores, 'total_cores', total_cores, 'lpc', load_per_cores

        # After rebalancing, all instances should have the same load-per-core
        # Hence the slowdown is proportional to the fraction lost: (ncores - int(ncores)) / total_cores
        frac_lost_and_j = [(1.0 * (ncores[j] - int(ncores[j])) / total_cores[j],j) for j in range(0,len(total_cores))]
        # print 'node', node, 'indices', indices, 'cores', ncores, 'total_cores', total_cores, 'frac_lost', frac_lost_and_j
        # print 'ncores_int', ncores_int

        extra_cores = C[instance] - sum(ncores_int)
        # print 'extra_cores', extra_cores
        for c in range(0,extra_cores):
            frac,j = frac_lost_and_j[c]
            ncores_int[j] = ncores_int[j] + 1

        for instance in range(0,ni):
            int_allocs[(instance,node)] = 0
        for j,instance in enumerate(indices):
            int_allocs[(instance,node)] = ncores_int[j]
    return int_allocs






ni, nn, ranks, allocs, loads = read_current_alloc()
# print 'ni =', ni
# print 'nn =', nn
# print 'allocs =', allocs
# print 'loads =', loads
Ll = []
Brows = []
for instance in range(0,ni):
    load = 0
    Brow = []
    for node in range(0,nn):
        if (instance,node) in loads:
            load += 1.0 * loads[(instance,node)]
            Brow.append(1.0)
        else:
            loads[(instance,node)] = 0.0
            allocs[(instance,node)] = 0.0
            Brow.append(0.0)
    Ll.append(load)
    Brows.append(Brow)

# Load vector
L = matrix(Ll)
# print(L)

# Topology matrix
B = matrix(Brows).trans()
# print(B)

# Available cores vector
C = matrix( [[48]] * nn).trans()
# print(C)

print 'Current allocation'
printout(ni,nn,ranks,allocs,L)


#print(sol['x'])
opt_allocs = optimize(ni, nn, L, B, C)
print 'Optimized allocation'
printout(ni,nn,ranks,opt_allocs,L)

integer_allocs = make_integer(ni, nn, opt_allocs, L, B, C)
print 'Integer allocation'
printout(ni,nn,ranks,integer_allocs,L)

write_new_alloc(ni, nn, B, integer_allocs)



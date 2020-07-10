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
            ranks[ (instance,node) ] = int(l[0])
            allocs[ (instance,node) ] = float(l[1])
            loads[ (instance,node) ] = float(l[2])
            f.close()
        
    return max_instance+1, max_node+1, ranks, allocs, loads
    
def write_new_alloc(ni, nn, B, opt_allocs):
    for instance in range(0,ni):
        f = open('.balance/alloc-%d' % instance, 'w')
        for node in range(0,nn):
            if B[instance,node] != 0:
                print >>f, opt_allocs[(instance,node)]
        f.close()

    
    

def printout(ni,nn,ranks,allocs,L):
    print '           ' + ' ' * 5 * nn + 'NUM CORES ' + ' ' * 5 * nn + 'Total_cores   Load   Load/Total_cores'
    print '           ' + ' ' * 5 * nn + '   node   ' + ' ' * 5 * nn + '    '
    print '           ' + ('%10d' * nn) % tuple(range(0,nn))

    for instance in range(0, ni):
        print 'Instance %2d' % instance,
        for node in range(0, nn):
            if (instance,node) in allocs:
                print '%10.2f' % allocs[(instance,node)],
            else:
                print '%10s' % '-',
        load = L[instance,0]
        total_c = 0
        for node in range(0, nn):
            total_c += allocs.get((instance,node),0)
        print '%10.2f' % total_c,
        print '%10.2f' % load,
        print '%10.3f' % (load/total_c)


    used = []
    for node in range(0,nn):
        u = 0
        for instance in range(0,ni):
            u += allocs.get((instance,node),0)
        used.append(u)
    print '           ' + ('%10s' % '---') * nn
    print '           ' + ('%10.2f' * nn) % tuple(used)

def optimize(ni, nn, ranks, L, B, C):
    # Minimise      -t                                        (i.e. maximise worst-case cores per load)
    #
    #                   N
    # subject to      Sum  a               <= C     all i       (available cores)           [1]
    #                 j=1   ij                 i
    #
    #                   N
    #               - Sum  a     +  L  t   <= 0     all j       (Sum cores >= load * t)     [2]
    #                 i=1   ij       j
    #
    #                      a               >= 0     all i,j                                 [3a]
    #                       ij
    #
    #                      a               >= 1     when j is the master node (as main cannot be migrated) [3b]

    # Variables are t, a_ij  (but only when B_ij = 1; otherwise will get bogus values for the other a_ij)

    #  For which (instance,node) do we need an a_ij
    entries = [(instance,node) for instance in range(0,ni) for node in range(0,nn) if B[instance,node] != 0]
    num_aij = len(entries)
    num_variables = 1 + num_aij   # including t

    # Objective function to minimise: t
    c = matrix( [[-1.0]] + [[0.0]] * num_aij).trans()

    # LHS of constraints: one per node and one per instance
    Ai = []
    bi = []
    # Constraints [1]
    for node in range(0,nn):
        row = [0.0] #  coefficient of t is zero in [3]
        for (instance,n) in entries:
            if n == node:
                row.append(1.0)   # This variable has coefficient 1 in [1]
            else:
                row.append(0)   # This variable has coefficient 0 in [1]
        Ai.append(row)          # LHS of [1]
        bi.append(C[node])      # RHS of [1]
    # Constraints [2]
    for instance in range(0,ni):
        # Sum up the nodes in this instance
        row = [L[instance]]   # Coefficient of t is Lj (j=instance)
        for (i,node) in entries:
            if i == instance:                 # Variable affects this instance
                row.append(-1.0)
            else:
                row.append(0.0)
        Ai.append(row)         # LHS of [2]
        bi.append(0.0)         # RHS of [2]
    # Constraints [3]
    for k,(instance,node) in enumerate(entries):
        row = [0.0] * num_variables
        row[k+1] = -1.0
        Ai.append(row)        # LHS of [3]
        if ranks[(instance,node)] == 0:
            bi.append(-1.0)          # RHS of [3b]
        else:
            bi.append(0.0)           # RHS of [3a]

    A = matrix(Ai).trans()
    # print 'A', A.size
    # print(A)
    b = matrix(bi)
    # print 'b', b.size
    # print(b)

    sol = solvers.lp( c, A, b)

    opt_allocs = {}
    for k,(instance,node) in enumerate(entries):
            opt_allocs[(instance,node)] = sol['x'][1 + k]
    return opt_allocs


def make_integer(ni, nn, allocs, L, B, C):
    int_allocs = {}
    for node in range(0,nn):
        indices = [instance for instance in range(0,ni) if B[(instance,node)]>0 and allocs[(instance,node)] > 0 ]

        if len(indices) > 0:
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
            frac_lost_and_j.sort()
            frac_lost_and_j.reverse()

            extra_cores = C[instance] - sum(ncores_int)
            # print 'extra_cores', extra_cores
            for c in range(0,extra_cores):
                # Sometimes it is not necessary to use all cores: in this case we just keep going
                # filling up the available cores anyway
                frac,j = frac_lost_and_j[c % len(frac_lost_and_j)]
                ncores_int[j] = ncores_int[j] + 1

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
            Brow.append(0.0)
    Ll.append(load)
    Brows.append(Brow)

# Load vector
L = matrix(Ll)

# Topology matrix
B = matrix(Brows).trans()

# Available cores vector
C = matrix( [[48]] * nn).trans()
print 'Current allocation'
printout(ni,nn,ranks,allocs,L)

if max(L) == 0.0:
    # Currently no work!!!
    opt_allocs = allocs
    print 'No work!'
else:
    opt_allocs = optimize(ni, nn, ranks, L, B, C)


print 'Optimized allocation'
printout(ni,nn,ranks,opt_allocs,L)

integer_allocs = make_integer(ni, nn, opt_allocs, L, B, C)
print 'Integer allocation'
printout(ni,nn,ranks,integer_allocs,L)

write_new_alloc(ni, nn, B, integer_allocs)



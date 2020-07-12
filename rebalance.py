#! /usr/bin/env python
import sys
import cvxopt
import re
import time
import getopt
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
    f = open('.map', 'r')
    line = f.readline()
    nn = line.strip().split(' ')
    ranktonode = [int(n) for n in nn]

    fnames = listdir('.balance')
    max_group = 0
    max_node = 0
    ranks = {}
    allocs = {}
    loads = {}
    for fname in fnames:
        m = re.match(r'load-([0-9*])-([0-9]*)', fname)
        if m:
            group = int(m.group(1))
            rank = int(m.group(2))
            node = ranktonode[rank]
            max_group = max(group, max_group)
            max_node = max(node, max_node)
            f = open('.balance/' + fname)
            l = f.readline().strip().split(' ')
            ranks[ (group,int(l[0])) ] = node
            allocs[ (group,node) ] = float(l[1])
            loads[ (group,node) ] = float(l[2])
            f.close()
        
    return max_group+1, max_node+1, ranks, allocs, loads
    
def write_new_alloc(ni, nn, ranks, B, opt_allocs):
    for group in range(0,ni):
        f = open('.balance/alloc-%d' % group, 'w')
        for rank in range(0,nn):
            if (group,rank) in ranks.keys():
                node = ranks[(group,rank)]
                print >>f, opt_allocs[(group,node)]
        f.close()

    
    

def printout(ni,nn,ranks,allocs,L):
    print '           ' + ' ' * 5 * nn + 'NUM CORES ' + ' ' * 5 * nn + 'Total_cores   Load   Load/Total_cores'
    print '           ' + ' ' * 5 * nn + '   node   ' + ' ' * 5 * nn + '    '
    print '           ' + ('%10d' * nn) % tuple(range(0,nn))

    for group in range(0, ni):
        print 'Instance %2d' % group,
        for node in range(0, nn):
            if (group,node) in allocs:
                print '%10.2f' % allocs[(group,node)],
            else:
                print '%10s' % '-',
        load = L[group,0]
        total_c = 0
        for node in range(0, nn):
            total_c += allocs.get((group,node),0)
        print '%10.2f' % total_c,
        print '%10.2f' % load,
        if total_c > 0:
            print '%10.3f' % (load/total_c)
        else:
            print '%10s' % 'inf'


    used = []
    for node in range(0,nn):
        u = 0
        for group in range(0,ni):
            u += allocs.get((group,node),0)
        used.append(u)
    print '           ' + ('%10s' % '---') * nn
    print '           ' + ('%10.2f' * nn) % tuple(used)

def optimize(ni, nn, ranks, L, B, C, policy):
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

    #  For which (group,node) do we need an a_ij

    # This depends on the policy
    if policy == 'optimized':
        # Optimized: have all a_ij
        entries = [(group,node) for group in range(0,ni)  
                                    for node in range(0,nn) 
                                        if B[group,node] != 0]
    elif policy == 'master':
        # Only have masters
        entries = [(group,node) for group in range(0,ni) 
                                    for node in range(0,nn) 
                                        if ranks[(group,0)] == node ]
    elif policy == 'slaves':
        # Only have slaves 
        entries = [(group,node) for group in range(0,ni)  
                                    for node in range(0,nn) 
                                        if B[group,node] != 0 and ranks[(group,0)] != node ]

    elif policy == 'slave1':
        # Only have slave 1
        entries = [(group,node) for group in range(0,ni)  
                                    for node in range(0,nn) 
                                        if ranks[(group,1)] == node ]
    else:
        assert False

    num_aij = len(entries)
    num_variables = 1 + num_aij   # including t

    # Objective function to minimise: t
    c = matrix( [[-1.0]] + [[0.0]] * num_aij).trans()

    # LHS of constraints: one per node and one per group
    Ai = []
    bi = []
    # Constraints [1]
    for node in range(0,nn):
        row = [0.0] #  coefficient of t is zero in [3]
        for (group,n) in entries:
            if n == node:
                if ranks[(group,0)] == node:
                    row.append(1.0)   # This variable has coefficient 1 in [1]
                else:
                    row.append(1.00001)   # This variable has coefficient 1 in [1]
            else:
                row.append(0)   # This variable has coefficient 0 in [1]
        Ai.append(row)          # LHS of [1]
        bi.append(C[node])      # RHS of [1]
    # Constraints [2]
    for group in range(0,ni):
        # Sum up the nodes in this group
        row = [L[group]]   # Coefficient of t is Lj (j=group)
        for (i,node) in entries:
            if i == group:                 # Variable affects this group
                row.append(-1.0)
            else:
                row.append(0.0)
        Ai.append(row)         # LHS of [2]
        bi.append(0.0)         # RHS of [2]
    # Constraints [3]
    for k,(group,node) in enumerate(entries):
        row = [0.0] * num_variables
        row[k+1] = -1.0
        Ai.append(row)        # LHS of [3]
        if ranks[(group,0)] == node:
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
    for k,(group,node) in enumerate(entries):
            opt_allocs[(group,node)] = sol['x'][1 + k]

    for group in range(0,ni):
        for node in range(0,nn):
            if B[group,node] != 0:
                if not (group,node) in opt_allocs:
                    assert policy != 'optimized' # only happens when possibilities set to zero
                    opt_allocs[(group,node)] = 0.0

    return opt_allocs


def make_integer(ni, nn, allocs, L, B, C):
    int_allocs = {}
    for node in range(0,nn):
        indices = [group for group in range(0,ni) if B[(group,node)]>0 and allocs[(group,node)] > 0 ]

        if len(indices) > 0:
            ncores      = [allocs[(group,node)] for group in indices]
            ncores_int  = [int(allocs[(group,node)]) for group in indices]
            total_cores = [sum([allocs[(group,n)] for n in range(0,nn) if B[(group,n)]>0]) for group in indices]
            # load_per_cores = [L[indices[j]] / total_cores[j] for j in range(0,len(total_cores))]
            # print 'node', node, 'indices', indices, 'cores', ncores, 'total_cores', total_cores, 'lpc', load_per_cores

            # After rebalancing, all groups should have the same load-per-core
            # Hence the slowdown is proportional to the fraction lost: (ncores - int(ncores)) / total_cores
            frac_lost_and_j = [(1.0 * (ncores[j] - int(ncores[j])) / total_cores[j],j) for j in range(0,len(total_cores))]
            # print 'node', node, 'indices', indices, 'cores', ncores, 'total_cores', total_cores, 'frac_lost', frac_lost_and_j
            # print 'ncores_int', ncores_int
            frac_lost_and_j.sort()
            frac_lost_and_j.reverse()

            extra_cores = C[group] - sum(ncores_int)
            # print 'extra_cores', extra_cores
            for c in range(0,extra_cores):
                # Sometimes it is not necessary to use all cores: in this case we just keep going
                # filling up the available cores anyway
                frac,j = frac_lost_and_j[c % len(frac_lost_and_j)]
                ncores_int[j] = ncores_int[j] + 1

        for j,group in enumerate(indices):
            int_allocs[(group,node)] = ncores_int[j]

    # Lost the zeros: put them back
    for group in range(0,ni):
        for node in range(0,nn):
            if B[(group,node)]>0:
                int_allocs[(group,node)] = int_allocs.get((group,node),0)
    return int_allocs




def run_policy(equal, policy):

    ni, nn, ranks, allocs, loads = read_current_alloc()

    # print 'ni =', ni
    # print 'nn =', nn
    # print 'allocs =', allocs
    # print 'loads =', loads
    Ll = []
    Brows = []
    for group in range(0,ni):
        load = 0
        Brow = []
        for node in range(0,nn):
            if (group,node) in loads:
                load += 1.0 * loads[(group,node)]
                Brow.append(1.0)
            else:
                Brow.append(0.0)
        Ll.append(load)
        Brows.append(Brow)

    # Modify problem depending on the policy
    if equal:
        # Ignore real loads, set all to 1.0
        Ll = [1.0] * ni


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
        opt_allocs = optimize(ni, nn, ranks, L, B, C, policy)
    print 'Optimized allocation'
    printout(ni,nn,ranks,opt_allocs,L)

    integer_allocs = make_integer(ni, nn, opt_allocs, L, B, C)
    print 'Integer allocation'
    printout(ni,nn,ranks,integer_allocs,L)

    write_new_alloc(ni, nn, ranks, B, integer_allocs)

def Usage(argv):
    print argv[0], '   opts <number-of-iterations>'
    print '   --help         Show this help'
    print '   --equal        Assume equal loads, not indicated loads'
    print '   --master       All work on the master'
    print '   --slaves       No work on the master'
    print '   --slave1       Work only on the first slave'

def main(argv):

    equal = False
    policy = 'optimized'
    try:
        opts, args = getopt.getopt( argv[1:], "h", ["help", "equal", "master", "slaves", "slave1"])
    except getopt.error, msg:
        print msg
        print "for help use --help"
        return 2
    for o,a in opts:
        if o in ('-h', '--help'):
            return Usage(argv)
        elif o == '--equal':
            equal = True
        elif o in ('--master', '--slaves', '--slave1'):
            policy = o[2:]
            print 'policy', policy

    niter = 1
    if len(args) == 1:
        niter = int(args[0])
    # print 'Number of iterations', niter
    for it in range(0,niter):
        if it > 0:
            time.sleep(2)
        run_policy(equal, policy)


if __name__ == '__main__':
    sys.exit(main(sys.argv))

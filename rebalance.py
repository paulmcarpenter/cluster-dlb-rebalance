#! /usr/bin/env python
import os
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

def isint(s):
	if re.match('^[0-9]*$', s):
		return True
	else:
		return False

def isfloat(s):
	if isint(s):
		return True
	if re.match('^[0-9]*\.[0-9]*$', s):
		return True
	else:
		return False

def read_current_alloc():
	# Read map files
	mapfiles = [f for f in os.listdir('.hybrid') if f.startswith('map')]
	ranks = {}
	extranktonode = []
	extranktogroup = []
	for j, mapfile in enumerate(sorted(mapfiles)):
		m = re.match('map([0-9]*)$', mapfile)
		if not m:
			return None
		extrank = int(m.group(1))
		if j != extrank:
			return
		f = open('.hybrid/' + mapfile, 'r')
		s = f.readline().split()  # externalRank r
		assert s[0] == 'externalRank'
		assert int(s[1]) == extrank
		s = f.readline().split()  # groupNum g
		assert s[0] == 'groupNum'
		groupNum = int(s[1])
		s = f.readline().split()  # internalRank i
		assert s[0] == 'internalRank'
		internalRank = int(s[1])
		s = f.readline().split()  # nodeNum n
		assert s[0] == 'nodeNum'
		nodeNum = int(s[1])

		extranktogroup.append(groupNum)
		extranktonode.append(nodeNum)
		ranks[ (groupNum,internalRank)] = nodeNum

	max_group = max(extranktogroup)
	max_node = max(extranktonode)

	fnames = os.listdir('.hybrid')
	allocs = {}
	loads = {}
	for fname in fnames:
		m = re.match(r'utilization([0-9*])', fname)
		if m:
			extrank = int(m.group(1))
			group = extranktogroup[extrank]
			node = extranktonode[extrank]

			f = open('.hybrid/' + fname)
			for line in f.readlines():
				pass
			l = line.strip().split(' ')
			if len(l) < 3:
				return None
			if not isfloat(l[0]) and isint(l[1]) and isfloat(l[2]):
				return None
			allocs[ (group,node) ] = float(l[1])
			loads[ (group,node) ] = float(l[2])
			f.close()
		
	return extranktonode, extranktogroup, max_group+1, max_node+1, ranks, allocs, loads
	
def write_new_alloc(ni, nn, ranks, B, opt_allocs):
	for group in range(0,ni):
		f = open('.hybrid/alloc%d' % group, 'w')
		for rank in range(0,nn):
			if (group,rank) in ranks.keys():
				node = ranks[(group,rank)]
				print >>f, opt_allocs[(group,node)]
			else:
				break
		f.close()

	
	

def printout(ni,nn,ranks,allocs,loads):
	print '        ' + ' ' * 5 * nn + 'NUM CORES ' + ' ' * 6 * nn + 'Total_cores   Load	 Load/Total_cores'
	print '        ' + ' ' * 5 * nn + '	 node	' + ' ' * 6 * nn + '	'
	print '           ' + ('%11d' * nn) % tuple(range(0,nn))

	for group in range(0, ni):
		print 'Instance %2d' % group,
		for node in range(0, nn):
			if (group,node) in allocs:
				print '%9.2f' % allocs[(group,node)],
				if ranks[(group,0)] == 0:
					print '#',
				else:
					print ' ',
			else:
				print '%11s' % '-',
		load = loads[group]
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
	print '			  ' + ('%10s' % '---') * nn
	print '			  ' + ('%10.2f' * nn) % tuple(used)

def optimize(ni, nn, ranks, L, B, C, policy, min_alloc):
	# Minimise		-t										  (i.e. maximise worst-case cores per load)
	#
	#					N
	# subject to	  Sum  a			   <= C		all i		(available cores)			[1]
	#				  j=1	ij				   i
	#
	#					N
	#				- Sum  a	 +	L  t   <= 0		all j		(Sum cores >= load * t)		[2]
	#				  i=1	ij		 j
	#
	#					   a			   >= 0		all i,j									[3a]
	#						ij
	#
	#					   a			   >= 1		when j is the master node (as main cannot be migrated) [3b]

	# Variables are t, a_ij  (but only when B_ij = 1; otherwise will get bogus values for the other a_ij)

	# All a_{ij}
	entries_all = [(group,node) for group in range(0,ni)  
								for node in range(0,nn) 
									if B[group,node] != 0]

	pos_all = range(0, len(entries_all))
	#  For which (group,node) are the a_{ij} variable?
	if policy == 'optimized':
		# Optimized: have all a_ij
		pos_var = pos_all
		
	elif policy == 'master':
		# Only have masters
		pos_var = [j for j,(group,node) in enumerate(entries_all) 
										if ranks[(group,0)] == node ]
	elif policy == 'slaves':
		# Only have slaves 
		pos_var = [j for j,(group,node) in enumerate(entries_all)
										if B[group,node] != 0 and ranks[(group,0)] != node ]

	elif policy == 'slave1':
		# Only have slave 1
		pos_var = [j for j,(group,node) in enumerate(entries_all)
										if ranks[(group,1)] == node ]
	else:
		assert False

	pos_fixed = [j for j in pos_all if not j in pos_var]

	num_aij_all = len(pos_all)
	num_aij_var = len(pos_var)
	num_variables = 1 + num_aij_all   # including t

	# Objective function to minimise: t
	c = matrix( [[-1.0]] + [[0.0]] * num_aij_all).trans()

	# LHS of constraints
	Ai = []
	bi = []
	# Constraints [1]: one per node
	for node in range(0,nn):
		row = [0.0] #  coefficient of t is zero in [3]
		for (group,n) in entries_all:
			if n == node:
				if ranks[(group,0)] == node:
					row.append(1.0)   # This variable has coefficient 1 in [1]
				else:
					row.append(1.00001)   # This variable has coefficient 1 in [1]
			else:
				row.append(0)	# This variable has coefficient 0 in [1]
		Ai.append(row)			# LHS of [1]
		bi.append(C[node])		# RHS of [1]
	# Constraints [2]: one per group
	for group in range(0,ni):
		# Sum up the nodes in this group
		row = [L[group]]   # Coefficient of t is Lj (j=group)
		for (i,node) in entries_all:
			if i == group:				   # Variable affects this group
				row.append(-1.0)
			else:
				row.append(0.0)
		Ai.append(row)		   # LHS of [2]
		bi.append(0.0)		   # RHS of [2]
	# Constraints [3]: for variable entries
	for k in pos_var:
		group,node = entries_all[k]
		row = [0.0] * num_variables
		row[k+1] = -1.0
		Ai.append(row)		  # LHS of [3]
		if ranks[(group,0)] == node:
			bi.append(-max(min_alloc,1.0))			# RHS of [3b]
		else:
			bi.append(-min_alloc)			# RHS of [3a]
	# Fixed entries (fixed to min_alloc)
	for k in pos_fixed:
		group,node = entries_all[k]
		print 'a_ij', group, node, 'must be', min_alloc
		# a_ij >= min_alloc
		row = [0.0] * num_variables
		row[k+1] = -1.0
		Ai.append(row)		  # LHS of [3]
		bi.append(-min_alloc)			# RHS of [3a]
		# a_ij <= min_alloc
		row = [0.0] * num_variables
		row[k+1] = 1.0
		Ai.append(row)		  # LHS of [3]
		bi.append(min_alloc)		   # RHS of [3a]

	A = matrix(Ai).trans()
	# print 'A', A.size
	# print(A)
	b = matrix(bi)
	# print 'b', b.size
	# print(b)

	solvers.options['show_progress'] = False
	sol = solvers.lp( c, A, b)

	opt_allocs = {}
	for k,(group,node) in enumerate(entries_all):
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
			ncores		= [allocs[(group,node)] for group in indices]
			ncores_int	= [int(allocs[(group,node)]) for group in indices]
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


def make_topology(ni, nn, nanosloads):
	Brows = []
	for group in range(0,ni):
		load = 0
		Brow = []
		for node in range(0,nn):
			if (group,node) in nanosloads:
				Brow.append(1.0)
			else:
				Brow.append(0.0)
		Brows.append(Brow)
	return Brows

def run_policy(ni, nn, ranks, allocs, topology, loads, policy, min_alloc):

	# print 'ni =', ni
	# print 'nn =', nn
	# print 'allocs =', allocs
	# print 'loads =', loads

	# Load vector
	L = matrix(loads)

	# Topology matrix
	B = matrix(topology).trans()

	# Available cores vector
	C = matrix( [[48]] * nn).trans()

	if max(L) == 0.0:
		# Currently no work!!!
		opt_allocs = allocs
		print 'No work!'
	else:
		opt_allocs = optimize(ni, nn, ranks, L, B, C, policy, min_alloc)
	integer_allocs = make_integer(ni, nn, opt_allocs, L, B, C)

	write_new_alloc(ni, nn, ranks, B, integer_allocs)

	return opt_allocs, integer_allocs

def Usage(argv):
	print argv[0], '   opts <number-of-iterations>'
	print '   --help		 Show this help'
	print '   --equal		 Assume equal loads, not indicated loads'
	print '   --master		 All work on the master'
	print '   --slaves		 No work on the master'
	print '   --slave1		 Work only on the first slave'
	print '   --loads l		 Specify the loads (comma-separated)'
	print '   --min m		 Set minimum allocation per instance to m cores'
	print '   --wait secs	 Wait time between rebalancings'

def main(argv):

	equal = False
	policy = 'optimized'
	cmdloads = None
	min_alloc = 1 # Every instance has at least one core
	wait_time = 2 # seconds
	try:
		opts, args = getopt.getopt( argv[1:], "h", ["help", "equal", "master", "slaves", "slave1", "loads=", "min=", 'wait='])
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
		elif o == '--loads':
			cmdloads = a
		elif o == '--min':
			min_alloc = int(a)
		elif o == '--wait':
			wait_time = int(a)

	niter = 1
	if len(args) == 1:
		niter = int(args[0])
	# print 'Number of iterations', niter
	did_rebalance = False
	for it in range(0,niter):
		if it > 0:
			if did_rebalance:
				time.sleep(wait_time)
			else:
				time.sleep(1)
			if os.path.exists('.kill'):
				print 'Rebalance.py killed by .kill file'
				os.system('rm .kill')
				break
		x = read_current_alloc()
		if x is None:
			# Happens if program changes: just wait and try again
			did_rebalance = False
			continue
		extranktonode, extranktogroup, ni, nn, ranks, allocs, nanosloads = x
		topology = make_topology(ni, nn, nanosloads)
		# print 'nanosloads', nanosloads

		# Modify problem depending on the policy
		if equal:
			# Ignore real loads, set all to 1.0
			loads = [1.0] * ni
		elif not cmdloads is None:
			# Use loads given on command line
			loads = [float(x) for x in cmdloads.split(',')]
		else:
			# Use loads provided by Nanos
			print 'ni=', ni
			loads = [0.0] * ni
			for (group,node),load in nanosloads.items():
				loads[group] += load

		print 'Current allocation'
		printout(ni,nn,ranks,allocs,loads)

		opt_allocs, integer_allocs = run_policy(ni, nn, ranks, allocs, topology, loads, policy, min_alloc)

		print 'Optimized allocation'
		printout(ni,nn,ranks,opt_allocs,loads)

		print 'Integer allocation'
		printout(ni,nn,ranks,integer_allocs,loads)
		did_rebalance = True



if __name__ == '__main__':
	sys.exit(main(sys.argv))

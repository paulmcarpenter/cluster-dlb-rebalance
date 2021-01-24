#! /usr/bin/env python
import os
import sys
import cvxopt
import re
import time
import getopt
from cvxopt import matrix, solvers

from os import listdir

monitor_time = 1.0

def average(l):
	return sum(l) / len(l)

# Check whether the string is an integer
def isint(s):
	if re.match('^[0-9]*$', s):
		return True
	else:
		return False

# Check whether the string is a float
def isfloat(s):
	if isint(s):
		return True
	if re.match('^[0-9]*\.[0-9]*$', s):
		return True
	else:
		return False


def read_map_entry(label, line):
	s = line.split()
	assert(s[0] == label)
	return int(s[1])

def read_current_alloc():
	global monitor_time
	if not os.path.exists('.hybrid'):
		return None
	# Read map files
	mapfiles = [f for f in os.listdir('.hybrid') if f.startswith('map')]

	# Get all the external ranks for which we have a map file
	extranks = []
	for mapfile in mapfiles:
		m = re.match('map([0-9]*)$', mapfile)
		if not m:
			return None
		extranks.append(int(m.group(1)))

	ranks = {}
	extranktonode = []
	extranktogroup = []
	for j, extrank in enumerate(sorted(extranks)):
		if j != extrank:
			return
		f = open('.hybrid/map%d' % extrank, 'r')

		extrank2 = read_map_entry('externalRank', f.readline())
		assert extrank2 == extrank
		groupNum = read_map_entry('groupNum', f.readline())
		internalRank = read_map_entry('internalRank', f.readline())
		nodeNum = read_map_entry('nodeNum', f.readline())

		extranktogroup.append(groupNum)
		extranktonode.append(nodeNum)
		ranks[ (groupNum,internalRank)] = nodeNum
		
	# print('ranks(group, internalrank): ', ranks)
	max_group = max(extranktogroup)
	max_node = max(extranktonode)

	fnames = os.listdir('.hybrid')
	allocs = {}
	loads = {}
	for fname in fnames:
		m = re.match(r'utilization([0-9]*)$', fname)
		if m:
			extrank = int(m.group(1))
			group = extranktogroup[extrank]
			node = extranktonode[extrank]

			f = open('.hybrid/' + fname)
			line = None
			log = []
			for line in f.readlines():
				s = line.strip().split(' ')
				if len(s) < 4:
					return None
				#           Timestamp     Global alloc   curr DLB alloc  Total busy cores
				log.append( (float(s[0]), int(s[1]),     int(s[2]),      float(s[3])) )
			if len(log) == 0:
				return None
			end_time = log[-1][0]
			recent_log = [ (t,alloc,dlballoc,load) for (t,alloc,dlballoc,load) in log  if t >= end_time - monitor_time ]
			# print('recent_log', recent_log)

			assert not (group,node) in allocs
			assert not (group,node) in loads
			allocs[ (group,node) ] = average( [ alloc for (t,alloc,dlballoc,load) in recent_log])
			loads[ (group,node) ] = average( [ load for (t,alloc,dlballoc,load) in recent_log])
			f.close()
		
	return extranktonode, extranktogroup, max_group+1, max_node+1, ranks, allocs, loads
	
def write_new_alloc(ni, nn, ranks, B, opt_allocs):
	for group in range(0,ni):
		f = open('.hybrid/alloc%d' % group, 'w')
		for rank in range(0,nn):
			if (group,rank) in ranks.keys():
				node = ranks[(group,rank)]
				print(opt_allocs[(group,node)], file=f)
			else:
				break
		f.close()

	
	

def printout(ni,nn,ranks,allocs,loads):
	print('     ' + ' ' * 5 * nn + 'NUM CORES ' + ' ' * 6 * nn + 'Total_cores   Load	 Load/Total_cores')
	print('     ' + ' ' * 5 * nn + '	 node	' + ' ' * 6 * nn + '	')
	print('       ' + ('%11d' * nn) % tuple(range(0,nn)))

	for group in range(0, ni):
		print('Group %2d ' % group, end='')
		for node in range(0, nn):
			if (group,node) in allocs:
				print('%9.2f ' % allocs[(group,node)], end='')
				if ranks[(group,0)] == node:
					print('# ', end='')
				else:
					print(' ', end='')
			else:
				print('%11s ' % '-', end='')
		load = loads[group]
		total_c = 0
		for node in range(0, nn):
			total_c += allocs.get((group,node),0)
		print('%10.2f ' % total_c, end='')
		print('%10.2f ' % load, end='')
		if total_c > 0:
			print('%10.3f ' % (load/total_c), end='')
		else:
			print('%10s' % 'inf', end='')
		print()


	used = []
	for node in range(0,nn):
		u = 0
		for group in range(0,ni):
			u += allocs.get((group,node),0)
		used.append(u)
	print('      ' + ('%10s' % '---') * nn)
	print('      ' + ('%10.2f' * nn) % tuple(used))

def optimize(ni, nn, ranks, L, B, C, policy, min_master, min_slave):
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
			bi.append(-min_master)			# RHS of [3b]
		else:
			bi.append(-min_slave)			# RHS of [3a]
	# Fixed entries (fixed to min_alloc)
	for k in pos_fixed:
		group,node = entries_all[k]
		if ranks[(group,0)] == node:
			min_alloc = min_master
		else:
			min_alloc = min_slave
		print('a_ij', group, node, 'must be', min_alloc)
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
	# print('A', A.size)
	# print(A)
	b = matrix(bi)
	# print('b', b.size)
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


def make_integer(ni, nn, allocs, L, B, C, fill_idle):
	int_allocs = {}
	for node in range(0,nn):
		indices = [group for group in range(0,ni) if B[(group,node)]>0 and allocs[(group,node)] > 0 ]

		if len(indices) > 0:
			ncores		= [allocs[(group,node)] for group in indices]
			ncores_int	= [int(allocs[(group,node)]) for group in indices]
			total_cores = [sum([allocs[(group,n)] for n in range(0,nn) if B[(group,n)]>0]) for group in indices]
			# load_per_cores = [L[indices[j]] / total_cores[j] for j in range(0,len(total_cores))]
			# print('node', node, 'indices', indices, 'cores', ncores, 'total_cores', total_cores, 'lpc', load_per_cores)

			# After rebalancing, all groups should have the same load-per-core
			# Hence the slowdown is proportional to the fraction lost: (ncores - int(ncores)) / total_cores
			frac_lost_and_j = [(1.0 * (ncores[j] - int(ncores[j])) / total_cores[j],j) for j in range(0,len(total_cores))]
			# print('node', node, 'indices', indices, 'cores', ncores, 'total_cores', total_cores, 'frac_lost', frac_lost_and_j)
			# print('ncores_int', ncores_int)
			frac_lost_and_j.sort()
			frac_lost_and_j.reverse()

			# Number of cores not used so far
			num_unused_cores = C[node] - sum(ncores_int)

			if fill_idle:
				# Use up all the unallocated cores on the node
				num_extra_cores = num_unused_cores
			else:
				# Not fill_idle: each truncated process can get at most one extra core
				num_extra_cores = min(num_unused_cores, len(frac_lost_and_j))
				# print 'extra_cores', extra_cores

			for c in range(0,num_extra_cores):
				# Sometimes it is not necessary to use all cores: in this case we just keep going
				# filling up the available cores anyway. The modulo operation just hands them
				# out round-robin.
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

def run_policy(ni, nn, ranks, allocs, topology, loads, policy, min_master, min_slave, fill_idle):

	# print('ni =', ni)
	# print('nn =', nn)
	# print('allocs =', allocs)
	# print('loads =', loads)

	# Load vector
	L = matrix(loads)

	# Topology matrix
	B = matrix(topology).trans()

	# Available cores vector
	C = matrix( [[48]] * nn).trans()

	if max(L) == 0.0:
		# Currently no work!!!
		opt_allocs = allocs
		print('No work!')
	else:
		opt_allocs = optimize(ni, nn, ranks, L, B, C, policy, min_master, min_slave)
	integer_allocs = make_integer(ni, nn, opt_allocs, L, B, C, fill_idle)

	write_new_alloc(ni, nn, ranks, B, integer_allocs)

	return opt_allocs, integer_allocs

def Usage(argv):
	print(argv[0], '   opts <number-of-iterations>')
	print('   --help		 Show this help')
	print('   --equal		 Assume equal loads, not indicated loads')
	print('   --master		 All work on the master')
	print('   --slaves		 No work on the master')
	print('   --slave1		 Work only on the first slave')
	print('   --loads l		 Specify the loads (comma-separated)')
	print('   --min m		 Set minimum allocation per instance to m cores')
	print('   --monitor secs Monitor load over past time period of given length')
	print('   --wait secs	 Wait time between rebalancings')
	print('   --no-fill      Do not use whole nodes if not necessary')

def main(argv):

	global monitor_time
	equal = False
	policy = 'optimized'
	cmdloads = None
	min_master = None
	min_slave = None
	wait_time = 2 # seconds
	fill_idle = True
	try:
		opts, args = getopt.getopt( argv[1:], "h", ["help", "equal", "master",
		                                            "slaves", "slave1", "loads=",
													"min=", 'min-master=', 'min-slave=',
													'wait=', 'monitor=', 'no-fill'])
	except getopt.error as msg:
		print(msg)
		print("for help use --help")
		return 2
	for o,a in opts:
		if o in ('-h', '--help'):
			return Usage(argv)
		elif o == '--equal':
			equal = True
		elif o in ('--master', '--slaves', '--slave1'):
			policy = o[2:]
			print('policy', policy)
		elif o == '--loads':
			cmdloads = a
		elif o == '--min-master':
			assert min_master is None
			min_master = int(a)
		elif o == '--min':
			assert min_master is None
			assert min_slave is None
			min_master = int(a)
			min_slave = int(a)
		elif o == '--min-slave':
			assert min_slave is None
			min_slave = int(a)
		elif o == '--wait':
			wait_time = int(a)
		elif o == '--monitor':
			monitor_time = float(a)
		elif o == '--no-fill':
			fill_idle = False

	if min_master is None:
		min_master = 1
	if min_slave is None:
		min_slave = 1

	niter = 1
	if len(args) == 1:
		niter = int(args[0])
	# print('Number of iterations', niter)
	did_rebalance = False
	for it in range(0,niter):
		if it > 0:
			sys.stdout.flush()
			if did_rebalance:
				time.sleep(wait_time)
			else:
				time.sleep(1)
			if os.path.exists('.kill'):
				print('Rebalance.py killed by .kill file')
				os.system('rm .kill')
				break
		x = read_current_alloc()
		if x is None:
			# Happens if program changes: just wait and try again
			did_rebalance = False
			continue
		extranktonode, extranktogroup, ni, nn, ranks, allocs, nanosloads = x
		topology = make_topology(ni, nn, nanosloads)
		# print('nanosloads', nanosloads)

		# Modify problem depending on the policy
		if equal:
			# Ignore real loads, set all to 1.0
			loads = [1.0] * ni
		elif not cmdloads is None:
			# Use loads given on command line
			loads = [float(x) for x in cmdloads.split(',')]
		else:
			# Use loads provided by Nanos
			loads = [0.0] * ni
			for (group,node),load in nanosloads.items():
				loads[group] += load

		print('Current allocation')
		printout(ni,nn,ranks,allocs,loads)

		opt_allocs, integer_allocs = run_policy(ni, nn, ranks, allocs, topology, loads, policy, min_master, min_slave, fill_idle)

		print('Optimized allocation')
		printout(ni,nn,ranks,opt_allocs,loads)

		print('Integer allocation')
		printout(ni,nn,ranks,integer_allocs,loads)
		did_rebalance = True



if __name__ == '__main__':
	sys.exit(main(sys.argv))

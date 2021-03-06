#! /usr/bin/env python
import os
import sys
import re
import time
import getopt

try:
	import networkx
	canImportNetworkX = True
except ImportError:
	canImportNetworkX = False

from os import listdir

num_vranks = 0
num_nodes = 0

def vrank(v):
        global num_vranks
        assert 0 <= v and v < num_vranks
        return 'Vrank%d' % v

def node(i):
        global num_nodes
        assert 0 <= i and i < num_nodes
        return 'Node%d' % i


monitor_time = 1.0
hybrid_directory = '.hybrid'

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
	global num_vranks
	global num_nodes
	if not os.path.exists(hybrid_directory):
		return None
	# Read map files
	mapfiles = [f for f in os.listdir(hybrid_directory) if f.startswith('map')]

	# Get all the external ranks for which we have a map file
	extranks = []
	for mapfile in mapfiles:
		m = re.match('map([0-9]*)$', mapfile)
		if not m:
			return None
		extranks.append(int(m.group(1)))

	ranks = {}
	extranktonode = []
	extranktoapprank = []
	extranktomaster = []
	extranksonnode = {}

	for j, extrank in enumerate(sorted(extranks)):
		if j != extrank:
			return
		f = open(hybrid_directory + '/map%d' % extrank, 'r')

		extrank2 = read_map_entry('externalRank', f.readline())
		assert extrank2 == extrank
		apprankNum = read_map_entry('apprankNum', f.readline())
		internalRank = read_map_entry('internalRank', f.readline())
		nodeNum = read_map_entry('nodeNum', f.readline())
		indexThisNode = read_map_entry('indexThisNode', f.readline())
		isMaster = (internalRank==0)

		extranktoapprank.append(apprankNum)
		extranktonode.append(nodeNum)
		extranktomaster.append(isMaster)
		ranks[ (apprankNum,internalRank)] = nodeNum
		if not nodeNum in extranksonnode:
			extranksonnode[nodeNum] = []
		extranksonnode[nodeNum].append(extrank)
		
	max_apprank = max(extranktoapprank)
	max_node = max(extranktonode)
	num_vranks = max_apprank+1
	num_nodes = max_node+1

	# Now calculate baseline allocation
	baseline_alloc = dict([((a,n),0) for a in range(0,num_vranks) for n in range(0,num_nodes)])
	for nodeNum in range(0,num_nodes):
		ee = extranksonnode[nodeNum]
		num_masters = len([e for e in ee if extranktomaster[e]])
		num_workers = len(ee) - num_masters
		print('Node', nodeNum, 'nm: ', num_masters, 'nw: ', num_workers)
		for e in ee:
			apprankNum = extranktoapprank[e]
			assert nodeNum == extranktonode[e]
			if extranktomaster[e]:
				num_cores = (48-num_workers) / num_masters   # TODO: remove fixed 48
			else:
				num_cores = 1
			baseline_alloc[(apprankNum, nodeNum)] = num_cores

	print('baseline_alloc', baseline_alloc)

	G = networkx.DiGraph()
	for nodeNum in range(0,num_nodes):
		G.add_node(node(nodeNum))
	for vrankNum in range(0,num_vranks):
		G.add_node(vrank(vrankNum))

	fnames = os.listdir(hybrid_directory)
	allocs = {}
	loads = {}
	for fname in fnames:
		m = re.match(r'utilization([0-9]*)$', fname)
		if m:
			extrank = int(m.group(1))
			apprank = extranktoapprank[extrank]
			nodeNum = extranktonode[extrank]

			f = open(hybrid_directory + '/' + fname)
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

			assert not (apprank,nodeNum) in allocs
			assert not (apprank,nodeNum) in loads
			delta_alloc = average( [ alloc for (t,alloc,dlballoc,load) in recent_log]) - baseline_alloc[(apprank,nodeNum)]
			if delta_alloc > 0:
				G.add_edge(vrank(apprank), node(nodeNum), capacity=delta_alloc)
			else:
				G.add_edge(node(nodeNum), vrank(apprank), capacity=-delta_alloc)

			allocs[ (apprank,nodeNum) ] = delta_alloc
			loads[ (apprank,nodeNum) ] = average( [ load for (t,alloc,dlballoc,load) in recent_log])
			f.close()
		
	return extranktonode, extranktoapprank, num_vranks, num_nodes, ranks, allocs, loads, G
	
def write_new_alloc(ni, nn, ranks, B, opt_allocs):
	for apprank in range(0,ni):
		f = open(hybrid_directory + '/alloc%d' % apprank, 'w')
		for rank in range(0,nn):
			if (apprank,rank) in ranks.keys():
				nodeNum = ranks[(apprank,rank)]
				print(opt_allocs[(apprank,nodeNum)], file=f)
			else:
				break
		f.close()

	
	

def printout(ni,nn,ranks,allocs):
	print('       ' + ' ' * 5 * nn + 'NUM CORES ' + ' ' * 6 * nn + 'Total_cores   Load	 Load/Total_cores')
	print('       ' + ' ' * 5 * nn + '	 node	' + ' ' * 6 * nn + '	')
	print('         ' + ('%11d' * nn) % tuple(range(0,nn)))

	for apprank in range(0, ni):
		print('Apprank %2d ' % apprank, end='')
		for nodeNum in range(0, nn):
			if (apprank,nodeNum) in allocs:
				print('%9.2f ' % allocs[(apprank,nodeNum)], end='')
				if ranks[(apprank,0)] == nodeNum:
					print('# ', end='')
				else:
					print(' ', end='')
			else:
				print('%11s ' % '-', end='')
		print()


	used = []
	for nodeNum in range(0,nn):
		u = 0
		for apprank in range(0,ni):
			u += allocs.get((apprank,nodeNum),0)
		used.append(u)
	print('        ' + ('%10s' % '---') * nn)
	print('        ' + ('%10.2f' * nn) % tuple(used))

#def optimize(ni, nn, ranks, L, B, C, policy, min_master, min_slave):
#	# Minimise		-t										  (i.e. maximise worst-case cores per load)
#	#
#	#					N
#	# subject to	  Sum  a			   <= C		all i		(available cores)			[1]
#	#				  j=1	ij				   i
#	#
#	#					N
#	#				- Sum  a	 +	L  t   <= 0		all j		(Sum cores >= load * t)		[2]
#	#				  i=1	ij		 j
#	#
#	#					   a			   >= 0		all i,j									[3a]
#	#						ij
#	#
#	#					   a			   >= 1		when j is the master node (as main cannot be migrated) [3b]
#
#	# Variables are t, a_ij  (but only when B_ij = 1; otherwise will get bogus values for the other a_ij)
#
#	# All a_{ij}
#	entries_all = [(apprank,node) for apprank in range(0,ni)  
#								for node in range(0,nn) 
#									if B[apprank,node] != 0]
#
#	pos_all = range(0, len(entries_all))
#	#  For which (apprank,node) are the a_{ij} variable?
#	if policy == 'optimized':
#		# Optimized: have all a_ij
#		pos_var = pos_all
#		
#	elif policy == 'master':
#		# Only have masters
#		pos_var = [j for j,(apprank,node) in enumerate(entries_all) 
#										if ranks[(apprank,0)] == node ]
#	elif policy == 'slaves':
#		# Only have slaves 
#		pos_var = [j for j,(apprank,node) in enumerate(entries_all)
#										if B[apprank,node] != 0 and ranks[(apprank,0)] != node ]
#
#	elif policy == 'slave1':
#		# Only have slave 1
#		pos_var = [j for j,(apprank,node) in enumerate(entries_all)
#										if ranks[(apprank,1)] == node ]
#	else:
#		assert False
#
#	pos_fixed = [j for j in pos_all if not j in pos_var]
#
#	num_aij_all = len(pos_all)
#	num_aij_var = len(pos_var)
#	num_variables = 1 + num_aij_all   # including t
#
#	# Objective function to minimise: t
#	c = matrix( [[-1.0]] + [[0.0]] * num_aij_all).trans()
#
#	# LHS of constraints
#	Ai = []
#	bi = []
#	# Constraints [1]: one per node
#	for node in range(0,nn):
#		row = [0.0] #  coefficient of t is zero in [3]
#		for (apprank,n) in entries_all:
#			if n == node:
#				if ranks[(apprank,0)] == node:
#					row.append(1.0)   # This variable has coefficient 1 in [1]
#				else:
#					row.append(1.00001)   # This variable has coefficient 1 in [1]
#			else:
#				row.append(0)	# This variable has coefficient 0 in [1]
#		Ai.append(row)			# LHS of [1]
#		bi.append(C[node])		# RHS of [1]
#	# Constraints [2]: one per apprank
#	for apprank in range(0,ni):
#		# Sum up the nodes in this apprank
#		row = [L[apprank]]   # Coefficient of t is Lj (j=apprank)
#		for (i,node) in entries_all:
#			if i == apprank:				   # Variable affects this apprank
#				row.append(-1.0)
#			else:
#				row.append(0.0)
#		Ai.append(row)		   # LHS of [2]
#		bi.append(0.0)		   # RHS of [2]
#	# Constraints [3]: for variable entries
#	for k in pos_var:
#		apprank,node = entries_all[k]
#		row = [0.0] * num_variables
#		row[k+1] = -1.0
#		Ai.append(row)		  # LHS of [3]
#		if ranks[(apprank,0)] == node:
#			bi.append(-min_master)			# RHS of [3b]
#		else:
#			bi.append(-min_slave)			# RHS of [3a]
#	# Fixed entries (fixed to min_alloc)
#	for k in pos_fixed:
#		apprank,node = entries_all[k]
#		if ranks[(apprank,0)] == node:
#			min_alloc = min_master
#		else:
#			min_alloc = min_slave
#		print('a_ij', apprank, node, 'must be', min_alloc)
#		# a_ij >= min_alloc
#		row = [0.0] * num_variables
#		row[k+1] = -1.0
#		Ai.append(row)		  # LHS of [3]
#		bi.append(-min_alloc)			# RHS of [3a]
#		# a_ij <= min_alloc
#		row = [0.0] * num_variables
#		row[k+1] = 1.0
#		Ai.append(row)		  # LHS of [3]
#		bi.append(min_alloc)		   # RHS of [3a]
#
#	A = matrix(Ai).trans()
#	# print('A', A.size)
#	# print(A)
#	b = matrix(bi)
#	# print('b', b.size)
#	# print(b)
#
#	solvers.options['show_progress'] = False
#	sol = solvers.lp( c, A, b)
#
#	opt_allocs = {}
#	for k,(apprank,node) in enumerate(entries_all):
#			opt_allocs[(apprank,node)] = sol['x'][1 + k]
#
#	for apprank in range(0,ni):
#		for node in range(0,nn):
#			if B[apprank,node] != 0:
#				if not (apprank,node) in opt_allocs:
#					assert policy != 'optimized' # only happens when possibilities set to zero
#					opt_allocs[(apprank,node)] = 0.0
#
#	return opt_allocs


# def make_integer(ni, nn, allocs, L, B, C, fill_idle):
# 	int_allocs = {}
# 	for node in range(0,nn):
# 		indices = [apprank for apprank in range(0,ni) if B[(apprank,node)]>0 and allocs[(apprank,node)] > 0 ]
# 
# 		if len(indices) > 0:
# 			ncores		= [allocs[(apprank,node)] for apprank in indices]
# 			ncores_int	= [int(allocs[(apprank,node)]) for apprank in indices]
# 			total_cores = [sum([allocs[(apprank,n)] for n in range(0,nn) if B[(apprank,n)]>0]) for apprank in indices]
# 			# load_per_cores = [L[indices[j]] / total_cores[j] for j in range(0,len(total_cores))]
# 			# print('node', node, 'indices', indices, 'cores', ncores, 'total_cores', total_cores, 'lpc', load_per_cores)
# 
# 			# After rebalancing, all appranks should have the same load-per-core
# 			# Hence the slowdown is proportional to the fraction lost: (ncores - int(ncores)) / total_cores
# 			frac_lost_and_j = [(1.0 * (ncores[j] - int(ncores[j])) / total_cores[j],j) for j in range(0,len(total_cores))]
# 			# print('node', node, 'indices', indices, 'cores', ncores, 'total_cores', total_cores, 'frac_lost', frac_lost_and_j)
# 			# print('ncores_int', ncores_int)
# 			frac_lost_and_j.sort()
# 			frac_lost_and_j.reverse()
# 
# 			# Number of cores not used so far
# 			num_unused_cores = C[node] - sum(ncores_int)
# 
# 			if fill_idle:
# 				# Use up all the unallocated cores on the node
# 				num_extra_cores = num_unused_cores
# 			else:
# 				# Not fill_idle: each truncated process can get at most one extra core
# 				num_extra_cores = min(num_unused_cores, len(frac_lost_and_j))
# 				# print 'extra_cores', extra_cores
# 
# 			for c in range(0,num_extra_cores):
# 				# Sometimes it is not necessary to use all cores: in this case we just keep going
# 				# filling up the available cores anyway. The modulo operation just hands them
# 				# out round-robin.
# 				frac,j = frac_lost_and_j[c % len(frac_lost_and_j)]
# 				ncores_int[j] = ncores_int[j] + 1
# 
# 		for j,apprank in enumerate(indices):
# 			int_allocs[(apprank,node)] = ncores_int[j]
# 
# 	# Lost the zeros: put them back
# 	for apprank in range(0,ni):
# 		for node in range(0,nn):
# 			if B[(apprank,node)]>0:
# 				int_allocs[(apprank,node)] = int_allocs.get((apprank,node),0)
# 	return int_allocs


#def make_topology(ni, nn, nanosloads):
#	Brows = []
#	for apprank in range(0,ni):
#		load = 0
#		Brow = []
#		for nodeNum in range(0,nn):
#			if (apprank,nodeNum) in nanosloads:
#				Brow.append(1.0)
#			else:
#				Brow.append(0.0)
#		Brows.append(Brow)
#	return Brows

# def run_policy(ni, nn, ranks, allocs, topology, loads, policy, min_master, min_slave, fill_idle):
# 
# 	# print('ni =', ni)
# 	# print('nn =', nn)
# 	# print('allocs =', allocs)
# 	# print('loads =', loads)
# 
# 	# Load vector
# 	L = matrix(loads)
# 
# 	# Topology matrix
# 	B = matrix(topology).trans()
# 
# 	# Available cores vector
# 	C = matrix( [[48]] * nn).trans()
# 
# 	if max(L) == 0.0:
# 		# Currently no work!!!
# 		opt_allocs = allocs
# 		print('No work!')
# 	else:
# 		opt_allocs = optimize(ni, nn, ranks, L, B, C, policy, min_master, min_slave)
# 	integer_allocs = make_integer(ni, nn, opt_allocs, L, B, C, fill_idle)
# 
# 	write_new_alloc(ni, nn, ranks, B, integer_allocs)
# 
# 	return opt_allocs, integer_allocs

def Usage(argv):
	print(argv[0], '   opts <number-of-iterations>')
	print('   --help		         Show this help')
	print('   --slaves		         No work on the master')
	print('   --slave1		         Work only on the first slave')
	print('   --min m		         Set minimum allocation per instance to m cores')
	print('   --monitor secs         Monitor load over past time period of given length')
	print('   --wait secs	         Wait time between rebalancings')
	print('   --no-fill              Do not use whole nodes if not necessary')
	print('   --hybrid-directory d   Path to hybrid-directory: default .hybrid')

def main(argv):

	if not canImportNetworkX:
		if len(argv) >= 2 and argv[1] == '--recurse':
			print('Error with recursive invocation')
			return 1
		ret = os.system('module load python/3.6.1; python ' + argv[0] + ' --recurse ' + ' '.join(argv[1:]))
		return ret


	print('hello')
	global monitor_time
	global hybrid_directory

	min_master = None
	min_slave = None
	wait_time = 2 # seconds
	fill_idle = True
	try:
		opts, args = getopt.getopt( argv[1:], "h", ["help", 
													"min=", 'min-master=', 'min-slave=',
													'wait=', 'monitor=', 'no-fill', 'recurse', 'hybrid-directory='])
	except getopt.error as msg:
		print(msg)
		print("for help use --help")
		return 2
	for o,a in opts:
		if o in ('-h', '--help'):
			return Usage(argv)
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
		elif o == '--recurse':
			pass
		elif o == '--hybrid-directory':
			hybrid_directory = a

	if min_master is None:
		min_master = 1
	if min_slave is None:
		min_slave = 1

	niter = 1
	if len(args) == 1:
		niter = int(args[0])
	print('Number of iterations', niter)
	did_rebalance = False
	for it in range(0,niter):
		if it > 0:
			sys.stdout.flush()
			if did_rebalance:
				time.sleep(wait_time)
			else:
				time.sleep(1)
			if os.path.exists('.kill'):
				print('killcycles.py killed by .kill file')
				os.system('rm .kill')
				break
		x = read_current_alloc()
		if x is None:
			# Happens if program changes: just wait and try again
			did_rebalance = False
			continue
		extranktonode, extranktoapprank, ni, nn, ranks, allocs, nanosloads, G = x
		#topology = make_topology(ni, nn, nanosloads)
		# print('nanosloads', nanosloads)

		print('Current allocation')
		printout(ni,nn,ranks,allocs)

		for vrankNum in range(0,ni):
			v = vrank(vrankNum)
			mf = networkx.maximum_flow(G, v, v)
			print('Apprank', vrankNum, 'flow', mf)

		# opt_allocs, integer_allocs = run_policy(ni, nn, ranks, allocs, topology, loads, policy, min_master, min_slave, fill_idle)

		# print('Optimized allocation')
		# printout(ni,nn,ranks,opt_allocs,loads)

		# print('Integer allocation')
		# printout(ni,nn,ranks,integer_allocs,loads)
		# did_rebalance = True



if __name__ == '__main__':
	sys.exit(main(sys.argv))

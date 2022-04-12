Readme
======

OmpSs-2@Cluster + DLB global allocation policy.

Usage
-----

	./rebalance.py    opts <number-of-iterations>
	   --help		         Show this help
	   --equal		         Assume equal loads, not indicated loads
	   --master		         All work on the master
	   --slaves		         No work on the master
	   --slave1		         Work only on the first slave
	   --loads l		     Specify the loads (comma-separated)
	   --min m		         Set minimum allocation per instance to m cores
	   --monitor secs        Monitor load over past time period of given length
	   --wait secs	         Wait time between rebalancings
	   --no-fill             Do not use whole nodes if not necessary
	   --hybrid-directory d  Path to hybrid-directory: default .hybrid
	   --local               Kill cycles when using local policy
	   --concise             Concise print (default if >1 iteration)
	   --verbose             Verbose print (default if 1 iteration)

It should not normally have to be executed directly by the user, as it would
normally be launched by the runhybrid.py script.

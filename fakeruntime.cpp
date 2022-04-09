#include <iostream>
#include <fstream>
#include <pthread.h>
#include <unistd.h>
#include <vector>

#include <stdarg.h>
#include <sys/stat.h>
#include <cassert>
#include <errno.h>

#define MAX_NODES 16
#define TOTAL_CORES_PER_NODE 48

int numNodes = 0;
int numAppranks = 0;

double loadPerCore[MAX_NODES];

class NanosInstance
{
    private:
        int _apprankNum;
        int _totalLoad;
        int _numIranks;
        pthread_t _thread;
        std::vector<int> _extranks; // for each i, which node is rank i on?
		std::vector<int> _nodes;
        std::vector<int> _cores; // for each i, how many cores does rank i have?

    public:

        NanosInstance(int apprankNum, int totalLoad, int numIranks, ...);

        void thread(void)
        {
			char filename[1024];
            std::cout << "Thread for " << _apprankNum << "\n";
			std::ofstream outfile[MAX_NODES];
			std::ifstream allocfile;
			for (int rank=0; rank<_numIranks; rank++) {
				sprintf(filename, ".hybrid/utilization%d", _extranks[rank]);
				outfile[rank].open(filename);
			}

            for (int i=0; i<_numIranks; i++)
            {
                std::cout << " " << _apprankNum << " " << _extranks[i] << "\n";
            }

			double timestamp = 0.0;
            for(int iter = 0; iter < 1000; iter++)
            {
				timestamp += 0.5;
                // Calculate total number of cores
                int totalCores = 0;
                for (int i=0; i<_numIranks; i++)
                {
                    totalCores += _cores[i];
                }
				loadPerCore[_apprankNum] = (float)_totalLoad / totalCores;
				std::cout << "loadPerCore " << _apprankNum << " is " << loadPerCore[_apprankNum] << "\n";
				double maxLoadPerCore = 0;
				for (int j=0; j< numAppranks; j++) {
					if (loadPerCore[j] > maxLoadPerCore) {
						maxLoadPerCore = loadPerCore[j];
					}
				}
				double fracBusy = loadPerCore[_apprankNum] / maxLoadPerCore;
                std::cout << _apprankNum << ": iter=" << iter << ", cores=" << totalCores << "\n";

                // Print theoretical load
                for (int rank=0; rank<_numIranks; rank++)
                {
                    double nodeLoad = (double)_totalLoad * _cores[rank] / totalCores;
					double busyCores = _cores[rank] * fracBusy;
                    outfile[rank] << timestamp << " " << _cores[rank] << " " << _cores[rank] << " " << busyCores << "\n";
					outfile[rank].flush();
                }

                // Get the allocation of cores
                char filename[1024];
                sprintf(filename, ".hybrid/alloc%d", _apprankNum);
                allocfile.open(filename);
                if (allocfile.is_open())
                {
                    for (int rank=0; rank<_numIranks; rank++)
                    {
                        int nodeAlloc;
                        allocfile >> nodeAlloc;
                        if (nodeAlloc != _cores[rank])
                        {
                            std::cout << "new allocation of " << nodeAlloc << " cores for instance " << _apprankNum << " rank " << rank << "\n";
                            _cores[rank] = nodeAlloc;
                        }
                    }
                    allocfile.close();
                }

                usleep(500000);
            }
			for (int rank=0; rank<_numIranks; rank++)
			{
				outfile[rank].close();
			}
        }

        void wait(void)
        {
            void *status;
            int rc = pthread_join(_thread, &status);
        }

        ~NanosInstance()
        {
            wait();
        }
};

static void *instance_thread(void *arg)
{
    NanosInstance *ni = static_cast<NanosInstance *>(arg);
    std::cout << "hi\n";
    ni->thread();
}


NanosInstance::NanosInstance(int apprankNum, int totalLoad, int numIranks, ...) 
                            : _apprankNum(apprankNum), _totalLoad(totalLoad), 
                              _numIranks(numIranks), _extranks(numIranks), _cores(numIranks)
{
    int rc;

    va_list args;
    char hostname[1024];
    char buffer[1024];
    gethostname(hostname, 1024);

	// Set up extrank
    va_start(args, numIranks);
	_extranks.resize(numIranks);
	_nodes.resize(numIranks);
    for (int i=0; i<numIranks; i++)
    {
        _extranks[i] = va_arg(args, int);
		_nodes[i] = _extranks[i] % numNodes;
		std::cout << "Apprank " << apprankNum << " irank " << i << " has extRank " << _extranks[i] << " on node " << _nodes[i] << "\n";
        _cores[i] = 0;
    }

    char filename[1024];
	// Create map files for the instance
    for (int i=0; i<numIranks; i++) {
		sprintf(filename, ".hybrid/map%d", _extranks[i]);
		std::ofstream myfile;
		myfile.open(filename);
		myfile << "externalRank " << _extranks[i] << "\n";
		myfile << "apprankNum " << apprankNum << "\n";
		myfile << "internalRank " << i << "\n";
		myfile << "nodeNum " << _nodes[i] << "\n";
		myfile << "indexThisNode " << _extranks[i] / numNodes << "\n";
		myfile << "cpusOnNode " << TOTAL_CORES_PER_NODE << "\n";
		myfile.close();
	}

    // Reset allocation for the instance
    sprintf(filename, ".hybrid/alloc%d", _apprankNum);
    std::ofstream myfile;
    myfile.open(filename);
    for (int i=0; i<numIranks; i++)
    {
        myfile << ((i==0)?48:1) << "\n";
    }
    myfile.close();

    _cores[0] = TOTAL_CORES_PER_NODE;
    va_end(args);

    rc = pthread_create(&_thread, NULL, instance_thread, (void *)this);

    if (rc) 
    {
        std::cout << "Unable to create thread\n";
        std::cout << "If it failed on login node, try ulimit -s 209715\n";
        assert(false);
    }
}


int main(void)
{
    system("rm -rf .hybrid/");
    system("mkdir -p .hybrid/");

    // Map mpi rank to node number
    std::vector<int>map {0,1,2,3, 0,1,2,3};
    std::ofstream myfile;
    myfile.open(".map");
    for (int node: map)
    {
        myfile << node << " ";
    }
    myfile << "\n";
    myfile.close();

	numNodes = 4;
	numAppranks = 4;
	assert(numNodes <= MAX_NODES);

    //                 apprankNum  totalLoad, numIranks   mpi ranks
    NanosInstance ins0(          0,        10,       2,  0,5);
    NanosInstance ins1(          1,        20,       2,  1,6);
    NanosInstance ins2(          2,         5,       2,  2,7);
    NanosInstance ins3(          3,        10,       2,  3,4);
#if 0
    ins0.wait();
    ins1.wait();
    ins2.wait();
    ins3.wait();
#endif
}

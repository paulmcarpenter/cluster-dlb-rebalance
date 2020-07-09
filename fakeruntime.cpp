#include <iostream>
#include <fstream>
#include <pthread.h>
#include <unistd.h>
#include <vector>

#include <stdarg.h>
#include <sys/stat.h>
#include <cassert>
#include <errno.h>

#define TOTAL_CORES_PER_NODE 48

class NanosInstance
{
    private:
        int _instanceNum;
        int _totalLoad;
        int _numRanks;
        pthread_t _thread;
        std::vector<int> _nodes; // for each i, which node is rank i on?
        std::vector<int> _cores; // for each i, how many cores does rank i have?

    public:

        NanosInstance(int instanceNum, int totalLoad, int numRanks, ...);

        void thread(void)
        {
            std::cout << "Thread for " << _instanceNum << "\n";
            for (int i=0; i<_numRanks; i++)
            {
                std::cout << " " << _instanceNum << " " << _nodes[i] << "\n";
            }

            for(int iter = 0; iter < 10; iter++)
            {
                // Calculate total number of cores
                int totalCores = 0;
                for (int i=0; i<_numRanks; i++)
                {
                    totalCores += _cores[i];
                }
                std::cout << _instanceNum << ": iter=" << iter << ", cores=" << totalCores << "\n";

                // Print theoretical load
                for (int rank=0; rank<_numRanks; rank++)
                {
                    double nodeLoad = (double)_totalLoad * _cores[rank] / totalCores;
                    char filename[1024];
                    sprintf(filename, ".balance/load-%d-%d", _instanceNum, rank);
                    std::ofstream myfile;
                    myfile.open(filename);
                    myfile << _cores[rank] << " " << nodeLoad;
                    myfile.close();
                }

                // Get the allocation of cores
                for (int rank=0; rank<_numRanks; rank++)
                {
                    char filename[1024];
                    sprintf(filename, ".balance/alloc-%d-%d", _instanceNum, rank);
                    std::ifstream myfile;
                    myfile.open(filename);
                    if (myfile.is_open())
                    {
                        int nodeAlloc;
                        myfile >> nodeAlloc;
                        if (nodeAlloc != _cores[rank])
                        {
                            std::cout << "new allocation of " << nodeAlloc << " cores for instance " << _instanceNum << " rank " << rank << "\n";
                            _cores[rank] = nodeAlloc;
                        }
                    }
                    myfile.close();
                }


                sleep(1);
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


NanosInstance::NanosInstance(int instanceNum, int totalLoad, int numRanks, ...) 
                            : _instanceNum(instanceNum), _totalLoad(totalLoad), 
                              _numRanks(numRanks), _nodes(numRanks), _cores(numRanks)
{
    int rc;

    va_list args;
    char hostname[1024];
    char buffer[1024];
    gethostname(hostname, 1024);

    va_start(args, numRanks);
    for (int i=0; i<numRanks; i++)
    {
        _nodes[i] = va_arg(args, int);
        _cores[i] = 0;
    }
    _cores[0] = TOTAL_CORES_PER_NODE;
    va_end(args);

    rc = pthread_create(&_thread, NULL, instance_thread, (void *)this);

    std::cout << "hello\n";
    if (rc) 
    {
        std::cout << "Unable to create thread " << rc << " (EAGAIN=" << EAGAIN << " EINVAL=" << EINVAL << " EPERM=" << EPERM << "\n";
        assert(false);
    }
}


int main(void)
{
    std::cout << "hello world\n";

    struct stat bal;
    if (stat(".balance/", &bal) == -1)
    {
        mkdir(".balance/", 0700);
    }


    //                 instanceNum  totalLoad, numRanks   nodes
    NanosInstance ins0(          0,        10,       2,  0,1);
    NanosInstance ins1(          1,        20,       2,  1,2);
    NanosInstance ins2(          2,         5,       2,  2,3);
    NanosInstance ins3(          3,        10,       2,  3,0);
#if 0
    ins0.wait();
    ins1.wait();
    ins2.wait();
    ins3.wait();
#endif
}

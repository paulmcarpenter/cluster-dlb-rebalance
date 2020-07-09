#include <iostream>
#include <pthread.h>
#include <unistd.h>

class NanosInstance
{
    private:
        int _instancenum;
        int _numranks;
        pthread_t _thread;

    public:

        NanosInstance(int instancenum, int numranks);

        void thread(void)
        {
            std::cout << "Thread for " << _instancenum << "\n";
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


NanosInstance::NanosInstance(int instancenum, int numranks) : _instancenum(instancenum), _numranks(numranks)
{
    int rc = pthread_create(&_thread, NULL, instance_thread, (void *)this);
    std::cout << "hello\n";
    if (rc) 
    {
        std::cout << "Unable to create thread\n";
    }
}


int main(void)
{
    std::cout << "hello world\n";
    NanosInstance ins0(0,3);
    NanosInstance ins1(1,2);
}

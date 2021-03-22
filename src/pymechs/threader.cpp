#include <pybind11/pybind11.h>
#include "threader.h"
#include "mechanism.h"

using namespace Mech;

Threader::Threader(int n_threads)
    :   n_threads_(n_threads),
        stopping_(false)
{
    if(n_threads <= 0)
        throw std::invalid_argument("Threader can only be initialized with a >0 number of threads");
    for(int i = 0;  i < n_threads_; i++) {
        threads_.emplace_back(&Threader::loop, this);
    }

}

Threader::~Threader()
{
    stopping_ = true;
    condition_.notify_all();

    for(std::thread& thread : threads_) {
        thread.join();
    }
}

void Threader::simulate(pybind11::list& mechanism_list, double t, int steps)
{
    for(auto mechanism : mechanism_list) {
        enqueue(std::bind(&Mechanism::simulate, mechanism.cast<Mechanism*>(), t, steps));
    }

    // Timed lock which waits for the queue to empty, prevents hangups when wait condition somehow fails to trigger
    std::unique_lock<std::mutex> lock(mutex_);
    while(!finishedSimulation(mechanism_list))  {
        finished_.wait_for(lock, std::chrono::milliseconds(1000));
    }
}

void Threader::evaluate(pybind11::list &mechanism_list, Objective& objective, double t, int steps)
{
    for(auto mechanism : mechanism_list) {
        enqueue(std::bind(&Mechanism::evaluate, mechanism.cast<Mechanism*>(), std::ref(objective), t, steps));
    }

    // Timed lock which waits for the queue to empty, prevents hangups when wait condition somehow fails to trigger
    std::unique_lock<std::mutex> lock(mutex_);
    while(!finishedEvaluation(mechanism_list))  {
        finished_.wait_for(lock, std::chrono::milliseconds(1000));
    }
}

bool Threader::finishedSimulation(const pybind11::list& mechanism_list) {

    bool finished = true;

    for(int i = 0; i < mechanism_list.size(); i++)
    {
        finished = finished && mechanism_list[mechanism_list.size()-i-1].cast<Mechanism*>()->isSimulated();
    }
    return finished;
}

bool Threader::finishedEvaluation(const pybind11::list& mechanism_list) {

    bool finished = true;

    for(int i = 0; i < mechanism_list.size(); i++)
    {
        finished = finished && mechanism_list[mechanism_list.size()-i-1].cast<Mechanism*>()->isEvaluated();
    }
    return finished;
}

void Threader::loop()
{
    std::function<void()> func;

    while(true) {
        {
            std::unique_lock<std::mutex> lock(mutex_);
            condition_.wait(lock, [this]() {return !queue_.empty() || stopping_;});
            if(stopping_ && queue_.empty()) {
                return;
            }
            func = queue_.front();
            queue_.pop();
        }
        func();
        if(!stopping_ && queue_.empty()) {
            finished_.notify_one();
        }
    }
}

void Threader::enqueue(std::function<void()> func)
{
    std::unique_lock<std::mutex> lock(mutex_);
    queue_.push(func);
    lock.unlock();
    condition_.notify_one();
}

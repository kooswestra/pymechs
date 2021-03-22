#ifndef MULTITHREAD_H
#define MULTITHREAD_H

#include <vector>
#include <thread>
#include <mutex>
#include <queue>
#include <condition_variable>
#include <pybind11/pytypes.h>
#include <Eigen/Dense>
#include "objective.h"


/** @brief Threader is a helper class to simulate mechanisms in a multithreaded manner
 *
 *  A simple multithreading implementation to provide significant performance gains on multicore
 *  CPU's. As entire populations of independent mechanisms need to be simulated threading
 *  easily provides a 3x + performance boost.
 */

class Threader
{
public:
    Threader(int n_threads);
    ~Threader();

    void simulate(pybind11::list& mechanism_list, double t, int steps);
    void evaluate(pybind11::list& mechanism_list, Mech::Objective& objective, double t, int steps);

private:
    int n_threads_;
    bool stopping_;
    std::vector<std::thread> threads_;
    std::queue<std::function<void()>> queue_;

    std::mutex mutex_;
    std::condition_variable condition_;
    std::condition_variable finished_;

    void loop();
    void enqueue(std::function<void()> func);

    bool finishedSimulation(const pybind11::list& mechanism_list);
    bool finishedEvaluation(const pybind11::list& mechanism_list);
};

#endif // MULTITHREAD_H

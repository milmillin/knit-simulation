#pragma once

#include <vector>
#include <thread>
#include <memory>
#include <mutex>
#include <queue>
#include <functional>
#include <type_traits>
#include <iostream>

namespace simulator {

class ParallelWorker
{
public:
  ParallelWorker(std::function<bool(void)> cancelled,
    size_t threads = std::thread::hardware_concurrency()) : m_numThreads(threads), m_cancelled(cancelled)
  {
    assert(threads != 0);
  }

  using Proc = std::function<void(void)>;

  void addWork(const Proc& f) noexcept
  {
    m_workQueue.push(f);
  }

  void run() {
    std::cout << "[ParallelWorker] Running " << m_workQueue.size() << " tasks with " << m_numThreads << " threads" << std::endl;
    reportInterval = m_workQueue.size() / 20;
    for (size_t i = 0; i < m_numThreads; ++i)
      m_threads.emplace_back(std::thread([this]() {
      while (!blockingQueueEmpty() && !m_cancelled())
      {
        Proc workItem = blockingPop();
        if (workItem) {
          workItem();
        }
        else {
          std::cout << "[ParallelWorker] Failed To Run" << std::endl;
        }
      }
        }));

    for (auto& thread : m_threads) {
      thread.join();
    }
  }

private:
  size_t m_numThreads;
  std::vector<std::thread> m_threads;
  std::queue<Proc> m_workQueue;
  mutable std::mutex m_queueLock;
  std::function<bool(void)> m_cancelled;
  int reportInterval;

  const Proc blockingPop() {
    std::lock_guard<std::mutex> lock(m_queueLock);
    if (m_workQueue.size() % reportInterval == 0) {
      std::cout << "[ParallelWorker] " << m_workQueue.size() << " tasks remain" << std::endl;
    }
    const Proc res = m_workQueue.front();
    m_workQueue.pop();
    return res;
  }

  bool blockingQueueEmpty() {
    std::lock_guard<std::mutex> lock(m_queueLock);
    return m_workQueue.empty();
  }
};

}  // namespace simulator
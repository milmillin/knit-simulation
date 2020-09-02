#pragma once

#include "spdlog/spdlog.h"

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
    SPDLOG_INFO("[ParallelWorker] Running {}  tasks with {} threads", m_workQueue.size(), m_numThreads);
    reportInterval = m_workQueue.size() / 20;
    for (size_t i = 0; i < m_numThreads; ++i)
      m_threads.emplace_back(std::thread([this]() {
      Proc workItem;
      while (blockingPop(workItem) && !m_cancelled())
      {
        if (workItem) {
          workItem();
        }
        else {
          SPDLOG_INFO("[ParallelWorker] Failed To Run");
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

  const bool blockingPop(Proc& res) {
    std::lock_guard<std::mutex> lock(m_queueLock);
    if (!m_workQueue.empty()) {
      if (m_workQueue.size() % reportInterval == 0) {
        SPDLOG_INFO("[ParallelWorker] {} tasks remain", m_workQueue.size());
      }

      res = m_workQueue.front();
      m_workQueue.pop();
      return true;
    }
    return false;
  }
};

}  // namespace simulator
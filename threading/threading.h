#pragma once

#include <functional>

#include "./ctpl_stl.h"

namespace threading {

typedef std::function<void(int, ctpl::thread_pool*)> TaskProducer;
typedef std::function<void(int, int, int)> SequencialTaskConsumer;
typedef std::function<void(int, size_t)> Task;

void submitProducerAndWait(ctpl::thread_pool& thread_pool,
  TaskProducer producer);
void runSequentialJob(ctpl::thread_pool& thread_pool,
  SequencialTaskConsumer consumer,
  int start, int end, int step);
void runSequentialJob(ctpl::thread_pool& thread_pool,
  SequencialTaskConsumer consumer,
  int start, int end);
void runParallelFor(ctpl::thread_pool& thread_pool, size_t start, size_t end,
  const Task& task);

}  // namespace threading

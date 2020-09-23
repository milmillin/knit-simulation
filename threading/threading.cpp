#include "./threading.h"

#include <functional>
#include <iostream>
#include <algorithm>
#include <cmath>

#include "./ctpl_stl.h"

namespace threading {

void submitProducerAndWait(ctpl::thread_pool &thread_pool,
                           TaskProducer producer) {
  auto producer_result = thread_pool.push(producer, &thread_pool);
  producer_result.get();
  thread_pool.wait();
}

static void sequencialJobProducer(int thread_id, ctpl::thread_pool *thread_pool,
    SequencialTaskConsumer consumer,
    int start, int end, int step) {
  for (int i = start; i < end; i += step) {
    thread_pool->push(consumer, i, std::min(end, i + step));
  }
}

void runSequentialJob(ctpl::thread_pool &thread_pool,
                      SequencialTaskConsumer consumer,
                      int start, int end, int step) {
  using namespace std::placeholders;

  auto producer = std::bind(sequencialJobProducer, _1, _2,
    consumer, start, end, step);
  submitProducerAndWait(thread_pool, producer);
}

void runSequentialJob(ctpl::thread_pool &thread_pool,
                      SequencialTaskConsumer consumer,
                      int start, int end) {
  runSequentialJob(thread_pool, consumer,
    start, end, 1 + ((end - start) / thread_pool.size()));
}

static void loopTask(int thread_id, size_t start, size_t end, const Task& task) {
  for (size_t i = start; i < end; i++) {
    task(thread_id, i);
  }
}

static void loopTaskProducer(int thread_id, ctpl::thread_pool* thread_pool, size_t start, size_t end, const Task& task) {
  size_t step = std::ceil((double)(end - start) / thread_pool->size());
  for (size_t i = start; i < end; i += step) {
    thread_pool->push(loopTask, i, std::min(end, i + step), std::ref(task));
  }
}

void runParallelFor(ctpl::thread_pool& thread_pool, size_t start, size_t end, const Task& task)
{
  using namespace std::placeholders;
  auto producer = std::bind(loopTaskProducer, _1, _2, start, end, std::ref(task));
  submitProducerAndWait(thread_pool, producer);
}

}  // namespace threading

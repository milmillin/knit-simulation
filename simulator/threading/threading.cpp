#include "./threading.h"

#include <functional>
#include <iostream>
#include <algorithm>

#include "./ctpl_stl.h"

namespace simulator {
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
  std::cout << "job submission done" << std::endl;
}

void runSequentialJob(ctpl::thread_pool &thread_pool,
                      SequencialTaskConsumer consumer,
                      int start, int end, int step) {
  using namespace std::placeholders;

  auto producer = std::bind(sequencialJobProducer, _1, _2,
    consumer, start, end, step);
  submitProducerAndWait(thread_pool, producer);
}

}  // namespace threading
}  // namespace simulator

#pragma once

#include <functional>

#include "./ctpl_stl.h"

namespace simulator {
namespace threading {

typedef std::function<void(int, ctpl::thread_pool*)> TaskProducer;
typedef std::function<void(int, int, int)> SequencialTaskConsumer;

void submitProducerAndWait(ctpl::thread_pool &thread_pool,
                           TaskProducer producer);

void runSequentialJob(ctpl::thread_pool &thread_pool,
                      SequencialTaskConsumer consumer,
                      int start, int end, int step);
void runSequentialJob(ctpl::thread_pool &thread_pool,
                      SequencialTaskConsumer consumer,
                      int start, int end);

}  // namespace threading
}  // namespace simulator

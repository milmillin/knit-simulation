#define USE_EASY_PROFILER

#ifdef USE_EASY_PROFILER
  #include <easy/profiler.h>
#else
  // Define empty macros for easy profiler when it's disabled
  #define EASY_PROFILER_ENABLE
  #define EASY_FUNCTION()
  #define EASY_BLOCK(name)
  #define EASY_END_BLOCK
#endif
#include <vector>
#include <unordered_map>
#include <mutex>

namespace simulator {

// A 2D lazy-initialized cache array
// T: cache type
template<class T>
class CacheArray2D {
 public:
  // sizeX: size in x-axis
  // sizeY: size in y-axis
  CacheArray2D(int sizeX, int sizeY) :
    cacheData(sizeX),
    locks(sizeY) {}

  // Return a reference to the cached element,
  // or a default-constructed element if the element doesn't exist.
  // `unlock` should be called when modifications to the element is done
  // x: x-axis position
  // y: y-axis position
  inline T &lock(int x, int y) {
    locks[x].lock();
    return cacheData[x][y];
  }

  // Release the element to other threads
  // x: x-axis position
  // y: y-axis position
  inline void unlock(int x, int y) {
    locks[x].unlock();
  }

 private:
  // `cacheData[x][y]` is the data
  std::vector<std::unordered_map<int, T> > cacheData;
  // `locks[x]` protects `cacheData[x]`
  std::vector<std::mutex> locks;
};

}

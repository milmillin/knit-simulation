#include <vector>
#include <unordered_map>
#include <mutex>

namespace simulator {

template<class T>
class Cache {
 private:
  std::vector<std::unordered_map<int, T> > cacheData;
  std::vector<std::mutex> locks;
 public:
  Cache(int sizeX, int sizeY) :
    cacheData(sizeX),
    locks(sizeY) {}

  T &lock(int x, int y) {
    locks[x].lock();
    return cacheData[x][y];
  }

  void unlock(int x, int y) {
    locks[x].unlock();
  }
};

}

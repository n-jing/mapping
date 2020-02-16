#ifndef HASH_KEY_JJ_H
#define HASH_KEY_JJ_H

#include <unordered_map>
#include <algorithm>


template <class T>
inline void hash_combine(std::size_t &seed, const T &v) {
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

template<typename T, size_t d>
struct HashFunc
{
  std::size_t operator()(const T &key) const
    {
      size_t value = 0;
      for (size_t i = 0; i < d; ++i)
        hash_combine(value, key[i]);

      return value;
    }
};

template<typename T, size_t d>
struct UnorderedHashFunc
{
  std::size_t operator()(const T &key) const
    {
      T k = key;
      std::sort(k.begin(), k.end());
  
      size_t value = 0;
      for (size_t i = 0; i < d; ++i)
        hash_combine(value, k[i]);
      return value;
    }
};


template<typename T, size_t d>
struct EqualKey
{
  bool operator()(const T &lhs, const T &rhs) const
  {
    for (size_t i = 0; i < d; ++i)
      if (lhs[i] != rhs[i])
        return false;

    return true;
  }
};


template<typename T, size_t d>
struct UnorderedEqualKey
{
  bool operator()(const T &lhs, const T &rhs) const
    {
      T lower = lhs;
      T upper = rhs;
      std::sort(lower.begin(), lower.end());
      std::sort(upper.begin(), upper.end());
  
      for (size_t i = 0; i < d; ++i)
      {
        if (lower[i] != upper[i])
          return false;
      }

      return true;
    }
};




#endif // HASH_KEY_JJ_H

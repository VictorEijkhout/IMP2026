#pragma once
#include "coordinate.hpp"

// The box of all p with P <= p < Q (componentwise).
template <std::size_t d>
class domain {
public:
  using coord = coordinate<d>;

  domain() = default;
  domain(const coord& P, const coord& Q) : lo_(P), hi_(Q) {}

  const coord& first() const { return lo_; }  // P (inclusive)
  const coord& last() const { return hi_; }   // Q (exclusive)

  bool contains(const coord& p) const { return p >= lo_ && p < hi_; }

  // Extent per dimension (0 where the box is empty in that dimension).
  coord extent() const {
    coord e;
    for (std::size_t i = 0; i < d; ++i)
      e[i] = hi_[i] > lo_[i] ? hi_[i] - lo_[i] : 0;
    return e;
  }
  long long size() const { return extent().product(); }
  bool empty() const { return size() == 0; }

  bool operator==(const domain& o) const {
    return (empty() && o.empty()) || (lo_ == o.lo_ && hi_ == o.hi_);
  }

  // Call f(coord) for every point in the domain, last index fastest.
  template <typename F>
  void for_each(F&& f) const {
    if (empty()) return;
    coord p = lo_;
    while (true) {
      f(p);
      std::size_t i = d;
      while (i-- > 0) {
        if (++p[i] < hi_[i]) break;
        p[i] = lo_[i];
        if (i == 0) return;
      }
    }
  }

private:
  coord lo_, hi_;
};

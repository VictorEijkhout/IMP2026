#pragma once
#include <stdexcept>
#include "domain.hpp"

// Block distribution of the index space [0,N) over a processor grid [0,P).
// In each dimension i, N[i] points are split over P[i] processors as evenly
// as possible: processor k owns [k*N/P, (k+1)*N/P) (integer division), so
// block sizes differ by at most one, and the blocks tile [0,N[i]) exactly.
template <std::size_t d>
class decomposition {
public:
  using coord = coordinate<d>;

  decomposition(const coord& N, const coord& P) : N_(N), P_(P) {
    for (std::size_t i = 0; i < d; ++i)
      if (P[i] < 1 || N[i] < 0)
        throw std::invalid_argument("need P>=1 and N>=0 in every dimension");
  }

  const coord& global_size() const { return N_; }
  const coord& processors() const { return P_; }
  domain<d> global_domain() const { return domain<d>(coord(0), N_); }

  // The sub domain of processor p (requires coord(0) <= p < P).
  domain<d> operator[](const coord& p) const {
    if (!(p >= coord(0) && p < P_))
      throw std::out_of_range("processor coordinate outside grid");
    coord lo, hi;
    for (std::size_t i = 0; i < d; ++i) {
      lo[i] = p[i] * N_[i] / P_[i];
      hi[i] = (p[i] + 1) * N_[i] / P_[i];
    }
    return domain<d>(lo, hi);
  }

  // The processor owning global point n (requires coord(0) <= n < N).
  coord owner(const coord& n) const {
    if (!global_domain().contains(n))
      throw std::out_of_range("point outside global domain");
    coord p;
    for (std::size_t i = 0; i < d; ++i) {
      // largest k with k*N/P <= n  <=>  k = floor(((n+1)*P - 1) / N)
      p[i] = ((n[i] + 1) * P_[i] - 1) / N_[i];
    }
    return p;
  }

  // Visit every processor coordinate.
  template <typename F>
  void for_each_processor(F&& f) const {
    domain<d>(coord(0), P_).for_each(f);
  }

private:
  coord N_, P_;
};

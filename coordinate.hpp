#pragma once
#include <array>
#include <cstddef>
#include <initializer_list>
#include <format>

// A point in d-dimensional integer space.
template <std::size_t d>
class coordinate {
public:
  using value_type = long long;

  coordinate() { data_.fill(0); }
  explicit coordinate(value_type v) { data_.fill(v); }  // all entries equal v
  coordinate(std::initializer_list<value_type> init) {
    data_.fill(0);
    std::size_t i = 0;
    for (auto v : init) if (i < d) data_[i++] = v;
  }

  static constexpr std::size_t dimension() { return d; }
  value_type& operator[](std::size_t i) { return data_[i]; }
  const value_type& operator[](std::size_t i) const { return data_[i]; }

  // Product of entries (number of points in a box of this extent).
  value_type product() const {
    value_type r = 1;
    for (auto v : data_) r *= v;
    return r;
  }

  bool operator==(const coordinate& o) const { return data_ == o.data_; }
  bool operator!=(const coordinate& o) const { return !(*this == o); }

  // Strict partial order: p < q iff EVERY component is strictly less.
  // Note: therefore !(p<q) does not imply p>=q.
  bool operator<(const coordinate& q) const {
    for (std::size_t i = 0; i < d; ++i)
      if (!(data_[i] < q.data_[i])) return false;
    return true;
  }
  bool operator>(const coordinate& q) const { return q < *this; }
  // Componentwise: all components >=.
  bool operator>=(const coordinate& q) const {
    for (std::size_t i = 0; i < d; ++i)
      if (!(data_[i] >= q.data_[i])) return false;
    return true;
  }
  bool operator<=(const coordinate& q) const { return q >= *this; }

private:
  std::array<value_type, d> data_;
};

// std::format / std::print support: format("{}", coordinate<2>{1,2}) -> "(1,2)".
// A standard integer format spec applies to each component: "{:03}" -> "(001,002)".
template <std::size_t d>
struct std::formatter<coordinate<d>> {
  std::formatter<long long> element;

  constexpr auto parse(std::format_parse_context& ctx) { return element.parse(ctx); }

  template <typename Out>
  auto format(const coordinate<d>& c, std::basic_format_context<Out, char>& ctx) const {
    auto out = ctx.out();
    *out++ = '(';
    for (std::size_t i = 0; i < d; ++i) {
      if (i) *out++ = ',';
      ctx.advance_to(out);
      out = element.format(c[i], ctx);
    }
    *out++ = ')';
    return out;
  }
};

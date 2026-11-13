// Catch2 (v3) unit tests, run for d = 1, 2, 3.
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <format>
#include <map>
#include <stdexcept>
#include "decomposition.hpp"

namespace {
template <std::size_t d>
using key_t = std::array<long long, d>;

template <std::size_t d>
key_t<d> key(const coordinate<d>& c) {
  key_t<d> k;
  for (std::size_t i = 0; i < d; ++i) k[i] = c[i];
  return k;
}

long long ipow(long long b, std::size_t e) {
  long long r = 1;
  while (e--) r *= b;
  return r;
}

// The defining property of a decomposition of [0,N) over processors [0,P):
// sub domains are disjoint, lie inside [0,N), and cover it exactly.
template <std::size_t d>
void check_decomposition(const coordinate<d>& N, const coordinate<d>& P) {
  using C = coordinate<d>;
  decomposition<d> dec(N, P);
  domain<d> global(C(0), N);
  std::map<key_t<d>, int> count;
  long long total = 0;

  dec.for_each_processor([&](const C& p) {
    domain<d> sub = dec[p];
    total += sub.size();
    sub.for_each([&](const C& n) {
      REQUIRE(global.contains(n));
      REQUIRE(dec.owner(n) == p);
      ++count[key(n)];
    });
    for (std::size_t i = 0; i < d; ++i) {  // balanced: floor or ceil of N/P
      long long lo = N[i] / P[i], hi = (N[i] + P[i] - 1) / P[i];
      REQUIRE(sub.extent()[i] >= lo);
      REQUIRE(sub.extent()[i] <= hi);
    }
  });

  REQUIRE(total == global.size());
  REQUIRE(static_cast<long long>(count.size()) == global.size());
  for (auto& kv : count) REQUIRE(kv.second == 1);
}
}  // namespace

template <std::size_t d>
using dim = std::integral_constant<std::size_t, d>;
using dim1 = dim<1>;
using dim2 = dim<2>;
using dim3 = dim<3>;
#define DIMS dim1, dim2, dim3

// ---- coordinate ----------------------------------------------------------
TEMPLATE_TEST_CASE("coordinate construction", "[coordinate]", DIMS) {
  constexpr std::size_t d = TestType::value;
  using C = coordinate<d>;
  CHECK(C::dimension() == d);
  CHECK(C() == C(0));
  CHECK(C(3)[d - 1] == 3);
  CHECK(C(5).product() == ipow(5, d));
  CHECK(C(0).product() == 0);

  C i{1, 2, 3};  // extra values ignored, missing ones are 0
  CHECK(i[0] == 1);
  if constexpr (d >= 3) CHECK(i[2] == 3);
}

TEMPLATE_TEST_CASE("coordinate ordering", "[coordinate]", DIMS) {
  constexpr std::size_t d = TestType::value;
  using C = coordinate<d>;
  C a(1), b(2);

  CHECK(a < b);
  CHECK_FALSE(b < a);
  CHECK_FALSE(a < a);  // strict
  CHECK(b > a);
  CHECK((a <= b && a >= a && a <= a));
  CHECK(a != b);

  if constexpr (d > 1) {
    C c(1);
    c[0] = 0;  // equal in all but one component
    CHECK_FALSE(c < a);
    CHECK_FALSE(a < c);
    CHECK(c <= a);

    C e = b;
    e[d - 1] = 0;  // smaller in one component, larger in the others
    CHECK_FALSE(e < a);
    CHECK_FALSE(a < e);
    CHECK_FALSE(e >= a);
  }
}

// ---- domain --------------------------------------------------------------
TEMPLATE_TEST_CASE("domain membership and size", "[domain]", DIMS) {
  constexpr std::size_t d = TestType::value;
  using C = coordinate<d>;
  domain<d> D(C(2), C(5));

  CHECK(D.first() == C(2));
  CHECK(D.last() == C(5));
  CHECK(D.extent() == C(3));
  CHECK(D.size() == ipow(3, d));
  CHECK_FALSE(D.empty());

  CHECK(D.contains(C(2)));  // lower bound inclusive
  CHECK(D.contains(C(4)));
  CHECK_FALSE(D.contains(C(5)));  // upper bound exclusive
  CHECK_FALSE(D.contains(C(1)));
  CHECK_FALSE(D.contains(C(7)));
  for (std::size_t i = 0; i < d; ++i) {  // outside in one dimension only
    C q(3);
    q[i] = 5;
    CHECK_FALSE(D.contains(q));
    q[i] = 1;
    CHECK_FALSE(D.contains(q));
  }
}

TEMPLATE_TEST_CASE("domain iteration", "[domain]", DIMS) {
  constexpr std::size_t d = TestType::value;
  using C = coordinate<d>;
  domain<d> D(C(2), C(5));

  std::map<key_t<d>, int> seen;
  D.for_each([&](const C& p) {
    CHECK(D.contains(p));
    ++seen[key(p)];
  });
  CHECK(static_cast<long long>(seen.size()) == D.size());
  for (auto& kv : seen) CHECK(kv.second == 1);
}

TEMPLATE_TEST_CASE("empty domains", "[domain]", DIMS) {
  constexpr std::size_t d = TestType::value;
  using C = coordinate<d>;
  domain<d> E(C(3), C(3));
  CHECK(E.empty());
  CHECK(E.size() == 0);
  CHECK_FALSE(E.contains(C(3)));

  int visits = 0;
  E.for_each([&](const C&) { ++visits; });
  CHECK(visits == 0);

  domain<d> inverted(C(5), C(2));
  CHECK(inverted.empty());
  CHECK(E == inverted);
}

// ---- decomposition -------------------------------------------------------
TEMPLATE_TEST_CASE("decomposition tiles the global domain", "[decomposition]", DIMS) {
  constexpr std::size_t d = TestType::value;
  using C = coordinate<d>;

  SECTION("even split") { check_decomposition<d>(C(12), C(3)); }
  SECTION("uneven split") { check_decomposition<d>(C(10), C(3)); }
  SECTION("one point per processor") { check_decomposition<d>(C(7), C(7)); }
  SECTION("single processor") { check_decomposition<d>(C(5), C(1)); }
  SECTION("more processors than points") { check_decomposition<d>(C(3), C(5)); }
  SECTION("empty global domain") { check_decomposition<d>(C(0), C(2)); }
  SECTION("anisotropic") {
    C N(8), P(2);
    N[0] = 11; P[0] = 4;
    N[d - 1] = 5; P[d - 1] = 3;  // for d=1 this overrides dimension 0
    check_decomposition<d>(N, P);
  }
}

TEMPLATE_TEST_CASE("decomposition concrete blocks", "[decomposition]", DIMS) {
  constexpr std::size_t d = TestType::value;
  using C = coordinate<d>;
  decomposition<d> dec(C(10), C(3));  // dimension 0: [0,3) [3,6) [6,10)

  C p0(0), p1(0), p2(0);
  p1[0] = 1;
  p2[0] = 2;
  CHECK(dec[p0].first()[0] == 0);
  CHECK(dec[p0].last()[0] == 3);
  CHECK(dec[p1].first()[0] == 3);
  CHECK(dec[p1].last()[0] == 6);
  CHECK(dec[p2].first()[0] == 6);
  CHECK(dec[p2].last()[0] == 10);
  CHECK(dec.global_domain() == domain<d>(C(0), C(10)));
  CHECK(dec.owner(C(9)) == C(2));
}

TEMPLATE_TEST_CASE("decomposition rejects bad input", "[decomposition]", DIMS) {
  constexpr std::size_t d = TestType::value;
  using C = coordinate<d>;
  decomposition<d> dec(C(10), C(3));

  CHECK_THROWS_AS(dec[C(3)], std::out_of_range);    // processor outside grid
  CHECK_THROWS_AS(dec.owner(C(10)), std::out_of_range);  // point outside domain
  CHECK_THROWS_AS(decomposition<d>(C(10), C(0)), std::invalid_argument);
  CHECK_THROWS_AS(decomposition<d>(C(-1), C(2)), std::invalid_argument);
}

TEMPLATE_TEST_CASE("coordinate formatting", "[coordinate]", DIMS) {
  constexpr std::size_t d = TestType::value;
  using C = coordinate<d>;
  std::string plain = "(", padded = "(";
  for (std::size_t i = 0; i < d; ++i) {
    plain += (i ? ",7" : "7");
    padded += (i ? ",007" : "007");
  }
  CHECK(std::format("{}", C(7)) == plain + ")");
  CHECK(std::format("{:03}", C(7)) == padded + ")");
}

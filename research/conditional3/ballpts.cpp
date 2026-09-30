// Max occupancy of projective-distance balls by single-qubit Clifford+T words of T-count <= tau.
//
// Words are enumerated in Matsumoto-Amano normal form (T|e)(HT|SHT)^* C (C one of the 24 Cliffords
// mod phase), so every element of Lambda_tau appears exactly once; |Lambda_tau| = 72*2^tau - 48.
// Each word is stored as a unit quaternion mod sign (SU(2) representative), and
//   dproj(A,B) = min_phase ||A - e^{i phi} B|| = min(|qA - qB|, |qA + qB|).
// Symmetry: Lambda_tau is invariant under right multiplication by the Clifford group G, so the
// ball-count function f(xi) = #(Lambda_tau cap B(xi,rho)) satisfies f(xi g) = f(xi). Hence
// max_xi f = max over xi in the Voronoi cell F of 1 for the right G-action. For each normal-form
// prefix P we keep the translates P*C that lie within a margin of F (one of them lies in F and is
// used as a centre) and count, for every centre, the kept points within dproj <= rho.
// The count at a data point is (for a Poisson process) distributed as 1 + Poisson(mu), where
// mu = N_tau * Haar(B(rho)).
//
// Usage: balls tau mu1 [mu2 ...] [--random seed]
//   --random: replace the words by Haar-random prefixes P (same number, same right-G symmetrization),
//             as a control with the same statistic and the same symmetry.
#include <bits/stdc++.h>
#include <omp.h>
#include <parallel/algorithm>
using namespace std;
typedef complex<double> cd;
struct SU { cd a, b; };  // [[a, b], [-conj b, conj a]]
static inline SU mul(const SU& p, const SU& q) { return {p.a * q.a - p.b * conj(q.b), p.a * q.b + p.b * conj(q.a)}; }
static inline SU inv(const SU& p) { return {conj(p.a), -p.b}; }
static inline double ip(const SU& p, const SU& q) {  // quaternion inner product = Re tr(p^dag q)/2
  return real(p.a) * real(q.a) + imag(p.a) * imag(q.a) + real(p.b) * real(q.b) + imag(p.b) * imag(q.b);
}
struct Pt { double x[4]; uint64_t key; uint64_t label; unsigned char centre; };

int TAU;
vector<SU> G;  // 24 Cliffords mod sign
SU A_, B_, Tg;
double MARGIN, CELL;

[[maybe_unused]] static inline void emit_unused(const SU& P, uint64_t label, vector<Pt>& out) {
  double v[24], vmax = -1; int arg = -1;
  for (int h = 0; h < 24; h++) { v[h] = fabs(ip(P, G[h])); if (v[h] > vmax) { vmax = v[h]; arg = h; } }
  for (int h = 0; h < 24; h++) {
    if (v[h] < vmax - MARGIN) continue;
    SU x = mul(P, inv(G[h]));  // x = P * C with C = G[h]^{-1}; <x,1> = <P,G[h]>
    Pt p;
    double s = (real(x.a) >= 0) ? 1.0 : -1.0;
    p.x[0] = s * real(x.a); p.x[1] = s * imag(x.a); p.x[2] = s * real(x.b); p.x[3] = s * imag(x.b);
    p.label = (label << 5) | (uint64_t)h;
    p.centre = (h == arg);
    out.push_back(p);
  }
}
SU XI; double RHO;
static inline void test(const SU& P, uint64_t label, int tc, vector<string>& out) {
  for (int h = 0; h < 24; h++) {
    SU x = mul(P, inv(G[h]));
    double d = sqrt(max(0.0, 2 - 2 * fabs(ip(x, XI))));
    if (d <= RHO) {
      double s = (ip(x, XI) >= 0) ? 1.0 : -1.0;
      char buf[256];
      snprintf(buf, sizeof buf, "%d %llu %d %.17g %.17g %.17g %.17g %.6g", tc, (unsigned long long)((label << 5) | h), h,
               s * real(x.a), s * imag(x.a), s * real(x.b), s * imag(x.b), d);
      out.push_back(buf);
    }
  }
}
static void dfs2(const SU& p, int d, uint64_t path, vector<string>& out) {
  test(p, path | ((uint64_t)d << 40), d, out);
  if (d + 1 <= TAU) {
    test(mul(Tg, p), path | ((uint64_t)d << 40) | (1ULL << 46), d + 1, out);
    dfs2(mul(A_, p), d + 1, path | (0ULL << d), out);
    dfs2(mul(B_, p), d + 1, path | (1ULL << d), out);
  }
}
int main(int argc, char** argv) {
  if (argc < 7) { fprintf(stderr, "usage: ballpts tau q0 q1 q2 q3 rho\n"); return 1; }
  TAU = atoi(argv[1]);
  double q[4]; for (int j = 0; j < 4; j++) q[j] = atof(argv[2 + j]);
  double n = sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);
  XI = {cd(q[0]/n, q[1]/n), cd(q[2]/n, q[3]/n)}; RHO = atof(argv[6]);
  const double r2 = sqrt(0.5);
  SU H = {cd(0, -r2), cd(0, -r2)};
  SU S = {polar(1.0, -M_PI / 4), cd(0, 0)};
  Tg = {polar(1.0, -M_PI / 8), cd(0, 0)};
  A_ = mul(H, Tg); B_ = mul(S, mul(H, Tg));
  {
    auto canon = [](const SU& x) {
      double q[4] = {real(x.a), imag(x.a), real(x.b), imag(x.b)};
      int i = 0; while (fabs(q[i]) < 1e-9) i++;
      double s = q[i] > 0 ? 1 : -1;
      array<long long, 4> k; for (int j = 0; j < 4; j++) k[j] = llround(s * q[j] * 1e6); return k;
    };
    map<array<long long, 4>, SU> seen; SU I = {cd(1, 0), cd(0, 0)};
    seen[canon(I)] = I; vector<SU> fr = {I};
    while (!fr.empty()) {
      vector<SU> nf;
      for (auto& m : fr) for (auto g : {H, S}) { SU nn = mul(g, m); auto k = canon(nn); if (!seen.count(k)) { seen[k] = nn; nf.push_back(nn); } }
      fr = nf;
    }
    for (auto& kv : seen) G.push_back(kv.second);
  }
  int split = min(TAU, 8);
  vector<pair<SU, pair<int, uint64_t>>> roots; vector<string> shallow;
  function<void(const SU&, int, uint64_t)> go = [&](const SU& p, int d, uint64_t path) {
    if (d == split) { roots.push_back({p, {d, path}}); return; }
    test(p, path | ((uint64_t)d << 40), d, shallow);
    if (d + 1 <= TAU) {
      test(mul(Tg, p), path | ((uint64_t)d << 40) | (1ULL << 46), d + 1, shallow);
      go(mul(A_, p), d + 1, path | (0ULL << d));
      go(mul(B_, p), d + 1, path | (1ULL << d));
    }
  };
  go(SU{cd(1, 0), cd(0, 0)}, 0, 0);
  vector<vector<string>> parts(roots.size());
#pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < roots.size(); i++) dfs2(roots[i].first, roots[i].second.first, roots[i].second.second, parts[i]);
  printf("# tcount label cliff q0 q1 q2 q3 dist\n");
  for (auto& l : shallow) printf("%s\n", l.c_str());
  for (auto& v : parts) for (auto& l : v) printf("%s\n", l.c_str());
  return 0;
}

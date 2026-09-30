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
//   --symrandom seed (review addition): Haar-random P, symmetrized under the FULL symmetry group of Lambda_tau,
//             i.e. left and right Cliffords and inversion (1152 = 24*24*2 images per P), same total number of points.
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

static inline void emit(const SU& P, uint64_t label, vector<Pt>& out) {
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
// label: bits [0..39] path, [40..45] depth, [46] leading T
static void dfs(const SU& p, int d, uint64_t path, vector<Pt>& out) {
  emit(p, path | ((uint64_t)d << 40), out);
  if (d + 1 <= TAU) {
    emit(mul(Tg, p), path | ((uint64_t)(d) << 40) | (1ULL << 46), out);
    dfs(mul(A_, p), d + 1, path | (0ULL << d), out);
    dfs(mul(B_, p), d + 1, path | (1ULL << d), out);
  }
}

int main(int argc, char** argv) {
  if (argc < 3) { fprintf(stderr, "usage: balls tau mu1 [mu2..] [--random seed]\n"); return 1; }
  TAU = atoi(argv[1]);
  vector<double> mus; bool rnd = false, symr = false; uint64_t seed = 1;
  for (int i = 2; i < argc; i++) {
    if (!strcmp(argv[i], "--random")) { rnd = true; seed = strtoull(argv[++i], 0, 10); }
    else if (!strcmp(argv[i], "--symrandom")) { rnd = true; symr = true; seed = strtoull(argv[++i], 0, 10); }
    else mus.push_back(atof(argv[i]));
  }
  if (mus.size() > 16) { fprintf(stderr, "at most 16 radii\n"); return 1;
  }
  const double r2 = sqrt(0.5);
  SU H = {cd(0, -r2), cd(0, -r2)};                       // -i H
  SU S = {polar(1.0, -M_PI / 4), cd(0, 0)};             // e^{-i pi/4} S
  Tg = {polar(1.0, -M_PI / 8), cd(0, 0)};               // e^{-i pi/8} T
  A_ = mul(H, Tg); B_ = mul(S, mul(H, Tg));
  // Clifford group mod sign by closure
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
      for (auto& m : fr) for (auto g : {H, S}) { SU n = mul(g, m); auto k = canon(n); if (!seen.count(k)) { seen[k] = n; nf.push_back(n); } }
      fr = nf;
    }
    for (auto& kv : seen) G.push_back(kv.second);
    if (G.size() != 24) { fprintf(stderr, "clifford count %zu\n", G.size()); return 1; }
  }
  const double Ntau = 72.0 * ldexp(1.0, TAU) - 48;
  // Haar measure of a dproj-ball of radius rho: two antipodal caps of angle th = 2 asin(rho/2) on S^3,
  // each of normalized volume (th - sin th cos th)/pi.
  auto haar = [](double rho) { double th = 2 * asin(rho / 2); return 2 * (th - sin(th) * cos(th)) / M_PI; };
  vector<double> rhos;
  for (double mu : mus) {  // solve Ntau * haar(rho) = mu
    double lo = 0, hi = 1.4;
    for (int it = 0; it < 200; it++) { double md = (lo + hi) / 2; if (Ntau * haar(md) < mu) lo = md; else hi = md; }
    rhos.push_back(hi);
  }
  double rmax = *max_element(rhos.begin(), rhos.end());
  MARGIN = 2 * rmax + 1e-12;
  CELL = rmax;
  double t0 = omp_get_wtime();
  // enumerate
  vector<vector<Pt>> parts;
  if (!rnd) {
    int split = min(TAU, 8);
    vector<pair<SU, pair<int, uint64_t>>> roots;  // nodes at depth `split` (plus shallower emits done here)
    vector<Pt> shallow;
    function<void(const SU&, int, uint64_t)> go = [&](const SU& p, int d, uint64_t path) {
      if (d == split) { roots.push_back({p, {d, path}}); return; }
      emit(p, path | ((uint64_t)d << 40), shallow);
      if (d + 1 <= TAU) {
        emit(mul(Tg, p), path | ((uint64_t)d << 40) | (1ULL << 46), shallow);
        go(mul(A_, p), d + 1, path | (0ULL << d));
        go(mul(B_, p), d + 1, path | (1ULL << d));
      }
    };
    go(SU{cd(1, 0), cd(0, 0)}, 0, 0);
    parts.resize(roots.size() + 1);
    parts.back() = shallow;
#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < roots.size(); i++) dfs(roots[i].first, roots[i].second.first, roots[i].second.second, parts[i]);
  } else {
    long long M = 3LL * (1LL << TAU) - 2;  // same number of prefixes as the words
    int nt = omp_get_max_threads(); parts.resize(nt);
#pragma omp parallel
    {
      int id = omp_get_thread_num();
      mt19937_64 rng(seed * 1000003ULL + id);
      normal_distribution<double> nd(0, 1);
      if (!symr) {
      for (long long i = id; i < M; i += nt) {
        double q[4]; double n = 0; for (int j = 0; j < 4; j++) { q[j] = nd(rng); n += q[j] * q[j]; }
        n = sqrt(n); SU P = {cd(q[0] / n, q[1] / n), cd(q[2] / n, q[3] / n)};
        emit(P, (uint64_t)i, parts[id]);
      }
      } else {
      long long M48 = (M + 47) / 48;  // each random P yields 48 left-Clifford/inverse prefixes
      for (long long i = id; i < M48; i += nt) {
        double q[4]; double n = 0; for (int j = 0; j < 4; j++) { q[j] = nd(rng); n += q[j] * q[j]; }
        n = sqrt(n); SU P = {cd(q[0] / n, q[1] / n), cd(q[2] / n, q[3] / n)};
        for (int iv = 0; iv < 2; iv++) { SU R = iv ? inv(P) : P;
          for (int h = 0; h < 24; h++) emit(mul(G[h], R), (uint64_t)((i * 2 + iv) * 24 + h), parts[id]); }
      }
      }
    }
  }
  size_t tot = 0; for (auto& v : parts) tot += v.size();
  vector<Pt> pts; pts.reserve(tot);
  for (auto& v : parts) { pts.insert(pts.end(), v.begin(), v.end()); vector<Pt>().swap(v); }
  size_t ncent = 0; for (auto& p : pts) ncent += p.centre;
  fprintf(stderr, "tau=%d kept=%zu centres=%zu (expect %.0f) %.1fs\n", TAU, pts.size(), ncent, Ntau / 24, omp_get_wtime() - t0);
  // 3D grid on (x1,x2,x3) with cell size CELL (projection to x' is 1-Lipschitz for the chordal metric)
  const long long OFF = 1 << 20;
  auto cellof = [&](double v) { return (long long)floor(v / CELL) + OFF; };
  auto mk = [&](long long i, long long j, long long k) { return (uint64_t)((i << 42) | (j << 21) | k); };
#pragma omp parallel for
  for (size_t i = 0; i < pts.size(); i++) pts[i].key = mk(cellof(pts[i].x[1]), cellof(pts[i].x[2]), cellof(pts[i].x[3]));
  __gnu_parallel::sort(pts.begin(), pts.end(), [](const Pt& a, const Pt& b) { return a.key < b.key; });
  vector<uint64_t> keys(pts.size()); for (size_t i = 0; i < pts.size(); i++) keys[i] = pts[i].key;
  fprintf(stderr, "sorted %.1fs\n", omp_get_wtime() - t0);
  int nm = rhos.size();
  vector<vector<long long>> hist(nm, vector<long long>(4096, 0));
  vector<int> best(nm, 0); vector<size_t> bestidx(nm, 0);
  int nt = omp_get_max_threads();
  vector<vector<vector<long long>>> lh(nt, vector<vector<long long>>(nm, vector<long long>(4096, 0)));
  vector<vector<int>> lb(nt, vector<int>(nm, 0)); vector<vector<size_t>> lbi(nt, vector<size_t>(nm, 0));
#pragma omp parallel for schedule(dynamic, 4096)
  for (size_t c = 0; c < pts.size(); c++) {
    if (!pts[c].centre) continue;
    int id = omp_get_thread_num();
    const double* x = pts[c].x;
    long long ci = cellof(x[1]), cj = cellof(x[2]), ck = cellof(x[3]);
    int cnt[16] = {0};
    for (int di = -1; di <= 1; di++) for (int dj = -1; dj <= 1; dj++) {
      uint64_t k0 = mk(ci + di, cj + dj, ck - 1), k1 = mk(ci + di, cj + dj, ck + 1);
      size_t s = lower_bound(keys.begin(), keys.end(), k0) - keys.begin();
      for (size_t u = s; u < keys.size() && keys[u] <= k1; u++) {
        const double* y = pts[u].x;
        double d2 = 0; for (int j = 0; j < 4; j++) { double e = x[j] - y[j]; d2 += e * e; }
        for (int m = 0; m < nm; m++) if (d2 <= rhos[m] * rhos[m]) cnt[m]++;
      }
    }
    for (int m = 0; m < nm; m++) {
      lh[id][m][min(cnt[m], 4095)]++;
      if (cnt[m] > lb[id][m]) { lb[id][m] = cnt[m]; lbi[id][m] = c; }
    }
  }
  for (int id = 0; id < nt; id++) for (int m = 0; m < nm; m++) {
    for (int k = 0; k < 4096; k++) hist[m][k] += lh[id][m][k];
    if (lb[id][m] > best[m]) { best[m] = lb[id][m]; bestidx[m] = lbi[id][m]; }
  }
  printf("# tau=%d N=%.0f %s time=%.1fs\n", TAU, Ntau, rnd ? (symr ? "SYMRANDOM" : "RANDOM") : "WORDS", omp_get_wtime() - t0);
  for (int m = 0; m < nm; m++) {
    const Pt& b = pts[bestidx[m]];
    printf("mu=%g rho=%.6g max=%d at q=(%.6f,%.6f,%.6f,%.6f) label=%llu\n", mus[m], rhos[m], best[m], b.x[0], b.x[1], b.x[2], b.x[3], (unsigned long long)b.label);
    printf("hist mu=%g:", mus[m]);
    for (int k = 0; k <= best[m]; k++) if (hist[m][k]) printf(" %d:%lld", k, hist[m][k]);
    printf("\n");
  }
  return 0;
}

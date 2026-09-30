// Moments of tube counts over Haar-random cosets G1 R_z G2.
//
// Hopf picture: W is in the eps-tube of G1 R_z G2 iff angle(R(W) x2, x1) <= phi(eps) = 4 asin(eps/2), where
// x1 = G1 z G1^-1, x2 = G2^-1 z G2 in S^2 and R(W) in SO(3) is the rotation of W.  Haar-random cosets correspond to
// independent uniform (x1, x2) on S^2 x S^2.  For each sampled x2 we compute the orbit {R(W) x2 : W in Lambda_tau}
// (72*2^tau points; words in Matsumoto-Amano normal form T^e (HT|SHT)^* C), grid it, and evaluate the cap counts
// N(x1, x2) at many uniform x1 for several radii, chosen so that the Haar mean E = N_tau (1 - cos phi)/2 takes the
// requested values.  Output: power sums of N (orders 1..10), maxima, histograms; the same for a control in which the
// orbit is replaced by N_tau i.i.d. uniform points.
//
// Usage: moments tau nx2 nx1 seed E1 [E2 ...]
#include <bits/stdc++.h>
#include <omp.h>
#include <parallel/algorithm>
using namespace std;
typedef complex<double> cd;
struct SU { cd a, b; };
static inline SU mul(const SU& p, const SU& q) { return {p.a * q.a - p.b * conj(q.b), p.a * q.b + p.b * conj(q.a)}; }
typedef array<double, 9> M3;
static M3 rot(const SU& m) {  // SO(3) matrix of X -> M X M^dag on x.sigma
  M3 R;
  for (int c = 0; c < 3; c++) {
    double e[3] = {0, 0, 0}; e[c] = 1;
    cd X00 = e[2], X01 = cd(e[0], -e[1]), X10 = cd(e[0], e[1]), X11 = -e[2];
    cd A = m.a, B = m.b, C = -conj(m.b), D = conj(m.a);  // M = [[A,B],[C,D]]
    // Y = M X M^dag
    cd T00 = A * X00 + B * X10, T01 = A * X01 + B * X11, T10 = C * X00 + D * X10, T11 = C * X01 + D * X11;
    cd Y00 = T00 * conj(A) + T01 * conj(B), Y10 = T10 * conj(A) + T11 * conj(B);
    R[0 * 3 + c] = real(Y10); R[1 * 3 + c] = imag(Y10); R[2 * 3 + c] = real(Y00);
  }
  return R;
}
static inline void app(const M3& R, const double* x, double* y) {
  for (int i = 0; i < 3; i++) y[i] = R[i * 3] * x[0] + R[i * 3 + 1] * x[1] + R[i * 3 + 2] * x[2];
}
int TAU;
M3 RA, RB, RT;
vector<M3> RC;  // 24 Clifford rotations
// DFS over prefixes P; vec holds the 24 vectors R(P) R(C) x2; emit R(P)R(C)x2 and R(T)R(P)R(C)x2
static void dfs(const double* vec, int d, vector<float>& out) {
  for (int c = 0; c < 24; c++) for (int i = 0; i < 3; i++) out.push_back((float)vec[3 * c + i]);
  if (d + 1 <= TAU) {
    double t[72], a[72], b[72];
    for (int c = 0; c < 24; c++) { app(RT, vec + 3 * c, t + 3 * c); app(RA, vec + 3 * c, a + 3 * c); app(RB, vec + 3 * c, b + 3 * c); }
    for (int c = 0; c < 24; c++) for (int i = 0; i < 3; i++) out.push_back((float)t[3 * c + i]);
    dfs(a, d + 1, out); dfs(b, d + 1, out);
  }
}
static void uniform_s2(mt19937_64& g, double* x) {
  normal_distribution<double> nd(0, 1);
  double n; do { x[0] = nd(g); x[1] = nd(g); x[2] = nd(g); n = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]); } while (n < 1e-12);
  x[0] /= n; x[1] /= n; x[2] /= n;
}
int main(int argc, char** argv) {
  if (argc < 6) { fprintf(stderr, "usage: moments tau nx2 nx1 seed E1 [E2..]\n"); return 1; }
  TAU = atoi(argv[1]); int nx2 = atoi(argv[2]); int nx1 = atoi(argv[3]); uint64_t seed = strtoull(argv[4], 0, 10);
  vector<double> Es; for (int i = 5; i < argc; i++) Es.push_back(atof(argv[i]));
  int nE = Es.size();
  const double r2 = sqrt(0.5);
  SU H = {cd(0, -r2), cd(0, -r2)}, S = {polar(1.0, -M_PI / 4), cd(0, 0)}, T = {polar(1.0, -M_PI / 8), cd(0, 0)};
  RA = rot(mul(H, T)); RB = rot(mul(S, mul(H, T))); RT = rot(T);
  {  // Clifford rotations by closure (as SO(3) matrices, 24 distinct)
    set<vector<long long>> seen; vector<M3> fr = {rot(SU{cd(1, 0), cd(0, 0)})};
    auto key = [](const M3& R) { vector<long long> k(9); for (int i = 0; i < 9; i++) k[i] = llround(R[i] * 1e6); return k; };
    seen.insert(key(fr[0])); RC.push_back(fr[0]);
    M3 GH = rot(H), GS = rot(S);
    while (!fr.empty()) {
      vector<M3> nf;
      for (auto& R : fr) for (auto& G : {GH, GS}) {
        M3 P; for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) { double s = 0; for (int k = 0; k < 3; k++) s += G[i * 3 + k] * R[k * 3 + j]; P[i * 3 + j] = s; }
        auto k = key(P); if (!seen.count(k)) { seen.insert(k); nf.push_back(P); RC.push_back(P); }
      }
      fr = nf;
    }
    if (RC.size() != 24) { fprintf(stderr, "clifford %zu\n", RC.size()); return 1; }
  }
  const double Ntau = 72.0 * ldexp(1.0, TAU) - 48;
  vector<double> chord(nE);
  for (int e = 0; e < nE; e++) { double c = 1 - 2 * Es[e] / Ntau; double phi = acos(max(-1.0, c)); chord[e] = 2 * sin(phi / 2); }
  double cmax = *max_element(chord.begin(), chord.end());
  const int MO = 10;
  // accumulators [mode 0 = words, 1 = control][E][order]
  vector<long double> ps(2 * nE * (MO + 1), 0); vector<long long> mx(2 * nE, 0);
  vector<vector<long long>> hist(2 * nE, vector<long long>(4096, 0));
  mt19937_64 g(seed);
  double t0 = omp_get_wtime();
  // optional: X2AXIS="ax,ay,az,psi" samples x2 uniformly in the cap of angular radius psi around the axis
  const char* ax = getenv("X2AXIS");
  double axv[3] = {0, 0, 1}, axpsi = -1;
  if (ax) { sscanf(ax, "%lf,%lf,%lf,%lf", &axv[0], &axv[1], &axv[2], &axpsi);
    double nn = sqrt(axv[0]*axv[0]+axv[1]*axv[1]+axv[2]*axv[2]); for (int i=0;i<3;i++) axv[i]/=nn; }
  for (int s2 = 0; s2 < nx2; s2++) {
    double x2[3]; uniform_s2(g, x2);
    if (axpsi > 0) {  // rejection-free: uniform in cap via cos-angle uniform in [cos psi, 1]
      uniform_real_distribution<double> ud(0, 1);
      double c = 1 - ud(g) * (1 - cos(axpsi)), sn = sqrt(max(0.0, 1 - c * c)), th = 2 * M_PI * ud(g);
      double e1[3], e2[3]; double t[3] = {1, 0, 0}; if (fabs(axv[0]) > 0.9) { t[0] = 0; t[1] = 1; }
      e1[0] = axv[1]*t[2]-axv[2]*t[1]; e1[1] = axv[2]*t[0]-axv[0]*t[2]; e1[2] = axv[0]*t[1]-axv[1]*t[0];
      double n1 = sqrt(e1[0]*e1[0]+e1[1]*e1[1]+e1[2]*e1[2]); for (int i=0;i<3;i++) e1[i]/=n1;
      e2[0] = axv[1]*e1[2]-axv[2]*e1[1]; e2[1] = axv[2]*e1[0]-axv[0]*e1[2]; e2[2] = axv[0]*e1[1]-axv[1]*e1[0];
      for (int i = 0; i < 3; i++) x2[i] = c * axv[i] + sn * (cos(th) * e1[i] + sin(th) * e2[i]);
    }
    for (int mode = 0; mode < 2; mode++) {
      vector<float> pts;
      if (mode == 0) {
        // start vectors R(C) x2
        double v0[72];
        for (int c = 0; c < 24; c++) app(RC[c], x2, v0 + 3 * c);
        // split the tree at depth 6 for parallelism
        int split = min(TAU, 6);
        vector<vector<double>> roots; vector<int> rdepth;
        vector<float> shallow;
        function<void(const double*, int)> go = [&](const double* vec, int d) {
          if (d == split) { roots.push_back(vector<double>(vec, vec + 72)); rdepth.push_back(d); return; }
          for (int c = 0; c < 24; c++) for (int i = 0; i < 3; i++) shallow.push_back((float)vec[3 * c + i]);
          if (d + 1 <= TAU) {
            double t[72], a[72], b[72];
            for (int c = 0; c < 24; c++) { app(RT, vec + 3 * c, t + 3 * c); app(RA, vec + 3 * c, a + 3 * c); app(RB, vec + 3 * c, b + 3 * c); }
            for (int c = 0; c < 24; c++) for (int i = 0; i < 3; i++) shallow.push_back((float)t[3 * c + i]);
            go(a, d + 1); go(b, d + 1);
          }
        };
        go(v0, 0);
        vector<vector<float>> parts(roots.size());
#pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < roots.size(); i++) dfs(roots[i].data(), rdepth[i], parts[i]);
        size_t tot = shallow.size(); for (auto& p : parts) tot += p.size();
        pts.reserve(tot); pts.insert(pts.end(), shallow.begin(), shallow.end());
        for (auto& p : parts) { pts.insert(pts.end(), p.begin(), p.end()); vector<float>().swap(p); }
      } else {
        size_t n = (size_t)Ntau; pts.resize(3 * n);
#pragma omp parallel
        {
          mt19937_64 gg(seed * 7919 + s2 * 104729 + omp_get_thread_num());
#pragma omp for
          for (long long i = 0; i < (long long)n; i++) { double x[3]; uniform_s2(gg, x); pts[3 * i] = x[0]; pts[3 * i + 1] = x[1]; pts[3 * i + 2] = x[2]; }
        }
      }
      size_t n = pts.size() / 3;
      // grid
      const long long OFF = 1 << 20;
      auto cell = [&](double v) { return (long long)floor(v / cmax) + OFF; };
      auto mk = [](long long i, long long j, long long k) { return (uint64_t)((i << 42) | (j << 21) | k); };
      vector<pair<uint64_t, uint32_t>> kv(n);
#pragma omp parallel for
      for (size_t i = 0; i < n; i++) kv[i] = {mk(cell(pts[3 * i]), cell(pts[3 * i + 1]), cell(pts[3 * i + 2])), (uint32_t)i};
      __gnu_parallel::sort(kv.begin(), kv.end());
      vector<uint64_t> keys(n); vector<float> sp(3 * n);
      for (size_t i = 0; i < n; i++) { keys[i] = kv[i].first; for (int c = 0; c < 3; c++) sp[3 * i + c] = pts[3 * kv[i].second + c]; }
      vector<pair<uint64_t, uint32_t>>().swap(kv); vector<float>().swap(pts);
      // queries
      int nt = omp_get_max_threads();
      vector<vector<long double>> lps(nt, vector<long double>(nE * (MO + 1), 0));
      vector<vector<long long>> lmx(nt, vector<long long>(nE, 0));
      vector<vector<vector<long long>>> lh(nt, vector<vector<long long>>(nE, vector<long long>(4096, 0)));
#pragma omp parallel
      {
        int id = omp_get_thread_num();
        mt19937_64 gq(seed * 31337 + s2 * 1009 + mode * 17 + id);
#pragma omp for schedule(static)
        for (int q = 0; q < nx1; q++) {
          double x1[3]; uniform_s2(gq, x1);
          long long ci = cell(x1[0]), cj = cell(x1[1]), ck = cell(x1[2]);
          long long cnt[16] = {0};
          for (int di = -1; di <= 1; di++) for (int dj = -1; dj <= 1; dj++) {
            uint64_t k0 = mk(ci + di, cj + dj, ck - 1), k1 = mk(ci + di, cj + dj, ck + 1);
            size_t s = lower_bound(keys.begin(), keys.end(), k0) - keys.begin();
            for (size_t u = s; u < n && keys[u] <= k1; u++) {
              double dx = sp[3 * u] - x1[0], dy = sp[3 * u + 1] - x1[1], dz = sp[3 * u + 2] - x1[2];
              double d = sqrt(dx * dx + dy * dy + dz * dz);
              for (int e = 0; e < nE; e++) if (d <= chord[e]) cnt[e]++;
            }
          }
          for (int e = 0; e < nE; e++) {
            long double p = 1;
            for (int o = 0; o <= MO; o++) { lps[id][e * (MO + 1) + o] += p; p *= cnt[e]; }
            lmx[id][e] = max(lmx[id][e], cnt[e]);
            lh[id][e][min(cnt[e], 4095LL)]++;
          }
        }
      }
      for (int id = 0; id < nt; id++) for (int e = 0; e < nE; e++) {
        for (int o = 0; o <= MO; o++) ps[(mode * nE + e) * (MO + 1) + o] += lps[id][e * (MO + 1) + o];
        mx[mode * nE + e] = max(mx[mode * nE + e], lmx[id][e]);
        for (int h = 0; h < 4096; h++) hist[mode * nE + e][h] += lh[id][e][h];
      }
    }
    if ((s2 + 1) % max(1, nx2 / 10) == 0) fprintf(stderr, "tau=%d x2 %d/%d %.0fs\n", TAU, s2 + 1, nx2, omp_get_wtime() - t0);
  }
  printf("# tau=%d N=%.0f nx2=%d nx1=%d seed=%llu\n", TAU, Ntau, nx2, nx1, (unsigned long long)seed);
  for (int mode = 0; mode < 2; mode++) for (int e = 0; e < nE; e++) {
    long double* p = &ps[(mode * nE + e) * (MO + 1)];
    printf("%s E=%g max=%lld moments:", mode ? "CONTROL" : "WORDS", Es[e], mx[mode * nE + e]);
    for (int o = 1; o <= MO; o++) printf(" %.6Le", p[o] / p[0]);
    printf("\n");
    printf("%s E=%g hist:", mode ? "CONTROL" : "WORDS", Es[e]);
    for (int h = 0; h < 4096; h++) if (hist[mode * nE + e][h]) printf(" %d:%lld", h, hist[mode * nE + e][h]);
    printf("\n");
  }
  return 0;
}

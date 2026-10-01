// orbit.c -- Hecke-ball tubes (rounds t >= 2): the orbit {W phi : W in Lambda_t} of a source state
// phi = V|0> (V a Clifford+T word of T-count h), binned in Bloch-sphere cells of side s = 2^{-beta t}.
// Lambda_t = all Clifford+T unitaries of T-count <= t (mod phase), in Matsumoto-Amano normal form
//   W = (T | e) (HT | SHT)^m C,  C in the 24-element Clifford group;  |Lambda_t| = 72*2^t - 48.
// Usage: ./orbit t beta h seed      (h = 0 gives phi = |0>, the two-round case)
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <omp.h>

typedef struct { double m[9]; } M3;
static M3 mul(M3 a, M3 b) {
  M3 c;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) c.m[3 * i + j] = a.m[3 * i] * b.m[j] + a.m[3 * i + 1] * b.m[3 + j] + a.m[3 * i + 2] * b.m[6 + j];
  return c;
}
static inline void app(const M3 *a, const double *p, double *q) {
  for (int i = 0; i < 3; i++) q[i] = a->m[3 * i] * p[0] + a->m[3 * i + 1] * p[1] + a->m[3 * i + 2] * p[2];
}

static M3 Tm, Hm, Sm, HT, SHT, CL[24];
static int T_MAX, nz, nphi;
static double s;
static uint32_t *cells;

static inline void bin(const double *p) {
  int iz = (int)((p[2] + 1.0) / s);
  if (iz >= nz) iz = nz - 1;
  if (iz < 0) iz = 0;
  double ph = atan2(p[1], p[0]) + M_PI;
  int ip = (int)(ph / s);
  if (ip >= nphi) ip = nphi - 1;
  __atomic_fetch_add(&cells[(size_t)iz * nphi + ip], 1u, __ATOMIC_RELAXED);
}

// node at depth d (d syllables applied): emit p (T-count d) and T p (T-count d+1), then recurse
static void dfs(const double *p, int d) {
  double q[3];
  bin(p);
  if (d + 1 <= T_MAX) {
    app(&Tm, p, q);
    bin(q);
    app(&HT, p, q);
    dfs(q, d + 1);
    app(&SHT, p, q);
    dfs(q, d + 1);
  }
}

int main(int argc, char **argv) {
  if (argc < 5) {
    fprintf(stderr, "usage: %s t beta h seed\n", argv[0]);
    return 1;
  }
  T_MAX = atoi(argv[1]);
  double beta = atof(argv[2]);
  int h = atoi(argv[3]);
  unsigned seed = (unsigned)atoi(argv[4]);
  double r = 1 / sqrt(2.0);
  Tm = (M3){{r, -r, 0, r, r, 0, 0, 0, 1}};
  Hm = (M3){{0, 0, 1, 0, -1, 0, 1, 0, 0}};
  Sm = (M3){{0, -1, 0, 1, 0, 0, 0, 0, 1}};
  HT = mul(Hm, Tm);
  SHT = mul(Sm, HT);
  // Clifford group: signed permutation matrices with det +1
  int nc = 0;
  int perm[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
  int psign[6] = {1, -1, -1, 1, 1, -1};
  for (int pi = 0; pi < 6; pi++)
    for (int sg = 0; sg < 8; sg++) {
      int s0 = (sg & 1) ? -1 : 1, s1 = (sg & 2) ? -1 : 1, s2 = (sg & 4) ? -1 : 1;
      if (psign[pi] * s0 * s1 * s2 != 1) continue;
      M3 c;
      memset(&c, 0, sizeof c);
      int sgn[3] = {s0, s1, s2};
      for (int i = 0; i < 3; i++) c.m[3 * i + perm[pi][i]] = sgn[i];
      CL[nc++] = c;
    }
  if (nc != 24) {
    fprintf(stderr, "clifford count %d\n", nc);
    return 1;
  }
  // source phi = V|0>, V a random MA word with h syllables (T-count h), V = (HT|SHT)^h C0
  srand(seed);
  double phi[3] = {0, 0, 1}, q[3];
  M3 V = CL[rand() % 24];
  for (int i = 0; i < h; i++) V = mul((rand() & 1) ? SHT : HT, V);
  app(&V, phi, q);
  memcpy(phi, q, sizeof q);

  s = pow(2.0, -beta * T_MAX);
  nz = (int)ceil(2.0 / s);
  nphi = (int)ceil(2 * M_PI / s);
  cells = calloc((size_t)nz * nphi, sizeof(uint32_t));
  double t0 = omp_get_wtime();

  // tasks: Clifford start x syllable prefixes of length D0 (inner syllables)
  int D0 = T_MAX < 8 ? T_MAX : 8;
  int ntask = 24 << D0;
  // nodes of depth < D0 are emitted here (serially), deeper ones inside the tasks
  for (int c = 0; c < 24; c++) {
    double p0[3];
    app(&CL[c], phi, p0);
    // breadth-first over depths 0..D0-1
    int cnt = 1;
    double *lev = malloc(3 * sizeof(double));
    memcpy(lev, p0, 3 * sizeof(double));
    for (int d = 0; d < D0; d++) {
      double *nxt = malloc(2 * cnt * 3 * sizeof(double));
      for (int i = 0; i < cnt; i++) {
        double *p = lev + 3 * i;
        bin(p);
        if (d + 1 <= T_MAX) {
          app(&Tm, p, q);
          bin(q);
        }
        app(&HT, p, nxt + 6 * i);
        app(&SHT, p, nxt + 6 * i + 3);
      }
      free(lev);
      lev = nxt;
      cnt *= 2;
    }
    free(lev);
  }
#pragma omp parallel for schedule(dynamic, 1)
  for (int tk = 0; tk < ntask; tk++) {
    int c = tk >> D0, bits = tk & ((1 << D0) - 1);
    double p[3], qq[3];
    app(&CL[c], phi, p);
    for (int i = 0; i < D0; i++) {
      app((bits >> i) & 1 ? &SHT : &HT, p, qq);
      memcpy(p, qq, sizeof p);
    }
    if (D0 <= T_MAX) dfs(p, D0);
  }
  double t1 = omp_get_wtime();

  // statistics on the bulk |z| <= 0.9, phi periodic
  int z0 = (int)((1 - 0.9) / s), z1 = nz - z0;
  double total = 72.0 * pow(2.0, T_MAX) - 48;
  double lam = total * s * s / (4 * M_PI);
  uint64_t sum = 0, sum2 = 0, nb = 0;
  uint32_t mx1 = 0;
  int bz = 0, bp = 0;
  uint64_t mx2 = 0, mx4 = 0;
  for (int iz = z0; iz < z1; iz++)
    for (int ip = 0; ip < nphi - 1; ip++) {
      uint32_t c = cells[(size_t)iz * nphi + ip];
      sum += c;
      sum2 += (uint64_t)c * c;
      nb++;
      if (c > mx1) {
        mx1 = c;
        bz = iz;
        bp = ip;
      }
    }
  for (int w = 2; w <= 4; w += 2) {
    uint64_t mx = 0;
    for (int iz = z0; iz + w <= z1; iz++)
      for (int ip = 0; ip < nphi - 1; ip++) {
        uint64_t a = 0;
        for (int dz = 0; dz < w; dz++)
          for (int dp = 0; dp < w; dp++) a += cells[(size_t)(iz + dz) * nphi + (ip + dp) % (nphi - 1)];
        if (a > mx) mx = a;
      }
    if (w == 2) mx2 = mx;
    else mx4 = mx;
  }
  uint64_t all = 0;
  for (size_t i = 0; i < (size_t)nz * nphi; i++) all += cells[i];
  if (all != (uint64_t)(72.0 * pow(2.0, T_MAX) - 48)) printf("COUNT MISMATCH %llu vs %.0f\n", (unsigned long long)all, 72.0 * pow(2.0, T_MAX) - 48);
  double mean = (double)sum / nb, var = (double)sum2 / nb - mean * mean;
  printf("t=%d beta=%.3f h=%d seed=%u phi=(%.4f,%.4f,%.4f) s=2^-%.2f lambda=%.2f mean=%.2f var/mean=%.3f "
         "max1=%u at (z=%.4f,phi=%.4f) max2=%llu max4=%llu  time %.1fs\n",
         T_MAX, beta, h, seed, phi[0], phi[1], phi[2], beta * T_MAX, lam, mean, var / mean, mx1, -1 + (bz + 0.5) * s,
         (bp + 0.5) * s - M_PI, (unsigned long long)mx2, (unsigned long long)mx4, t1 - t0);
  return 0;
}

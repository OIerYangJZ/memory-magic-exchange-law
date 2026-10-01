// sl.c -- numerics for research/rate3/NOTES.md §10(A) (band decomposition).
//
// Clifford+T states: (u,t) in Z[w]^2, w = e^{i pi/4}, |u|^2 + |t|^2 = 2^k in Z[sqrt2] (both embeddings).
// u = a + b w + c w^2 + d w^3:  |u|^2 = A + B sqrt2,  A = a^2+b^2+c^2+d^2,  B = ab+bc+cd-da,
// sigma1(u) = (a + (b-d)/sqrt2) + i (c + (b+d)/sqrt2).
// Orbit representatives under u -> w u: arg sigma1(u) in [0, pi/4).
//
// Test A: S_l(W_j) = sum_{m in W_j} a_l(m) a_l(2^k - m),  a_l(m) = sum_{|u|^2 = m} cos(l arg sigma1 u),
//         l = 8n, n = 0..N (N = nb/8), W_j = { m : sigma1(m)/2^k in [j/nb, (j+1)/nb) }.
//         S_0(W_j) = number of pairs (u,t) in the band.  a_l is real by conjugation symmetry.
// Test B: counts of orbit pairs (u-orbit, t-orbit) in cells of the folded Bloch sphere
//         z = (|u|^2-|t|^2)/2^k in [-1,1], phi = arg t - arg u mod pi/4; cell side 2/nb in z and phi.
//
// Usage: ./sl k beta outprefix      (nb = round(2^{beta k}))
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <omp.h>

static const double SQ2 = 1.41421356237309504880;
static int64_t P;
static int64_t *off;
static int32_t *bm;

static inline int64_t midx(int64_t A, int64_t B) { return off[A] + B + bm[A]; }

static inline int isrep(int64_t a, int64_t b, int64_t c, int64_t d) {
  if (a == 0 && b == d) return 0;  // x == 0 exactly
  double x = a + (b - d) / SQ2, y = c + (b + d) / SQ2;
  if (!(x > 0)) return 0;
  if (!(c == 0 && b == -d) && !(y > 0)) return 0;  // y >= 0 (exact zero handled)
  if (a == c && d == 0) return 0;                  // y == x exactly
  return y < x;
}

static void cheb(const double *t, uint32_t n, int N, double *v) {
  for (int j = 0; j <= N; j++) v[j] = 0.0;
  for (uint32_t i = 0; i < n; i++) {
    double x = cos(8.0 * t[i]), x2 = 2.0 * x, t0 = 1.0, t1 = x;
    v[0] += 1.0;
    if (N >= 1) v[1] += x;
    for (int j = 2; j <= N; j++) {
      double t2 = x2 * t1 - t0;
      v[j] += t2;
      t0 = t1;
      t1 = t2;
    }
  }
}

int main(int argc, char **argv) {
  if (argc < 4) {
    fprintf(stderr, "usage: %s k beta outprefix\n", argv[0]);
    return 1;
  }
  int K = atoi(argv[1]);
  double beta = atof(argv[2]);
  const char *pre = argv[3];
  P = (int64_t)1 << K;
  int nb = (int)llround(pow(2.0, beta * K));
  int N = nb / 8;
  int nphi = (int)ceil(M_PI * nb / 8.0);
  double t0 = omp_get_wtime();

  // m = A + B sqrt2 with 0 <= sigma_{1,2}(m) <= 2^k  <=>  |B| sqrt2 <= min(A, P-A)
  off = malloc((P + 2) * sizeof(int64_t));
  bm = malloc((P + 1) * sizeof(int32_t));
  int64_t acc = 0;
  for (int64_t A = 0; A <= P; A++) {
    int64_t M = A < P - A ? A : P - A;
    int64_t b = (int64_t)floor(M / SQ2);
    while (2 * (b + 1) * (b + 1) <= M * M) b++;
    while (b > 0 && 2 * b * b > M * M) b--;
    bm[A] = (int32_t)b;
    off[A] = acc;
    acc += 2 * b + 1;
  }
  off[P + 1] = acc;
  int64_t NM = acc;
  uint32_t *cnt = calloc(NM, sizeof(uint32_t));
  uint32_t *beg = malloc((NM + 1) * sizeof(uint32_t));
  int64_t R = (int64_t)floor(sqrt((double)P));
  while ((R + 1) * (R + 1) <= P) R++;

  // two passes over u: count representatives per m, then fill their angles
  int64_t nu = 0, nrep = 0;
  double *th = NULL;
  for (int pass = 0; pass < 2; pass++) {
    int64_t nu_ = 0, nrep_ = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : nu_, nrep_)
    for (int64_t a = -R; a <= R; a++)
      for (int64_t b = -R; b <= R; b++) {
        int64_t s2 = a * a + b * b;
        if (s2 > P) continue;
        for (int64_t c = -R; c <= R; c++) {
          int64_t s3 = s2 + c * c;
          if (s3 > P) continue;
          int64_t dm = (int64_t)floor(sqrt((double)(P - s3)));
          while ((dm + 1) * (dm + 1) <= P - s3) dm++;
          while (dm * dm > P - s3) dm--;
          for (int64_t d = -dm; d <= dm; d++) {
            int64_t A = s3 + d * d, B = a * b + b * c + c * d - d * a, Q = P - A;
            if (2 * B * B > Q * Q) continue;
            nu_++;
            if (!isrep(a, b, c, d)) continue;
            nrep_++;
            int64_t im = midx(A, B);
            if (pass == 0) {
              __atomic_fetch_add(&cnt[im], 1u, __ATOMIC_RELAXED);
            } else {
              uint32_t k = __atomic_sub_fetch(&cnt[im], 1u, __ATOMIC_RELAXED);
              double x = a + (b - d) / SQ2, y = c + (b + d) / SQ2;
              th[beg[im] + k] = atan2(y, x);
            }
          }
        }
      }
    if (pass == 0) {
      nu = nu_;
      nrep = nrep_;
      if (nrep >= 4294967295LL) {
        fprintf(stderr, "too many reps\n");
        return 1;
      }
      beg[0] = 0;
      for (int64_t i = 0; i < NM; i++) beg[i + 1] = beg[i] + cnt[i];
      th = malloc(nrep * sizeof(double));
    }
  }
  free(cnt);
  double t1 = omp_get_wtime();

  // Test A
  int nth = omp_get_max_threads();
  size_t W = (size_t)nb * (N + 1);
  double *bs_all = calloc((size_t)nth * W, sizeof(double));
  int64_t *Lc_all = calloc((size_t)nth * nb, sizeof(int64_t));
#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    double *bs = bs_all + (size_t)tid * W;
    int64_t *Lc = Lc_all + (size_t)tid * nb;
    double *va = malloc((N + 1) * sizeof(double)), *vb = malloc((N + 1) * sizeof(double));
#pragma omp for schedule(dynamic, 8)
    for (int64_t A = 0; A <= P / 2; A++)
      for (int64_t B = -bm[A]; B <= bm[A]; B++) {
        int64_t Ap = P - A, Bp = -B;
        if (A == Ap && B > 0) continue;
        int64_t im = midx(A, B), ip = midx(Ap, Bp);
        uint32_t na = beg[im + 1] - beg[im], np = beg[ip + 1] - beg[ip];
        if (!na || !np) continue;
        cheb(th + beg[im], na, N, va);
        cheb(th + beg[ip], np, N, vb);
        int jm = (int)((A + B * SQ2) / P * nb), jp = (int)((Ap + Bp * SQ2) / P * nb);
        if (jm >= nb) jm = nb - 1;
        if (jp >= nb) jp = nb - 1;
        double *r1 = bs + (size_t)jm * (N + 1), *r2 = bs + (size_t)jp * (N + 1);
        int self = (im == ip);
        for (int n = 0; n <= N; n++) {
          double F = 64.0 * va[n] * vb[n];
          r1[n] += F;
          if (!self) r2[n] += F;
        }
        Lc[jm]++;
        if (!self) Lc[jp]++;
      }
    free(va);
    free(vb);
  }
  for (int t = 1; t < nth; t++) {
    for (size_t i = 0; i < W; i++) bs_all[i] += bs_all[(size_t)t * W + i];
    for (int j = 0; j < nb; j++) Lc_all[j] += Lc_all[(size_t)t * nb + j];
  }
  double t2 = omp_get_wtime();

  // Test B
  uint32_t *cells = calloc((size_t)nb * nphi, sizeof(uint32_t));
  int64_t npairs = 0;
#pragma omp parallel for schedule(dynamic, 8) reduction(+ : npairs)
  for (int64_t A = 0; A <= P; A++)
    for (int64_t B = -bm[A]; B <= bm[A]; B++) {
      int64_t im = midx(A, B), ip = midx(P - A, -B);
      uint32_t na = beg[im + 1] - beg[im], np = beg[ip + 1] - beg[ip];
      if (!na || !np) continue;
      int iz = (int)((A + B * SQ2) / P * nb);
      if (iz >= nb) iz = nb - 1;
      for (uint32_t i = 0; i < na; i++)
        for (uint32_t j = 0; j < np; j++) {
          double ph = th[beg[ip] + j] - th[beg[im] + i];
          if (ph < 0) ph += M_PI / 4;
          int ip_ = (int)(ph * nb / 2.0);
          if (ip_ >= nphi) ip_ = nphi - 1;
          __atomic_fetch_add(&cells[(size_t)iz * nphi + ip_], 1u, __ATOMIC_RELAXED);
          npairs++;
        }
    }
  double t3 = omp_get_wtime();

  char fn[512];
  snprintf(fn, sizeof fn, "%s_A.bin", pre);
  FILE *f = fopen(fn, "wb");
  fwrite(bs_all, sizeof(double), W, f);
  fclose(f);
  snprintf(fn, sizeof fn, "%s_L.bin", pre);
  f = fopen(fn, "wb");
  fwrite(Lc_all, sizeof(int64_t), nb, f);
  fclose(f);
  snprintf(fn, sizeof fn, "%s_B.bin", pre);
  f = fopen(fn, "wb");
  fwrite(cells, sizeof(uint32_t), (size_t)nb * nphi, f);
  fclose(f);
  snprintf(fn, sizeof fn, "%s.txt", pre);
  f = fopen(fn, "w");
  fprintf(f, "k %d\nbeta %.6f\nnb %d\nN %d\nnphi %d\nNM %lld\nnu %lld\nnrep %lld\nnpairs %lld\n", K, beta, nb, N,
          nphi, (long long)NM, (long long)nu, (long long)nrep, (long long)npairs);
  fprintf(f, "time_enum %.1f\ntime_A %.1f\ntime_B %.1f\n", t1 - t0, t2 - t1, t3 - t2);
  fclose(f);
  printf("k=%d nb=%d N=%d NM=%lld nu=%lld nrep=%lld npairs=%lld  t=%.1f/%.1f/%.1f s\n", K, nb, N, (long long)NM,
         (long long)nu, (long long)nrep, (long long)npairs, t1 - t0, t2 - t1, t3 - t2);
  return 0;
}

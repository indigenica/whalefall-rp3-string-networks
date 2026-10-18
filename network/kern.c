/* U2 lattice kernel: RP3 Q-model (10 comps, 4x4 symmetric) and U(1) control (2 comps) in radiation domination.
   Conformal time, a = s*tau, physical strings.  pi = a^2 phi'.  EOM: pi' = a^2 (lap phi - a^2 dV/dphi), phi' = pi / a^2.
   Fields stored component-major: f[c*V + idx], idx = (z*N + y)*N + x, float32.  OpenMP over z.
   Build: gcc -O3 -march=native -fopenmp -shared -fPIC kern.c -o kern.so -lm                                      */
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <omp.h>

#define LAM_Z  0.8    /* lambda / Z_Q   (320/400) */
#define LAMP_Z 0.48   /* lambda_P / Z_Q (192/400) */
#define LAM_U  0.24   /* U(1): V = LAM_U/4 (phi^2 - V2)^2 */
#define V2_U   0.5
static const double ISQ2 = 0.70710678118654752440;

/* --- RP3 helpers: q[10] -> Q[4][4]; offdiag order 01,02,03,12,13,23 --- */
static const int OI[6] = {0, 0, 0, 1, 1, 2}, OJ[6] = {1, 2, 3, 2, 3, 3};
static inline void q2Q(const double *q, double Q[4][4]) {
    for (int i = 0; i < 4; i++) Q[i][i] = q[i];
    for (int k = 0; k < 6; k++) { double v = q[4 + k] * ISQ2; Q[OI[k]][OJ[k]] = v; Q[OJ[k]][OI[k]] = v; }
}
static inline void mm(double A[4][4], double B[4][4], double C[4][4]) {
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) { double s = 0; for (int k = 0; k < 4; k++) s += A[i][k] * B[k][j]; C[i][j] = s; }
}
/* U19 (ESPEC_U19): explicit breaking (N_DW = 1).  RP3: V -> V - EPS_G Q_11 (m_G^2 = EPS_G); U(1): V -> V - EPS_G phi_1
   (m_G^2 = EPS_G / v).  Off (0) by default: then nothing is added and the evolution is bit-identical to U2-U16. */
static double EPS_G = 0.0;
void set_tilt(double eps) { EPS_G = eps; }
/* gradient dV/dq (10 comps) and potential value */
static inline double rp3_grad(const double *q, double *g) {
    double Q[4][4], Q2[4][4], M[4][4], T[4][4], W[4][4];
    q2Q(q, Q); mm(Q, Q, Q2);
    double tr = Q[0][0] + Q[1][1] + Q[2][2] + Q[3][3];
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) { M[i][j] = Q2[i][j] - Q[i][j]; W[i][j] = 2 * Q[i][j] - (i == j); }
    mm(M, W, T);
    double trM2 = 0; for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) trM2 += M[i][j] * M[j][i];
    for (int i = 0; i < 4; i++) g[i] = 0.5 * LAM_Z * (tr - 1) + 0.5 * LAMP_Z * T[i][i];
    for (int k = 0; k < 6; k++) g[4 + k] = 0.5 * LAMP_Z * 0.5 * (T[OI[k]][OJ[k]] + T[OJ[k]][OI[k]]) * 2 * ISQ2;
    double V = 0.25 * LAM_Z * (tr - 1) * (tr - 1) + 0.25 * LAMP_Z * trM2;
    if (EPS_G != 0.0) { g[0] -= EPS_G; V -= EPS_G * q[0]; }
    return V;
}
static inline double u1_grad(const double *p, double *g) {
    double r2 = p[0] * p[0] + p[1] * p[1], d = r2 - V2_U;
    g[0] = LAM_U * d * p[0]; g[1] = LAM_U * d * p[1];
    double V = 0.25 * LAM_U * d * d;
    if (EPS_G != 0.0) { g[0] -= EPS_G; V -= EPS_G * p[0]; }
    return V;
}

/* force F = a^2 (lap phi - a^2 dV) into out; returns nothing */
/* U3b: core-growth preparation. For tau < TAU_C the potential uses a frozen scale factor a_p = s*TAU_C (constant comoving
   core width); afterwards a_p = a (physical strings, identical to U2). */
static double TAU_C = 0.0;
void set_core(double tau_c) { TAU_C = tau_c; }
/* U19b (ESPEC_U19b): late frozen core.  For tau > TAU_F the symmetric potential uses a_p = s*TAU_F (constant comoving core
   width), while the explicit tilt keeps the physical a^2.  Off (0) by default: then a_p = a and nothing changes. */
static double TAU_F = 0.0;
void set_core_late(double tau_f) { TAU_F = tau_f; }
static void force(int model, int N, const float *f, float *out, double a, double ap) {
    const long V = (long)N * N * N; const int C = model ? 10 : 2; const double a2 = a * a, ap2 = ap * ap;
    #pragma omp parallel for schedule(static)
    for (int z = 0; z < N; z++) {
        int zp = (z + 1) % N, zm = (z + N - 1) % N;
        double q[10], g[10];
        for (int y = 0; y < N; y++) {
            int yp = (y + 1) % N, ym = (y + N - 1) % N;
            for (int x = 0; x < N; x++) {
                int xp = (x + 1) % N, xm = (x + N - 1) % N;
                long i = ((long)z * N + y) * N + x;
                for (int c = 0; c < C; c++) q[c] = f[c * V + i];
                if (model) rp3_grad(q, g); else u1_grad(q, g);
                if (EPS_G != 0.0 && ap2 != a2) g[0] += EPS_G * (1.0 - a2 / ap2);   /* U19b: tilt with a^2, not a_p^2 */
                for (int c = 0; c < C; c++) {
                    const float *fc = f + c * V;
                    double lap = fc[((long)z * N + y) * N + xp] + fc[((long)z * N + y) * N + xm] + fc[((long)z * N + yp) * N + x]
                               + fc[((long)z * N + ym) * N + x] + fc[((long)zp * N + y) * N + x] + fc[((long)zm * N + y) * N + x] - 6.0 * q[c];
                    out[c * V + i] = (float)(a2 * (lap - ap2 * g[c]));
                }
            }
        }
    }
}

/* U10a (ESPEC_U10a): drain of the massive part of the velocity for tau < DAMP_TAU, at physical rate DAMP_G, only where
   the order parameter is near the vacuum (RP3: gap lambda1 - lambda2 > 0.6; U(1): |phi|/v > 0.6). Off when DAMP_TAU <= 0. */
static double DAMP_TAU = 0.0, DAMP_G = 0.0;
#define DAMP_EVERY 4
void set_damp(double tau_d, double gamma) { DAMP_TAU = tau_d; DAMP_G = gamma; }
static void jacobi4(double A[4][4], double w[4], double Vv[4][4]);
void damp_massive(int model, int N, const float *f, float *p, double fac) {
    const long V = (long)N * N * N;
    #pragma omp parallel for schedule(static)
    for (long i = 0; i < V; i++) {
        if (model) {
            double q[10], Q[4][4], w[4], Vv[4][4], P[4][4], n[4], u[4], wv[4];
            for (int c = 0; c < 10; c++) q[c] = f[c * V + i];
            q2Q(q, Q); jacobi4(Q, w, Vv);
            int b = 0; for (int k = 1; k < 4; k++) if (w[k] > w[b]) b = k;
            double s2 = -1e300; for (int k = 0; k < 4; k++) if (k != b && w[k] > s2) s2 = w[k];
            if (w[b] - s2 <= 0.6) continue;
            for (int k = 0; k < 4; k++) n[k] = Vv[k][b];
            double pc[10]; for (int c = 0; c < 10; c++) pc[c] = p[c * V + i];
            q2Q(pc, P);
            double nPn = 0;
            for (int r = 0; r < 4; r++) { u[r] = 0; for (int c = 0; c < 4; c++) u[r] += P[r][c] * n[c]; nPn += n[r] * u[r]; }
            for (int r = 0; r < 4; r++) wv[r] = u[r] - nPn * n[r];
            double Pn[4][4];
            for (int r = 0; r < 4; r++) for (int c = 0; c < 4; c++) {
                double pt = wv[r] * n[c] + n[r] * wv[c];
                Pn[r][c] = pt + fac * (P[r][c] - pt);
            }
            for (int k = 0; k < 4; k++) p[k * V + i] = (float)Pn[k][k];
            for (int k = 0; k < 6; k++) p[(4 + k) * V + i] = (float)(Pn[OI[k]][OJ[k]] * 2 * ISQ2);
        } else {
            double f0 = f[i], f1 = f[V + i], r2 = f0 * f0 + f1 * f1;
            if (r2 <= 0.36 * V2_U) continue;
            double rad = (f0 * p[i] + f1 * p[V + i]) / r2;
            p[i] -= (float)((1 - fac) * rad * f0); p[V + i] -= (float)((1 - fac) * rad * f1);
        }
    }
}

/* evolve nsteps of KDK leapfrog from tau with step dtau; scratch has C*V floats. Minkowski if s<=0 (a=1). */
void evolve(int model, int N, float *f, float *p, float *scratch, double tau, double dtau, double s, int nsteps) {
    const long V = (long)N * N * N; const int C = model ? 10 : 2; const long n = C * V;
    #define AOF(t) (s > 0 ? s * (t) : 1.0)
    #define APF(t) ((s > 0 && TAU_C > 0 && (t) < TAU_C) ? s * TAU_C : (s > 0 && TAU_F > 0 && (t) > TAU_F) ? s * TAU_F : AOF(t))
    force(model, N, f, scratch, AOF(tau), APF(tau));
    for (int k = 0; k < nsteps; k++) {
        double t0 = tau + k * dtau, th = t0 + 0.5 * dtau, t1 = t0 + dtau;
        double ih = 1.0 / (AOF(th) * AOF(th));
        #pragma omp parallel for schedule(static)
        for (long i = 0; i < n; i++) { p[i] += (float)(0.5 * dtau * scratch[i]); f[i] += (float)(dtau * ih * p[i]); }
        force(model, N, f, scratch, AOF(t1), APF(t1));
        #pragma omp parallel for schedule(static)
        for (long i = 0; i < n; i++) p[i] += (float)(0.5 * dtau * scratch[i]);
        /* drain applied every DAMP_EVERY steps with the accumulated factor (same rate; step index from t1, so it does not
           depend on how the evolution is chunked) */
        if (DAMP_TAU > 0 && t1 < DAMP_TAU && ((long)llround(t1 / dtau)) % DAMP_EVERY == 0)
            damp_massive(model, N, f, p, exp(-DAMP_G * AOF(t1) * dtau * DAMP_EVERY));
    }
}

/* total energy (physical density * a^4, i.e. conformal-frame energy  sum[ p^2/(2a^4)*a^2 ... ]):
   returns comoving energy E = sum [ (pi^2/a^2)/2 + (grad phi)^2/2 * a^2/a^2... ]  -- we return physical mean densities:
   out[0]=kinetic, out[1]=gradient, out[2]=potential (mean over sites, physical units)                                */
void energy(int model, int N, const float *f, const float *p, double a, double *out) {
    const long V = (long)N * N * N; const int C = model ? 10 : 2;
    double K = 0, G = 0, P = 0;
    #pragma omp parallel for reduction(+:K,G,P) schedule(static)
    for (int z = 0; z < N; z++) {
        int zp = (z + 1) % N; double q[10], g[10];
        for (int y = 0; y < N; y++) {
            int yp = (y + 1) % N;
            for (int x = 0; x < N; x++) {
                int xp = (x + 1) % N; long i = ((long)z * N + y) * N + x;
                for (int c = 0; c < C; c++) {
                    const float *fc = f + c * V; double v = fc[i]; q[c] = v;
                    double pc = p[c * V + i]; K += 0.5 * pc * pc;
                    double dx = fc[((long)z * N + y) * N + xp] - v, dy = fc[((long)z * N + yp) * N + x] - v, dz = fc[((long)zp * N + y) * N + x] - v;
                    G += 0.5 * (dx * dx + dy * dy + dz * dz);
                }
                P += model ? rp3_grad(q, g) : u1_grad(q, g);
            }
        }
    }
    /* phi' = pi/a^2 ; physical kinetic = phi'^2/(2 a^2) = pi^2/(2 a^6); gradient phys = (grad)^2/(2a^2) */
    out[0] = K / V / pow(a, 6); out[1] = G / V / (a * a); out[2] = P / V;
}

/* Jacobi eigen-decomposition of symmetric 4x4; returns eigenvector of largest eigenvalue in v[4], eigenvalues in w */
static void jacobi4(double A[4][4], double w[4], double Vv[4][4]) {
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) Vv[i][j] = (i == j);
    for (int sweep = 0; sweep < 30; sweep++) {
        double off = 0; for (int i = 0; i < 4; i++) for (int j = i + 1; j < 4; j++) off += A[i][j] * A[i][j];
        if (off < 1e-26) break;
        for (int pp = 0; pp < 3; pp++) for (int qq = pp + 1; qq < 4; qq++) {
            if (fabs(A[pp][qq]) < 1e-300) continue;
            double th = 0.5 * (A[qq][qq] - A[pp][pp]) / A[pp][qq];
            double t = (th >= 0 ? 1.0 : -1.0) / (fabs(th) + sqrt(th * th + 1)), c = 1 / sqrt(t * t + 1), sn = t * c;
            for (int k = 0; k < 4; k++) { double akp = A[k][pp], akq = A[k][qq]; A[k][pp] = c * akp - sn * akq; A[k][qq] = sn * akp + c * akq; }
            for (int k = 0; k < 4; k++) { double apk = A[pp][k], aqk = A[qq][k]; A[pp][k] = c * apk - sn * aqk; A[qq][k] = sn * apk + c * aqk; }
            for (int k = 0; k < 4; k++) { double vkp = Vv[k][pp], vkq = Vv[k][qq]; Vv[k][pp] = c * vkp - sn * vkq; Vv[k][qq] = sn * vkp + c * vkq; }
        }
    }
    for (int i = 0; i < 4; i++) w[i] = A[i][i];
}

/* principal direction n (4 comps) of Q at every site; also gap = lambda1 - lambda2 */
void principal(int N, const float *f, float *nvec, float *gap) {
    const long V = (long)N * N * N;
    #pragma omp parallel for schedule(static)
    for (long i = 0; i < V; i++) {
        double q[10], Q[4][4], w[4], Vv[4][4];
        for (int c = 0; c < 10; c++) q[c] = f[c * V + i];
        q2Q(q, Q); jacobi4(Q, w, Vv);
        int b = 0; for (int k = 1; k < 4; k++) if (w[k] > w[b]) b = k;
        double s2 = -1e300; for (int k = 0; k < 4; k++) if (k != b && w[k] > s2) s2 = w[k];
        for (int k = 0; k < 4; k++) nvec[k * V + i] = (float)Vv[k][b];
        if (gap) gap[i] = (float)(w[b] - s2);
    }
}

/* string plaquettes. RP3: Z2 holonomy prod sign(n_i.n_j) = -1.  U(1): winding of phase != 0.
   seed[i] |= 1 for the 4 corners of every pierced plaquette.  returns number of pierced plaquettes.            */
static inline double dotn(const float *n, long V, long i, long j) {
    return (double)n[i] * n[j] + (double)n[V + i] * n[V + j] + (double)n[2 * V + i] * n[2 * V + j] + (double)n[3 * V + i] * n[3 * V + j];
}
static inline double wrap(double d) { while (d > M_PI) d -= 2 * M_PI; while (d < -M_PI) d += 2 * M_PI; return d; }
long plaquettes(int model, int N, const float *n_or_phi, uint8_t *seed) {
    const long V = (long)N * N * N; long cnt = 0;
    memset(seed, 0, V);
    double *th = NULL;                     /* U(1): phase per site, computed once (U4b speed-up; identical values) */
    if (!model) {
        th = (double *)malloc(sizeof(double) * V);
        #pragma omp parallel for schedule(static)
        for (long i = 0; i < V; i++) th[i] = atan2(n_or_phi[V + i], n_or_phi[i]);
    }
    #pragma omp parallel for reduction(+:cnt) schedule(static)
    for (int z = 0; z < N; z++) {
        int zp = (z + 1) % N;
        for (int y = 0; y < N; y++) {
            int yp = (y + 1) % N;
            for (int x = 0; x < N; x++) {
                int xp = (x + 1) % N;
                long i0 = ((long)z * N + y) * N + x;
                long ix = ((long)z * N + y) * N + xp, iy = ((long)z * N + yp) * N + x, iz = ((long)zp * N + y) * N + x;
                long ixy = ((long)z * N + yp) * N + xp, iyz = ((long)zp * N + yp) * N + x, izx = ((long)zp * N + y) * N + xp;
                long sq[3][4] = {{i0, ix, ixy, iy}, {i0, iy, iyz, iz}, {i0, iz, izx, ix}};
                for (int o = 0; o < 3; o++) {
                    int pierced = 0;
                    if (model) {
                        double sg = 1;
                        for (int k = 0; k < 4; k++) sg *= (dotn(n_or_phi, V, sq[o][k], sq[o][(k + 1) % 4]) >= 0 ? 1 : -1);
                        pierced = sg < 0;
                    } else {
                        double w = 0;
                        for (int k = 0; k < 4; k++) {
                            long a = sq[o][k], b = sq[o][(k + 1) % 4];
                            w += wrap(th[b] - th[a]);
                        }
                        pierced = fabs(w) > M_PI;
                    }
                    if (pierced) {
                        cnt++;
                        for (int k = 0; k < 4; k++) seed[sq[o][k]] = 1;   /* benign race: all writers write 1 */
                    }
                }
            }
        }
    }
    free(th);
    return cnt;
}

/* Goldstone and massive kinetic density per site (physical units): KG = |w|^2 with w = (I-P)Q'n  (=1/2 Tr(Q'_G^2)),
   KM = 1/2 Tr(Q'^2) - KG.  For U(1): KG = (phi1 phi2' - phi2 phi1')^2 / (2|phi|^2), KM = (phi.phi')^2/(2|phi|^2).
   phi' = pi / a^2; physical kinetic = phi'^2 / (2 a^2).                                                           */
void kinetic_split(int model, int N, const float *f, const float *p, const float *nvec, double a, float *KG, float *KM) {
    const long V = (long)N * N * N; const double a6 = pow(a, 6);
    #pragma omp parallel for schedule(static)
    for (long i = 0; i < V; i++) {
        if (model) {
            double qd[10], Qd[4][4], nn[4], u[4], tot = 0;
            for (int c = 0; c < 10; c++) { qd[c] = p[c * V + i]; tot += 0.5 * qd[c] * qd[c]; }
            q2Q(qd, Qd);
            for (int k = 0; k < 4; k++) nn[k] = nvec[k * V + i];
            double nQn = 0;
            for (int r = 0; r < 4; r++) { u[r] = 0; for (int c = 0; c < 4; c++) u[r] += Qd[r][c] * nn[c]; nQn += nn[r] * u[r]; }
            double uu = u[0] * u[0] + u[1] * u[1] + u[2] * u[2] + u[3] * u[3];
            double kg = uu - nQn * nQn;
            KG[i] = (float)(kg / a6); KM[i] = (float)((tot - kg) / a6);
        } else {
            double f0 = f[i], f1 = f[V + i], p0 = p[i], p1 = p[V + i], r2 = f0 * f0 + f1 * f1 + 1e-30;
            double ang = f0 * p1 - f1 * p0, rad = f0 * p0 + f1 * p1;
            KG[i] = (float)(0.5 * ang * ang / r2 / a6); KM[i] = (float)(0.5 * rad * rad / r2 / a6);
        }
    }
}

/* tangent (Goldstone) velocity components for spectra: RP3 -> 10 comps of Q'_G = w n^T + n w^T (Frobenius basis);
   U(1) -> 1 comp  (phi1 phi2' - phi2 phi1')/|phi|.  Physical normalisation: divided by a^3 (so |out|^2/2 = KG-like). */
void tangent_velocity(int model, int N, const float *f, const float *p, const float *nvec, double a, float *out) {
    const long V = (long)N * N * N; const double a3 = a * a * a;
    #pragma omp parallel for schedule(static)
    for (long i = 0; i < V; i++) {
        if (model) {
            double qd[10], Qd[4][4], nn[4], u[4], w[4], G[4][4];
            for (int c = 0; c < 10; c++) qd[c] = p[c * V + i];
            q2Q(qd, Qd);
            for (int k = 0; k < 4; k++) nn[k] = nvec[k * V + i];
            double nQn = 0;
            for (int r = 0; r < 4; r++) { u[r] = 0; for (int c = 0; c < 4; c++) u[r] += Qd[r][c] * nn[c]; nQn += nn[r] * u[r]; }
            for (int r = 0; r < 4; r++) w[r] = u[r] - nQn * nn[r];
            for (int r = 0; r < 4; r++) for (int c = 0; c < 4; c++) G[r][c] = w[r] * nn[c] + nn[r] * w[c];
            for (int k = 0; k < 4; k++) out[k * V + i] = (float)(G[k][k] / a3);
            for (int k = 0; k < 6; k++) out[(4 + k) * V + i] = (float)(G[OI[k]][OJ[k]] * 2 * ISQ2 / a3);
        } else {
            double f0 = f[i], f1 = f[V + i], r = sqrt(f0 * f0 + f1 * f1 + 1e-30);
            out[i] = (float)((f0 * p[V + i] - f1 * p[i]) / r / a3);
        }
    }
}

/* Exact periodic Euclidean distance to the nearest seed site (seed != 0), separable (Felzenszwalb & Huttenlocher 2012),
   OpenMP over lines.  Squared distances are integers < 2^24 for N <= 2300, hence exact in float32; the output is
   (float)sqrt((double)d2), the same rounding as scipy's distance_transform_edt cast to float32.  U4b speed-up. */
#define EDT_INF 1e30f
static void edt_line(const float *f, int N, float *out, int *v, double *zb) {
    /* lower envelope of parabolas (q - j)^2 + f[j mod N] over the periodic images j in [-N, 2N); evaluated at q in [0, N) */
    int k = -1;
    for (int j = -N; j < 2 * N; j++) {
        float fj = f[((j % N) + N) % N];
        if (fj >= EDT_INF) continue;
        double fjd = fj;
        while (k >= 0) {
            int vk = v[k]; double fvk = f[((vk % N) + N) % N];
            double sx = ((fjd + (double)j * j) - (fvk + (double)vk * vk)) / (2.0 * (j - vk));
            if (sx <= zb[k]) k--; else { k++; v[k] = j; zb[k] = sx; break; }
        }
        if (k < 0) { k = 0; v[0] = j; zb[0] = -1e300; }
    }
    if (k < 0) { for (int q = 0; q < N; q++) out[q] = EDT_INF; return; }
    zb[k + 1] = 1e300;
    int r = 0;
    for (int q = 0; q < N; q++) {
        while (zb[r + 1] < q) r++;
        double d = (double)(q - v[r]); out[q] = (float)(d * d + f[((v[r] % N) + N) % N]);
    }
}
void edt_periodic(int N, const uint8_t *seed, float *d) {
    const long N2 = (long)N * N, V = N2 * N;
    #pragma omp parallel
    {
        float *f = (float *)malloc(sizeof(float) * N), *o = (float *)malloc(sizeof(float) * N);
        int *v = (int *)malloc(sizeof(int) * 3 * N + 4); double *zb = (double *)malloc(sizeof(double) * (3 * N + 4));
        #pragma omp for schedule(static)
        for (long l = 0; l < N2; l++) {                 /* x lines: init 0 / INF */
            const uint8_t *s = seed + l * N; float *dl = d + l * N;
            for (int x = 0; x < N; x++) f[x] = s[x] ? 0.0f : EDT_INF;
            edt_line(f, N, dl, v, zb);
        }
        #pragma omp for schedule(static)
        for (long l = 0; l < N2; l++) {                 /* y lines: (z, x) */
            long z = l / N, x = l % N, base = z * N2 + x;
            for (int y = 0; y < N; y++) f[y] = d[base + (long)y * N];
            edt_line(f, N, o, v, zb);
            for (int y = 0; y < N; y++) d[base + (long)y * N] = o[y];
        }
        #pragma omp for schedule(static)
        for (long l = 0; l < N2; l++) {                 /* z lines: (y, x) */
            for (int z = 0; z < N; z++) f[z] = d[l + (long)z * N2];
            edt_line(f, N, o, v, zb);
            for (int z = 0; z < N; z++) d[l + (long)z * N2] = o[z];
        }
        #pragma omp for schedule(static)
        for (long i = 0; i < V; i++) d[i] = d[i] >= EDT_INF ? 1e9f : (float)sqrt((double)d[i]);
        free(f); free(o); free(v); free(zb);
    }
}

/* U10a (ESPEC_U10a): local massive -> Goldstone transfer X and Goldstone energy (kinetic, gradient), physical units, summed
   over sites with keep[i] != 0 (all sites if keep == NULL).  Central differences.  out = {sum X, sum e_kin, sum e_grad, count}.
   U(1): e_G = r^2 (thdot^2 + |grad th|^2)/2,  X = 2 (rdot/r)(e_grad - e_kin)            (exact identity, minus a divergence)
   RP3 : per plane (e1, ek), e_G = sum_k u_k^2 + sum_i g_ki^2,  u_k = ek.Qdot.e1,  g_ki = ek.d_iQ.e1,
         X = sum_k 2 (ldot_1 - ldot_k)/(l_1 - l_k) (sum_i g_ki^2 - u_k^2)               (abelian analogue; checked by V5) */
static inline void loadQ(const float *f, long V, long j, double Q[4][4]) { double q[10]; for (int c = 0; c < 10; c++) q[c] = f[c * V + j]; q2Q(q, Q); }
void exchange(int model, int N, const float *f, const float *p, double a, const uint8_t *keep, double *out) {
    const long V = (long)N * N * N; const double a3 = a * a * a;
    double SX = 0, SK = 0, SG = 0, CNT = 0;
    #pragma omp parallel for reduction(+:SX,SK,SG,CNT) schedule(static)
    for (int z = 0; z < N; z++) {
        int zp = (z + 1) % N, zm = (z + N - 1) % N;
        for (int y = 0; y < N; y++) {
            int yp = (y + 1) % N, ym = (y + N - 1) % N;
            for (int x = 0; x < N; x++) {
                long i = ((long)z * N + y) * N + x;
                if (keep && !keep[i]) continue;
                int xp = (x + 1) % N, xm = (x + N - 1) % N;
                long nb[3][2] = {{((long)z * N + y) * N + xp, ((long)z * N + y) * N + xm}, {((long)z * N + yp) * N + x, ((long)z * N + ym) * N + x},
                                 {((long)zp * N + y) * N + x, ((long)zm * N + y) * N + x}};
                if (model) {
                    double Q[4][4], A[4][4], w[4], Vv[4][4], Qd[4][4], dQ[3][4][4], Qp[4][4], Qm[4][4], pc[10];
                    loadQ(f, V, i, Q);
                    for (int r = 0; r < 4; r++) for (int c = 0; c < 4; c++) A[r][c] = Q[r][c];
                    jacobi4(A, w, Vv);
                    for (int c = 0; c < 10; c++) pc[c] = p[c * V + i] / a3;
                    q2Q(pc, Qd);
                    for (int d = 0; d < 3; d++) {
                        loadQ(f, V, nb[d][0], Qp); loadQ(f, V, nb[d][1], Qm);
                        for (int r = 0; r < 4; r++) for (int c = 0; c < 4; c++) dQ[d][r][c] = 0.5 * (Qp[r][c] - Qm[r][c]) / a;
                    }
                    int b = 0; for (int k = 1; k < 4; k++) if (w[k] > w[b]) b = k;
                    double e1[4]; for (int r = 0; r < 4; r++) e1[r] = Vv[r][b];
                    double ld1 = 0; for (int r = 0; r < 4; r++) for (int c = 0; c < 4; c++) ld1 += e1[r] * Qd[r][c] * e1[c];
                    double X = 0, K = 0, G = 0;
                    for (int k = 0; k < 4; k++) {
                        if (k == b) continue;
                        double ek[4]; for (int r = 0; r < 4; r++) ek[r] = Vv[r][k];
                        double uk = 0, ldk = 0, gg = 0;
                        for (int r = 0; r < 4; r++) for (int c = 0; c < 4; c++) { uk += ek[r] * Qd[r][c] * e1[c]; ldk += ek[r] * Qd[r][c] * ek[c]; }
                        for (int d = 0; d < 3; d++) { double g = 0; for (int r = 0; r < 4; r++) for (int c = 0; c < 4; c++) g += ek[r] * dQ[d][r][c] * e1[c]; gg += g * g; }
                        K += uk * uk; G += gg;
                        double den = w[b] - w[k];
                        if (den > 1e-3) X += 2.0 * (ld1 - ldk) / den * (gg - uk * uk);
                    }
                    SX += X; SK += K; SG += G; CNT += 1;
                } else {
                    double f0 = f[i], f1 = f[V + i], r2 = f0 * f0 + f1 * f1 + 1e-30;
                    double p0 = p[i] / a3, p1 = p[V + i] / a3;
                    double ang = f0 * p1 - f1 * p0, K = 0.5 * ang * ang / r2, G = 0;
                    for (int d = 0; d < 3; d++) {
                        double d0 = 0.5 * (f[nb[d][0]] - f[nb[d][1]]) / a, d1 = 0.5 * (f[V + nb[d][0]] - f[V + nb[d][1]]) / a;
                        double cr = f0 * d1 - f1 * d0; G += 0.5 * cr * cr / r2;
                    }
                    double rdr = (f0 * p0 + f1 * p1) / r2;
                    SX += 2.0 * rdr * (G - K); SK += K; SG += G; CNT += 1;
                }
            }
        }
    }
    out[0] = SX; out[1] = SK; out[2] = SG; out[3] = CNT;
}

/* U10a: massive (normal) velocity components for spectra, physical normalisation /a^3 (|out|^2/2 = KM-like).
   RP3 -> 10 comps of Q'_M = Q' - (w n^T + n w^T);  U(1) -> 1 comp (phi.phi')/|phi|. */
void normal_velocity(int model, int N, const float *f, const float *p, const float *nvec, double a, float *out) {
    const long V = (long)N * N * N; const double a3 = a * a * a;
    #pragma omp parallel for schedule(static)
    for (long i = 0; i < V; i++) {
        if (model) {
            double qd[10], Qd[4][4], nn[4], u[4], w[4];
            for (int c = 0; c < 10; c++) qd[c] = p[c * V + i];
            q2Q(qd, Qd);
            for (int k = 0; k < 4; k++) nn[k] = nvec[k * V + i];
            double nQn = 0;
            for (int r = 0; r < 4; r++) { u[r] = 0; for (int c = 0; c < 4; c++) u[r] += Qd[r][c] * nn[c]; nQn += nn[r] * u[r]; }
            for (int r = 0; r < 4; r++) w[r] = u[r] - nQn * nn[r];
            double M[4][4];
            for (int r = 0; r < 4; r++) for (int c = 0; c < 4; c++) M[r][c] = Qd[r][c] - (w[r] * nn[c] + nn[r] * w[c]);
            for (int k = 0; k < 4; k++) out[k * V + i] = (float)(M[k][k] / a3);
            for (int k = 0; k < 6; k++) out[(4 + k) * V + i] = (float)(M[OI[k]][OJ[k]] * 2 * ISQ2 / a3);
        } else {
            double f0 = f[i], f1 = f[V + i], r = sqrt(f0 * f0 + f1 * f1 + 1e-30);
            out[i] = (float)((f0 * p[i] + f1 * p[V + i]) / r / a3);
        }
    }
}

int nthreads(void) { return omp_get_max_threads(); }

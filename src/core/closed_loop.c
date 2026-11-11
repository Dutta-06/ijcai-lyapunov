/* The closed loop of Proposition D‴.49 in C, called from closed_loop.py through ctypes.
 *
 * The scalar plant ẋ = f(x) + u, u = ẋ_m − k_e e − f̂(x̃), e = x̃ − x_m, x̃ = x + n(t) the measured state (n a sum of tones, zero without
 * noise), x_m = x0 + a sin(ωt) + a₂ sin(ω₂t); f = Σ_k A_k exp(−(x − c*_k)²/(2s_k²)); f̂ = Σ_i w_iφ_i, φ_i = exp(−(x − c_i)²/(2σ²)). The
 * kernels follow the paradigm's free law on their curves from α = 0 with the masses b0_i — β_i = b0_i + 2c sinh²(α_i/a_max), so that
 * Q_i = β_i − c cosh(2α_i/a_max) is kept exactly —, w_i = β_ia_i/m, a_i = a_max tanh(α_i/a_max):
 *     α̇_i = Gβ_ia′_i e φ_i(x_r),   ċ_i = v_i e ∂_cφ_i(x_r),   v_i = Gβ_ia_i (metric 0, the paradigm's), Ga_i (1, Wasserstein–Fisher–Rao),
 *     0 (2, fixed centres); x_r = x̃ (regx = 0) or x_m (regx = 1) the regressor's evaluation point.
 * Proposition D‴.46's leak toward the start (η_p = 0): the drive e φ_i(x_r) − λη_i, η_i = (α_i + a_maxκ_i arctan(κ_i tanh(α_i/a_max)))/(2Gc)
 * the mirror coordinate on the curve from α = 0 (κ_i² = (2c − b0_i)/b0_i); Q is kept (λ = 0: no leak).
 *
 * Proposition D‴.57's dead zone (gate > 0): a filtered error energy ż = (e² − z)/T_f joins the state, and the adaptation's drive (not the
 * leak) of kernel i is multiplied by g_i: g = S((z − z_lo)/(z_hi − z_lo)) for every kernel (gate 1, the energy dead zone), or
 * g_i = s_i + (1 − s_i)g with s_i = S((|w_i| − θ/2)/(θ/2)) — the formed kernels driven, the others gated (gate 2, the formation dead zone);
 * S(u) = u²(3 − 2u) on [0, 1], 0 below, 1 above (z_lo, z_hi the squared thresholds).
 *
 * Two options for referee round 50's tests of the mechanism: etrue = 1 drives the adaptation (not the control, not the gate) with the
 * true error x − x_m instead of the measured one; cbox < ∞ bounds the centres to [−cbox, cbox] (an outward velocity at the bound is
 * zeroed).
 *
 * Dormand–Prince 5(4) with the usual error control; the state is hit exactly at each output time. Returns the number of accepted steps,
 * or −1 if max_steps was exceeded, −2 if the step fell below 1e-14. Also returns the time each kernel first exceeds |w| > th
 * (linear interpolation within the step; NaN if never) and, in emax: max|e| over [0, per] and [T − per, T] (the measured error), the
 * true error x − x_m's max and its integrated square over [T − per, T] (trapezoid over the accepted steps), and the state's least and
 * largest value over the run (at the accepted steps).
 */
#include <math.h>
#include <stddef.h>

typedef struct {
    int d, nb, metric, regx, nn, gate, etrue;
    double m, G, cc, amax, sig, ke, lam, Tf, zlo, zhi, theta, cbox;
    double x0, amp, om, a2, om2;
    const double *bc, *ba, *bs, *b0, *nw, *nph, *nam;
} Par;

static double eta_of(const Par *p, double al, double b0)
{
    /* the mirror coordinate η(α) = ∫₀^α dα'/(G(b0 sech² + 2c tanh²)) in closed form (b0 < 2c) */
    double kp = sqrt((2.0 * p->cc - b0) / b0);
    return (al + p->amax * kp * atan(kp * tanh(al / p->amax))) / (2.0 * p->G * p->cc);
}

static double smooth01(double u)
{
    return u <= 0.0 ? 0.0 : (u >= 1.0 ? 1.0 : u * u * (3.0 - 2.0 * u));
}

static void rhs(const Par *p, double t, const double *y, double *dy, double *e_out)
{
    int i, k, d = p->d;
    double g = 1.0;
    double x = y[0];
    double xm = p->x0 + p->amp * sin(p->om * t) + p->a2 * sin(p->om2 * t);
    double dxm = p->amp * p->om * cos(p->om * t) + p->a2 * p->om2 * cos(p->om2 * t);
    double n = 0.0, fs = 0.0, fh = 0.0, xt, e, xr, s2 = p->sig * p->sig;
    for (k = 0; k < p->nn; k++) n += p->nam[k] * cos(p->nw[k] * t + p->nph[k]);
    xt = x + n;
    e = xt - xm;
    double ea = p->etrue ? x - xm : e;                       /* the adaptation's drive */
    xr = p->regx ? xm : xt;
    for (k = 0; k < p->nb; k++) {
        double u = x - p->bc[k];
        fs += p->ba[k] * exp(-u * u / (2.0 * p->bs[k] * p->bs[k]));
    }
    if (p->gate) {
        double z = y[1 + 2 * d];
        g = smooth01((z - p->zlo) / (p->zhi - p->zlo));
        dy[1 + 2 * d] = (e * e - z) / p->Tf;
    }
    for (i = 0; i < d; i++) {
        double al = y[1 + i], c = y[1 + d + i];
        double th = tanh(al / p->amax), sh = sinh(al / p->amax);
        double beta = p->b0[i] + 2.0 * p->cc * sh * sh, a = p->amax * th, da = 1.0 - th * th;
        double w = beta * a / p->m;
        double uc = xt - c, ur = xr - c;
        double phc = exp(-uc * uc / (2.0 * s2));
        double phr = p->regx ? exp(-ur * ur / (2.0 * s2)) : phc;
        double v = (p->metric == 0) ? p->G * beta * a : ((p->metric == 1) ? p->G * a : 0.0);
        double gi = 1.0;
        if (p->gate == 1) gi = g;
        else if (p->gate == 2) {
            double si = smooth01((fabs(w) - 0.5 * p->theta) / (0.5 * p->theta));
            gi = si + (1.0 - si) * g;
        }
        fh += w * phc;
        dy[1 + i] = gi * p->G * beta * da * ea * phr;
        if (p->lam > 0.0) {
            /* the leak −λη on the mirror coordinate: α̇ −= λ η (dα/dη), dα/dη = G(b0 sech² + 2c tanh²) */
            double sch = 1.0 / cosh(al / p->amax);
            double mob = p->G * (p->b0[i] * sch * sch + 2.0 * p->cc * th * th);
            dy[1 + i] -= p->lam * eta_of(p, al, p->b0[i]) * mob;
        }
        {
            double dc = gi * v * ea * phr * ur / s2;
            if ((c >= p->cbox && dc > 0.0) || (c <= -p->cbox && dc < 0.0)) dc = 0.0;
            dy[1 + d + i] = dc;
        }
    }
    dy[0] = fs - fh + dxm - p->ke * e;
    if (e_out) *e_out = e;
}

static double wabs(const Par *p, int i, double al)
{
    double th = tanh(al / p->amax), sh = sinh(al / p->amax);
    return fabs((p->b0[i] + 2.0 * p->cc * sh * sh) * p->amax * th / p->m);
}

/* ip = [d, nb, metric, regx, nn, nout, gate, etrue]; dp = [m, G, cc, amax, sig, ke, lam, x0, amp, om, a2, om2, rtol, atol, hmax, th,
 * per, max_steps, Tf, z_lo, z_hi, theta, cbox]; the state (x, α, c) and, with a gate, z */
long simulate(const int *ip, const double *dp, const double *bc, const double *ba, const double *bs, const double *b0, const double *nw,
              const double *nph, const double *nam, const double *y0, const double *tout, double *yout, double *ent, double *emax,
              double *work)
{
    Par p;
    int d = ip[0], nout = ip[5], n = 1 + 2 * d + (ip[6] > 0), i, j, io = 0;
    double rtol = dp[12], atol = dp[13], hmax = dp[14], th = dp[15], per = dp[16];
    long max_steps = (long)dp[17], steps = 0;
    double *y = work, *yn = work + n, *k1 = work + 2 * n, *k2 = work + 3 * n, *k3 = work + 4 * n, *k4 = work + 5 * n, *k5 = work + 6 * n,
           *k6 = work + 7 * n, *k7 = work + 8 * n, *yt = work + 9 * n;
    double t, h, T = tout[nout - 1], e;
    static const double a21 = 1.0 / 5, a31 = 3.0 / 40, a32 = 9.0 / 40, a41 = 44.0 / 45, a42 = -56.0 / 15, a43 = 32.0 / 9,
        a51 = 19372.0 / 6561, a52 = -25360.0 / 2187, a53 = 64448.0 / 6561, a54 = -212.0 / 729, a61 = 9017.0 / 3168, a62 = -355.0 / 33,
        a63 = 46732.0 / 5247, a64 = 49.0 / 176, a65 = -5103.0 / 18656, b1 = 35.0 / 384, b3 = 500.0 / 1113, b4 = 125.0 / 192,
        b5 = -2187.0 / 6784, b6 = 11.0 / 84, e1 = 71.0 / 57600, e3 = -71.0 / 16695, e4 = 71.0 / 1920, e5 = -17253.0 / 339200,
        e6 = 22.0 / 525, e7 = -1.0 / 40;
    p.d = d; p.nb = ip[1]; p.metric = ip[2]; p.regx = ip[3]; p.nn = ip[4]; p.gate = ip[6]; p.etrue = ip[7];
    p.Tf = dp[18]; p.zlo = dp[19]; p.zhi = dp[20]; p.theta = dp[21]; p.cbox = dp[22];
    p.m = dp[0]; p.G = dp[1]; p.cc = dp[2]; p.amax = dp[3]; p.sig = dp[4]; p.ke = dp[5]; p.lam = dp[6];
    p.x0 = dp[7]; p.amp = dp[8]; p.om = dp[9]; p.a2 = dp[10]; p.om2 = dp[11];
    p.bc = bc; p.ba = ba; p.bs = bs; p.b0 = b0; p.nw = nw; p.nph = nph; p.nam = nam;
    for (i = 0; i < n; i++) y[i] = y0[i];
    for (i = 0; i < d; i++) ent[i] = NAN;
    emax[0] = 0.0; emax[1] = 0.0; emax[2] = 0.0; emax[3] = 0.0; emax[4] = y[0]; emax[5] = y[0];
    t = tout[0];
    while (io < nout && tout[io] <= t) {
        for (i = 0; i < n; i++) yout[io * n + i] = y[i];
        io++;
    }
    rhs(&p, t, y, k1, &e);
    if (t <= per) emax[0] = fabs(e);
    if (t >= T - per) emax[1] = fabs(e);
    h = hmax > 0 ? fmin(hmax, 1e-3) : 1e-3;
    while (io < nout) {
        double tn, err = 0.0, fac;
        int hit = 0;
        if (steps >= max_steps) return -1;
        if (hmax > 0 && h > hmax) h = hmax;
        if (t + h >= tout[io]) { h = tout[io] - t; hit = 1; }
        for (i = 0; i < n; i++) yt[i] = y[i] + h * a21 * k1[i];
        rhs(&p, t + h / 5, yt, k2, NULL);
        for (i = 0; i < n; i++) yt[i] = y[i] + h * (a31 * k1[i] + a32 * k2[i]);
        rhs(&p, t + 3 * h / 10, yt, k3, NULL);
        for (i = 0; i < n; i++) yt[i] = y[i] + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
        rhs(&p, t + 4 * h / 5, yt, k4, NULL);
        for (i = 0; i < n; i++) yt[i] = y[i] + h * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
        rhs(&p, t + 8 * h / 9, yt, k5, NULL);
        for (i = 0; i < n; i++) yt[i] = y[i] + h * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
        rhs(&p, t + h, yt, k6, NULL);
        for (i = 0; i < n; i++) yn[i] = y[i] + h * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]);
        tn = hit ? tout[io] : t + h;
        rhs(&p, tn, yn, k7, &e);
        for (i = 0; i < n; i++) {
            double sc = atol + rtol * fmax(fabs(y[i]), fabs(yn[i]));
            double ei = h * (e1 * k1[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * k7[i]) / sc;
            err += ei * ei;
        }
        err = sqrt(err / n);
        if (!(err <= 1.0)) {
            fac = (err != err) ? 0.2 : fmax(0.2, 0.9 * pow(err, -0.2));
            h *= fac;
            if (h < 1e-14) return -2;
            continue;
        }
        /* accepted */
        steps++;
        for (j = 0; j < d; j++) {
            if (ent[j] != ent[j]) {
                double wo = wabs(&p, j, y[1 + j]), wn = wabs(&p, j, yn[1 + j]);
                if (wn > th) ent[j] = (wn > wo) ? t + (tn - t) * (th - wo) / (wn - wo) : tn;
            }
        }
        if (tn <= per && fabs(e) > emax[0]) emax[0] = fabs(e);
        if (tn >= T - per && fabs(e) > emax[1]) emax[1] = fabs(e);
        {
            double xmo = p.x0 + p.amp * sin(p.om * t) + p.a2 * sin(p.om2 * t);
            double xmn = p.x0 + p.amp * sin(p.om * tn) + p.a2 * sin(p.om2 * tn);
            double eo = y[0] - xmo, en = yn[0] - xmn;
            if (tn >= T - per) {
                double ta = fmax(t, T - per);
                if (fabs(en) > emax[2]) emax[2] = fabs(en);
                if (tn > ta) emax[3] += 0.5 * (tn - ta) * (eo * eo + en * en);
            }
            if (yn[0] < emax[4]) emax[4] = yn[0];
            if (yn[0] > emax[5]) emax[5] = yn[0];
        }
        for (i = 0; i < n; i++) { y[i] = yn[i]; k1[i] = k7[i]; }
        fac = (err > 0) ? fmin(5.0, fmax(0.2, 0.9 * pow(err, -0.2))) : 5.0;
        if (hit) {
            for (i = 0; i < n; i++) yout[io * n + i] = y[i];
            io++;
            t = tn;
            /* keep the step that was cut short to hit the output time from shrinking the next one */
            h = fmax(h, h * fac);
        } else {
            t = tn;
            h *= fac;
        }
    }
    return steps;
}

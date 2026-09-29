# push2: attempts to go beyond alpha ~2.23 (fork report, 2026-09-29)

## Best rigorous result: alpha >= 223/100 (exact certificate, same lemma set as App. app:elem)
`fastcert.py 179/100 179/100 737/500 400/223` (exact Fractions; same hypotheses as
research/cert/certificate.py minus its superfluous 1-2a <= XS; analytic best-XL search):
- a0 = a = 179/100, b = 737/500: 2a+3b = 8.002 > 8, a >= b, a0 + a = 3.58 >= 2b = 2.948
- E1 = 2(a-a0)+b = 1.474, max_g min N = 0.317 (worst at g = a = 1.79), E(a0) = 1.791
- kappa = E/4 = 0.44775 < 100/223 = 0.44843  => rate 223/100
Paper changes needed (no new lemma): thm:elem for t <= (400/179)L with exponent 0.4478;
thm:rate115 with lambda = 100/223, tau_* = ceil(2.23 log2 K) <= (400/179)L for L >= L0 = 534.
Numerical ceiling of this lemma set: E(a0) = a0 at a0 ~ 1.7905, i.e. alpha ~ 2.234.

## What was tried and why it does not go further (worst sphere r ~ R'^{-0.1}, tangent to the core)
1. sigma_2 radius r2 / height h(nu) of H: refined wedge threshold h(nu) with N(m) = h^2 r^2 r2^2 gives only
   h >= 1/(r r2), trivial in the worst case. The worst configuration is arithmetically realisable:
   low-height H (e.g. x1 = c) with sigma_1(m) tiny, sigma_2(m) ~ R'^2, carrying ~ r r2 ~ R'^{0.9} points;
   proving its rim patch has fewer points is a small-cap equidistribution problem for ternary
   Z[sqrt2]-spheres (Linnik type), beyond elementary methods.
2. Counting bad boxes: per-box worst case R'^{0.3} is consistent with the spectral total R'^2, so no
   contradiction is available; would need the same ternary small-cap input.
3. Better level 2 on round spheres: at the worst g the cells are already 1-D strips (x_S = e_u), so
   high-degree determinants on Y (method A) + Bezout level 3 give no gain; adding a total-turning
   argument for the level-3 curve gains <= 0.011 in the exponent and needs an unproved lemma.
   The 4-point volume Theta(x_L^3 delta / r) in a rod of cross-section delta is sharp.
4. Level 1: sigma_2-localisation costs more than it saves (net +0.33 per unit at the worst case);
   degree-2 auxiliary surfaces (12a + 16b > 44) fail at level 2 because thin tori around the core
   destroy flatness; rigid core-adapted families lose integrality; a < b boxes are worse.
Estimated gain of unproved directions: ~0.01 (turning lemma); ceilings 5/2, 8/3 need new ideas on
the tangent small-sphere case.

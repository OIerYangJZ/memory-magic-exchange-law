"""Zero case (r > R'/2, first g-interval [0, w]) is not part of the z3 claim C.  Check exactly:
E1 + cost(S0 = 1-b, U0 = 1-a, R0 = 1-w, y = 0) <= E1 for all mu in (0, 1/10], w in [0, 1/100]
(cell X_L = 1-b, X_S = 1-a), i.e. no split is needed.  Also: HIGH(g, 0) = E1 for g <= 2/5 - mu."""
import z3, indep as I, z3_indep as Z
mu, w, g = z3.Reals('mu w g')
a0 = z3.Q(8, 5) + mu; a = a0; b = (8 - 2 * a0) / 3
s = z3.Solver()
s.add(mu > 0, mu <= z3.Q(1, 10), w >= 0, w <= z3.Q(1, 100))
s.add(I.cost(Z.O, 1 - b, 1 - a, 1 - w, 0 * mu, 'cand') > 0)
print('zero case: exists (mu, w) with cost > 0 ?', s.check(), '(unsat = zero case closes at E1)')
s = z3.Solver()
M = I.model(Z.O, a0, b + mu / 2, g, 0 * mu, 'cand')
s.add(mu > 0, mu <= z3.Q(1, 10), g > 0, g <= z3.Q(2, 5) - mu, M['HIGH'] > b)
print('HIGH(g, 0) > E1 for some g <= 2/5 - mu ?', s.check(), '(unsat = no split needed there)')

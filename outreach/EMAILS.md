# Email drafts

Attachments for both:
1. `single_norm_tubes.pdf` (4 pp.): the question, what is proven, where the methods stop.
2. `Yang-Li-Deng_memory-magic-exchange-law_v2.pdf` (39 pp.): the paper as it will appear as arXiv:2609.37368v2.
   The note cites its theorem numbers (Thm. 8, Apps. E–H), which are the v2 numbers. v1 on arXiv has different
   numbers and the 17/7 result.

---

## 1. Raphael S. Steiner

**To:** raphael.steiner.academic@gmail.com
**Subject:** A single-norm count of quaternions near a torus (from a quantum compilation problem)

Dear Dr. Steiner,

I work on lower bounds for quantum circuit compilation. In arXiv:2609.37368, joint with Yangyang Li and Xiu-Hao Deng,
we reduce an exchange rate between classical memory and magic (T) gates to a lattice-point count. The count is in the
maximal order of the quaternion algebra (−1,−1) over Q(√2). It counts the elements of one fixed norm, a power of
√2, whose image in the first real embedding lies within ε of a given great circle of S³. There is no condition in the
second embedding.

The Ramanujan bound stops at the square-root barrier. A determinant method, with a case analysis checked exactly by
an SMT solver, takes us exactly to the scale ε ≈ 2^{−2t/5}, where t is the T-count. At that scale it gives one point
per ε-cell along the circle.
- Any power saving there would push our result past 5/2. That means showing that a power fraction of the cells
  along *every* great circle is empty.
- The conjectured volume law would give the optimal value 3.

In the sup-norm literature, the analogous single-norm counts of Hecke returns near a torus seem to be handled by
geometry of numbers. The improvements come from averaging over the norm, which is unavailable here because all
norms are powers of a single prime.

Given your work on the sup-norm problem on S³, and the theta-function method with Khayutin and Nelson, I wondered:
- Do you know of any single-norm result of this kind?
- Or do you see a reason it should be out of reach at present?

The attached four-page note states the question precisely (Section 3), with what we can prove and where our methods
stop. The paper is also attached. Any pointer would be very welcome; I certainly do not expect you to work on it.

With best regards,
Jinze Yang
School of Physics, Xidian University, Xi'an, China

---

## 2. Simon Marshall

**To:** simon.marshall@unimelb.edu.au
**Subject:** Single-norm Hecke returns near a torus — a question from quantum compilation

Dear Professor Marshall,

I work on lower bounds for quantum circuit compilation. In arXiv:2609.37368, joint with Yangyang Li and Xiu-Hao Deng,
an exchange rate between classical memory and magic (T) gates comes down to a count of Hecke returns near a torus,
at a single norm. Concretely, we count the elements of the maximal order of (−1,−1) over Q(√2) with one fixed norm,
a power of √2, whose first-embedding image lies within ε of a given great circle of S³. There is no condition in the
second embedding.

Your paper on geodesic restrictions was very helpful in seeing where we stand.
- Our elementary count is a determinant method in the spirit of the single-norm bound of your Lemma 3.2, with an
  exactly checked case analysis.
- It stops exactly at the scale ε ≈ 2^{−2t/5}, where it gives one point per ε-cell.
- The spectral averaging of your Prop. 5.2 is what would give the volume law. But here every norm is a power of the
  one prime above 2: the local factor has its poles on |X| = 2^{−1/2}, so there is nothing to average over.
- Any power saving at that scale would improve our rate beyond 5/2. That means showing that a power fraction of the
  cells along every great circle is empty. The volume law would give 3.

I wanted to ask:
- Do you know of any result, in any arithmetic setting, that beats geometry of numbers for such single-norm counts
  near a *generic* torus?
- Is there a reason to think this is currently out of reach?
- Or is there a known example where such counts are genuinely large?

The attached four-page note states the question precisely (Section 3), together with what we can prove and a
toy model over Z⁴. The paper is also attached. Any pointer would be greatly appreciated; I certainly do not expect
you to work on it.

With best regards,
Jinze Yang
School of Physics, Xidian University, Xi'an, China

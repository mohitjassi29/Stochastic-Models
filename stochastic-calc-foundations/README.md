# Stochastic Calculus Foundations

From-scratch derivations of Itô's Lemma, the GBM log-transform, and Girsanov's Theorem, with Monte Carlo validation.

- `derivations.md` — full derivations, each step justified rather than asserted
- `validate.py` — numerical checks:
  - GBM log-drift correction (`mu - 0.5*sigma^2` vs naive `mu`)
  - Girsanov martingale property `E^P[Z_T] = 1`
  - Measure-change drift check: `E^Q[S_T] = S0 * e^(rT)`, computed as `E^P[Z_T * S_T]`

Run: `python validate.py`

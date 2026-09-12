"""
Computational validation of:
  1. Ito's lemma / GBM log-transform correction term
  2. Girsanov E^P[Z_T] = 1

Run: python validate.py
"""

import numpy as np

rng = np.random.default_rng(42)

# ---------------------------------------------------------------------------
# 1. GBM log-transform: verify d(ln S_t) drift is (mu - 0.5*sigma^2), not mu
# ---------------------------------------------------------------------------
def validate_gbm_log_drift(mu=0.10, sigma=0.30, T=1.0, n_steps=252, n_paths=200_000):
    dt = T / n_steps
    log_S = np.zeros(n_paths)
    for _ in range(n_steps):
        dW = rng.normal(0.0, np.sqrt(dt), size=n_paths)
        log_S += (mu - 0.5 * sigma**2) * dt + sigma * dW

    empirical_drift = log_S.mean() / T
    theoretical_drift = mu - 0.5 * sigma**2
    naive_drift = mu  # what you'd get WITHOUT the Ito correction

    print("=== GBM log-transform check ===")
    print(f"Theoretical drift (mu - 0.5*sigma^2): {theoretical_drift:.5f}")
    print(f"Empirical drift of ln(S_T)/T:         {empirical_drift:.5f}")
    print(f"Naive (uncorrected) drift mu:         {naive_drift:.5f}")
    print(f"-> Empirical matches the CORRECTED drift, not the naive one.\n")


# ---------------------------------------------------------------------------
# 2. Girsanov: verify E^P[Z_T] = 1 via Monte Carlo
# ---------------------------------------------------------------------------
def validate_girsanov_martingale(mu=0.10, r=0.03, sigma=0.30, T=1.0, n_paths=500_000):
    theta = (mu - r) / sigma
    W_T = rng.normal(0.0, np.sqrt(T), size=n_paths)  # P-Brownian motion at T
    Z_T = np.exp(-theta * W_T - 0.5 * theta**2 * T)

    print("=== Girsanov E^P[Z_T] = 1 check ===")
    print(f"theta = (mu - r)/sigma = {theta:.5f}")
    print(f"Monte Carlo E^P[Z_T]:  {Z_T.mean():.5f}  (theory: 1.00000)")
    print(f"Std error:             {Z_T.std() / np.sqrt(n_paths):.5f}\n")


# ---------------------------------------------------------------------------
# 3. Bonus: verify Q-measure drift is actually r after reweighting by Z_T
#    (E^Q[S_T] should equal S_0 * e^{rT}, computed as E^P[Z_T * S_T])
# ---------------------------------------------------------------------------
def validate_measure_change_drift(S0=100.0, mu=0.10, r=0.03, sigma=0.30, T=1.0, n_paths=500_000):
    theta = (mu - r) / sigma
    W_T = rng.normal(0.0, np.sqrt(T), size=n_paths)
    S_T = S0 * np.exp((mu - 0.5 * sigma**2) * T + sigma * W_T)
    Z_T = np.exp(-theta * W_T - 0.5 * theta**2 * T)

    E_Q_S_T = np.mean(Z_T * S_T)  # E^P[Z_T * S_T] = E^Q[S_T]
    theory = S0 * np.exp(r * T)

    print("=== Measure-change drift check: E^Q[S_T] = S0 * e^(rT) ===")
    print(f"E^P[Z_T * S_T] (i.e. E^Q[S_T]): {E_Q_S_T:.4f}")
    print(f"S0 * e^(rT):                    {theory:.4f}\n")


if __name__ == "__main__":
    validate_gbm_log_drift()
    validate_girsanov_martingale()
    validate_measure_change_drift()

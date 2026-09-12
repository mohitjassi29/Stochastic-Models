# Stochastic Calculus Foundations: Itô's Lemma, GBM Log-Transform, and Girsanov's Theorem

From-scratch derivations underpinning risk-neutral pricing, with computational validation in `validate.py`.

---

## 1. Itô's Lemma

**Setup.** Let $X_t$ follow an Itô process:
$$dX_t = \mu(X_t,t)\,dt + \sigma(X_t,t)\,dW_t$$

For a smooth $f(x,t)$, Taylor-expand $f(X_t + dX_t,\, t+dt)$ to second order:

$$df = f_t\,dt + f_x\,dX_t + \tfrac12\left(f_{tt}(dt)^2 + 2f_{tx}\,dt\,dX_t + f_{xx}(dX_t)^2\right) + \dots$$

**Order-in-$dt$ argument.** Using $dW_t \sim \sqrt{dt}$: terms $f_{tt}(dt)^2$ and $f_{tx}\,dt\,dX_t$ are $O(dt^2)$ or smaller and vanish. The term $f_{xx}(dX_t)^2$ does **not** vanish, because

$$(dX_t)^2 = \mu^2(dt)^2 + 2\mu\sigma\,dt\,dW_t + \sigma^2(dW_t)^2 \;\longrightarrow\; \sigma^2\,dt$$

using the Itô multiplication rule $(dW_t)^2 = dt$, $dt\cdot dW_t = 0$, $(dt)^2=0$ (quadratic variation of Brownian motion).

**Result:**
$$\boxed{df = \left(f_t + \mu f_x + \tfrac12\sigma^2 f_{xx}\right)dt + \sigma f_x\,dW_t}$$

---

## 2. GBM Log-Transform: $d(\ln S_t)$

Let $dS_t = \mu S_t\,dt + \sigma S_t\,dW_t$ and $f(S_t)=\ln S_t$. Then $f_t=0$, $f_x = 1/S_t$, $f_{xx} = -1/S_t^2$.

Applying Itô's lemma with $\mu \to \mu S_t$, $\sigma \to \sigma S_t$:

$$d(\ln S_t) = \left(\mu - \tfrac12\sigma^2\right)dt + \sigma\,dW_t$$

The $-\tfrac12\sigma^2$ Itô correction is why $\mathbb{E}[\ln S_T] \ne \ln \mathbb{E}[S_T]$ — the log and level processes have different drifts.

---

## 3. Girsanov's Theorem

**Goal.** Under $\mathbb{P}$, $dS_t = \mu S_t\,dt + \sigma S_t\,dW_t^{\mathbb{P}}$. Find a measure $\mathbb{Q}$ under which the drift becomes $r$.

**Step 1 — choose the shift.** Define $\theta_t$ via $dW_t^{\mathbb{Q}} = dW_t^{\mathbb{P}} + \theta_t\,dt$. Substituting into $dS_t$ and matching to $rS_t\,dt$:
$$\theta_t = \frac{\mu - r}{\sigma}$$

**Step 2 — find the Radon-Nikodym derivative $Z_t = \frac{d\mathbb{Q}}{d\mathbb{P}}\big|_{\mathcal F_t}$.**

Real goal: $W_t^{\mathbb{Q}}$ is a $\mathbb{Q}$-Brownian motion (Lévy's characterization: $\mathbb{Q}$-martingale + quadratic variation $t$).

Abstract Bayes rule converts this to: $Z_tW_t^{\mathbb{Q}}$ must be a $\mathbb{P}$-martingale.

$Z_t$ itself is a $\mathbb{P}$-martingale automatically (it's $\mathbb{E}^{\mathbb{P}}[Z_T\mid\mathcal F_t]$), so by the Martingale Representation Theorem $dZ_t = \gamma_t\,dW_t^{\mathbb{P}}$ for some $\gamma_t$.

Itô product rule on $Z_tW_t^{\mathbb{Q}}$, using $dW_t^{\mathbb{Q}} = dW_t^{\mathbb{P}}+\theta_t dt$:
$$d(Z_tW_t^{\mathbb{Q}}) = (Z_t\theta_t + \gamma_t)\,dt + (\dots)\,dW_t^{\mathbb{P}}$$

Zero-drift requirement $\Rightarrow \gamma_t = -Z_t\theta_t$, i.e. $dZ_t = -Z_t\theta_t\,dW_t^{\mathbb{P}}$, $Z_0=1$.

**Step 3 — solve the SDE.** Since $Z_t>0$ always, write $Z_t=e^{X_t}$. Itô's lemma on $e^{X_t}$ with $dX_t = a_t\,dt+b_t\,dW_t^{\mathbb{P}}$ gives drift $a_t+\tfrac12 b_t^2$ and diffusion $b_t$. Matching to the target SDE: $b_t=-\theta_t$, and zero drift forces $a_t=-\tfrac12\theta_t^2$. Hence:

$$\boxed{Z_t = \exp\left(-\int_0^t\theta_s\,dW_s^{\mathbb{P}} - \tfrac12\int_0^t\theta_s^2\,ds\right)}$$

**Sanity check: $\mathbb{E}^{\mathbb{P}}[Z_T]=1$.** For constant $\theta$, $Z_T = \exp(-\theta W_T^{\mathbb{P}} - \tfrac12\theta^2T)$. Since $-\theta W_T^{\mathbb{P}}\sim\mathcal N(0,\theta^2T)$, using the Gaussian MGF $\mathbb{E}[e^X]=e^{v/2}$ for $X\sim\mathcal N(0,v)$:
$$\mathbb{E}^{\mathbb{P}}[Z_T] = e^{-\theta^2T/2}\cdot e^{\theta^2T/2} = 1$$

Validated computationally below via Monte Carlo.

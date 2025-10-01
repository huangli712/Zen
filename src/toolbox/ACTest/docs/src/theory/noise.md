Assuming that we already construct the spectral function ``A(\omega)``, it is then straightforward to calculate ``G(\tau)`` or ``G(i\omega_n)`` via Laplace transformation. At this point, the synthetic Green's function ``G`` is exact, containing no numerical noise. We name it ``G_{\text{exact}}``. However, the Green's functions obtained from quantum many-body calculations are often noisy. This is especially the case in finite-temperature quantum Monte Carlo simulations, where numerical noise is inevitable. To make things worse, when the fermionic sign problem is severe, the noise correspondingly increases. To simulate this scenario, the ACTest toolkit can introduce artificial noise into the synthetic Green's function as follows:
```math
\begin{align}
G_{\text{noisy}} = G_{\text{exact}}[1 + \delta \mathcal{N}_{C}(0,1)]
\end{align}
```
where ``\mathcal{N}_{C}(0,1)`` represents complex-valued Gaussian noise with zero mean and unit variance, and the parameter ``\delta`` is used to control the noise level (``0 \le \delta \le 1``).

In Eq.(1), the noise is uncorrelated. Shao *et al.*[^1] proposed a new method to generate correlated noise for imaginary time Green's function:
```math
\begin{align}
G_{\text{noisy}}(\tau_i) = G_{\text{exact}}(\tau_i) +
    \frac{\sum_j e^{-|\tau_j-\tau_i|/\xi}R_j}
         {\sqrt{\sum_j e^{-2|\tau_j-\tau_i|/\xi}}},
\end{align}
```
where the sum is performed assuming periodic boundary conditions, ``\xi`` denotes the correlation length, and ``R_j \sim \mathcal{N}(0,\delta)``. Note that in the case that a normal distribution is used it is possible for ``G_{\text{noisy}}(\tau_i)`` to have a different sign to ``G(\tau_i)``. The ACTest toolkit also supports this feature.

[^1]: See Phys. Rev. X 7, 041072 (2017) for more details.

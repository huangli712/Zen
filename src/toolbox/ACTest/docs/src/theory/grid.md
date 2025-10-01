As mentioned before, Green's functions in quantum many-body physics are usually defined on imaginary time or Matsubara frequency axes. Imaginary time Green's functions ``G(\tau)`` and Matsubara Green's functions ``G(i\omega_n)`` are related via Fourier transformation:

```math
\begin{align}
G(\tau) = \frac{1}{\beta} \sum_n e^{i\omega_n \tau} G(i\omega_n),
\end{align}
```

and

```math
\begin{align}
G(i\omega_n) = \int^{\beta}_0 d\tau\ e^{-i\omega_n \tau} G(\tau),
\end{align}
```

where ``\beta`` denotes inverse temperature (``\beta \equiv 1/T``). The grids for Green's functions are linear. Specifically, ``\tau_i = i \beta/N_{\tau}``, where ``N_{\tau}`` denotes number of time slices and ``i \in [0,N_{\tau}]``. ``\omega_n = (2n+1)\pi/\beta`` for fermions and ``2n\pi/\beta`` for bosons, where ``n \in [0, N]`` and ``N`` means number of Matsubara frequency points.

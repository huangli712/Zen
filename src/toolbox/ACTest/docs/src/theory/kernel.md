## Kernel Function

Just as stated above, the spectral function and the imaginary time or Matsubara Green's function are related with each other by the Laplace transformation:
```math
\begin{align}
G(x) = \int d\omega~K(x,\omega) A(\omega).
\end{align}
```
Here ``K(x,\omega)`` is the so-called kernel function. It plays a key role in this equation. In this section, we would like to introduce the kernels that have been implemented in ACTest.

## Fermionic Kernels

For fermionic Green's function, we have
```math
\begin{align}
G(\tau) = \int^{+\infty}_{-\infty} d\omega
          \frac{e^{-\tau\omega}}{1 + e^{-\beta\omega}} A(\omega),
\end{align}
```
and
```math
\begin{align}
G(i\omega_n) = \int^{+\infty}_{-\infty} d\omega
               \frac{1}{i\omega_n - \omega} A(\omega).
\end{align}
```
The kernels are defined as
```math
\begin{align}
K(\tau,\omega) = \frac{e^{-\tau\omega}}{1 + e^{-\beta\omega}},
\end{align}
```
and
```math
\begin{align}
K(\omega_n,\omega) = \frac{1}{i\omega_n - \omega}.
\end{align}
```
For fermionic systems, ``A(\omega)`` is defined on ``(-\infty,\infty)``. It is causal, i.e., ``A(\omega) \ge 0``.

## Bosonic Kernels

For bosonic system, the spectral function obeys the following constraint:
```math
\begin{align}
\text{sign}(\omega) A(\omega) \ge 0.
\end{align}
```
It is quite convenient to introduce a new variable ``\tilde{A}(\omega)``:
```math
\begin{align}
\tilde{A}(\omega) = \frac{A(\omega)}{\omega}.
\end{align}
```
Clearly, ``\tilde{A}(\omega) \ge 0``. It means that ``\tilde{A}(\omega)`` is positive definite. So, we have
```math
\begin{align}
G(\tau)
= \int^{+\infty}_{-\infty} d\omega
          \frac{e^{-\tau\omega}}{1 - e^{-\beta\omega}}
          A(\omega)
=\int^{+\infty}_{-\infty} d\omega
          \frac{\omega e^{-\tau\omega}}{1 - e^{-\beta\omega}}
          \tilde{A}(\omega),
\end{align}
```
and
```math
\begin{align}
G(i\omega_n) = \int^{+\infty}_{-\infty} d\omega
               \frac{1}{i\omega_n - \omega} A(\omega)
                 = \int^{+\infty}_{-\infty} d\omega
               \frac{\omega}{i\omega_n - \omega} \tilde{A}(\omega).
\end{align}
```
The corresponding kernels read:
```math
\begin{align}
K(\tau,\omega) = \frac{\omega e^{-\tau\omega}}{1 - e^{-\beta\omega}},
\end{align}
```
and
```math
\begin{align}
K(\omega_n,\omega) = \frac{\omega}{i\omega_n - \omega}.
\end{align}
```
Especially,
```math
\begin{align}
K(\tau,\omega = 0) \equiv \frac{1}{\beta},
\end{align}
```
and
```math
\begin{align}
K(\omega_n = 0,\omega = 0) \equiv -1.
\end{align}
```

## Symmetric Bosonic Kernels

This is a special case for the bosonic Green's function with Hermitian bosonic operators. Here, the spectral function ``A(\omega)`` is an odd function. Let us introduce ``\tilde{A}(\omega) = A(\omega)/\omega`` again. Since ``\tilde{A}(\omega)`` is an even function, we can restrict it in ``(0,\infty)``. Now we have
```math
\begin{align}
G(\tau)
= \int^{\infty}_{0} d\omega
              \frac{\omega [e^{-\tau\omega} + e^{-(\beta - \tau)\omega}]}
                   {1 - e^{-\beta\omega}}
              \tilde{A}(\omega),
\end{align}
```
and
```math
\begin{align}
G(i\omega_n) = \int^{\infty}_{0} d\omega
                   \frac{-2\omega^2}{\omega_n^2 + \omega^2} \tilde{A}(\omega).
\end{align}
```
The corresponding kernel functions read:
```math
\begin{align}
K(\tau,\omega) =
    \frac{\omega [e^{-\tau\omega} + e^{-(\beta - \tau)\omega}]}
    {1 - e^{-\beta\omega}},
\end{align}
```
and
```math
\begin{align}
K(\omega_n, \omega) = \frac{-2\omega^2}{\omega_n^2 + \omega^2}.
\end{align}
```
There are two special cases:
```math
\begin{align}
K(\tau,\omega = 0) = \frac{2}{\beta},
\end{align}
```
and
```math
\begin{align}
K(\omega_n = 0,\omega = 0) = -2.
\end{align}
```

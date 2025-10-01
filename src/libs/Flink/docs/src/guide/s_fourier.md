It supports forward and backward fourier transformation. We usually use the following subroutines to perform fourier transformation for the Green's function.

## Type

subroutines

## Source

`src/s_fourier.f90`

## Usage

**(1)** Forward FFT, from $G(\tau) \to G(i\omega_n)$.

```fortran
subroutine s_fft_forward(ntime, tmesh, ftau, mfreq, rmesh, fmat)
```

**(2)** Backward FFT, from $G(i\omega_n) \to G(\tau)$.

```fortran
subroutine s_fft_backward(mfreq, rmesh, fmat, ntime, tmesh, ftau, beta)
```

Here, `ntime` is number of imaginary time points, `tmesh` means $\tau$, `ftau` means $G(\tau)$, `mfreq` is number of Matsubara frequency points, `rmesh` means $\omega_n$, `fmat` means $G(i\omega_n)$, `beta` means inverse temperature of system $\beta$.

!!! note

    The `s_fft_tails()` subroutine is called by the `s_fft_backward()` subroutine internally. DO NOT call it directly!

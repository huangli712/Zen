!!! info

    In the **actest/test** directory, there are twelve typical test cases. Users can modify them to meet their requirements. This section will use an independent example to demonstrate the basic usage of the ACTest toolkit.

Now let us test the performance of the maximum entropy method in the ACFlow toolkit. We just consider four typical scenarios: (1) Fermionic Green's functions, ``A(\omega) > 0``; (2) Fermionic Green's functions, ``A(\omega)`` is non-positive definite; (3) Bosonic Green's functions, ``A(\omega) > 0``; (4) Bosonic Green's functions, ``A(\omega)`` is non-positive definite. For each scenario, we apply the ACTest toolkit to randomly generate 100 spectral functions and corresponding Green's functions. The spectral functions are continuum. They are constructed with Gaussian peaks. The number of possible peaks in each ``A(\omega)`` ranges from 1 to 6. The synthetic Green's functions are on the Matsubara frequency axis. The number of Matsubara frequency points is 10. The noise level is ``10^{-6}``. The **act.toml** file for scenario (3) is shown below.
```toml
[Test]
solver  = "MaxEnt" # Analytic continuation solver in the ACFlow toolkit
ptype   = "gauss"  # Type of peaks
ktype   = "boson"  # Type of kernels
grid    = "bfreq"  # Type of grids
mesh    = "linear" # Type of meshes
ngrid   = 10       # Number of grid points
nmesh   = 801      # Number of mesh points
ntest   = 100      # Number of tests
nbins   = 1        # Number of data bins per test
wmax    = 8.0      # Right boundary of frequency mesh
wmin    = -8.0     # Left boundary of frequency mesh
pmax    = 4.0      # Right boundary of peaks
pmin    = -4.0     # Left boundary of peaks
beta    = 20.0     # Inverse temperature
noise   = 1.0e-6   # Noise level
lcorr   = 0.5      # Correlation length
tcorr   = false    # Is noise correlated in imaginary time axis
fnpd    = false    # Whether the spectrum is non-positive definite
fpbc    = false    # Whether the physical boundary condition is applied
lpeak   = [1,2,3,4,5,6] # Possible number of peaks
```
Once the **act.toml** file is prepared, the following command should be executed in the terminal:
```shell
$ actest/util/acgen.jl act.toml
```
Then, the ACTest toolkit will generate the required data in the present directory. Now there are 100 **image.data.i** and **green.data.i** files, where ``i`` ranges from 1 to 100. These correspond to ``A(\omega)`` and ``G(i\omega_n)`` for bosonic systems. We can further verify the data to make sure ``A(\omega) > 0``. Finally, we have to copy the **act.toml** file to another directory, change the *ktype*, *grid*, and *fnpd* parameters in it, and then execute the above command again to generate testing datasets for the other scenarios.

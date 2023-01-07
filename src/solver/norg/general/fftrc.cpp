#include "fftrc.h"

MatCmplx fftrc(const MatReal& x)
{
    const Int N1 = x.nrows();
    const Int N2 = x.ncols();
    MKL_LONG status = 0;    // execution status
    DFTI_DESCRIPTOR_HANDLE hand = 0;    // DFTI descriptor

    // create DFTI descriptor
    {
        MKL_LONG N[2];  N[0] = N1;  N[1] = N2;
        status = DftiCreateDescriptor(&hand, DFTI_DOUBLE, DFTI_REAL, 2, N);
        if (status != DFTI_NO_ERROR) { DftiFreeDescriptor(&hand); ERR(NAV3(status, N1, N2)); }
    }

    // set configuration: out-of-place
    status = DftiSetValue(hand, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    if (status != DFTI_NO_ERROR) { DftiFreeDescriptor(&hand); ERR(NAV3(status, N1, N2)); }

    // set configuration: CCE storage
    status = DftiSetValue(hand, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    if (status != DFTI_NO_ERROR) { DftiFreeDescriptor(&hand); ERR(NAV3(status, N1, N2)); }

    /* this is not needed for DFTI_COMPLEX_COMPLEX storage */
    /* status = DftiSetValue(hand, DFTI_PACKED_FORMAT, DFTI_CCE_FORMAT); */
    /* if (status != DFTI_NO_ERROR) { DftiFreeDescriptor(&hand); ERR(NAV3(status, N1, N2)); } */

    // consecutive index n = rs[0] + rs[1] * n1 + rs[2] * n2 + rs[3] * n3 + ...;
    // set input strides
    {
        MKL_LONG rs[3]; rs[0] = 0; rs[1] = N2; rs[2] = 1;
        status = DftiSetValue(hand, DFTI_INPUT_STRIDES, rs);
        if (status != DFTI_NO_ERROR) { DftiFreeDescriptor(&hand); ERR(NAV3(status, N1, N2)); }
    }

    // consecutive index k = cs[0] + cs[1] * k1 + cs[2] * k2 + cs[3] * k3 + ...;
    // set output strides
    {
        MKL_LONG cs[3]; cs[0] = 0; cs[1] = N2 / 2 + 1; cs[2] = 1;
        status = DftiSetValue(hand, DFTI_OUTPUT_STRIDES, cs);
        if (status != DFTI_NO_ERROR) { DftiFreeDescriptor(&hand); ERR(NAV3(status, N1, N2)); }
    }

    // commit the descriptor
    status = DftiCommitDescriptor(hand);
    if (status != DFTI_NO_ERROR) { DftiFreeDescriptor(&hand); ERR(NAV3(status, N1, N2)); }

    // z is the result of the forward fourier tranform of x
    MatCmplx z_cce(N1, N2 / 2 + 1);
    // intialize data arrays
    const Real* x_real = x.p();
    Cmplx* z_cmplx = z_cce.p();

    // compute real-to-complex transform
    status = DftiComputeForward(hand, (void*)x_real, z_cmplx);
    if (status != DFTI_NO_ERROR) { DftiFreeDescriptor(&hand); ERR(NAV3(status, N1, N2)); }

    // free DFTI descriptor
    DftiFreeDescriptor(&hand);

    MatCmplx z(N1, N2);
    for_Int(n1, 0, N1) for_Int(n2, 0, N2 / 2 + 1) {
        z[n1][n2] = z_cce[n1][n2];
    }
    for_Int(n1, 0, N1) for_Int(n2, N2 / 2 + 1, N2) {
        z[n1][n2] = cnjg(z_cce[(N1 - n1) % N1][(N2 - n2) % N2]);
    }
    return z;
}

VecCmplx fftrc(const VecReal& x)
{
    const Int N1 = x.size();
    MKL_LONG status = 0;    // execution status
    DFTI_DESCRIPTOR_HANDLE hand = 0;    // DFTI descriptor

    // create DFTI descriptor
    {
        status = DftiCreateDescriptor(&hand, DFTI_DOUBLE, DFTI_REAL, 1, (MKL_LONG)N1);
        if (status != DFTI_NO_ERROR) { DftiFreeDescriptor(&hand); ERR(NAV2(status, N1)); }
    }

    // set configuration: out-of-place
    status = DftiSetValue(hand, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    if (status != DFTI_NO_ERROR) { DftiFreeDescriptor(&hand); ERR(NAV2(status, N1)); }

    // set configuration: CCE storage
    status = DftiSetValue(hand, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    if (status != DFTI_NO_ERROR) { DftiFreeDescriptor(&hand); ERR(NAV2(status, N1)); }

    /* this is not needed for DFTI_COMPLEX_COMPLEX storage */
    /* status = DftiSetValue(hand, DFTI_PACKED_FORMAT, DFTI_CCE_FORMAT); */
    /* if (status != DFTI_NO_ERROR) { DftiFreeDescriptor(&hand); ERR(NAV2(status, N1)); } */

    // commit the descriptor
    status = DftiCommitDescriptor(hand);
    if (status != DFTI_NO_ERROR) { DftiFreeDescriptor(&hand); ERR(NAV2(status, N1)); }

    // z is the result of the forward fourier tranform of x
    VecCmplx z_cce(N1 / 2 + 1);
    // intialize data arrays
    const Real* x_real = x.p();
    Cmplx* z_cmplx = z_cce.p();

    // compute real-to-complex transform
    status = DftiComputeForward(hand, (void*)x_real, z_cmplx);
    if (status != DFTI_NO_ERROR) { DftiFreeDescriptor(&hand); ERR(NAV2(status, N1)); }

    // free DFTI descriptor
    DftiFreeDescriptor(&hand);

    VecCmplx z(N1);
    for_Int(n1, 0, N1 / 2 + 1) {
        z[n1] = z_cce[n1];
    }
    for_Int(n1, N1 / 2 + 1, N1) {
        z[n1] = cnjg(z_cce[(N1 - n1) % N1]);
    }
    return z;
}
